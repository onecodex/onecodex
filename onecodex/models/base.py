import inspect
import json


from typing import ClassVar, Optional, List, Type, TypedDict
from pprint import pformat
from html import escape

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict, Field, model_validator
from pydantic._internal._model_construction import ModelMetaclass

from onecodex.exceptions import MethodNotSupported, OneCodexException
from onecodex.models.helpers import (
    generate_potion_sort_clause,
    generate_potion_keyword_where,
)


# This Pydantic magic callable needs to be removed
_EXCLUDE_PYDANTIC_CALLABLES = {"model_post_init"}

DEFAULT_PAGE_SIZE = 200


class AllowedMethods(TypedDict, total=False):
    update: Type["OneCodexBase"]
    delete: Type["OneCodexBase"]
    save: Type["OneCodexBase"]


class ApiRef(PydanticBaseModel):
    """Represents an unresolved API reference."""

    model_config = ConfigDict(populate_by_alias=True, extra="forbid")

    ref: str = Field(..., alias="$ref")

    @model_validator(mode="before")
    @classmethod
    def validate_ref(cls, data):
        # Special case wherein we autoconvert another model to a $ref and ignore the other fields
        if "$uri" in data and "$ref" not in data:
            return {"$ref": data["$uri"]}
        if "field_uri" in data and "$ref" not in data:
            return {"$ref": data["field_uri"]}
        return data

    def __init__(self, **data):
        super().__init__(**data)
        self._resolved_cache = None

    @property
    def id(self) -> str:
        """Extract the ID from the reference URI."""
        return self.ref.split("/")[-1]

    @property
    def resource_type(self) -> str:
        """Extract the resource type from the reference URI."""
        parts = self.ref.split("/")
        if len(parts) >= 4 and parts[1] == "api" and parts[2] == "v1":
            return parts[3]
        raise OneCodexException(f"Unexpected reference with no resource type: {self.ref}")

    def __repr__(self):
        return f"{self.resource_type.title()}('{self.id}')"

    def _resolve(self):
        """Resolve this reference to an actual object."""
        # Make API call to resolve
        target_class = self._get_target_class()
        if not target_class:
            raise ValueError(f"Cannot determine target class for {self.ref}")

        resp = target_class._client.get(
            f"{target_class._api._base_url}{target_class._resource_path}/{self.id}"
        )
        return target_class.model_validate(resp.json())

    def _get_target_class(self):
        """Determine the target class from the reference URI."""
        from onecodex.models import get_model_class

        if not hasattr(self, "_target_class"):
            self._target_class = get_model_class(self.ref)

        return self._target_class


def _get_dir_fields(
    cls, include_fields=True, include_classmethods=True, include_staticmethods=True
):
    if include_fields:
        fields = list(cls.model_fields.keys()) + list(cls.model_computed_fields.keys()) + ["id"]
    else:
        fields = []
    methods = set()
    mro = inspect.getmro(cls)
    try:
        onecodex_base_index = mro.index(OneCodexBase)
        relevant_classes = mro[: onecodex_base_index + 1]
    except ValueError:
        relevant_classes = [cls]

    for class_ in relevant_classes:
        # Use __dict__ directly to avoid any dir() calls
        for name in class_.__dict__:
            if not name.startswith("_") and name not in _EXCLUDE_PYDANTIC_CALLABLES:
                attr = class_.__dict__[name]

                # Check if it's callable (handles regular methods, staticmethods, classmethods, etc.)
                if callable(attr):
                    methods.add(name)
                if include_classmethods and isinstance(attr, classmethod):
                    methods.add(name)
                if include_staticmethods and isinstance(attr, staticmethod):
                    methods.add(name)

    return fields + sorted(list(methods))


class _DirMeta(ModelMetaclass):
    # A metaclass that overrides the `__dir__` method to return a fixed list of attributes.
    def __dir__(cls):
        return _get_dir_fields(cls, include_fields=False)


class OneCodexBase(PydanticBaseModel, metaclass=_DirMeta):
    _api: ClassVar[Optional["Api"]] = (  # noqa: F821
        None  # Forward reference to avoid circular imports
    )
    _client: ClassVar[Optional["HTTPClient"]] = None  # noqa: F821
    _resource_path: ClassVar[str]  # Default resource path, subclasses override
    _allowed_methods: ClassVar[AllowedMethods] = {}

    model_config = ConfigDict(
        populate_by_alias=True,
        extra="allow",
        arbitrary_types_allowed=True,
    )

    @property
    def id(self) -> Optional[str]:
        if self.field_uri is None:
            return None
        return self.field_uri.split("/")[-1]

    def __init__(self, **data):
        super().__init__(**data)
        self._resolved_cache = {}

    @classmethod
    def _convert_id_to_uri(cls, uuid):
        if not uuid.startswith(cls._resource_path):
            uuid = "{}/{}".format(cls._resource_path, uuid)
        return uuid

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.id}>"

    def _repr_html_(self):
        fields = []
        for k, v in self.__class__.model_fields.items():
            # For backwards compatibility, exclude `field_uri` manually from rich printing
            if k == "field_uri":
                continue
            if v.exclude:
                continue
            fields.append(
                "<tr><td>{}</td><td><code>{}</code></td>".format(
                    escape(k), escape(pformat(self.__dict__[k]))
                )
            )
        fields = "\n".join(fields)
        return """<table>
        <thead>
            <tr>
                <th colspan="2"><code>{cls}({id})</code></th>
            </tr>
        </thead>
        <tbody>{fields}</tbody>
        </table>""".format(
            cls=self.__class__.__name__,
            id=escape(repr(self.id)),
            fields=fields,
        )

    def __dir__(self):
        return _get_dir_fields(self.__class__, include_classmethods=False)

    def __setattr__(self, key, value):
        # Note: I'm not sure this is *really* desirable behavior, but this preserves backwards compatibility for <=0.18.0.
        #       Consider deprecating and removing this behavior and only validating on save.
        if key in self.__class__.model_fields and (
            "update" not in self._allowed_methods
            or key not in self._allowed_methods["update"].model_fields
        ):
            raise MethodNotSupported(f"Cannot set attribute {key} on {self.__class__.__name__}")
        super().__setattr__(key, value)

    def _resolve_and_cache(self, api_ref, cache_key):
        """Resolve an ApiRef and cache the result."""
        if cache_key in self._resolved_cache:
            return self._resolved_cache[cache_key]

        # Check for simple circular reference (a.b.a pattern)
        if hasattr(self, "_immediate_parent") and self._immediate_parent.id == api_ref.id:
            return self._immediate_parent

        # Make API call to resolve, add the parent, cache it
        resolved = api_ref._resolve()
        resolved._immediate_parent = self
        self._resolved_cache[cache_key] = resolved
        return resolved

    def __getattribute__(self, name):
        """Automatically resolve ApiRef objects when accessed."""
        value = super().__getattribute__(name)

        # Only process field attributes - get model_fields from class to avoid deprecation warning
        try:
            # Simple access - just use the class directly
            model_fields = type(self).model_fields
        except AttributeError:
            return value

        if name not in model_fields:
            return value

        # Handle OneCodexBase objects (`?expand` data) - set up parent context
        if isinstance(value, OneCodexBase):
            value._immediate_parent = self
            return value

        # Handle ApiRef objects
        if isinstance(value, ApiRef):
            return self._resolve_and_cache(value, f"{name}:{value.ref}")

        # Handle lists of ApiRef objects
        if isinstance(value, list) and all(isinstance(item, ApiRef) for item in value):
            return [self._resolve_and_cache(item, f"{name}:{item.ref}") for item in value]

        return value

    @classmethod
    def all(cls, sort=None, limit=None) -> List["OneCodexBase"]:
        return cls.where(sort=sort, limit=limit)

    @classmethod
    def get(cls, id: str) -> Optional["OneCodexBase"]:
        resp = cls._client.get(f"{cls._api._base_url}{cls._resource_path}/{id}?expand=all")
        if resp.status_code == 200:
            return cls.model_validate(resp.json())
        return None

    @classmethod
    def where(cls, *filters, **keyword_filters) -> List["OneCodexBase"]:
        """Filter model records of this type from the One Codex server.

        This method works for all OneCodex model types including Samples,
        Classifications, Projects, and Panels.

        Parameters
        ----------
        filters : `object`
            Advanced filters to use (not implemented)
        sort : `str` or `list`, optional
            Sort the results by this field (or list of fields). By default in descending order,
            but if any of the fields start with the special character ^, sort in ascending order.
            For example, sort=['size', '^filename'] will sort by size from largest to smallest and
            filename from A-Z for items with the same size.
        limit : `int`, optional
            Number of records to return. For smaller searches, this can reduce the number of
            network requests made.
        keyword_filters : `str` or `object`
            Filter the results by specific keywords (or filter objects, in advanced usage)

        Examples
        --------
        You can filter objects that are returned locally using a lambda function

        All models implement the same method. To see a complete list of
        models which can be retrieved by ID, see the `models page
        <models.html>`_.
        ::

            # Filter samples by filename
            my_samples = Samples.where(filter=lambda s: s.filename.endswith('.gz'))

            # Filter classifications by completion status
            completed = Classifications.where(filter=lambda c: c.complete == True)

            # Filter projects by name
            my_projects = Projects.where(filter=lambda p: 'test' in p.name)

        Returns
        -------
        `list`
            A list of all objects matching these filters. If no filters are passed, this
            matches all objects.
        """
        # Special case a `filter` kwarg, which is a local function
        filter_func = keyword_filters.pop("filter", None)

        public = keyword_filters.pop("public", False)
        if public is True and "instances_public" not in cls._allowed_methods:
            if "public" in cls.model_fields:
                public = False
                keyword_filters["public"] = True
            else:
                raise MethodNotSupported(f"Cannot search public objects for {cls.__name__}")

        instances_route = keyword_filters.pop(
            "_instances", "instances" if not public else "instances_public"
        )
        where_schema = {}

        sort = generate_potion_sort_clause(keyword_filters.pop("sort", None), None)
        limit = keyword_filters.pop("limit", None if not public else 1000)
        where = {}

        # we're filtering by fancy objects (like SQLAlchemy's filter)
        if len(filters) > 0:
            if len(filters) == 1 and isinstance(filters[0], dict):
                where = filters[0]
            elif all(isinstance(f, str) for f in filters):
                # if it's a list of strings, treat it as multiple "GET" request
                where = {"$uri": {"$in": [cls._convert_id_to_uri(f) for f in filters]}}
            else:
                raise NotImplementedError("Advanced filtering hasn't been implemented yet")

        # we're filtering by keyword arguments (like SQLAlchemy's filter_by)
        if len(keyword_filters) > 0:
            for k, v in generate_potion_keyword_where(keyword_filters, where_schema, cls).items():
                if k in where:
                    raise AttributeError("Multiple definitions for same field {}".format(k))
                where[k] = v

        # the potion-client method returns an iterator (which lazily fetchs the records
        # using `per_page` instances per request) so for limiting we only want to fetch the first
        # n (and not instantiate all the available which is what would happen if we just sliced)
        instances = []

        page = 1
        while True:
            # TODO: This feels like a hack and should be moved up into the model/allowed_methods attribute
            instances_suffix = ""
            if instances_route == "instances_public":
                instances_suffix = "/public"
            elif instances_route == "instances_organization":
                instances_suffix = "/organization"

            resp = cls._client.get(
                f"{cls._api._base_url}{cls._resource_path}{instances_suffix}",
                params={
                    "where": json.dumps(where),
                    "sort": json.dumps(sort),
                    "page": page,
                    "per_page": DEFAULT_PAGE_SIZE,
                    "expand": "all",
                },
            )
            resp.raise_for_status()
            instances.extend(cls.model_validate(r) for r in resp.json())
            n_instances = len(instances)
            if limit is not None and n_instances >= limit:
                break
            if n_instances >= int(resp.headers.get("X-Total-Count", 0)):
                break
            page += 1

        # finally, apply local filtering function on objects before returning
        if filter_func:
            if callable(filter_func):
                instances = [obj for obj in instances if filter_func(obj) is True]
            else:
                raise OneCodexException(
                    "Expected callable for filter, got: {}".format(type(filter_func).__name__)
                )
        return instances

    @classmethod
    def create(cls, **kwargs):
        if "create" not in cls._allowed_methods:
            raise MethodNotSupported(f"Cannot create {cls.__name__} objects")
        resp = cls._client.post(f"{cls._api._base_url}{cls._resource_path}", json=kwargs)
        return cls.model_validate(resp.json())

    def update(self, **kwargs):
        if "update" not in self._allowed_methods:
            raise MethodNotSupported("Cannot update {self.__class__.__name__} objects")

        if not kwargs:
            # TODO: Consider using `pydantic_changedetect` to only send the fields that have changed
            kwargs = self.model_dump(exclude_unset=True, by_alias=True)

        patch_set = self._allowed_methods["update"].model_validate(kwargs)
        resp = self._client.patch(
            f"{self._api._base_url}{self._resource_path}/{self.id}",
            json=patch_set.model_dump(exclude_unset=True, by_alias=True),
        )
        # TODO: Nicely format Pydantic error messages here
        if resp.status_code >= 400:
            raise OneCodexException(resp.json().get("message", f"Unknown error updating {self}"))

        # Finally, we update the model in-place with the validated server-side response.
        # Note that this is lazy and does *not* trigger `getattr` and cache resolution
        # whereas `setattr(self, field, getattr(new_pydantic_obj, field))` does.
        updated_data = resp.json()
        new_pydantic_obj = self.__class__.model_validate(updated_data)
        for field, value in new_pydantic_obj:
            try:
                setattr(self, field, value)
            except MethodNotSupported:
                pass

    def delete(self) -> bool:
        if "delete" not in self._allowed_methods:
            raise MethodNotSupported(f"Cannot delete {self.__class__.__name__} objects")
        if self.id is None:
            raise MethodNotSupported("Cannot delete an unsaved object")
        resp = self._client.delete(f"{self._api._base_url}{self._resource_path}/{self.id}")
        return resp.status_code == 204

    def save(self):
        if self.id is not None:
            self.update()
        else:
            self.create()
