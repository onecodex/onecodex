import inspect
import json


from typing import ClassVar, Optional, List, Type, TypedDict
from pprint import pformat
from html import escape

from pydantic import BaseModel as PydanticBaseModel
from pydantic import ConfigDict, Field

from onecodex.exceptions import MethodNotSupported, OneCodexException
from onecodex.models.helpers import (
    generate_potion_sort_clause,
    generate_potion_keyword_where,
)


# TODO: Unclear why this is necessary, but it is. See if we can remove.
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
        return "unknown"

    def __repr__(self):
        return f"<ApiRef {self.resource_type}:{self.id}>"

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

        return get_model_class(self.ref)


class OneCodexBase(PydanticBaseModel):
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

    field_uri: str = Field(..., alias="$uri")

    def __init__(self, **data):
        super().__init__(**data)
        self._resolved_cache = {}

    @property
    def id(self) -> Optional[str]:
        if self.field_uri is None:
            return None
        return self.field_uri.split("/")[-1]

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.id}>"

    def _repr_html_(self):
        return """<table>
        <thead>
            <tr>
                <th colspan="2"><code>{cls}({id})</code></th>
            </tr>
        </thead>
        <tbody>{properties}</tbody>
        </table>""".format(
            cls=self.__class__.__name__,
            id=escape(repr(self.id)),
            properties="\n".join(
                "<tr><td>{}</td><td><code>{}</code></td>".format(escape(k), escape(pformat(v)))
                for k, v in self._properties.items()
                if not k.startswith("$")
            ),
        )

    def __dir__(self):
        # We expose all the defined model_fields, any computed model fields, plus
        # the special `id` property for backwards compatibility (but this is not a computed field
        # as then it would be included in serialized JSON output)
        fields = (
            list(self.__class__.model_fields.keys())
            + list(self.__class__.model_computed_fields.keys())
            + ["id"]
        )

        # Get the MRO and find where to stop
        mro = inspect.getmro(self.__class__)

        # Find OneCodexBase in the MRO
        try:
            onecodex_base_index = mro.index(OneCodexBase)
            relevant_classes = mro[: onecodex_base_index + 1]
        except ValueError:
            relevant_classes = [self.__class__]

        methods = set()
        is_instance = self != self.__class__

        for cls in relevant_classes:
            for name in dir(cls):
                if not name.startswith("_") and name not in _EXCLUDE_PYDANTIC_CALLABLES:
                    attr = getattr(cls, name, None)
                    if callable(attr):
                        # Check if it's a classmethod by looking at __self__
                        is_classmethod = hasattr(attr, "__self__") and attr.__self__ is cls

                        # Skip classmethods if we're on an instance
                        if is_instance and is_classmethod:
                            continue

                        # Only include if method is defined in this class
                        if name in cls.__dict__:
                            methods.add(name)

        return fields + sorted(methods)

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
    def get(cls, id: str) -> "OneCodexBase":
        resp = cls._client.get(f"{cls._api._base_url}{cls._resource_path}/{id}?expand=all")
        return cls.model_validate(resp.json())

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
                    "per_page": DEFAULT_PAGE_SIZE,
                    "expand": "all",
                },
            )
            instances.extend([cls.model_validate(r) for r in resp.json()])
            if limit is not None and len(instances) >= limit:
                break
            if resp.links.get("next") is None:
                break

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
            kwargs = self.model_dump(exclude_unset=True, by_alias=True)

        # if not kwargs:
        #     # noop
        #     return self

        patch_set = self._allowed_methods["update"].model_validate(kwargs)
        resp = self._client.patch(
            f"{self._api._base_url}{self._resource_path}/{self.id}",
            data=patch_set.model_dump_json(),
        )
        updated_data = resp.json()
        copied_cache = self._resolved_cache.copy()
        self = self.__class__.model_validate(updated_data)
        self._resolved_cache = copied_cache

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
