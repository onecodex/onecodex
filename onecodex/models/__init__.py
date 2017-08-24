from datetime import datetime
import inspect
import itertools
import sys

from dateutil.parser import parse
from requests.exceptions import HTTPError
from potion_client.resource import Resource
import six

from onecodex.exceptions import MethodNotSupported, PermissionDenied, ServerError
from onecodex.models.helpers import (check_bind, generate_potion_sort_clause,
                                     generate_potion_keyword_where)


DEFAULT_PAGE_SIZE = 200


class OneCodexBase(object):
    """
    A parent object for all the One Codex objects that wraps the Potion-Client API and makes
    access and usage easier.
    """
    def __init__(self, _resource=None, **kwargs):
        # FIXME: allow setting properties via kwargs?
        # FIXME: get a resource from somewhere instead of setting to None (lots of stuff assumes
        # non-None) if we have a class.resource?
        if _resource is not None:
            if not isinstance(_resource, Resource):
                raise TypeError('Use the .get() method to fetch an individual resource.')
            self._resource = _resource
        elif hasattr(self.__class__, '_resource'):
            for key, val in kwargs.items():
                # This modifies kwargs in place to be the underlying
                # Potion-Client resource
                if isinstance(val, OneCodexBase):
                    kwargs[key] = val._resource
            self._resource = self.__class__._resource(**kwargs)

    def __repr__(self):
        return '<{} {}>'.format(self.__class__.__name__, self.id)

    def _repr_html_(self):
        return self._resource._repr_html_()

    def __dir__(self):
        # this only gets called on instances, so we're okay to add all the properties because
        # this won't appear when you call, e.g. dir(ocx.Samples)

        fields = [str(f) if f != '$uri' else 'id' for f in
                  self.__class__._resource._schema['properties']]

        # this might be a little too clever, but we mask out class methods/fxns from the instances
        base_object_names = []
        for name, obj in inspect.getmembers(self.__class__):
            if inspect.isfunction(obj):  # .save() and .delete() are functions in Py3
                base_object_names.append(name)
            if inspect.ismethod(obj) and obj.__self__ is not self.__class__:
                base_object_names.append(name)

        return fields + base_object_names

    def __getattr__(self, key):
        if hasattr(self, '_resource') and hasattr(self.__class__, '_resource'):
            schema_key = key if key != 'id' else '$uri'
            schema = self.__class__._resource._schema['properties'].get(schema_key)
            if schema is not None:
                value = getattr(self._resource, key)
                if isinstance(value, Resource):
                    # convert potion resources into wrapped ones
                    resource_path = value._uri.rsplit('/', 1)[0]
                    return _model_lookup[resource_path](_resource=value)
                elif isinstance(value, list) and len(value) > 0 and isinstance(value[0], Resource):
                    # convert lists of potion resources into wrapped ones
                    resource_path = value[0]._uri.rsplit('/', 1)[0]
                    return [_model_lookup[resource_path](_resource=o) for o in value]
                else:
                    if key == 'id':
                        # undo the bad coercion from potion_client/resource.py#L111
                        return str(value)
                    if schema.get('format') == 'date-time' and value is not None:
                        return parse(value)
                    return value
        elif key == 'id' or key in self.__class__._resource._schema['properties']:
            # make fields appear blank if there's no _resource bound to me
            return None

        raise AttributeError('\'{}\' object has no attribute \'{}\''.format(
            self.__class__.__name__, key
        ))

    def __setattr__(self, key, value):
        if key.startswith("_"):  # Allow directly setting _attributes, incl. _resource
            # these are any fields that have to be settable normally
            super(OneCodexBase, self).__setattr__(key, value)
            return
        elif key == 'id':
            raise AttributeError('can\'t set attribute')
        elif isinstance(value, OneCodexBase):
            self._resource[key] = value._resource
            return
        elif isinstance(value, (list, tuple)):
            # convert any fancy items into their underlying resources
            new_value = []
            for v in value:
                new_value.append(v._resource if isinstance(v, OneCodexBase) else v)
            # coerce back to the value passed in
            self._resource[key] = type(value)(new_value)
            return
        elif hasattr(self, '_resource') and hasattr(self.__class__, '_resource'):
            schema = self.__class__._resource._schema['properties'].get(key)

            if schema is not None:
                # do some type checking against the schema
                if not self.__class__._has_schema_method('update'):
                    raise MethodNotSupported('{} do not support editing.'.format(
                        self.__class__.__name__
                    ))
                if schema.get('readOnly', False):
                    raise MethodNotSupported('{} is a read-only field'.format(key))

                if schema.get('format') == 'date-time':
                    if isinstance(value, datetime):
                        if value.tzinfo is None:
                            value = value.isoformat() + 'Z'
                        else:
                            value = value.isoformat()

                # changes on this model also change the potion resource
                self._resource[key] = value
                return

        raise AttributeError('\'{}\' object has no attribute \'{}\''.format(
            self.__class__.__name__, key
        ))

    def __delattr__(self, key):
        if not self.__class__._has_schema_method('update'):
            raise MethodNotSupported('{} do not support editing.'.format(self.__class__.__name__))

        if hasattr(self, '_resource') and key in self._resource.keys():
            # changes on this model also change the potion resource
            del self._resource[key]

    def __eq__(self, other):
        # TODO: We should potentially check that both resources are up-to-date
        return self._resource._uri == other._resource._uri

    @classmethod
    def _convert_id_to_uri(cls, uuid):
        base_uri = cls._resource._schema._uri.replace('/schema#', '')
        if not uuid.startswith(base_uri):
            uuid = '{}/{}'.format(base_uri, uuid)
        return uuid

    @classmethod
    def _has_schema_method(cls, method_name):
        # potion-client is too stupid to check the schema before allowing certain operations
        # so we manually check it before allowing some instance methods

        # FIXME: this doesn't actually work though, because potion creates these routes for all
        # items :/
        method_links = cls._resource._schema['links']
        return any(True for l in method_links if l['rel'] == method_name)

    @classmethod
    def all(cls, sort=None, limit=None):
        """
        Returns all of the {classname}. Alias for {classname}.where() (without filter arguments).

        See `{classname}.where` for documentation on the `sort` and `limit` parameters.
        """.format(classname=cls.__name__)
        return cls.where(sort=sort, limit=limit)

    @classmethod
    def where(cls, *filters, **keyword_filters):
        """
        Retrieves {classname} from the One Codex server.

        Parameters
        ----------
        filters : objects
            Advanced filters to use (not implemented)
        sort : string | list, optional
            Sort the results by this field (or list of fields). By default in descending order,
            but if any of the fields start with the special character ^, sort in ascending order.
            For example, sort=['size', '^filename'] will sort by size from largest to smallest and
            filename from A-Z for items with the same size.
        limit : integer, optional
            Number of records to return. For smaller searches, this can reduce the number of
            network requests made.
        keyword_filters : strings | objects
            Filter the results by specific keywords (or filter objects, in advanced usage)

        Returns
        -------
        list
            A list of all {classname} matching these filters. If no filters are passed, this
            matches all {classname}.
        """.format(classname=cls.__name__)
        check_bind(cls)

        public = False
        if any(x['rel'] == 'instances_public' for x in cls._resource._schema['links']):
            public = keyword_filters.pop('public', False)

        instances_route = keyword_filters.pop('_instances',
                                              'instances' if not public else 'instances_public')

        schema = next(l for l in cls._resource._schema['links'] if l['rel'] == instances_route)
        sort_schema = schema['schema']['properties']['sort']['properties']
        where_schema = schema['schema']['properties']['where']['properties']

        sort = generate_potion_sort_clause(keyword_filters.pop('sort', None), sort_schema)
        limit = keyword_filters.pop('limit', None if not public else 1000)
        where = {}

        # we're filtering by fancy objects (like SQLAlchemy's filter)
        if len(filters) > 0:
            if len(filters) == 1 and isinstance(filters[0], dict):
                where = filters[0]
            elif all(isinstance(f, six.string_types) for f in filters):
                # if it's a list of strings, treat it as an multiple "get" request
                where = {'$uri': {'$in': [cls._convert_id_to_uri(f) for f in filters]}}
            else:
                # we're doing some more advanced filtering
                raise NotImplementedError('Advanced filtering hasn\'t been implemented yet')

        # we're filtering by keyword arguments (like SQLAlchemy's filter_by)
        if len(keyword_filters) > 0:
            for k, v in generate_potion_keyword_where(keyword_filters, where_schema, cls).items():
                if k in where:
                    raise AttributeError('Multiple definitions for same field {}'.format(k))
                where[k] = v

        # the potion-client method returns an iterator (which lazily fetchs the records
        # using `per_page` instances per request) so for limiting we only want to fetch the first
        # n (and not instantiate all the available which is what would happen if we just sliced)
        cursor = getattr(cls._resource, instances_route)(where=where, sort=sort, per_page=DEFAULT_PAGE_SIZE)
        if limit is not None:
            cursor = itertools.islice(cursor, limit)
        return [cls(_resource=r) for r in cursor]

    @classmethod
    def get(cls, uuid):
        """
        Retrieve one specific {classname} object from the server by its UUID
        (unique 16-character id). UUIDs can be found in the web browser's address bar while
        viewing analyses and other objects.

        Parameters
        ----------
        uuid : string
            UUID of the {classname} object to retrieve.

        Returns
        -------
        OneCodexBase | None
            The {classname} object with that UUID or None if no {classname} object could be found.

        Examples
        --------
        >>> api.Samples.get('xxxxxxxxxxxxxxxx')
        <Sample xxxxxxxxxxxxxxxx>
        """.format(classname=cls.__name__)
        check_bind(cls)

        # we're just retrieving one object from its uuid
        try:
            resource = cls._resource.fetch(uuid)
            if isinstance(resource, list):
                # TODO: Investigate why potion .fetch()
                #       method is occassionally returning a list here...
                if len(resource) == 1:
                    resource = resource[0]
                else:
                    raise TypeError("Potion-Client error in fetching resource")
        except HTTPError as e:
            # 404 error means this doesn't exist
            if e.response.status_code == 404:
                return None
            else:
                raise e
        return cls(_resource=resource)

    def delete(self):
        """
        Delete this {classname} object off the One Codex server.
        """.format(classname=self.__class__.__name__)
        check_bind(self)
        if self.id is None:
            raise ServerError('{} object does not exist yet'.format(self.__class__.name))
        elif not self.__class__._has_schema_method('destroy'):
            raise MethodNotSupported('{} do not support deletion.'.format(self.__class__.__name__))

        try:
            self._resource.delete()
        except HTTPError as e:
            if e.response.status_code == 403:
                raise PermissionDenied('')  # FIXME: is this right?
            else:
                raise e

    def save(self):
        """
        Either create or persist changes on this {classname} object back to the One Codex server.
        """.format(classname=self.__class__.__name__)
        check_bind(self)

        creating = self.id is None
        if creating and not self.__class__._has_schema_method('create'):
            raise MethodNotSupported('{} do not support creating.'.format(self.__class__.__name__))
        if not creating and not self.__class__._has_schema_method('update'):
            raise MethodNotSupported('{} do not support updating.'.format(self.__class__.__name__))

        try:
            self._resource.save()
        except HTTPError as e:
            if e.response.status_code == 400:
                err_json = e.response.json().get('errors', [])
                msg = pretty_print_error(err_json)
                raise ServerError(msg)
            elif e.response.status_code == 404:
                action = 'creating' if creating else 'updating'
                raise MethodNotSupported('{} do not support {}.'.format(self.__class__.__name__,
                                                                        action))
            elif e.response.status_code == 409:
                raise ServerError('This {} object already exists'.format(self.__class__.__name__))
            else:
                raise e


from onecodex.models.analysis import Analyses, Classifications, Alignments, Panels  # noqa
from onecodex.models.misc import Jobs, Projects, Tags, Users  # noqa
from onecodex.models.sample import Samples, Metadata  # noqa


__all__ = ['Samples', 'Classifications', 'Alignments', 'Panels', 'Jobs', 'Projects', 'Tags',
           'Users', 'Metadata']


def pretty_print_error(err_json):
    """Pretty print Flask-Potion error messages for the user
    """
    # Special case validation errors
    if len(err_json) == 1 and 'validationOf' in err_json[0]:
        required_fields = ', '.join(err_json[0]['validationOf']['required'])
        return 'Validation error. Requires properties: {}.'.format(required_fields)

    # General error handling
    msg = '; '.join(err.get('message', '') for err in err_json)

    # Fallback
    if not msg:
        msg = 'Bad request.'
    return msg


# go through all the models and generate a lookup table (to use in binding in the API and elsewhere)
def is_oc_class(cls):
    return inspect.isclass(cls) and issubclass(cls, OneCodexBase)


_model_lookup = {}
for name, obj in inspect.getmembers(sys.modules[__name__], is_oc_class):
    if hasattr(obj, '_resource_path'):
        _model_lookup[obj._resource_path] = obj
