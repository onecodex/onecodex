import collections
from pprint import pformat

import six

from .exceptions import ItemNotFound, MultipleItemsFound
from .utils import escape


def uri_for(reference):
    """
    Utility function for accessing the URI or a reference.

    :param reference:
    :return:
    """
    return reference._uri


class Reference(collections.Mapping):
    """

    This implementation makes the assumption that a {$ref} object always points to an object, never an array or
    any of the other types.
    """
    _client = None
    __properties = None

    def __init__(self, uri, client=None):
        self._status = None
        self._uri = uri
        self.__properties = {'$uri': uri}
        if client is not None:
            self._client = client

    @classmethod
    def _resolve(self, client, uri):
        return client.fetch(uri, uri_to_instance=False)

    @property
    def _properties(self):
        if self._uri and self._status is None:
            self.__properties = self._resolve(self._client, self._uri)
            self._status = 200
        return self.__properties

    @_properties.setter
    def _properties(self, value):
        self.__properties = value
        self._status = 200

    def __contains__(self, item):
        return item in self._properties

    def __getitem__(self, item):
        return self._properties[item]

    def __iter__(self):
        return iter(self._properties)

    def __len__(self):
        return len(self._properties)

    def __repr__(self):
        return '{cls}({uri})'.format(cls=self.__class__.__name__,
                                     uri=repr(self._uri))


class Resource(Reference):
    _client = None
    _links = None
    _self = None
    _instances = None
    _create = None
    _destroy = None
    _update = None

    def __new__(cls, uri=None, **kwargs):
        instance = None
        if uri is not None:
            if not (isinstance(uri, six.string_types) and uri.startswith('/')) and cls._self is not None:
                uri = cls._self.href.format(id=uri)

            instances = cls._client._instances
            instance = instances.get(uri, None)

        if instance is None:
            instance = super(Resource, cls).__new__(cls)
            super(Resource, instance).__init__(uri)
            instance._properties = {'$uri': uri}
            if not kwargs:
                instance._status = None
            else:
                instance._status = 200
                instance._properties.update(kwargs)

            if uri is not None:
                instances[uri] = instance

        # NOTE ensures that there is a single instance of a Resource with a given URL unless one creates an item
        # without URL and creates an item with the URL the first item is going to have, before saving the first item.
        return instance

    def __init__(self, uri=None, **kwargs):
        pass  # Must be blank. See __new__()

    @property
    def id(self):
        if self._uri is not None:
            id_ = self._uri[self._uri.rfind("/") + 1:]
            return id_
        return None

    # TODO cache this property
    @property
    def _validator(self):
        return None

    def __delitem__(self, item):
        del self._properties[item]

    def __setitem__(self, item, value):
        self._properties[item] = value

    def update(self, *args, **kwargs):
        self._properties.update(*args, **kwargs)
        self.save()

    @classmethod
    def first(cls, **params):
        try:
            return cls._instances(per_page=1, **params)[0]
        except IndexError:
            raise ItemNotFound("No '{}' item found matching: {}".format(cls.__name__, repr(params)))

    @classmethod
    def one(cls, **params):
        matching = cls._instances(per_page=1, **params)
        if len(matching) > 1:
            raise MultipleItemsFound("Multiple items found matching: {}".format(repr(params)))
        try:
            return matching[0]
        except IndexError:
            raise ItemNotFound("No '{}' item found matching: {}".format(cls.__name__, repr(params)))

    @classmethod
    def fetch(cls, id):
        return cls._self(id=id)

    def check(self):
        pass

    def save(self):
        # TODO only save the appropriate properties defined in the schemas
        # That logic should live in the links themselves.

        if self._uri is None:
            return self._create(**self)
        else:
            return self._update(**self)

    def delete(self):
        return self._destroy(id=self.id)

    def _repr_html_(self):
        return '''<table>
        <thead>
            <tr>
                <th colspan="2"><code>{cls}({id})</code></th>
            </tr>
        </thead>
        <tbody>{properties}</tbody>
        </table>'''.format(
            cls=self.__class__.__name__,
            id=escape(repr(self.id)),
            properties='\n'.join('<tr><td>{}</td><td><code>{}</code></td>'.format(escape(k),
                                                                                  escape(pformat(v)))
                                 for k, v in self._properties.items() if not k.startswith('$')))

    def __repr__(self):
        # if self._status == 200:
        #     parts = ['{}={}'.format(k, repr(v)) for k, v in self._properties.items() if not k.startswith('$')]
        #     if parts:
        #         return '{cls}({id}, {properties})'.format(cls=self.__class__.__name__,
        #                                                   id=repr(self.id),
        #                                                   properties=', '.join(parts))
        # TODO some way to define a good key for display
        return '{cls}({id})'.format(cls=self.__class__.__name__, id=repr(self.id))
