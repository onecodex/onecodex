# flake8: noqa
from collections.abc import MutableMapping
from functools import partial
from operator import getitem, delitem, setitem
from six.moves.urllib.parse import urlparse, urljoin
from weakref import WeakValueDictionary
import requests

from .converter import PotionJSONDecoder, PotionJSONSchemaDecoder
from .resource import Reference, Resource, uri_for
from .links import Link
from .utils import upper_camel_case, snake_case


class Client(object):
    # TODO optional HTTP/2 support: this makes multiple queries simultaneously.

    def __init__(self, api_root_url, schema_path='/schema', fetch_schema=True, **session_kwargs):
        self._instances = WeakValueDictionary()
        self._resources = {}

        self.session = session = requests.Session()
        for key, value in session_kwargs.items():
            setattr(session, key, value)

        parse_result = urlparse(api_root_url)
        self._root_url = '{}://{}'.format(parse_result.scheme, parse_result.netloc)
        self._api_root_url = api_root_url  # '{}://{}'.format(parse_result.scheme, parse_result.netloc)
        self._root_path = parse_result.path
        self._schema_url = api_root_url + schema_path

        if fetch_schema:
            self._fetch_schema()

    def _fetch_schema(self):
        schema = self.session \
            .get(self._schema_url) \
            .json(cls=PotionJSONSchemaDecoder,
                  referrer=self._schema_url,
                  client=self)

        # NOTE these should perhaps be definitions in Flask-Potion
        for name, resource_schema in schema['properties'].items():
            resource = self.resource_factory(name, resource_schema)
            setattr(self, upper_camel_case(name), resource)

    def instance(self, uri, cls=None, default=None, **kwargs):
        instance = self._instances.get(uri, None)

        if instance is None:
            if cls is None:
                try:
                    cls = self._resources[uri[:uri.rfind('/')]]
                except KeyError:
                    cls = Reference

            if isinstance(default, Resource) and default._uri is None:
                default._status = 200
                default._uri = uri
                instance = default
            else:
                instance = cls(uri=uri, **kwargs)
            self._instances[uri] = instance
        return instance

    def fetch(self, uri, cls=PotionJSONDecoder, **kwargs):
        # TODO handle URL fragments (#properties/id etc.)
        response = self.session \
            .get(urljoin(self._root_url, uri, True))

        response.raise_for_status()

        return response.json(cls=cls,
                             client=self,
                             referrer=uri,
                             **kwargs)

    def resource_factory(self, name, schema, resource_cls=None):
        """
        Registers a new resource with a given schema. The schema must not have any unresolved references
        (such as `{"$ref": "#"}` for self-references, or otherwise). A subclass of :class:`Resource`
        may be provided to add specific functionality to the resulting :class:`Resource`.

        :param str name:
        :param dict schema:
        :param Resource resource_cls: a subclass of :class:`Resource` or None
        :return: The new :class:`Resource`.
        """
        cls = type(str(upper_camel_case(name)), (resource_cls or Resource, MutableMapping), {
            '__doc__': schema.get('description', '')
        })

        cls._schema = schema
        cls._client = self
        cls._links = links = {}

        for link_schema in schema['links']:
            link = Link(self,
                        rel=link_schema['rel'],
                        href=link_schema['href'],
                        method=link_schema['method'],
                        schema=link_schema.get('schema', None),
                        target_schema=link_schema.get('targetSchema', None))

            # Set Resource._self, etc. for the special methods as they are managed by the Resource class
            if link.rel in ('self', 'instances', 'create', 'update', 'destroy'):
                setattr(cls, '_{}'.format(link.rel), link)
            links[link.rel] = link

            if link.rel != 'update':  # 'update' is a special case because of MutableMapping.update()
                setattr(cls, snake_case(link.rel), link)

        # TODO routes (instance & non-instance)

        for property_name, property_schema in schema.get('properties', {}).items():
            # skip $uri and $id as these are already implemented in Resource and overriding them causes unnecessary
            # fetches.
            if property_name.startswith('$'):
                continue

            if property_schema.get('readOnly', False):
                # TODO better error message. Raises AttributeError("can't set attribute")
                setattr(cls,
                        property_name,
                        property(fget=partial((lambda name, obj: getitem(obj, name)), property_name),
                                 doc=property_schema.get('description', None)))
            else:
                setattr(cls,
                        property_name,
                        property(fget=partial((lambda name, obj: getitem(obj, name)), property_name),
                                 fset=partial((lambda name, obj, value: setitem(obj, name, value)), property_name),
                                 fdel=partial((lambda name, obj: delitem(obj, name)), property_name),
                                 doc=property_schema.get('description', None)))

        root = None
        if 'instances' in links:
            root = cls._instances.href
        elif 'self' in links:
            root = cls._self.href[:cls._self.href.rfind('/')]
        else:
            root = self._root_path + '/' + name.replace('_', '-')

        self._resources[root] = cls
        return cls


ASC = ASCENDING = False
DESC = DESCENDING = True
