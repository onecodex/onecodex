import calendar
from functools import partial
try:
    from simplejson import JSONEncoder, JSONDecoder
except ImportError:
    from json import JSONEncoder, JSONDecoder
from datetime import date, datetime
from six.moves.urllib.parse import urljoin
import six

from .resource import Reference

try:
    from datetime import timezone
except ImportError:
    from datetime import tzinfo, timedelta

    # pragma: no cover
    class timezone(tzinfo):
        def __init__(self, utcoffset, name=None):
            self._utcoffset = utcoffset
            self._name = name

        def utcoffset(self, dt):
            return self._utcoffset

        def tzname(self, dt):
            return self._name

        def dst(self, dt):
            return timedelta(0)

    timezone.utc = timezone(timedelta(0), 'UTC')


try:
    import pandas as pd
except ImportError:
    pd = None


class PotionJSONEncoder(JSONEncoder):
    def encode(self, o):
        root_id = id(o)
        if self.check_circular:
            markers = {}
        else:
            markers = None

        def _encode(o):
            if isinstance(o, (list, tuple, dict)):
                if markers is not None:
                    marker_id = id(o)
                    if marker_id in markers:
                        if marker_id == root_id:
                            return {"$ref": "#"}
                        raise ValueError("Circular reference detected")
                    markers[marker_id] = o
                try:
                    if isinstance(o, dict):
                        return {k: _encode(v) for k, v in o.items()}
                    else:
                        return [_encode(v) for v in o]
                finally:
                    if markers is not None:
                        del markers[marker_id]

            if isinstance(o, date):
                return {"$date": int(calendar.timegm(o.timetuple()) * 1000)}
            if isinstance(o, datetime):
                return {"$date": int(calendar.timegm(o.utctimetuple()) * 1000)}

            if isinstance(o, Reference):
                # FIXME if reference is not saved, save it first here
                return {"$ref": o._uri}

            try:
                if pd and pd.api.types.is_integer(o):
                    return int(o)
            except AttributeError:
                # some versions of pandas may not have this API
                pass

            try:
                if pd and pd.api.types.is_float(o):
                    return float(o)
            except AttributeError:
                pass

            return o

        return JSONEncoder.encode(self, _encode(o))


class PotionJSONDecoder(JSONDecoder):
    def __init__(self, client, referrer=None, uri_to_instance=True, default_instance=None, *args, **kwargs):
        self.client = client
        self.referrer = referrer
        self.uri_to_instance = uri_to_instance
        self.default_instance = default_instance
        JSONDecoder.__init__(self, *args, **kwargs)

    def _decode(self, o, depth=0):
        if isinstance(o, dict):
            if len(o) == 1:
                if "$date" in o:
                    return datetime.fromtimestamp(o["$date"] / 1000.0, timezone.utc)
                if "$ref" in o and isinstance(o["$ref"], six.string_types):
                    reference = o["$ref"]
                    if reference.startswith("#"):
                        reference = urljoin(self.referrer, reference, True)
                    return self.client.instance(reference)
            elif self.uri_to_instance and "$uri" in o and isinstance(o["$uri"], six.string_types):
                # TODO handle or ("$id" in o and "$type" in o)
                if depth == 0:
                    instance = self.client.instance(o['$uri'], default=self.default_instance)
                else:
                    instance = self.client.instance(o['$uri'])

                instance._status = 200
                instance._properties.update({k: self._decode(v, depth + 1) for k, v in o.items()})
                return instance

            return {k: self._decode(v, depth + 1) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [self._decode(v, depth + 1) for v in o]
        return o

    def decode(self, s, *args, **kwargs):
        o = JSONDecoder.decode(self, s, *args, **kwargs)
        return self._decode(o)


class JSONSchemaReference(Reference):
    @classmethod
    def _resolve(self, client, uri):
        return client.fetch(uri, cls=PotionJSONSchemaDecoder)


class PotionJSONSchemaDecoder(JSONDecoder):
    def __init__(self, client, referrer=None, *args, **kwargs):
        self.client = client
        self.referrer = referrer
        JSONDecoder.__init__(self, *args, **kwargs)

    def decode(self, s, *args, **kwargs):
        o = JSONDecoder.decode(self, s, *args, **kwargs)
        return schema_resolve_refs(o, partial(self.client.instance,
                                              cls=JSONSchemaReference,
                                              client=self.client))


def schema_resolve_refs(schema, ref_resolver=None, root=None):
    """
    Helper method for decoding references. Self-references are resolved automatically; other references are
    resolved using a callback function.

    :param object schema:
    :param callable ref_resolver:
    :param None root:
    :return:
    """
    # FIXME more stable implementation that only attempts to resolve {"$ref"} objects where they are allowed.
    if isinstance(schema, dict):
        if len(schema) == 1 and "$ref" in schema and isinstance(schema["$ref"], six.string_types):
            reference = schema["$ref"]
            if reference.startswith("#"):
                # TODO should also resolve any paths within the reference, which would need to be deferred.
                return root
            return ref_resolver(reference)

        resolved = {}
        for k, v in schema.items():
            resolved[k] = schema_resolve_refs(v,
                                              ref_resolver=ref_resolver,
                                              root=root if root is not None else resolved)
        return resolved
    if isinstance(schema, (list, tuple)):
        return [schema_resolve_refs(v, ref_resolver=ref_resolver, root=root) for v in schema]
    return schema
