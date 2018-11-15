import collections
import re


class Schema(collections.Mapping):
    def __init__(self, schema):
        if isinstance(schema, Schema):
            schema = schema._schema
        self._schema = schema or {}

    @property
    def type(self):
        try:
            type = self._schema['type']

            if isinstance(type, (list, tuple)):
                return tuple(type)

            return (type,)
        except KeyError:
            return None

    @property
    def readonly_properties(self):
        if self.type != 'object':
            return ()

        properties = []
        for name, schema in self._schema.get('properties', {}).items():
            if schema.get('readOnly', False):
                properties.append(name)

        return tuple(properties)

    @property
    def required_properties(self):
        if 'object' not in self.type:
            return ()

        return tuple(self._schema.get('required', []))

    def can_include_property(self, name):
        # empty schema "{}" allows all properties
        if not self._schema:
            return True

        # only objects can have properties
        if 'object' not in self.type:
            return False

        try:
            property_schema = self._schema['properties'][name]

            if property_schema.get('readOnly', False):
                return False

            return True
        except KeyError:
            if self._schema.get('additionalProperties', True):
                return True

            for pattern in self._schema.get('patternProperties', {}).values():
                if re.match(pattern, name):
                    return True

            return False

    def __contains__(self, item):
        return item in self._schema

    def __getitem__(self, item):
        return self._schema[item]

    def __iter__(self):
        return iter(self._schema)

    def __len__(self):
        return len(self._schema)

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr(self._schema))
