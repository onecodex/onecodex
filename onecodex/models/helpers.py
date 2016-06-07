import inspect

from onecodex.exceptions import UnboundObject


def as_uri(uuid, base_class):
    return base_class._resource._schema._uri.rstrip('#') + '/' + uuid


def coerce_search_value(search_value, field_name, base_class):
    from onecodex.models import OneCodexBase  # in here to prevent circular import
    if field_name == '$uri':
        return as_uri(field_name, base_class)
    elif isinstance(search_value, OneCodexBase):
        return {'$ref': search_value._resource._uri}
    return search_value


def check_bind(self_or_cls):
    if not hasattr(self_or_cls, '_resource'):
        name = 'class' if inspect.isclass(self_or_cls) else 'instance'
        raise UnboundObject('This {} is not associated with an API binding.'.format(name))


def generate_potion_sort_clause(sort_items, sort_schema):
    if sort_items is None:
        return {}

    if not isinstance(sort_items, list):
        sort_items = [sort_items]

    sort_clause = {}
    for item in sort_items:
        if item.lstrip('^') not in sort_schema:
            raise AttributeError('Attribute {} can not be sorted on'.format(item.lstrip('^')))
        if item.startswith('^'):
            sort_clause[item[1:]] = False
        else:
            sort_clause[item] = True
    return sort_clause


def generate_potion_keyword_where(keyword_filters, where_schema, base_class):
    where = {}
    for keyword in keyword_filters:
        search_value = keyword_filters[keyword]

        if keyword == 'id':
            keyword = '$uri'

        if keyword not in where_schema:
            raise AttributeError('{} can not be searched on {}'.format(
                base_class.__name__, keyword
            ))

        avail_searches = [v['required'] for v in where_schema[keyword]['anyOf']
                          if 'required' in v]
        # flatten the list
        avail_searches = [item for sublist in avail_searches for item in sublist]

        # TODO: do schema type checking here too?
        if '$eq' not in avail_searches and '$containsall' in avail_searches:
            if not isinstance(search_value, list):
                search_value = [search_value]
            where[keyword] = {'$containsall': [coerce_search_value(v, keyword, base_class)
                                               for v in search_value]}
        elif isinstance(search_value, list):
            where[keyword] = {'$in': [coerce_search_value(v, keyword, base_class)
                                      for v in search_value]}
        else:
            where[keyword] = coerce_search_value(search_value, keyword, base_class)
    return where


def truncate_string(s, length=24):
    if len(s) < length - 3:
        return s
    else:
        s = s[0:(length - 3)]
        if s[-1] == '.':
            s = s + '..'
        else:
            s = s + '...'
        return s
