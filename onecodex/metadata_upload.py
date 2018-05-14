import datetime


def validate_appendables(appendables, api):
    appendables['valid_tags'] = []
    appendables['custom_metadata'] = {}
    appendables['valid_metadata'] = {}
    validate_tags(appendables)
    validate_metadata(appendables, api)
    return appendables


def validate_tags(appendables):
    if 'tags' not in appendables:
        return
    tag_array = appendables['tags']
    for tag in tag_array:
        # FIXME - change this number and specify the exception
        if len(tag) > 500:
            raise Exception

        appendables['valid_tags'].append(tag)


def validate_metadata(appendables, api):
    if 'metadata' not in appendables:
        return
    schema_props = metadata_properties(api)
    for key, value in appendables['metadata'].items():
        if is_blacklisted(key):
            # FIXME - Specify this Exception
            raise Exception

        if key in schema_props.keys():
            settable_value = validate_metadata_against_schema(schema_props, key, value)
            appendables['valid_metadata'][key] = settable_value
        else:
            appendables['custom_metadata'][key] = value


def validate_metadata_against_schema(schema_props, key, value):
    schema_rules = schema_props[key]
    if 'enum' in schema_rules:
        return validate_enum(value, schema_rules)
    elif 'number' in schema_rules['type']:
        return validate_number(value, schema_rules)
    elif 'boolean' in schema_rules['type']:
        return validate_boolean(value)
    elif 'format' in schema_rules and 'date-time' in schema_rules['format']:
        return validate_datetime(value)
    else:
        return value


def validate_enum(value, schema_rules):
    if value not in schema_rules['enum']:
        # FIXME - Make this more explicit
        raise Exception
    return value


def validate_number(value, schema_rules):
    num_value = float(value)
    if 'minimum' in schema_rules and schema_rules['minimum'] > num_value:
        # FIXME - Make this more explicit
        raise Exception
    if 'maximum' in value and schema_rules['maximum'] < num_value:
        # FIXME - Make this more explicit
        raise Exception
    return num_value


def validate_boolean(value):
    if value.lower() in truthy_values():
        return True
    elif value.lower() in falsy_values():
        return False
    else:
        # FIXME - Make this more explicit
        raise Exception


def validate_datetime(value):
    datetime_vals = value.split(',')
    datetime_ints = list(map(int, datetime_vals))
    try:
        return datetime.datetime(*datetime_ints)
    except TypeError:
        # FIXME - Make this more explicit
        raise Exception


def is_blacklisted(key):
    key in ['$uri', 'custom']


def truthy_values():
    return ['true', '1', 't', 'y', 'yes']


def falsy_values():
    return ['false', 0, 'f', 'n', 'no']


def set_valid_appendables(api, sample_uuids, appendables):
    # TODO - Change this to something like
    # `samples = Sample.where('sample_uuid in ?', sample_uuids)
    # for sample in samples:

    for sample_uuid in sample_uuids:
        sample = api.Samples.get(sample_uuid)
        new_tag_array = []
        for tag in appendables['valid_tags']:
            existing_tags = api.Tags.where(name=tag)
            if existing_tags:
                unsaved_tag = existing_tags[0]
                new_tag_array.append(unsaved_tag)
            else:
                unsaved_tag = api.Tags(name=tag, sample=sample)
                unsaved_tag.save()
                # TODO - this following line should not be necessary. Consider eventually update the Potion client to not require this.
                new_tag_array.append(unsaved_tag)

        if new_tag_array:
            potential_tags = sample.tags
            potential_tags.extend(new_tag_array)
            sample.tags = potential_tags

        for k, v in appendables['valid_metadata'].items():
            setattr(sample.metadata, k, v)

        for k, v in appendables['custom_metadata'].items():
            sample.metadata.custom[k] = v
        sample.metadata.save()
        sample.save()


def metadata_properties(api):
    return api.Metadata._resource._schema['properties']
