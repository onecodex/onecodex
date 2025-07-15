import re

import pydantic

from onecodex.exceptions import ValidationError
from onecodex.models.sample import _MetadataPatchSchema
from onecodex.models.misc import Tags


def validate_appendables(appendables, api):
    appendables["valid_tags"] = []
    appendables["valid_metadata"] = {"custom": {}}
    validate_tags(appendables, api)
    validate_metadata(appendables, api)
    return appendables


def validate_tags(appendables, api):
    if "tags" not in appendables:
        return

    tag_array = appendables["tags"]
    try:
        validated_tags = [Tags.model_validate({"name": t}) for t in tag_array]
        appendables["valid_tags"] = [{"name": t.name} for t in validated_tags]
    except pydantic.ValidationError as e:
        raise ValidationError(f"Invalid tag: {e}")


def validate_metadata(appendables, api):
    if "metadata" not in appendables:
        return

    for key, value in appendables["metadata"].items():
        if is_blacklisted(key):
            raise ValidationError("{} cannot be manually updated".format(key))

        if key in _MetadataPatchSchema.model_fields:
            settable_value = validate_metadata_against_schema(key, value)
            appendables["valid_metadata"][key] = settable_value
        else:
            coerced_value = coerce_custom_value(value)
            appendables["valid_metadata"]["custom"][key] = coerced_value


def validate_metadata_against_schema(key, value):
    try:
        value = getattr(_MetadataPatchSchema.model_validate({key: value}), key)
    except pydantic.ValidationError as exc_info:
        # Get all Pydantic validation errors for the given key
        msg = []
        for error in exc_info.errors():
            if error["loc"][0] == key:
                msg.append(error["msg"])

        msg = ", ".join(msg) if len(msg) > 0 else str(exc_info)
        raise ValidationError(f"Invalid metadata: {msg}")
    return value


def validate_enum(value, schema_rules):
    if value not in schema_rules["enum"]:
        # This is gross, but is necessary for string comparison between Python2 ([u'Illumina HiSeq']) and Python3 (['Illumina Hiseq']). On the plus side, it makes the error message more Human Readable.
        error_array = []
        for rule in schema_rules["enum"]:
            if rule:
                error_array.append(str(rule))
            else:
                error_array.append(rule)
        raise ValidationError(
            "{} is not a valid value for this key. Value must be one of the following options: {}".format(
                value, error_array
            )
        )
    return value


def validate_number(value, schema_rules):
    num_value = value
    try:
        num_value = float(value)
    except ValueError:
        raise ValidationError("{} must be a number".format(value))
    if "minimum" in schema_rules and num_value <= schema_rules["minimum"]:
        raise ValidationError(
            "{} must be larger than the minimum value: {}".format(value, schema_rules["minimum"])
        )
    if "maximum" in schema_rules and num_value >= schema_rules["maximum"]:
        raise ValidationError(
            "{} must be smaller than the maximum value: {}".format(value, schema_rules["maximum"])
        )
    return num_value


def validate_boolean(value):
    if value.lower() in truthy_values():
        return True
    elif value.lower() in falsy_values():
        return False
    else:
        raise ValidationError('{} must be either "true" or "false"'.format(value))


def validate_datetime(value):
    if not is_iso_8601_compliant(value):
        raise ValidationError(
            '"{}" must be formatted in iso8601 compliant date format. Example: "2018-05-15T16:21:36+00:00"'.format(
                value
            )
        )

    return value


def is_blacklisted(key):
    return key in ["$uri", "custom"]


def truthy_values():
    return ["true", "1", "t", "y", "yes"]


def falsy_values():
    return ["false", "0", "f", "n", "no"]


def coerce_custom_value(value):
    try:
        return float(value)
    except ValueError:
        pass

    if value.lower() in truthy_values():
        return True

    if value.lower() in falsy_values():
        return False

    # String case, includes dates in JSON
    return value


def is_iso_8601_compliant(value):
    iso8601 = re.compile(
        r"^(?P<full>((?P<year>\d{4})([/-]?(?P<mon>(0[1-9])|(1[012]))([/-]?(?P<mday>(0[1-9])|([12]\d)|(3[01])))?)?(?:T(?P<hour>([01][0-9])|(?:2[0123]))(\:?(?P<min>[0-5][0-9])(\:?(?P<sec>[0-5][0-9]([\,\.]\d{1,10})?))?)?(?:Z|([\-+](?:([01][0-9])|(?:2[0123]))(\:?(?:[0-5][0-9]))?))?)?))$"
    )
    return iso8601.match(value)
