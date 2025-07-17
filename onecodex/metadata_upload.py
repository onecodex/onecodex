import pydantic

from onecodex.exceptions import ValidationError
from onecodex.models.schemas.sample import MetadataPatchSchema
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

        if key in MetadataPatchSchema.model_fields:
            settable_value = validate_metadata_against_schema(key, value)
            appendables["valid_metadata"][key] = settable_value
        else:
            coerced_value = coerce_custom_value(value)
            appendables["valid_metadata"]["custom"][key] = coerced_value


def validate_metadata_against_schema(key, value):
    try:
        value = getattr(MetadataPatchSchema.model_validate({key: value}), key)
    except pydantic.ValidationError as exc_info:
        # Get all Pydantic validation errors for the given key
        msg = []
        for error in exc_info.errors():
            if error["loc"][0] == key:
                msg.append(error["msg"])

        msg = ", ".join(msg) if len(msg) > 0 else str(exc_info)
        raise ValidationError(f"Invalid metadata: {msg}")
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
