class OneCodexException(Exception):
    pass


class MethodNotSupported(OneCodexException):
    """
    The object does not support this operation.
    """
    pass


class PermissionDenied(OneCodexException):
    pass


class ServerError(OneCodexException):
    pass


class UnboundObject(OneCodexException):
    """
    To use against the One Codex server, all classes must be derived from an Api instance.
    """
    pass


class ValidationError(OneCodexException):
    pass


class ValidationWarning(Warning):
    pass


class UploadException(Exception):
    """
    An exception for when things go wrong with uploading
    """
    pass


def process_api_error(error):
    error_code = error['status']
    error_message = ''
    if error_code == 402:
        error_message = api_insufficient_funds_error_message()
    else:
        error_message = api_default_error_message()
    raise UploadException(error_message)


def api_default_error_message():
    return "The attempt to initiate your upload failed. Please make sure you are logged in (`onecodex login`) and try again. If you continue to experience problems, contact us at help@onecodex.com for assistance."


def api_insufficient_funds_error_message():
    return "Please add a payment method to upload more samples. If you continue to experience problems, contact us at help@onecodex.com for assistance."
