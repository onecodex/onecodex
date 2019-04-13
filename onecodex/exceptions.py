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


class UploadException(Exception):
    """
    An exception for when things go wrong with uploading
    """

    pass


class RetryableUploadException(UploadException):
    pass


def raise_api_error(resp, state=None):
    """Raise an exception with a pretty message in various states of upload"""
    # TODO: Refactor into an Exception class
    error_code = resp.status_code

    if error_code == 402:
        error_message = (
            "Please add a payment method to upload more samples. If you continue to "
            "experience problems, contact us at help@onecodex.com for assistance."
        )
    elif error_code == 403:
        error_message = "Please login to your One Codex account or pass the appropriate API key."
    else:
        try:
            error_json = resp.json()
        except ValueError:
            error_json = {}

        if "msg" in error_json:
            error_message = error_json["msg"].rstrip(".")
        elif "message" in error_json:
            error_message = error_json["message"].rstrip(".")
        else:
            error_message = None

        if state == "init" and not error_message:
            error_message = (
                "Could not initialize upload. Are you logged in? If this problem "
                "continues, please contact help@onecodex.com for assistance."
            )
        elif state == "upload" and not error_message:
            error_message = (
                "File could not be uploaded. If this problem continues, please contact "
                "help@onecodex.com for assistance."
            )
        elif state == "callback" and not error_message:
            error_message = (
                "Callback could not be completed. If this problem continues, please "
                "contact help@onecodex.com for assistance."
            )

    if error_message is None:
        error_message = "Upload failed. Please contact help@onecodex.com for assistance."

    raise UploadException(error_message)


def raise_connectivity_error(file_name):
    # TODO: This is really a general NonRetryableUploadError
    #       with a default msg. Refactor.
    raise UploadException(
        "The command line client is experiencing connectivity issues and "
        "cannot complete the upload of {} at this time. Please try again "
        "later. If the problem persists, contact us at help@onecodex.com "
        "or assistance.".format(file_name)
    )
