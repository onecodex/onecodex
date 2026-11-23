import datetime
import os
import socket
import time

import filelock
import requests
from pydantic import BaseModel, ValidationError, field_serializer, field_validator, model_validator

# Trigger a warning if the token expires in less than the following
TOKEN_RENEWAL_WINDOW = datetime.timedelta(hours=24)

DEVICE_NAME_HEADER = "X-OneCodex-Device-Name"

MAX_POLL_INTERVAL_SECONDS = 60

# Logging out shouldn't hang on an unreachable server
LOGOUT_TIMEOUT_SECONDS = 5

# Format of the timestamps in the credentials file
DATE_FORMAT = "%Y-%m-%d %H:%M"


class DeviceLoginError(Exception):
    """Browser-based login was denied, expired, or failed."""


class InvalidCredentialsError(Exception):
    """The credentials file is not valid JSON, or its contents are not valid credentials."""


class BearerTokenAuth(requests.auth.AuthBase):
    """Attache Bearer Auth headers to the given Request object."""

    def __init__(self, token):
        self.token = token

    def __call__(self, request):
        request.headers["Authorization"] = "Bearer {}".format(self.token)
        return request


class DeviceAuthorization(BaseModel):
    """A pending browser-based login, returned by `start_device_login`."""

    device_code: str
    user_code: str
    verification_uri: str
    verification_uri_complete: str
    expires_in: int
    interval: int


class DeviceToken(BaseModel):
    """An approved browser-based login, returned by `poll_device_token`."""

    access_token: str
    token_type: str
    expires_in: int
    email: str


class OneCodexCredentials(BaseModel):
    """The contents of the credentials file (~/.onecodex).

    Holds either a bearer token from browser-based login or an `api_key`.
    """

    email: str | None = None
    api_key: str | None = None
    bearer_token: str | None = None
    expires_at: datetime.datetime | None = None
    # Local time
    saved_at: datetime.datetime | None = None
    updated_at: datetime.datetime | None = None

    @field_validator("expires_at")
    @classmethod
    def assume_utc(cls, value: datetime.datetime | None) -> datetime.datetime | None:
        if value is not None and value.tzinfo is None:
            return value.replace(tzinfo=datetime.timezone.utc)
        return value

    @model_validator(mode="after")
    def require_token_expiry(self) -> "OneCodexCredentials":
        if self.bearer_token and self.expires_at is None:
            raise ValueError("expires_at is required with bearer token")
        return self

    @field_serializer("expires_at")
    def serialize_expires_at(self, value: datetime.datetime | None) -> str | None:
        return None if value is None else value.isoformat(timespec="seconds")

    @field_serializer("saved_at", "updated_at")
    def serialize_timestamp(self, value: datetime.datetime | None) -> str | None:
        return None if value is None else value.strftime(DATE_FORMAT)

    @property
    def time_left(self) -> datetime.timedelta | None:
        """Return the time left until the login token expires, or None for no token."""
        if not self.bearer_token or self.expires_at is None:
            return None
        return self.expires_at - datetime.datetime.now(datetime.timezone.utc)

    @classmethod
    def load(cls, path: str) -> "OneCodexCredentials | None":
        """Read the credentials file at `path`, or return None if there isn't one.

        Raises InvalidCredentialsError if the file doesn't hold valid credentials.
        """
        if not os.path.exists(path):
            return None

        with filelock.FileLock("{}.lock".format(path)):
            with open(path, "r") as fp:
                content = fp.read()

        try:
            return cls.model_validate_json(content)
        except ValidationError as e:
            raise InvalidCredentialsError("Invalid credentials file: {}".format(path)) from e

    def save(self, path: str) -> None:
        """Write the credentials file at `path`, readable only by the current user."""
        with filelock.FileLock("{}.lock".format(path)):
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
            with os.fdopen(fd, "w") as fp:
                fp.write(self.model_dump_json(exclude_none=True))


def _device_name() -> str:
    """Return this machine's hostname, e.g. "MacBook-Pro-2", or "" if it has none."""
    # Strip ".local" added by macOS
    return socket.gethostname().removesuffix(".local").encode("ascii", "ignore").decode()


def start_device_login(server_url: str) -> DeviceAuthorization:
    """Start a browser-based login (OAuth 2.0 device auth grant)."""
    device_name = _device_name()
    headers = {DEVICE_NAME_HEADER: device_name} if device_name else {}
    resp = requests.post(server_url + "api/cli/device/code", headers=headers, timeout=30)
    resp.raise_for_status()
    return DeviceAuthorization.model_validate(resp.json())


def poll_device_token(server_url: str, authorization: DeviceAuthorization) -> DeviceToken:
    """Wait for the user to approve a browser-based login in their browser.

    Raises DeviceLoginError if the user denies the request, it expires, or the server returns an error.
    """
    interval = authorization.interval
    deadline = time.monotonic() + authorization.expires_in
    while time.monotonic() < deadline:
        time.sleep(interval)
        try:
            resp = requests.post(
                server_url + "api/cli/device/token",
                json={"device_code": authorization.device_code},
                timeout=15,
            )
        except requests.exceptions.Timeout:
            # Back off exponentially
            interval = min(interval * 2, MAX_POLL_INTERVAL_SECONDS)
            continue
        except requests.exceptions.ConnectionError:
            continue  # keep waiting through network hiccups

        if resp.status_code == 200:
            return DeviceToken.model_validate(resp.json())

        try:
            error = resp.json().get("error")
        except ValueError:
            error = None

        if error == "authorization_pending":
            continue
        if error == "slow_down" or resp.status_code == 429:
            interval += 5
            continue
        if error == "access_denied":
            raise DeviceLoginError("The login request was denied.")
        if error == "expired_token":
            break
        raise DeviceLoginError(
            "Unexpected response from One Codex ({}): {}".format(resp.status_code, error)
        )

    raise DeviceLoginError("The login request expired. Please run `onecodex login` again.")


def revoke_token(server_url: str, bearer_token: str) -> bool:
    """Ask One Codex to revoke the bearer token.

    Returns False if the server can't be reached or returns an error.
    """
    try:
        resp = requests.post(
            server_url + "api/cli/logout",
            auth=BearerTokenAuth(bearer_token),
            timeout=LOGOUT_TIMEOUT_SECONDS,
        )
    except requests.exceptions.RequestException:
        return False
    return resp.ok


def check_version(version, server):
    """Check if the current CLI version is supported by the One Codex backend.

    Parameters
    ----------
    version : `string`
        Current version of the One Codex client library
    server : `string`
        Complete URL to One Codex server, e.g., https://app.onecodex.com

    Returns
    -------
    `tuple` containing two values:
        - True if the user *must* upgrade their software, otherwise False
        - An error message if the user should upgrade, otherwise None.
    """

    def version_inadequate(client_version, server_version):
        """Check for version inequality.

        Could use python package `semver` if we need more precise checks in edge cases, but this
        generally works for now.
        """
        client_version = tuple([int(x) for x in client_version.split("-")[0].split(".")])
        server_version = tuple([int(x) for x in server_version.split(".")])
        return client_version < server_version

    # this will probably live on /api/v0 forever for compat with older CLI versions
    data = requests.post(server + "api/v0/check_for_cli_update", data={"version": version})

    if data.status_code != 200:
        return False, "Error connecting to server"

    data = data.json()
    latest_version = data["latest_version"]

    if version_inadequate(version, latest_version):
        return (
            True,
            (
                "Please upgrade your client to the latest version (v{}) using the command "
                "`pip install --upgrade onecodex`".format(latest_version)
            ),
        )

    return False, None
