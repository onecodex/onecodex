from __future__ import print_function
import click
import datetime
import errno
from functools import wraps
import json
import logging
import filelock
import os
import sys
from typing import Any
import webbrowser

import requests

from onecodex.api import Api
from onecodex.lib.auth import (
    TOKEN_RENEWAL_WINDOW,
    DeviceLoginError,
    check_version,
    poll_device_token,
    start_device_login,
    token_time_left,
)
from onecodex.utils import collapse_user
from onecodex.version import __version__


log = logging.getLogger("onecodex")
DATE_FORMAT = "%Y-%m-%d %H:%M"
API_KEY_LEN = 32


def _default_creds_file() -> str:
    return os.path.expanduser("~/.onecodex")


def _server_url(server: str) -> str:
    # fetch_api_key and check_version expect server to end in /
    return server if server.endswith("/") else server + "/"


def _format_time_left(time_left: datetime.timedelta) -> str:
    hours = int(time_left.total_seconds() // 3600)
    if hours >= 1:
        return "{} hour{}".format(hours, "" if hours == 1 else "s")
    minutes = max(int(time_left.total_seconds() // 60), 1)
    return "{} minute{}".format(minutes, "" if minutes == 1 else "s")


def _store_creds(creds_file: str, creds: dict[str, Any]) -> None:
    try:
        fd = os.open(creds_file, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(fd, "w") as fp:
            json.dump(creds, fp)
    except Exception as e:
        if getattr(e, "errno", None) == errno.EACCES:
            click.echo(
                "Please check the permissions on {}".format(collapse_user(creds_file)), err=True
            )
            sys.exit(1)
        raise


def _read_creds(server: str, creds_file: str | None = None) -> dict[str, Any] | None:
    """Read the credentials file, or return None if there isn't one.

    Also checks once a day whether this client version needs an upgrade.
    """
    server = _server_url(server)
    if creds_file is None:
        creds_file = _default_creds_file()

    if not os.path.exists(creds_file):
        return None

    if not os.access(creds_file, os.R_OK):
        click.echo("Please check the permissions on {}".format(collapse_user(creds_file)), err=True)
        sys.exit(1)

    with filelock.FileLock("{}.lock".format(creds_file)):
        with open(creds_file, "r") as fp:
            try:
                creds = json.load(fp)
            except ValueError:
                click.echo(
                    "Your ~/.onecodex credentials file appears to be corrupted. "  # noqa
                    "Please delete it and re-authorize.",
                    err=True,
                )
                sys.exit(1)

        # check for updates if logged in more than one day ago
        last_update = creds.get("updated_at") or creds.get("saved_at")
        last_update = last_update if last_update else datetime.datetime.now().strftime(DATE_FORMAT)
        diff = datetime.datetime.now() - datetime.datetime.strptime(last_update, DATE_FORMAT)

        if diff.days >= 1:
            # if creds_file is old, check for updates
            upgrade_required, msg = check_version(__version__, server)
            creds["updated_at"] = datetime.datetime.now().strftime(DATE_FORMAT)
            _store_creds(creds_file, creds)

            if upgrade_required:
                click.echo("\nWARNING: {}\n".format(msg), err=True)

    return creds


def _can_open_browser() -> bool:
    # Attempt to determine whether it's a remote session.
    if os.environ.get("SSH_CONNECTION") or os.environ.get("SSH_TTY"):
        return False
    if sys.platform.startswith("linux") and not (
        os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")
    ):
        return False
    return True


def _browser_login(server: str, open_browser: bool = True) -> dict[str, Any]:
    """Log in by approving this CLI in a browser. Returns the new credentials."""
    try:
        authorization = start_device_login(server)
    except requests.exceptions.RequestException as e:
        click.echo("Could not start login with One Codex: {}".format(e), err=True)
        sys.exit(1)

    click.echo(
        "\nTo log in, open this page in your browser and confirm the code {}:\n\n    {}\n".format(
            authorization.user_code, authorization.verification_uri_complete
        ),
        err=True,
    )
    if open_browser and _can_open_browser():
        try:
            webbrowser.open(authorization.verification_uri_complete)
        except webbrowser.Error:
            pass

    click.echo("Waiting for you to log in in your browser...", err=True)
    try:
        token = poll_device_token(server, authorization)
    except DeviceLoginError as e:
        click.echo(str(e), err=True)
        sys.exit(1)

    expires_at = datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(
        seconds=token.expires_in
    )
    return {
        "access_token": token.access_token,
        "expires_at": expires_at.isoformat(timespec="seconds"),
        "email": token.email,
    }


def _login(
    server: str,
    creds_file: str | None = None,
    api_key: str | None = None,
    open_browser: bool = True,
) -> str | None:
    """Login main function. Returns the logged-in user's email."""
    server = _server_url(server)
    if creds_file is None:
        creds_file = _default_creds_file()

    creds = _read_creds(server, creds_file)
    if creds is not None:
        # A login token that has expired, or soon will, can be replaced without logging out.
        # Otherwise, the user needs to log out first
        time_left = token_time_left(creds)
        if time_left is None or time_left > TOKEN_RENEWAL_WINDOW:
            click.echo(
                "Credentials file already exists ({}). Logout first.".format(
                    collapse_user(creds_file)
                ),
                err=True,
            )
            return creds.get("email", None)

    if api_key is not None:
        email = click.prompt("Please enter your One Codex email")
        new_creds = {"api_key": api_key, "email": email}
    else:
        new_creds = _browser_login(server, open_browser=open_browser)

    new_creds.update(
        {
            "saved_at": datetime.datetime.now().strftime(DATE_FORMAT),
            "updated_at": None,
        }
    )

    with filelock.FileLock("{}.lock".format(creds_file)):
        _store_creds(creds_file, new_creds)

    click.echo("Your ~/.onecodex credentials file was successfully created.", err=True)

    return new_creds["email"]


def _remove_creds(creds_file=None):
    """Remove ~/.onecodex file, returning True if successul or False if the file didn't exist."""
    if creds_file is None:
        creds_file = _default_creds_file()

    try:
        os.remove(creds_file)
    except Exception as e:
        if e.errno == errno.ENOENT:
            return False
        elif e.errno == errno.EACCES:
            click.echo(
                "Please check the permissions on {}".format(collapse_user(creds_file)), err=True
            )
            sys.exit(1)
        else:
            raise

    try:
        os.remove("{}.lock".format(creds_file))
    except (IOError, OSError):
        pass

    return True


def _logout(creds_file=None):
    """Logout main function, just rm ~/.onecodex more or less."""
    if _remove_creds(creds_file=creds_file):
        click.echo("Successfully removed One Codex credentials.", err=True)
        sys.exit(0)
    else:
        click.echo("No One Codex API keys found.", err=True)
        sys.exit(1)


def login_required_with_args(experimental_api=False):
    def decorator(fn):
        """Require login before proceeding, but does not prompt the user to login.

        Decorator should be used only on Click CLI commands.

        Notes
        -----
        Different means of authentication will be attempted in this order:
            1. An API key present in the Click context object from a previous successful authentication.
            2. A bearer token (ONE_CODEX_BEARER_TOKEN) in the environment.
            3. An API key (ONE_CODEX_API_KEY) in the environment.
            4. A login token or API key in the credentials file (~/.onecodex).
        """

        @wraps(fn)
        def login_wrapper(ctx, *args, **kwargs):
            base_url = os.environ.get("ONE_CODEX_API_BASE", "https://app.onecodex.com")

            api_kwargs = {"telemetry": ctx.obj["TELEMETRY"], "load_extensions": False}

            api_key_prior_login = ctx.obj.get("API_KEY")
            bearer_token_env = os.environ.get("ONE_CODEX_BEARER_TOKEN")
            api_key_env = os.environ.get("ONE_CODEX_API_KEY")
            creds = _read_creds(base_url) or {}

            if api_key_prior_login is not None:
                api_kwargs["api_key"] = api_key_prior_login
            elif bearer_token_env is not None:
                api_kwargs["bearer_token"] = bearer_token_env
            elif api_key_env is not None:
                api_kwargs["api_key"] = api_key_env
            elif creds.get("access_token"):
                time_left = token_time_left(creds)
                if time_left <= datetime.timedelta(0):
                    click.echo(
                        "Your One Codex login has expired. Please run `onecodex login` again.",
                        err=True,
                    )
                    ctx.exit(1)
                elif time_left < TOKEN_RENEWAL_WINDOW:
                    click.echo(
                        "Your One Codex login expires in {}. Run `onecodex login` to renew "
                        "it.".format(_format_time_left(time_left)),
                        err=True,
                    )
                api_kwargs["bearer_token"] = creds["access_token"]
            elif creds.get("api_key"):
                api_kwargs["api_key"] = creds["api_key"]
            else:
                click.echo(
                    "The command you specified requires authentication. Please login first.\n",
                    err=True,
                )
                ctx.exit(1)

            if experimental_api:
                api_kwargs["experimental"] = True

            ctx.obj["API"] = Api(**api_kwargs)

            return fn(ctx, *args, **kwargs)

        return login_wrapper

    return decorator


login_required = login_required_with_args(experimental_api=False)
login_required_experimental_api = login_required_with_args(experimental_api=True)
