from __future__ import print_function
import click
import datetime
import errno
from functools import wraps
import logging
import os
import sys
import webbrowser

import requests

from onecodex.api import Api
from onecodex.lib.auth import (
    TOKEN_RENEWAL_WINDOW,
    DeviceLoginError,
    InvalidCredentialsError,
    OneCodexCredentials,
    check_version,
    poll_device_token,
    revoke_token,
    start_device_login,
)
from onecodex.utils import collapse_user
from onecodex.version import __version__


log = logging.getLogger("onecodex")
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


def _store_creds(creds_file: str, creds: OneCodexCredentials) -> None:
    try:
        creds.save(creds_file)
    except Exception as e:
        if getattr(e, "errno", None) == errno.EACCES:
            click.echo(
                "Please check the permissions on {}".format(collapse_user(creds_file)), err=True
            )
            sys.exit(1)
        raise


def _read_creds(server: str, creds_file: str | None = None) -> OneCodexCredentials | None:
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

    try:
        creds = OneCodexCredentials.load(creds_file)
    except InvalidCredentialsError:
        click.echo(
            "Your ~/.onecodex credentials file appears to be corrupted. "  # noqa
            "Please delete it and re-authorize.",
            err=True,
        )
        sys.exit(1)

    if creds is None:
        return None

    # check for updates if logged in more than one day ago
    now = datetime.datetime.now()
    last_update = creds.updated_at or creds.saved_at or now

    if (now - last_update).days >= 1:
        # if creds_file is old, check for updates
        upgrade_required, msg = check_version(__version__, server)
        creds.updated_at = now
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


def _browser_login(server: str, open_browser: bool = True) -> OneCodexCredentials:
    """Log in by approving this CLI in a browser. Returns the new credentials."""
    try:
        authorization = start_device_login(server)
    except requests.exceptions.RequestException as e:
        click.echo("Could not start login with One Codex: {}".format(e), err=True)
        sys.exit(1)

    user_code = click.style(authorization.user_code, fg="green", bold=True)
    verification_uri = click.style(authorization.verification_uri_complete, fg="cyan")
    if open_browser and _can_open_browser():
        try:
            webbrowser.open(authorization.verification_uri_complete)
            click.echo(
                f"Opening your browser to authorize the login...\n\nIf it doesn't open, visit:\n    {verification_uri}\n",
                err=True,
            )
        except webbrowser.Error:
            pass
    else:
        click.echo(f"Please open: {verification_uri}\n", err=True)

    click.echo(
        f"Confirm that the browser shows: {user_code}\n\nWaiting for you to authorize this login in your browser...",
        err=True,
    )
    try:
        token = poll_device_token(server, authorization)
    except DeviceLoginError as e:
        click.echo(str(e), err=True)
        sys.exit(1)

    expires_at = datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(
        seconds=token.expires_in
    )
    return OneCodexCredentials(
        bearer_token=token.access_token,
        expires_at=expires_at,
        email=token.email,
        saved_at=datetime.datetime.now(),
    )


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
        time_left = creds.time_left
        if time_left is None or time_left > TOKEN_RENEWAL_WINDOW:
            click.echo(
                "Credentials file already exists ({}). Logout first.".format(
                    collapse_user(creds_file)
                ),
                err=True,
            )
            return creds.email

    if api_key is not None:
        email = click.prompt("Please enter your One Codex email")
        new_creds = OneCodexCredentials(
            api_key=api_key, email=email, saved_at=datetime.datetime.now()
        )
    else:
        new_creds = _browser_login(server, open_browser=open_browser)

    _store_creds(creds_file, new_creds)

    click.echo("Your ~/.onecodex credentials file was successfully created.", err=True)

    return new_creds.email


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


def _logout(server: str, creds_file: str | None = None) -> None:
    """Logout main function: revoke the bearer token and rm ~/.onecodex."""
    server = _server_url(server)
    if creds_file is None:
        creds_file = _default_creds_file()

    try:
        creds = OneCodexCredentials.load(creds_file)
    except (InvalidCredentialsError, OSError):
        creds = None  # a corrupted or unreadable file is still removed below

    if creds is not None and creds.bearer_token and creds.time_left > datetime.timedelta(0):
        # Skip revoking expired tokens
        if not revoke_token(server, creds.bearer_token):
            expires_at = creds.expires_at.astimezone(datetime.timezone.utc)
            click.echo(
                "Your login token was not revoked on the One Codex server. It remains valid "
                "until {}.".format(expires_at.strftime("%Y-%m-%d %H:%M UTC")),
                err=True,
            )

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
            creds = _read_creds(base_url) or OneCodexCredentials()

            if api_key_prior_login is not None:
                api_kwargs["api_key"] = api_key_prior_login
            elif bearer_token_env is not None:
                api_kwargs["bearer_token"] = bearer_token_env
            elif api_key_env is not None:
                api_kwargs["api_key"] = api_key_env
            elif creds.bearer_token:
                time_left = creds.time_left
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
                api_kwargs["bearer_token"] = creds.bearer_token
            elif creds.api_key:
                api_kwargs["api_key"] = creds.api_key
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
