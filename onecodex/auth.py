from __future__ import print_function
import click
import datetime
import errno
from functools import wraps
import json
import logging
import os
import sys

from onecodex.api import Api
from onecodex.lib.auth import check_version, fetch_api_key_from_uname
from onecodex.utils import collapse_user
from onecodex.version import __version__


log = logging.getLogger(__name__)
DATE_FORMAT = "%Y-%m-%d %H:%M"
API_KEY_LEN = 32


def login_uname_pwd(server, api_key=None):
    """
    Prompts user for username and password, gets API key from server
    if not provided.
    """
    username = click.prompt("Please enter your One Codex (email)")
    if api_key is not None:
        return username, api_key

    password = click.prompt("Please enter your password (typing will be hidden)", hide_input=True)

    # now get the API key
    api_key = fetch_api_key_from_uname(username, password, server)
    return username, api_key


def _login(server, creds_file=None, api_key=None, silent=False):
    """
    Login main function
    """
    # fetch_api_key and check_version expect server to end in /
    if server[-1] != "/":
        server = server + "/"

    # creds file path setup
    if creds_file is None:
        creds_file = os.path.expanduser("~/.onecodex")

    # check if the creds file exists and is readable
    if not os.path.exists(creds_file):
        if silent:
            return None

        creds = {}
    elif not os.access(creds_file, os.R_OK):
        click.echo("Please check the permissions on {}".format(collapse_user(creds_file)), err=True)
        sys.exit(1)
    else:
        # it is, so let's read it!
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

            try:
                json.dump(creds, open(creds_file, "w"))
            except Exception as e:
                if e.errno == errno.EACCES:
                    click.echo(
                        "Please check the permissions on {}".format(collapse_user(creds_file)),
                        err=True,
                    )
                    sys.exit(1)
                else:
                    raise

            if upgrade_required:
                click.echo("\nWARNING: {}\n".format(msg), err=True)

        # finally, give the user back what they want (whether silent or not)
        if silent:
            return creds.get("api_key", None)

        click.echo(
            "Credentials file already exists ({}). Logout first.".format(collapse_user(creds_file)),
            err=True,
        )
        return creds.get("email", None)

    # creds_file was not found and we're not silent, so prompt user to login
    email, api_key = login_uname_pwd(server, api_key=api_key)

    if api_key is None:
        click.echo(
            "We could not verify your credentials. Either you mistyped your email "
            "or password, or your account does not currently have API access. "
            "Please contact help@onecodex.com if you continue to experience problems."
        )
        sys.exit(1)

    creds.update(
        {
            "api_key": api_key,
            "saved_at": datetime.datetime.now().strftime(DATE_FORMAT),
            "updated_at": None,
            "email": email,
        }
    )

    try:
        json.dump(creds, open(creds_file, "w"))
    except Exception as e:
        if e.errno == errno.EACCES:
            click.echo("Please check the permissions on {}".format(creds_file), err=True)
            sys.exit(1)
        else:
            raise

    click.echo("Your ~/.onecodex credentials file was successfully created.", err=True)

    return email


def _remove_creds(creds_file=None):
    """
    Remove ~/.onecodex file, returning True if successul or False if the file didn't exist
    """
    if creds_file is None:
        creds_file = os.path.expanduser("~/.onecodex")

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

    return True


def _logout(creds_file=None):
    """
    Logout main function, just rm ~/.onecodex more or less
    """
    if _remove_creds(creds_file=creds_file):
        click.echo("Successfully removed One Codex credentials.", err=True)
        sys.exit(0)
    else:
        click.echo("No One Codex API keys found.", err=True)
        sys.exit(1)


def login_required(fn):
    """Requires login before proceeding, but does not prompt the user to login. Decorator should
    be used only on Click CLI commands.

    Notes
    -----
    Different means of authentication will be attempted in this order:
        1. An API key present in the Click context object from a previous successful authentication.
        2. A bearer token (ONE_CODEX_BEARER_TOKEN) in the environment.
        3. An API key (ONE_CODEX_API_KEY) in the environment.
        4. An API key in the credentials file (~/.onecodex).
    """

    @wraps(fn)
    def login_wrapper(ctx, *args, **kwargs):
        base_url = os.environ.get("ONE_CODEX_API_BASE", "https://app.onecodex.com")

        api_kwargs = {"telemetry": ctx.obj["TELEMETRY"]}

        api_key_prior_login = ctx.obj.get("API_KEY")
        bearer_token_env = os.environ.get("ONE_CODEX_BEARER_TOKEN")
        api_key_env = os.environ.get("ONE_CODEX_API_KEY")
        api_key_creds_file = _login(base_url, silent=True)

        if api_key_prior_login is not None:
            api_kwargs["api_key"] = api_key_prior_login
        elif bearer_token_env is not None:
            api_kwargs["bearer_token"] = bearer_token_env
        elif api_key_env is not None:
            api_kwargs["api_key"] = api_key_env
        elif api_key_creds_file is not None:
            api_kwargs["api_key"] = api_key_creds_file
        else:
            click.echo(
                "The command you specified requires authentication. Please login first.\n", err=True
            )
            ctx.exit()

        ctx.obj["API"] = Api(**api_kwargs)

        return fn(ctx, *args, **kwargs)

    return login_wrapper
