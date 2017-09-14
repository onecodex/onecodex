"""
auth.py
author: @mbiokyle29 / @boydgreenfield
- Adapted from v0 auth.py
"""
from __future__ import print_function
import datetime
import json
import logging
import os
import sys

import click

from onecodex.utils import collapse_user
from onecodex.version import __version__
from onecodex.lib.auth import check_version, fetch_api_key_from_uname


log = logging.getLogger(__name__)
DATE_FORMAT = "%Y-%m-%d %H:%M"
API_KEY_LEN = 32


def login_uname_pwd(server):
    """
    Prompts user for username and password, gets API key from server
    """
    username = click.prompt("Please enter your One Codex  (email)")
    password = click.prompt("Please enter your password (typing will be hidden)",
                            hide_input=True)

    # fetech_api_key expects server to end in /
    if server[-1] != "/":
        server = server + "/"

    # now get the API key
    api_key = fetch_api_key_from_uname(username, password, server)
    return username, api_key


def _login(server, creds_file=None):
    """
    Login main function
    """
    # creds file path setup
    if creds_file is None:
        fp = os.path.expanduser("~/.onecodex")
    else:
        fp = creds_file

    # check if the creds file exists and has an api_key in it
    creds = {}
    if os.path.exists(fp):
        try:
            with open(fp, mode='r') as f:
                creds = json.load(f)
                if 'email' in creds:
                    click.echo('Credentials file already exists ({})'.format(collapse_user(fp)),
                               err=True)
                    return
        except ValueError:
            click.echo("Your ~/.onecodex credentials file appears to be corrupted."  # noqa
                       "Please delete it and re-authorize.", err=True)
            sys.exit(1)

    # else make it
    email, api_key = login_uname_pwd(server)

    if api_key is None:
        click.echo("We could not verify your credentials. Either you mistyped your email "
                   "or password, or your account does not currently have API access. "
                   "Please contact help@onecodex.com if you continue to experience problems.")
        sys.exit(1)

    now = datetime.datetime.now().strftime(DATE_FORMAT)
    creds.update({
        'api_key': api_key,
        'saved_at': now,
        'updated_at': None,
        'email': email,
    })
    with open(fp, mode='w') as f:
        json.dump(creds, f)
    click.echo("Your ~/.onecodex credentials file was successfully created.", err=True)


def _logout(creds_file=None):
    """
    Logout main function, just rm ~/.onecodex more or less
    """

    # creds file path setup
    if creds_file is None:
        fp = os.path.expanduser("~/.onecodex")
    else:
        fp = creds_file

    if os.path.exists(fp):
        # we might not want to do this if there's there are cached schema in it?
        os.remove(fp)
        click.echo("Successfully removed One Codex credentials.", err=True)
        sys.exit(0)
    else:
        click.echo("No One Codex API keys found.", err=True)
        sys.exit(1)


def _silent_login(check_for_update=True):
    """
    This attempts to get an API key when user doesn't pass it. Used only in the CLI.
    """
    fp = os.path.expanduser('~/.onecodex')
    if not os.path.exists(fp):
        return None

    try:
        with open(fp, mode='r') as f:
            creds = json.load(f)

            # Temporary format bridge: Warn user to log out and log back in
            # if they need an email but have an API key. Can remove on next
            # major version bump (to v0.3.0 or later)
            if creds.get('api_key') is not None and creds.get('email') is None:
                click.echo('Your login has expired. Please login again using `onecodex login`')
                return

            if creds.get('email') is None or creds.get('api_key') is None:
                return

        if check_for_update:
            last_update = creds['updated_at'] if creds.get('updated_at') else creds['saved_at']
            diff = datetime.datetime.now() - datetime.datetime.strptime(last_update,
                                                                        DATE_FORMAT)
            if diff.days >= 1:
                # Check and print warning. TODO: Consider moving this to login command as well
                upgrade_required, msg = check_version(__version__, 'https://app.onecodex.com/')
                creds['updated_at'] = datetime.datetime.now().strftime(DATE_FORMAT)
                json.dump(creds, open(fp, mode='w'))
                if upgrade_required:
                    click.echo('\nWARNING: {}\n'.format(msg), err=True)

        return creds['api_key']
    except (KeyError, ValueError):
        click.echo("Your ~/.onecodex credentials "
                   "file appears to be corrupted. "
                   "Please delete it and re-authorize.", err=True)
        sys.exit(1)
