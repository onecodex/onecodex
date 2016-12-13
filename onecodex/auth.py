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
from onecodex.lib.auth import fetch_api_key_from_uname

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
    return api_key


def _login(server, check_for_update=True, creds_file=None):
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
                if 'api_key' in creds:
                    click.echo('Credentials file already exists ({})'.format(collapse_user(fp)),
                               err=True)
                    return
        except ValueError:
            click.echo("Your ~/.onecodex credentials file appears to be corrupted."  # noqa
                       "Please delete it and re-authorize.", err=True)
            sys.exit(1)

    # else make it
    api_key = login_uname_pwd(server)

    if api_key is None:
        click.echo("We could not verify your credentials. Either you mistyped your email "
                   "or password, or your account does not currently have API access. "
                   "Please contact help@onecodex.com if you continue to experience problems.")
        sys.exit(1)

    now = datetime.datetime.now().strftime(DATE_FORMAT)
    creds.update({'api_key': api_key, 'saved_at': now, 'updated_at': None})
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


def _silent_login():
    """
    This attempts to get an API key when user doesn't pass it
    """
    fp = os.path.expanduser('~/.onecodex')

    if os.path.exists(fp):
        try:
            with open(fp, mode='r') as f:
                creds = json.load(f)
                return creds.get('api_key')
        except ValueError:
            # TODO: NEED A LOGGER CONFIGURED FOR THIS FILE...
            log.error("Your ~/.onecodex credentials "
                      "file appears to be corrupted. "
                      "Please delete it and re-authorize.")
            sys.exit(1)
