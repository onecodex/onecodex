from __future__ import print_function
import datetime
import getpass
import json
import os
import sys
from onecodex.api_v0 import get_update_message
from onecodex.helpers import stderr


DATE_FORMAT = "%Y-%m-%d %H:%M"


def get_api_key():
    while True:
        api_key = getpass.getpass("API Key (typing will be hidden): ")
        if not api_key:
            continue
        if len(api_key.strip()) != 32:
            print("Your API appears to be too short (it should be 32 characters). Try again.")
            continue
        return api_key.strip()


class OneCodexAuth(object):
    """
    Class to store and delete API keys
    at ~/.onecodex. Tested on Linux/OS X only.

    Note this adds a 'credentials' field
    to the parsed arguments.
    """
    def __init__(self, args, check_for_update=True, creds_file=None):
        args.credentials = {}
        if creds_file is None:
            fp = os.path.expanduser('~/.onecodex')
        else:
            fp = creds_file
        if args.api_key is not None:
            if len(args.api_key) != 32:
                stderr("Invalid API key length (should be 32 characters)")
                sys.exit(1)

            args.credentials["api_key"] = args.api_key
            args.credentials["saved_at"] = None
            args.credentials["updated_at"] = None
            self._check_for_update(args, fp)
            return  # temp login; don't save credentials

        if args.which == 'login':
            if os.path.exists(fp):
                stderr("Credentials file already exists (~/.onecodex)")
                sys.exit(1)

        if args.which == 'logout':
            if os.path.exists(fp):
                os.remove(fp)
                print("Successfully removed One Codex credentials.")
                sys.exit(0)
            else:
                stderr("No One Codex API keys found.")
                sys.exit(1)

        if os.path.exists(fp):
            try:
                args.credentials = json.load(open(fp, mode='r'))
            except ValueError:
                stderr("Your ~/.onecodex credentials file appears to be corrupted. "
                       "Please delete it and re-authorize.")
                sys.exit(1)
        else:
            args.credentials["api_key"] = get_api_key()
            args.credentials["saved_at"] = datetime.datetime.now().strftime(DATE_FORMAT)
            args.credentials["updated_at"] = None
            json.dump(args.credentials, open(fp, mode='w'))

        # Finally perform a version check as needed
        if check_for_update:
            self._check_for_update(args, fp)

    def _check_for_update(self, args, fp):
        time_diff = None
        if args.credentials["updated_at"] is not None:
            last_update = datetime.datetime.strptime(args.credentials["updated_at"],
                                                     DATE_FORMAT)
            time_diff = datetime.datetime.now() - last_update

        if time_diff is None or time_diff.days >= 1:
            msg = get_update_message()
            if msg:
                stderr(msg)

            if args.api_key is None:
                args.credentials["updated_at"] = datetime.datetime.now().strftime(DATE_FORMAT)
                json.dump(args.credentials, open(fp, mode='w'))
