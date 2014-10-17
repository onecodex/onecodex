import datetime
import getpass
import json
import os
import sys


def get_api_key():
    while True:
        api_key = getpass.getpass("API Key (typing will be hidden): ")
        if not api_key:
            continue
        if len(api_key.strip()) != 32:
            print "Your API appears to be too short (it should be 32 characters). Try again."
            continue
        return api_key.strip()


class OneCodexAuth(object):
    """
    Class to store and delete API keys
    at ~/.onecodex. Tested on Linux/OS X only.

    Note this adds a 'credentials' field
    to the parsed arguments.
    """
    def __init__(self, args):
        args.credentials = {}
        fp = os.path.expanduser('~/.onecodex')
        if args.api_key is not None:
            if len(args.api_key) != 32:
                print "Invalid API key length (should be 32 characters)"
                sys.exit(1)

            args.credentials["api_key"] = args.api_key
            args.credentials["saved_at"] = None
            return  # temp login; don't save credentials

        if args.which == 'login':
            if os.path.exists(fp):
                print "Credentials file already exists (~/.onecodex)"
                sys.exit(1)

        if args.which == 'logout':
            if os.path.exists(fp):
                os.remove(fp)
                print "Successfully removed One Codex credentials."
                sys.exit(0)
            else:
                print "No One Codex API keys found."
                sys.exit(1)

        if os.path.exists(fp):
            try:
                args.credentials = json.load(open(fp, mode='r'))
            except ValueError:
                print ("Your ~/.onecodex credentials file appears to be corrupted. "
                       "Please delete it and re-authorize.")
                sys.exit(1)
        else:
            args.credentials["api_key"] = get_api_key()
            args.credentials["saved_at"] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
            json.dump(args.credentials, open(fp, mode='w'))
