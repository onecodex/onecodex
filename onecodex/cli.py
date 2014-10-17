#!/usr/bin/env python
import argparse
import sys
from onecodex.auth import OneCodexAuth
from onecodex import api_v0 as api
from onecodex import version


def hello_world(*args):
    print "Hello world."


class OneCodexArgParser(argparse.ArgumentParser):
    """Parser for the One Codex CLI"""
    HELP = {
        'upload': 'Upload one or more files to the One Codex platform',
        'samples': 'Retrieve uploaded samples',
        'analyses': 'Retrieve performed analyses',
        'references': 'Describe available Reference databses',
        'login': 'Add an API key (saved in ~/.onecodex)',
        'logout': 'Delete your API key (saved in ~/.onecodex)',
        'version': '%(prog)s {} (API {}).\nDocs available at: {}'.format(
            version.VERSION,
            version.API_VERSION,
            version.API_LINK),
        'api_key': 'Manually provide a One Codex Beta API key',
        'threads': 'Use multiple background threads to upload files',
        'file': 'One or more FASTA or FASTQ files to upload. Optionally gzip-compressed.',
        'samples': 'One or more Samples to lookup. If absent returns all Samples.',
        'analyses': 'One or more Analyses to lookup. If absent returns all Analyses.',
        'references': 'One or more References to lookup. If absent returns all current References.',
    }

    def __init__(self, *args, **kwargs):
        super(OneCodexArgParser, self).__init__(*args, **kwargs)
        self._positionals.title = 'One Codex Commands'
        self._optionals.title = 'One Codex Options'
        self.add_argument('--no-pretty-print', dest='pprint',
                          action="store_false", help=self.HELP['api_key'])
        self.add_argument('--threads', action='store_true', help=self.HELP['threads'])
        self.add_argument('--api-key', help=self.HELP['api_key'])
        self.add_argument('--version', action='version',
                          version=self.HELP['version'])

    def _init_subparsers(self):
        """
        Needs to be outside of __init__ to avoid infinite recursion.
        """
        self.subparsers = self.add_subparsers()
        self.upload_parser = self.subparsers.add_parser("upload", help=self.HELP['upload'])
        self.upload_parser.add_argument("file", help=self.HELP['file'], nargs="+")
        self.upload_parser.set_defaults(which="upload")
        self.upload_parser.set_defaults(run=api.upload)

        self.samples_parser = self.subparsers.add_parser("samples", help=self.HELP['samples'])
        self.samples_parser.add_argument("samples", help=self.HELP['samples'], nargs="*")
        self.samples_parser.set_defaults(which="samples")
        self.samples_parser.set_defaults(run=api.samples)

        self.analyses_parser = self.subparsers.add_parser("analyses", help=self.HELP['analyses'])
        self.analyses_parser.add_argument("analyses", help=self.HELP['analyses'], nargs="*")
        self.analyses_parser.set_defaults(which="analyses")
        self.analyses_parser.set_defaults(run=api.analyses)

        self.references_parser = self.subparsers.add_parser("references",
                                                            help=self.HELP['references'])
        self.references_parser.set_defaults(which="references")
        self.references_parser.set_defaults(run=api.references)

        self.logout_parser = self.subparsers.add_parser("logout", help=self.HELP['logout'])
        self.logout_parser.set_defaults(which="logout")

        self.login_parser = self.subparsers.add_parser("login", help=self.HELP['login'])
        self.login_parser.set_defaults(which="login")

    def parse_args(self, *args, **kwargs):
        self._init_subparsers()
        return super(OneCodexArgParser, self).parse_args(*args, **kwargs)


def main(argv=sys.argv[1:]):
    parser = OneCodexArgParser()
    args = parser.parse_args(argv)

    OneCodexAuth(args)  # Check and add credentials
    try:
        args.run(args)
    except AttributeError:
        pass  # For sub-commands w/o run
