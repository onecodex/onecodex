#!/usr/bin/env python
import argparse
import sys
from onecodex.auth import OneCodexAuth
from onecodex import api_v0 as api
from onecodex.api_v0 import DEFAULT_THREADS
from onecodex import version


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
        'pprint': 'Do not pretty-print JSON responses',
        'threads': 'Do not use multiple background threads to upload files',
        'max_threads': 'Specify a different max # of N upload threads (defaults to 4)',
    }

    HELP_SUB = {
        'file': 'One or more FASTA or FASTQ files to upload. Optionally gzip-compressed.',
        'samples': 'One or more Sample IDs to lookup. If absent returns all Samples.',
        'analyses': 'One or more Analysis IDs to lookup. If absent returns all Analyses.',
        'references': ('One or more Reference IDs to lookup. '
                       'If absent returns all current References.'),
        'table': 'Get a JSON array of the Analysis results table.',
        'raw': ('Output path or directory for a .tsv file with '
                'the raw read-level results. Defaults to original filename '
                'in the current working directory. Note if not specifying a '
                'download path this argument must come last, e.g.: '
                '`onecodex analyses <uuid> --raw`.')
    }

    def __init__(self, *args, **kwargs):
        kwargs['prog'] = 'onecodex'
        super(OneCodexArgParser, self).__init__(*args, **kwargs)
        self._positionals.title = 'One Codex Commands'
        self._optionals.title = 'One Codex Options'
        self.add_argument('--no-pretty-print', dest='pprint',
                          action="store_false", help=self.HELP['pprint'])
        self.add_argument('--no-threads', dest='threads',
                          action='store_false', help=self.HELP['threads'])
        self.add_argument('--max-threads', default=DEFAULT_THREADS,
                          type=int, help=self.HELP['max_threads'],
                          metavar="N")
        self.add_argument('--api-key', help=self.HELP['api_key'])
        self.add_argument('--version', action='version',
                          version=self.HELP['version'])

    def _init_subparsers(self):
        """
        Needs to be outside of __init__ to avoid infinite recursion.
        """
        self.subparsers = self.add_subparsers()
        self.upload_parser = self.subparsers.add_parser("upload", help=self.HELP['upload'])
        self.upload_parser.add_argument("file", help=self.HELP_SUB['file'], nargs="+")
        self.upload_parser.set_defaults(which="upload")
        self.upload_parser.set_defaults(run=api.upload)

        self.samples_parser = self.subparsers.add_parser("samples", help=self.HELP['samples'])
        self.samples_parser.add_argument("samples", help=self.HELP_SUB['samples'], nargs="*")
        self.samples_parser.set_defaults(which="samples")
        self.samples_parser.set_defaults(run=api.samples)

        self.analyses_parser = self.subparsers.add_parser("analyses", help=self.HELP['analyses'])
        self.analyses_parser.add_argument("analyses", help=self.HELP_SUB['analyses'], nargs="*")
        self.analyses_parser.add_argument("--table", help=self.HELP_SUB['table'],
                                          action="store_true")
        self.analyses_parser.add_argument("--raw", help=self.HELP_SUB['raw'],
                                          metavar="RAW_DL_PATH", default=None,
                                          const=".", type=str, nargs="?")
        self.analyses_parser.set_defaults(which="analyses")
        self.analyses_parser.set_defaults(run=api.analyses)

        self.references_parser = self.subparsers.add_parser("references",
                                                            help=self.HELP['references'])
        self.references_parser.add_argument("references",
                                            help=self.HELP_SUB['references'], nargs="*")
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


if __name__ == "__main__":
    main()
