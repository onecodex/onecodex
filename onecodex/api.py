"""
api.py
author: @mbiokyle29

One Codex API
"""
from __future__ import print_function
from datetime import datetime
import json
import logging
import os
import warnings

from potion_client import Client as PotionClient
from potion_client.converter import PotionJSONSchemaDecoder, PotionJSONDecoder, PotionJSONEncoder
from potion_client.utils import upper_camel_case
from requests.auth import HTTPBasicAuth

from onecodex.lib.auth import BearerTokenAuth
from onecodex.models import _model_lookup
from onecodex.utils import ModuleAlias, get_raven_client

log = logging.getLogger(__name__)


class Api(object):
    """
    This is the base One Codex Api object class. It instantiates a Potion-Client
        object under the hood for making requests.
    """

    def __init__(self, api_key=None,
                 bearer_token=None, cache_schema=False,
                 base_url=None, telemetry=None,
                 schema_path='/api/v1/schema'):

        if base_url is None:
            base_url = os.environ.get('ONE_CODEX_API_BASE', 'https://app.onecodex.com')
            if base_url != 'https://app.onecodex.com':
                warnings.warn("Using base API URL: {}".format(base_url))

        self._req_args = {}
        self._base_url = base_url
        self._schema_path = schema_path

        # Attempt to automatically fetch API key from
        # ~/.onecodex file, API key, or bearer token environment vars
        # *if and only if* no auth is explicitly passed to Api
        #
        # TODO: Consider only doing this if an add'l env var like
        #       'ONE_CODEX_AUTO_LOGIN' or similar is set.
        if api_key is None and bearer_token is None:
            try:
                api_key = json.load(open(os.path.expanduser('~/.onecodex')))['api_key']
            except Exception:
                pass
            if api_key is None:
                api_key = os.environ.get('ONE_CODEX_API_KEY')
            if bearer_token is None:
                bearer_token = os.environ.get('ONE_CODEX_BEARER_TOKEN')

        if bearer_token:  # prefer bearer token where available
            self._req_args['auth'] = BearerTokenAuth(bearer_token)
        elif api_key:
            self._req_args['auth'] = HTTPBasicAuth(api_key, '')

        # Create client instance
        self._client = ExtendedPotionClient(self._base_url, schema_path=self._schema_path,
                                            fetch_schema=False, **self._req_args)
        self._client._fetch_schema(cache_schema=cache_schema)
        self._session = self._client.session
        self._copy_resources()

        # Optionally configure Raven
        if telemetry is True or (telemetry is None and os.environ.get('ONE_CODEX_AUTO_TELEMETRY', False)):
            self._raven_client = get_raven_client(user_context={'email': self._fetch_account_email()})
            self._telemetry = True
        else:
            self._raven_client = None
            self._telemetry = False

        # Try to import and copy key modules onto the Api object
        for module_name in ['onecodex.viz', 'onecodex.distance']:
            module = ModuleAlias(module_name)
            if module._imported:
                setattr(self, module._name, module)

    def _fetch_account_email(self):
        creds_fp = os.path.expanduser('~/.onecodex')
        if os.path.exists(creds_fp):
            with open(creds_fp) as f:
                creds = json.load(f)
                if 'email' in creds:
                    return creds['email']
        return os.environ.get('ONE_CODEX_USER_EMAIL', os.environ.get('ONE_CODEX_USER_UUID'))

    def _copy_resources(self):
        """
        Copy all of the resources over to the toplevel client

        -return: populates self with a pointer to each ._client.Resource
        """

        for resource in self._client._resources:
            # set the name param, the keys now have / in them
            potion_resource = self._client._resources[resource]

            try:
                oc_cls = _model_lookup[resource]
                oc_cls._api = self
                oc_cls._resource = potion_resource
                setattr(self, oc_cls.__name__, oc_cls)
            except KeyError:  # Ignore resources we don't explicitly model
                pass


class ExtendedPotionClient(PotionClient):
    """
    An extention of the PotionClient that caches schema
    """
    DATE_FORMAT = "%Y-%m-%d %H:%M"
    SCHEMA_SAVE_DURATION = 1  # day

    def fetch(self, uri, cls=PotionJSONDecoder, **kwargs):
        if uri in self._cached_schema:
            return self._cached_schema[uri]
        return super(ExtendedPotionClient, self).fetch(uri, cls=cls, **kwargs)

    def _fetch_schema(self, cache_schema=False, creds_file=None):
        self._cached_schema = {}
        creds_fp = os.path.expanduser('~/.onecodex') if creds_file is None else creds_file

        if os.path.exists(creds_fp):
            creds = json.load(open(creds_fp, 'r'))
        else:
            creds = {}

        schema = None
        serialized_schema = None
        if cache_schema:
            # Determine if we need to update
            schema_update_needed = True
            last_update = creds.get('schema_saved_at')
            if last_update is not None:
                last_update = datetime.strptime(last_update, self.DATE_FORMAT)
                time_diff = datetime.now() - last_update
                schema_update_needed = time_diff.days > self.SCHEMA_SAVE_DURATION

            if not schema_update_needed:
                # get the schema from the credentials file (as a string)
                serialized_schema = creds.get('schema')

        if serialized_schema is not None:
            # Catch schema caching issues and fall back to remote URL
            try:
                base_schema = serialized_schema.pop(self._schema_url)
                schema = json.loads(base_schema, cls=PotionJSONSchemaDecoder,
                                    referrer=self._schema_url, client=self)

                for route, route_schema in serialized_schema.items():
                    object_schema = json.loads(route_schema, cls=PotionJSONSchemaDecoder,
                                               referrer=self._schema_url, client=self)
                    self._cached_schema[route] = object_schema
            except KeyError:  # Caches issue with schema_url not existing
                pass

        if schema is None:
            # if the schema wasn't cached or if it was expired, get it anew
            schema = self.session.get(self._schema_url).json(cls=PotionJSONSchemaDecoder,
                                                             referrer=self._schema_url,
                                                             client=self)
            if cache_schema:
                # serialize the schemas back out
                creds['schema_saved_at'] = datetime.strftime(datetime.now(), self.DATE_FORMAT)

                # serialize the main schema
                serialized_schema = {}
                serialized_schema[self._schema_url] = json.dumps(schema, cls=PotionJSONEncoder)

                # serialize the object schemas
                for schema_ref in schema['properties'].values():
                    serialized_schema[schema_ref._uri] = json.dumps(schema_ref._properties,
                                                                    cls=PotionJSONEncoder)

                creds['schema'] = serialized_schema
            else:
                if 'schema_saved_at' in creds:
                    del creds['schema_saved_at']
                if 'schema' in creds:
                    del creds['schema']

            # always resave the creds (to make sure we're removing schema if we need to be or
            # saving if we need to do that instead)
            json.dump(creds, open(creds_fp, mode='w'))

        for name, resource_schema in schema['properties'].items():
            class_name = upper_camel_case(name)
            setattr(self, class_name, self.resource_factory(name, resource_schema))
