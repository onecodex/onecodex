"""
test_api.py
author: @mbiokyle29
"""
import json
from pkg_resources import resource_string
import responses
import unittest

from onecodex import Api


# TODO: Refactor this all into test_api_models.py
class TestVanillaApi(unittest.TestCase):

    @responses.activate
    def setUp(self):
        self.schema_url = "http://localhost:3000/api/v1/schema"
        self.schema_file = resource_string(__name__, 'data/schema.json').decode('utf-8')
        responses.add(responses.GET, self.schema_url,
                      json=json.loads(self.schema_file))
        self.api = Api(cache_schema=False, base_url="http://localhost:3000",
                       api_key='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')

    def test_attrs_copied(self):

        attrs = {"Samples", "Tags", "Analyses", "Classifications",
                 "Markerpanels", "Metadata", "Users"}

        not_yet_implemented = {"Markerpanels", "Metadata"}

        for attr in (attrs - not_yet_implemented):
            self.assertTrue(hasattr(self.api, attr))


class TestExtendedApi(unittest.TestCase):

    @responses.activate
    def setUp(self):

        self.schema_url = "http://localhost:3000/api/v1/schema"
        self.schema_file = resource_string(__name__, 'data/schema.json').decode('utf-8')
        responses.add(responses.GET, self.schema_url,
                      json=json.loads(self.schema_file))
        self.api = Api(cache_schema=False, base_url="http://localhost:3000",
                       api_key='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')

        self.uuid_one = "4a668ac6daf74364"
        self.uuid_two = "d61f459d077d4cd0"

        self.classification_one = resource_string(__name__, "data/api/{}.json".format(self.uuid_one)).decode('utf-8')  # noqa
        self.table_one = resource_string(__name__, "data/api/{}_table.json".format(self.uuid_one)).decode('utf-8')  # noqa

        self.classification_two = resource_string(__name__, "data/api/{}.json".format(self.uuid_two)).decode('utf-8')  # noqa
        self.table_two = resource_string(__name__, "data/api/{}_table.json".format(self.uuid_two)).decode('utf-8')  # noqa

        self.classification_one_url = "http://localhost:3000/api/v1/classifications/4a668ac6daf74364"
        self.classification_one_table_url = "http://localhost:3000/api/v1/classifications/4a668ac6daf74364/results"

        self.classification_two_url = "http://localhost:3000/api/v1/classifications/d61f459d077d4cd0"
        self.classification_two_table_url = "http://localhost:3000/api/v1/classifications/d61f459d077d4cd0/results"

    @responses.activate
    def test_classification_abundance_data(self):

        responses.add(responses.GET, self.classification_one_url,
                      json=json.loads(self.classification_one))
        responses.add(responses.GET, self.classification_one_table_url,
                      json=json.loads(self.table_one))

        # make sure the potion decorated method is there
        classification = self.api.Classifications.get(self.uuid_one)
        self.assertTrue(hasattr(classification, "table"))
        self.assertTrue(hasattr(classification.table, '__call__'))
        table = classification.abundances()

        for index_col in ['readcount', 'name', 'rank', 'tax_id']:
            self.assertTrue(index_col in table.columns.values)

        # from the json flat file
        expected_rows = 79
        self.assertEqual(len(table.index), expected_rows)

        # make sure we can query with ids
        test_ids = [1214158, 1357, 1397276, 1288394]
        subset_table = classification.abundances(ids=test_ids)
        self.assertEqual(len(subset_table.index), len(test_ids))
        self.assertEqual(test_ids, list(subset_table['tax_id']))


class TestQueryInterface(unittest.TestCase):

    @responses.activate
    def setUp(self):
        self.schema_url = "http://localhost:3000/api/v1/schema"
        self.schema_file = resource_string(__name__, 'data/schema.json').decode('utf-8')
        responses.add(responses.GET, self.schema_url,
                      json=json.loads(self.schema_file))
        self.api = Api(cache_schema=False, base_url="http://localhost:3000",
                       api_key='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')

        self.num_complete = 2
        self.num_incomplete = 1

        self.complete_analyses = resource_string(
            __name__, 'data/api/complete_analyses.json').decode('utf-8')

        self.incomplete_analyses = resource_string(
            __name__, 'data/api/incomplete_analyses.json').decode('utf-8')

        self.sample = resource_string(__name__, 'data/api/sample.json').decode('utf-8')
        self.metadata = resource_string(__name__, 'data/api/metadata.json').decode('utf-8')

        self.uuid_one = '1'
        self.analysis_one = resource_string(
            __name__, 'data/api/complete_analysis_one.json').decode('utf-8')

        self.uuid_two = '2'
        self.analysis_two = resource_string(
            __name__, 'data/api/complete_analysis_two.json').decode('utf-8')

        self.uuid_three = '3'
        self.analysis_three = resource_string(
            __name__, 'data/api/incomplete_analysis_one.json').decode('utf-8')

        # URLS for responses
        self.analysis_one_url = 'http://localhost:3000/api/v1/analyses/1'
        self.analysis_two_url = 'http://localhost:3000/api/v1/analyses/2'
        self.analysis_three_url = 'http://localhost:3000/api/v1/analyses/3'
        self.complete_url = 'http://localhost:3000/api/v1/analyses?per_page=200&where=%7B%22complete%22%3A+true%7D&page=1'
        self.complete_url_new = 'http://localhost:3000/api/v1/analyses?sort=%7B%7D&per_page=200&where=%7B%22complete%22%3A+true%7D&page=1'
        self.incomplete_url = 'http://localhost:3000/api/v1/analyses?per_page=200&where=%7B%22complete%22%3A+false%7D&page=1'
        self.incomplete_url_new = 'http://localhost:3000/api/v1/analyses?sort=%7B%7D&per_page=200&where=%7B%22complete%22%3A+false%7D&page=1'
        self.sample_url = 'http://localhost:3000/api/v1/samples/d66e901ea9854f1f'
        self.metadata_url = 'http://localhost:3000/api/v1/metadata/50fd2f0513d34372'

    @responses.activate
    def test_get(self):

        self.assertTrue(hasattr(self.api.Analyses, 'get'))
        self.assertTrue(hasattr(self.api.Analyses.get, '__call__'))

        responses.add(responses.GET, self.analysis_one_url,
                      json=json.loads(self.analysis_one))

        # old way, no longer supported
        # with pytest.raises(AttributeError):
        #     old_analysis_one = self.api.Analyses(self.uuid_one)
        old_analysis_one = self.api.Analyses.get(self.uuid_one)
        self.assertTrue(isinstance(old_analysis_one, self.api.Analyses))
        self.assertTrue(hasattr(old_analysis_one, 'complete'))
        self.assertTrue(old_analysis_one.complete)

        # new way
        new_analysis_one = self.api.Analyses.get(self.uuid_one)
        self.assertTrue(isinstance(new_analysis_one, self.api.Analyses))
        self.assertTrue(hasattr(new_analysis_one, 'complete'))
        self.assertTrue(new_analysis_one.complete)

        self.assertEqual(old_analysis_one, new_analysis_one)

    @responses.activate
    def test_where(self):

        self.assertTrue(hasattr(self.api.Analyses, 'where'))
        self.assertTrue(hasattr(self.api.Analyses.where, '__call__'))

        responses.add(responses.GET, self.complete_url,
                      json=json.loads(self.complete_analyses),
                      match_querystring=True)

        # old way
        old_complete_analyses = self.api.Analyses._resource.instances(where={'complete': True}, per_page=200)  # noqa
        self.assertEqual(len(old_complete_analyses), self.num_complete)

        responses.add(responses.GET, self.complete_url_new,
                      json=json.loads(self.complete_analyses),
                      match_querystring=True)

        # new way
        new_complete_analyses = self.api.Analyses.where(complete=True)
        self.assertEqual(len(new_complete_analyses), self.num_complete)

        responses.add(responses.GET, self.incomplete_url,
                      json=json.loads(self.incomplete_analyses),
                      match_querystring=True)

        # old way
        old_incomplete_analyses = self.api.Analyses._resource.instances(where={'complete': False}, per_page=200)  # noqa
        self.assertEqual(len(old_incomplete_analyses), self.num_incomplete)

        responses.add(responses.GET, self.incomplete_url_new,
                      json=json.loads(self.incomplete_analyses),
                      match_querystring=True)

        # new way
        new_incomplete_analyses = self.api.Analyses.where(complete=False)
        self.assertEqual(len(new_incomplete_analyses), self.num_incomplete)
