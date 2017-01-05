import datetime
from collections import defaultdict, OrderedDict

from onecodex.models import OneCodexBase


class Analyses(OneCodexBase):
    _resource_path = '/api/v1/analyses'
    _cached_result = None

    def results(self, json=True):
        """
        Fetch the results for an Analyses resource.

        Parameters
        ----------
        json : bool, optional
            Return a JSON result (raw API result)? Default True.

        Returns
        -------
        Return type varies by Analyses resource sub-type. See, e.g.,
        Classifications or Panels for documentation.
        """
        if json is True:
            return self._results()
        else:
            raise NotImplementedError('No non-JSON result format implemented.')

    def _results(self):
        try:
            if not self._cached_result:
                self._cached_result = self._resource.results()
            return self._cached_result
        except AttributeError:
            raise NotImplementedError('.results() not implemented for this Analyses resource.')


class Alignments(Analyses):
    _resource_path = '/api/v1/alignments'


class Classifications(Analyses):
    _resource_path = '/api/v1/classifications'
    _cached_table = None

    @staticmethod
    def to_otu(classifications):
        """
        Converts a list of classifications into a dictionary resembling
        an OTU table.

        Parameters
        ----------
        classifications : list of Classifications
            List of Classifications, i.e., ocx.Classifications.where(job=job).
            If 1 Classification is passed, it will be coerced into a list.

        Returns
        -------
        otu_table : OrderedDcit
            A BIOM OTU table, returned as a Python OrderedDict (can be dumped to JSON)
        """
        otu_format = 'Biological Observation Matrix 0.9.1-dev'
        otu_url = 'http://biom-format.org/documentation/format_versions/biom-1.0.html'  # noqa

        otu = OrderedDict({'format': otu_format,
                           'format_url': otu_url,
                           'type': 'OTU table',
                           'generated_by': 'One Codex API V1',
                           'date': datetime.datetime.now().isoformat(),
                           'rows': [],
                           'columns': [],
                           'matrix_type': 'sparse',
                           'matrix_element_type': 'int'})

        rows = defaultdict(dict)
        if not isinstance(classifications, list):
            classifications = [classifications]

        for classification in classifications:
            col_id = len(otu['columns'])  # 0 index
            columns_entry = {'id': str(classification.id)}
            otu['columns'].append(columns_entry)
            sample_df = classification.table()

            for tax_id in sample_df['tax_id']:
                rows[tax_id][col_id] = int(
                    sample_df[sample_df['tax_id'] == tax_id]['readcount'])

        num_rows = len(rows)
        num_cols = len(otu['columns'])

        otu['shape'] = [num_rows, num_cols]
        otu['data'] = []

        for present_taxa in sorted(rows):
            # add the row entry
            row_id = len(otu['rows'])
            otu['rows'].append({'id': present_taxa})

            for sample_with_hit in rows[present_taxa]:
                counts = rows[present_taxa][sample_with_hit]
                otu['data'].append([row_id, sample_with_hit, counts])

        return otu

    def results(self, json=True):
        """
        Returns the complete results table for the classification.

        Parameters
        ----------
        json : bool, optional
            Return result as JSON? Default True.

        Returns
        -------
        table : dict | DataFrame
            Return a JSON object with the classification results or a Pandas DataFrame
            if json=False.
        """
        if json is True:
            return self._results()
        else:
            return self._table()

    def readlevel(self):
        return self._resource.readlevel()

    def _table(self):
        import pandas as pd
        if self._cached_table is None:
            self._cached_table = pd.DataFrame(self._results()['table'])
        return self._cached_table

    def table(self):
        """
        Returns the complete results table for the classification.

        Returns
        -------
        table : DataFrame
            A Pandas DataFrame of the classification results.
        """
        return self.results(json=False)

    def abundances(self, ids=None):
        """
        Query the results table to get abundance data for all or some tax ids
        """
        # TODO: Consider removing this method... since it's kind of trivial
        #       May want to replace with something that actually gets genome-size adjusted
        #       abundances from the results table
        if ids is None:
            # get the data frame
            return self.table()

        else:
            res = self.table()
            return res[res['tax_id'].isin(ids)]


class Panels(Analyses):
    _resource_path = '/api/v1/panels'

    def results(self, json=True):
        if json is True:
            return self._results()
        else:
            raise NotImplementedError('Panel results only available as JSON at this time.')
