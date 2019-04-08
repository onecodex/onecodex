from onecodex.models import OneCodexBase


class Analyses(OneCodexBase):
    _resource_path = "/api/v1/analyses"
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
            raise NotImplementedError("No non-JSON result format implemented.")

    def _results(self):
        try:
            if not getattr(self._resource, "_cached_result", None):
                self._resource._cached_result = self._resource.results()
            return self._resource._cached_result
        except AttributeError:
            raise NotImplementedError(".results() not implemented for this Analyses resource.")


class Alignments(Analyses):
    _resource_path = "/api/v1/alignments"


class Classifications(Analyses):
    _resource_path = "/api/v1/classifications"
    _cached_table = None

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

    def _readlevel(self):
        return self._resource.readlevel()

    def _table(self):
        import pandas as pd

        if self._cached_table is None:
            self._cached_table = pd.DataFrame(self._results()["table"])
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
            return res[res["tax_id"].isin(ids)]

    @classmethod
    def where(cls, *filters, **keyword_filters):
        from onecodex.models.collection import SampleCollection

        wrapped = super(Classifications, cls).where(*filters, **keyword_filters)
        return SampleCollection([w._resource for w in wrapped], Classifications)


class Panels(Analyses):
    _resource_path = "/api/v1/panels"

    def results(self, json=True):
        if json is True:
            return self._results()
        else:
            raise NotImplementedError("Panel results only available as JSON at this time.")
