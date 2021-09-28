import warnings

from onecodex.lib.enums import Metric
from onecodex.models import OneCodexBase


class Analyses(OneCodexBase):
    _resource_path = "/api/v1/analyses"
    _cached_result = None

    def results(self, json=True):
        """Fetch the results of an Analyses resource.

        Parameters
        ----------
        json : bool, optional
            Return a JSON result (raw API result)? Default True.

        Returns
        -------
        Return type varies by Analyses resource sub-type. See, e.g., Classifications or Panels for
        documentation.
        """
        if json is True:
            return self._results()
        else:
            raise NotImplementedError("No non-JSON result format implemented.")

    @classmethod
    def _transform_api_results(cls, results):
        """Transform API results before caching if required."""
        return results

    def _results(self):
        try:
            if not getattr(self._resource, "_cached_result", None):
                self._resource._cached_result = self._transform_api_results(
                    self._resource.results()
                )
            return self._resource._cached_result
        except AttributeError:
            raise NotImplementedError(".results() not implemented for this Analyses resource.")


class Alignments(Analyses):
    _resource_path = "/api/v1/alignments"


class Classifications(Analyses):
    _resource_path = "/api/v1/classifications"
    _cached_table = None

    @classmethod
    def _append_abundance_rollups(cls, results):
        """Append cumulative abundances into the results returned from the API."""
        table = {t["tax_id"]: t.copy() for t in results}
        for tax_id in table:
            table[tax_id][Metric.AbundanceWChildren.value] = 0.0

        # Older results might have a bug where the parent is missing,
        # so we need to zero-out those results and renormalize the data.
        renormalize = False
        for tax_id, result in table.items():
            table[tax_id] = result
            if result[Metric.Abundance] is not None and result["parent_tax_id"] not in table:
                renormalize = True
                table[tax_id][Metric.Abundance] = None

        if renormalize:
            warnings.warn(
                "Taxa with an abundance metric but no assigned reads have been removed. In order to avoid this, re-run your samples on the latest One Codex Database."
            )
            abundance_sum = sum([t.get(Metric.Abundance, 0) or 0 for tax_id, t in table.items()])
            for tax_id, result in table.items():
                if Metric.Abundance in result and result[Metric.Abundance] is not None:
                    result[Metric.Abundance] = result[Metric.Abundance] / abundance_sum

        # Roll-up abundances to parent taxa
        for tax_id, result in table.items():
            if result["parent_tax_id"] not in table:
                continue

            parent = table[result["parent_tax_id"]]
            result[Metric.AbundanceWChildren.value] += result[Metric.Abundance] or 0.0

            while parent:
                parent[Metric.AbundanceWChildren.value] += result[Metric.Abundance] or 0.0
                parent = table.get(parent["parent_tax_id"])

        return list(table.values())

    @classmethod
    def _transform_api_results(cls, results):
        # For forward-compatibility with adding the rollups server-side
        if all(x.get(Metric.AbundanceWChildren) is not None for x in results["table"]):
            return results

        results["table"] = cls._append_abundance_rollups(results["table"])
        return results

    def results(self, json=True):
        """Return the complete results table for a classification.

        Parameters
        ----------
        json : `bool`, optional
            Return result as JSON? Default True.

        Returns
        -------
        table : `dict` or `pd.DataFrame`
            Return a JSON object with the classification results or a `pd.DataFrame` if json=False.
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
        """Return the complete results table for the classification.

        Returns
        -------
        table : `pd.DataFrame`
            A Pandas DataFrame of the classification results.
        """
        return self.results(json=False)

    def abundances(self, ids=None):
        """Query the results table to get abundance data for all or some tax ids."""
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
