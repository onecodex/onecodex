from typing import Optional

from onecodex.lib.enums import FunctionalAnnotations, FunctionalAnnotationsMetric
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

    @classmethod
    def where(cls, *filters, **keyword_filters):
        from onecodex.models.collection import SampleCollection

        wrapped = super(Classifications, cls).where(*filters, **keyword_filters)
        return SampleCollection([w._resource for w in wrapped], Classifications)


class FunctionalProfiles(Analyses):
    _resource_path = "/api/v1/functional_profiles"

    def _filtered_results(
        self,
        annotation: FunctionalAnnotations,
        metric: FunctionalAnnotationsMetric,
        taxa_stratified: bool,
    ):
        return self._resource.filtered_results(
            functional_group=annotation, metric=metric, taxa_stratified=taxa_stratified
        )

    def results(self, json: bool = True):
        """Return the complete results table for a functional analysis.

        Parameters
        ----------
        json : `bool`, optional
            Return result as JSON? Default True.

        Returns
        -------
        table : `dict` or `pd.DataFrame`
            Return a JSON object with the functional analysis results or a `pd.DataFrame` if json=False.
        """
        if json is True:
            return self._results()
        else:
            return self.table()

    def table(
        self,
        annotation: Optional[FunctionalAnnotations] = None,
        taxa_stratified: bool = True,
    ):
        """Return a results table for the functional analysis.

        Parameters
        ----------
        annotation : {None, onecodex.lib.enum.FunctionalAnnotation}, optional
            If None, return a table with all annotations, otherwise filter to
            one of `onecodex.lib.enum.FunctionalAnnotation`
        taxa_stratified : bool, optional
            If False, return data only by annotation ID, ignoring taxonomic stratification

        Returns
        -------
        results_df : pd.DataFrame
            A Pandas DataFrame of the functional results.
        """
        import pandas as pd

        result_json = self._results()
        if not result_json["table"]:
            return pd.DataFrame(
                {
                    "group_name": pd.Series(dtype="str"),
                    "id": pd.Series(dtype="str"),
                    "name": pd.Series(dtype="str"),
                    "metric": pd.Series(dtype="str"),
                    "value": pd.Series(dtype="float"),
                    "taxa_stratified": pd.Series(dtype="bool"),
                    "taxon_id": pd.Series(dtype="str"),
                    "taxon_name": pd.Series(dtype="str"),
                }
            )
        results_df = pd.DataFrame(result_json["table"])

        if annotation is None:
            return (
                results_df[results_df["taxa_stratified"]]
                if taxa_stratified
                else results_df[~results_df["taxa_stratified"]]
            )
        else:
            # Validate functional annotation
            FunctionalAnnotations(annotation)

            return (
                results_df[
                    (results_df["group_name"] == annotation) & (results_df["taxa_stratified"])
                ]
                if taxa_stratified
                else results_df[
                    (results_df["group_name"] == annotation) & ~results_df["taxa_stratified"]
                ]
            )

    def filtered_table(
        self,
        annotation: FunctionalAnnotations,
        metric: FunctionalAnnotationsMetric,
        taxa_stratified: bool = True,
    ):
        """Return a results table for the functional analysis.

        Parameters
        ----------
        annotation : onecodex.lib.enum.FunctionalAnnotation, required
            Return a table for a given `onecodex.lib.enum.FunctionalAnnotation`
        metric : onecodex.lib.enum.FunctionalAnnotationMetric, required
            Return a table for a given `onecodex.lib.enum.FunctionalAnnotationMetric`
        taxa_stratified : bool, optional
            If False, return data only by annotation ID, ignoring taxonomic stratification

        Returns
        -------
        results_df : pd.DataFrame
            A Pandas DataFrame of the functional results.
        """
        import pandas as pd

        result_json = self._filtered_results(
            annotation=annotation, metric=metric, taxa_stratified=taxa_stratified
        )
        if not result_json["table"]:
            return pd.DataFrame(
                {
                    "id": pd.Series(dtype="str"),
                    "name": pd.Series(dtype="str"),
                    "value": pd.Series(dtype="float"),
                }
            )
        return pd.DataFrame(result_json["table"])


class Panels(Analyses):
    _resource_path = "/api/v1/panels"

    def results(self, json=True):
        if json is True:
            return self._results()
        else:
            raise NotImplementedError("Panel results only available as JSON at this time.")
