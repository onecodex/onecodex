from typing import Optional, Union, Dict

from onecodex.exceptions import OneCodexException
from onecodex.models.base import OneCodexBase, ApiRef
from onecodex.lib.enums import FunctionalAnnotations, FunctionalAnnotationsMetric

from onecodex.models.schemas.analysis import AnalysisSchema
from onecodex.models.schemas.analysis import AlignmentSchema
from onecodex.models.schemas.analysis import ClassificationSchema
from onecodex.models.schemas.analysis import FunctionalRunSchema
from onecodex.models.schemas.analysis import PanelSchema


class _AnalysesBase(OneCodexBase):
    _allowed_methods = {
        "instances_public": None,
    }
    _cached_result: Dict = {}
    sample: Union["Samples", ApiRef]  # noqa: F821
    job: Union["Jobs", ApiRef]  # noqa: F821
    error_msg: Optional[str] = None

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
        if getattr(self, "_cached_result", None):
            return self._cached_result

        resp = self._client.get(f"{self._api._base_url}{self.field_uri}/results")
        self._cached_result = resp.json()
        return self._cached_result


class Analyses(_AnalysesBase, AnalysisSchema):
    _resource_path = "/api/v1/analyses"
    _cached_result = None


class Alignments(_AnalysesBase, AlignmentSchema):
    _resource_path = "/api/v1/alignments"


class Classifications(_AnalysesBase, ClassificationSchema):
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
        resp = self._client.get(f"{self._api._base_url}{self.field_uri}/readlevel")
        return resp.json()

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

        classifications = super(Classifications, cls).where(*filters, **keyword_filters)
        return SampleCollection(classifications, Classifications)


class FunctionalProfiles(_AnalysesBase, FunctionalRunSchema):
    _resource_path = "/api/v1/functional_profiles"

    def _filtered_results(
        self,
        annotation: FunctionalAnnotations,
        metric: FunctionalAnnotationsMetric,
        taxa_stratified: bool,
    ):
        resp = self._client.get(
            f"{self._api._base_url}{self.field_uri}/filtered_results",
            params={
                "functional_group": annotation,
                "metric": metric,
                "taxa_stratified": taxa_stratified,
            },
        )
        if resp.status_code != 200:
            raise OneCodexException(resp.json()["message"])
        return resp.json()

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


class Panels(_AnalysesBase, PanelSchema):
    _resource_path = "/api/v1/panels"

    def results(self, json=True):
        if json is True:
            return self._results()
        else:
            raise NotImplementedError("Panel results only available as JSON at this time.")
