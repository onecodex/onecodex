"""
Analysis model classes extending the generated models.
"""

from typing import Any, Dict, List, Optional, Union

from onecodex.lib.enums import FunctionalAnnotations, FunctionalAnnotationsMetric
from onecodex.models.base import OneCodexModel


class Analyses(OneCodexModel):
    """Base analysis model."""
    
    _resource_path = "/api/v1/analyses"
    _cached_result = None
    
    # Core analysis fields
    complete: Optional[bool] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None
    sample: Optional["Samples"] = None
    
    def results(self, json: bool = True) -> Union[Dict[str, Any], Any]:
        """Fetch the results of an analysis.
        
        Args:
            json: Return JSON result (raw API result)
            
        Returns:
            Analysis results in requested format
        """
        if json:
            return self._results()
        else:
            raise NotImplementedError("No non-JSON result format implemented.")
    
    def _results(self) -> Dict[str, Any]:
        """Fetch raw results from API."""
        if not hasattr(self, '_cached_result') or self._cached_result is None:
            # Make API call to fetch results
            self._check_bound()
            endpoint = f"{self._resource_path}/{self.id}/results"
            self._cached_result = self._client.http_client.get(endpoint)
        return self._cached_result


class Classifications(Analyses):
    """Classification analysis model."""
    
    _resource_path = "/api/v1/classifications"
    _cached_table = None
    
    def results(self, json: bool = True) -> Union[Dict[str, Any], "pd.DataFrame"]:
        """Return classification results.
        
        Args:
            json: Return as JSON (True) or DataFrame (False)
            
        Returns:
            Classification results
        """
        if json:
            return self._results()
        else:
            return self._table()
    
    def _table(self) -> "pd.DataFrame":
        """Convert results to pandas DataFrame."""
        import pandas as pd
        
        if self._cached_table is None:
            results = self._results()
            self._cached_table = pd.DataFrame(results.get("table", []))
        return self._cached_table
    
    def table(self) -> "pd.DataFrame":
        """Return results as pandas DataFrame."""
        return self.results(json=False)
    
    def _readlevel(self) -> Dict[str, Any]:
        """Get read-level results."""
        self._check_bound()
        endpoint = f"{self._resource_path}/{self.id}/readlevel"
        return self._client.http_client.get(endpoint)
    
    @classmethod
    def where(cls, *args, **kwargs):
        """Override where to return SampleCollection."""
        from onecodex.models.collection import SampleCollection
        
        models = super().where(*args, **kwargs)
        return SampleCollection([m for m in models], cls)


class Alignments(Analyses):
    """Alignment analysis model."""
    
    _resource_path = "/api/v1/alignments"


class FunctionalProfiles(Analyses):
    """Functional profile analysis model."""
    
    _resource_path = "/api/v1/functional_profiles"
    
    def _filtered_results(
        self,
        annotation: FunctionalAnnotations,
        metric: FunctionalAnnotationsMetric,
        taxa_stratified: bool,
    ) -> Dict[str, Any]:
        """Get filtered functional results."""
        self._check_bound()
        endpoint = f"{self._resource_path}/{self.id}/filtered_results"
        params = {
            "functional_group": annotation,
            "metric": metric,
            "taxa_stratified": taxa_stratified
        }
        return self._client.http_client.get(endpoint, params)
    
    def results(self, json: bool = True) -> Union[Dict[str, Any], "pd.DataFrame"]:
        """Return functional analysis results.
        
        Args:
            json: Return as JSON (True) or DataFrame (False)
            
        Returns:
            Functional analysis results
        """
        if json:
            return self._results()
        else:
            return self.table()
    
    def table(
        self,
        annotation: Optional[FunctionalAnnotations] = None,
        taxa_stratified: bool = True,
    ) -> "pd.DataFrame":
        """Return functional results as DataFrame.
        
        Args:
            annotation: Filter to specific annotation
            taxa_stratified: Include taxonomic stratification
            
        Returns:
            Functional results DataFrame
        """
        import pandas as pd
        
        result_json = self._results()
        if not result_json.get("table"):
            # Return empty DataFrame with expected columns
            return pd.DataFrame({
                "group_name": pd.Series(dtype="str"),
                "id": pd.Series(dtype="str"),
                "name": pd.Series(dtype="str"),
                "metric": pd.Series(dtype="str"),
                "value": pd.Series(dtype="float"),
                "taxa_stratified": pd.Series(dtype="bool"),
                "taxon_id": pd.Series(dtype="str"),
                "taxon_name": pd.Series(dtype="str"),
            })
        
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
    ) -> "pd.DataFrame":
        """Return filtered functional results as DataFrame.
        
        Args:
            annotation: Functional annotation type
            metric: Functional annotation metric
            taxa_stratified: Include taxonomic stratification
            
        Returns:
            Filtered functional results DataFrame
        """
        import pandas as pd
        
        result_json = self._filtered_results(
            annotation=annotation, metric=metric, taxa_stratified=taxa_stratified
        )
        
        if not result_json.get("table"):
            return pd.DataFrame({
                "id": pd.Series(dtype="str"),
                "name": pd.Series(dtype="str"),
                "value": pd.Series(dtype="float"),
            })
        
        return pd.DataFrame(result_json["table"])


class Panels(Analyses):
    """Panel analysis model."""
    
    _resource_path = "/api/v1/panels"
    
    def results(self, json: bool = True) -> Dict[str, Any]:
        """Return panel results (JSON only)."""
        if json:
            return self._results()
        else:
            raise NotImplementedError("Panel results only available as JSON at this time.")


# Forward reference resolution for type hints
from onecodex.models.samples import Samples
Classifications.model_rebuild()
Alignments.model_rebuild()
FunctionalProfiles.model_rebuild()
Panels.model_rebuild()