from __future__ import annotations

import os.path
import time
from datetime import datetime
from functools import lru_cache
from typing import IO, TYPE_CHECKING, Any, ClassVar, List, Optional, Union

import click
import requests
from typing_extensions import Self

from onecodex.exceptions import OneCodexException
from onecodex.lib.enums import FunctionalAnnotations, FunctionalAnnotationsMetric
from onecodex.models.base import UNSET, ApiRef, OneCodexBase
from onecodex.models.filters import (
    BoolFilter,
    DatetimeFilter,
    RefFilter,
    StrFilter,
)
from onecodex.models.schemas.analysis import (
    AlignmentSchema,
    AnalysisSchema,
    ClassificationSchema,
    FunctionalRunSchema,
    MlstSchema,
    PanelSchema,
    WorkflowSchema,
)
from onecodex.models.schemas.misc import FileDetailSchema

if TYPE_CHECKING:
    import pandas as pd

    from onecodex.models.collection import SampleCollection
    from onecodex.models.misc import Jobs
    from onecodex.models.sample import Samples


def _decompress(data: bytes) -> bytes:
    """Decompress bytes, detecting format from magic bytes. Supports zstd, gzip, and plain."""
    if data[:4] == b"\x28\xb5\x2f\xfd":  # zstd magic
        import zstandard

        return zstandard.decompress(data)
    elif data[:2] == b"\x1f\x8b":  # gzip magic
        import gzip

        return gzip.decompress(data)
    return data


@lru_cache(maxsize=1)
def _get_s3_session():
    """Return a requsets session that can be used for retrieving results from S3 using presigned URLs."""

    from onecodex.utils import get_requests_session

    return get_requests_session()


def _load_results_uri(uri: str) -> dict:
    """Retrieve and load a results.json file from a presigned S3 URL or a local file."""

    import orjson

    if uri.startswith(("http://", "https://")):
        resp = _get_s3_session().get(uri, timeout=30)
        resp.raise_for_status()

        data = resp.content
    else:
        with open(uri, "rb") as f:
            data = f.read()

    return orjson.loads(_decompress(data))


def _rehydrate_functional_results(
    condensed_results: dict,
    annotation_filter: Optional[str] = None,
    metric_filter: Optional[str] = None,
    taxa_stratified_filter: Optional[bool] = None,
) -> dict:
    """Rehydrate condensed functional results into the public API results format.

    the annotation, metric, and taxa stratified filters can be used to control what gets
    rehydrated. If any of these are set, all other annotations/metrics/stratifications will
    be removed prior to returning the rehydrated results.
    """

    # maps tax IDs -> names
    taxa_map = {node["id"]: node.get("name") for node in condensed_results["taxonomy"]["nodes"]}
    taxa_map["0"] = "unclassified"

    table = []
    pathway_table = []

    def add_row(
        group_name: str,
        feature_id: str,
        feature_name: Optional[str],
        metric: str,
        value: float,
        taxa_stratified: bool = False,
        taxon_id: Optional[str] = None,
        destination: Optional[list] = None,
    ) -> None:
        """Add a row formatted for results['table']."""
        # no need to rehydrate rows we don't need
        if (
            (annotation_filter is not None and group_name != annotation_filter)
            or (metric_filter is not None and metric != metric_filter)
            or (taxa_stratified_filter is not None and taxa_stratified != taxa_stratified_filter)
        ):
            return

        destination = table if destination is None else destination
        destination.append(
            {
                "group_name": group_name,
                "id": feature_id,
                "name": feature_name,
                "metric": metric,
                "value": value,
                "taxa_stratified": taxa_stratified,
                "taxon_id": taxon_id,
                "taxon_name": taxa_map.get(taxon_id) if taxa_stratified else None,
            }
        )

    def add_species_rows(
        group_name: str,
        feature_id: str,
        feature_name: Optional[str],
        taxon_id: str,
        metric_values: tuple[tuple[str, float], ...],
        destination: Optional[list] = None,
    ) -> None:
        """Add a species-level (taxa stratified) row formatted for results['table'].

        For all functional groups, except pathways, a taxa stratified row will have one of
        two metrics: cpm or rpk. For pathways, these metrics are abundance (synonymous with RPK)
        and coverage.
        """
        for metric, value in metric_values:
            add_row(
                group_name,
                feature_id,
                feature_name,
                metric,
                value,
                taxa_stratified=True,
                taxon_id=taxon_id,
                destination=destination,
            )

    # standard functional groups use the following format:
    # [id, name, total_cpm, total_rpk, [[taxid, cpm, rpk], ...]]
    for group_name, features in condensed_results["results"].items():
        # we'll do pathways separately because we neeed to split pathways into
        # the metacyc functional group, has cpm/rpk metrics, and the pathways
        # group, which has abundance/coverage metrics
        if group_name == "pathways":
            continue

        for feature_id, encoded_name, community_cpm, total_rpk, contributions in features:
            # Missing metadata is encoded by repeating the feature ID.
            feature_name = None if encoded_name == feature_id else encoded_name

            if taxa_stratified_filter is not True:
                add_row(group_name, feature_id, feature_name, "cpm", community_cpm)
                add_row(group_name, feature_id, feature_name, "rpk", total_rpk)

            if taxa_stratified_filter is not False:
                for taxon_id, cpm, rpk in contributions:
                    add_species_rows(
                        group_name,
                        feature_id,
                        feature_name,
                        taxon_id,
                        (("rpk", rpk), ("cpm", cpm)),
                    )

    pathways = condensed_results["results"].get("pathways", [])

    # metacyc values are folded into each condensed pathway:
    # [id, name, community_abundance, community_coverage, metacyc_cpm,
    #  [[taxid, species_abundance, species_coverage, metacyc_cpm], ...]]
    for (
        pathway_id,
        pathway_name,
        community_abundance,
        community_coverage,
        community_cpm,
        contributions,
    ) in pathways:
        # these are only added when not stratifying by taxa
        if taxa_stratified_filter is not True:
            add_row("metacyc", pathway_id, pathway_name, "cpm", community_cpm)
            add_row("metacyc", pathway_id, pathway_name, "rpk", community_abundance)

            add_row(
                "pathways",
                pathway_id,
                pathway_name,
                "coverage",
                community_coverage,
                destination=pathway_table,
            )
            add_row(
                "pathways",
                pathway_id,
                pathway_name,
                "abundance",
                community_abundance,
                destination=pathway_table,
            )

        if taxa_stratified_filter is not False:
            for taxon_id, abundance, coverage, cpm in contributions:
                # metacyc ids/names are synonymous with pathway ids/names
                add_species_rows(
                    "metacyc",
                    pathway_id,
                    pathway_name,
                    taxon_id,
                    (("rpk", abundance), ("cpm", cpm)),
                )
                # add these to pathways separately so they aren't interleaved with metacyc
                add_species_rows(
                    "pathways",
                    pathway_id,
                    pathway_name,
                    taxon_id,
                    (("coverage", coverage), ("abundance", abundance)),
                    destination=pathway_table,
                )

    table.extend(pathway_table)

    return {
        "table": table,
        "n_reads": condensed_results["n_reads"],
        "n_mapped": condensed_results["n_mapped"],
    }


def _rehydrate_filtered_functional_results(
    condensed_results: dict,
    annotation: FunctionalAnnotations,
    metric: FunctionalAnnotationsMetric,
    taxa_stratified: bool,
) -> dict:
    """Rehydrate condensed functional results into the public filtered API results format."""

    _SKIP_FUNCTIONAL_IDS = {"UNMAPPED", "UNGROUPED", "UNINTEGRATED"}

    annotation = FunctionalAnnotations.from_value(annotation)
    metric = FunctionalAnnotationsMetric.from_value(metric)

    if metric not in FunctionalAnnotationsMetric.metrics_for_annotation(annotation):
        raise OneCodexException(
            f"metric {metric} cannot be retrieved for functional group {annotation}"
        )

    effective_metric = metric

    if annotation in (FunctionalAnnotations.Pathways, FunctionalAnnotations.MetaCyc):
        selected = [
            pathway
            for pathway in condensed_results["results"].get("pathways", [])
            if pathway[0] not in _SKIP_FUNCTIONAL_IDS
        ]

        # complete abundance means annotations where the community level coverage is 1.0
        if (
            annotation == FunctionalAnnotations.Pathways
            and metric == FunctionalAnnotationsMetric.CompleteAbundance
        ):
            # pathway[3] == community_coverage
            selected = [pathway for pathway in selected if pathway[3] == 1.0]
            effective_metric = FunctionalAnnotationsMetric.Abundance

        selected_results = {"pathways": selected}
    else:
        selected_results = {
            annotation.value: [
                feature
                for feature in condensed_results["results"].get(annotation.value, [])
                if feature[0] not in _SKIP_FUNCTIONAL_IDS
            ]
        }

    # shouldn't mutate the original otherwise we change what @lru_cache._condensed_results()
    # stores
    filtered_condensed_results = {
        **condensed_results,
        "results": selected_results,
    }

    hydrated = _rehydrate_functional_results(
        filtered_condensed_results,
        annotation_filter=annotation.value,
        metric_filter=effective_metric.value,
        taxa_stratified_filter=taxa_stratified,
    )

    table = []

    for row in hydrated["table"]:
        filtered_row = {
            "id": row["id"],
            "name": row["name"],
            "value": row["value"],
        }

        if taxa_stratified:
            filtered_row["taxon_id"] = row["taxon_id"]
            filtered_row["taxon_name"] = row["taxon_name"]

        table.append(filtered_row)

    return {
        "table": table,
        "n_reads": condensed_results["n_reads"],
        "n_mapped": condensed_results["n_mapped"],
    }


class _AnalysesBase(OneCodexBase):
    _allowed_methods = {
        "instances_public": None,
    }
    sample: Union["Samples", ApiRef]  # noqa: F821
    job: Union["Jobs", ApiRef]  # noqa: F821
    error_msg: Optional[str] = None

    def __hash__(self) -> int:
        return hash(self.id)

    def __eq__(self, other: object) -> bool:
        # Required alongside __hash__: lru_cache uses __eq__ to check for cache hits, and
        # Pydantic's default __eq__ compares __dict__ which recurses infinitely on nested models.
        if not isinstance(other, _AnalysesBase):
            return NotImplemented
        return self.id == other.id

    @classmethod
    def where(
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        created_at: datetime | DatetimeFilter = UNSET,
        updated_at: datetime | DatetimeFilter = UNSET,
        complete: bool | BoolFilter = UNSET,
        draft: bool | BoolFilter = UNSET,
        success: bool | BoolFilter = UNSET,
        error_msg: str | StrFilter | None = UNSET,
        analysis_type: str | StrFilter = UNSET,
        job: Jobs | str | RefFilter = UNSET,
        sample: Samples | str | RefFilter = UNSET,
        **keyword_filters: Any,
    ) -> list[Self]:
        """Query analyses.

        An analysis is the result of running a job on a sample. Use this
        when you want results across analysis types (classifications,
        panels, functional profiles, etc.); for a single type, prefer the
        dedicated subclass (e.g. :class:`Classifications`).

        Examples
        --------
        Find completed analyses since yesterday::

            from datetime import datetime, timedelta, timezone
            since = (datetime.now(timezone.utc) - timedelta(days=1)).isoformat()
            ocx.Analyses.where(complete=True, created_at={"$gte": since})

        Find failed analyses (terminal but unsuccessful)::

            ocx.Analyses.where(complete=True, success=False)

        Find all analyses for a given sample::

            ocx.Analyses.where(sample=sample)

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        return super().where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            created_at=created_at,
            updated_at=updated_at,
            complete=complete,
            draft=draft,
            success=success,
            error_msg=error_msg,
            analysis_type=analysis_type,
            job=job,
            sample=sample,
            **keyword_filters,
        )

    def refresh(self) -> None:
        """Fetch the current state from the API and update this object's state fields in-place."""
        resp = self._client.get(f"{self._api._base_url}{self._resource_path}/{self.id}?expand=all")
        resp.raise_for_status()
        new_obj = self.__class__.model_validate(resp.json())
        for field, value in new_obj:
            setattr(self, field, value)
        self._snapshot = self._object_snapshot()

    MIN_POLL_INTERVAL_SECONDS: ClassVar[int] = 5

    def await_completion(
        self,
        timeout: Optional[float] = None,
        initial_interval: int = 5,
        max_interval: int = 120,
        backoff: float = 1.5,
    ) -> "_AnalysesBase":
        """Poll the API until the analysis reaches a terminal state.

        An analysis is considered terminal when ``complete`` is True. The method refreshes
        ``self`` in place and returns it.

        Parameters
        ----------
        timeout : float, optional
            Maximum number of seconds to wait. ``None`` (default) waits indefinitely.
        initial_interval : int, optional
            Seconds to wait between the first polls. Must be at least
            ``MIN_POLL_INTERVAL_SECONDS`` (5). Defaults to 5.
        max_interval : int, optional
            Upper bound on the polling interval as it backs off. Must be at least
            ``MIN_POLL_INTERVAL_SECONDS`` (5). Defaults to 120 seconds.
        backoff : float, optional
            Multiplier applied to the interval after each poll. Defaults to 1.5.

        Returns
        -------
        self : the analysis, refreshed to its terminal state.

        Raises
        ------
        TimeoutError
            If ``timeout`` elapses before the analysis completes.
        """
        if (
            initial_interval < self.MIN_POLL_INTERVAL_SECONDS
            or max_interval < self.MIN_POLL_INTERVAL_SECONDS
        ):
            raise ValueError(
                f"Polling intervals must be >= {self.MIN_POLL_INTERVAL_SECONDS} seconds."
            )
        if backoff < 1:
            raise ValueError("backoff must be >= 1.")

        deadline = None if timeout is None else time.monotonic() + timeout
        interval = initial_interval

        self.refresh()
        while not self.complete:
            if deadline is not None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(
                        f"Analysis {self.id} did not complete within {timeout} seconds."
                    )
                sleep_for = min(interval, remaining)
            else:
                sleep_for = interval
            time.sleep(sleep_for)
            interval = min(interval * backoff, max_interval)
            self.refresh()

        return self

    def logs(self, tail: Optional[int] = None) -> str:
        """Fetch the job run logs for this analysis.

        Parameters
        ----------
        tail : int, optional
            If set, return only the last ``tail`` lines of the logs. Must be >= 1.

        Returns
        -------
        str
            The job run logs as a plain-text string. Empty if the analysis has not
            produced any logs.
        """
        params = {}
        if tail is not None:
            if tail < 1:
                raise ValueError("tail must be >= 1.")
            params["tail"] = tail
        resp = self._client.get(f"{self._api._base_url}{self.field_uri}/logs", params=params)
        try:
            resp.raise_for_status()
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                raise OneCodexException(f"Logs not found for analysis {self.id}.")
            raise
        return resp.text

    def results(self, json: bool = True):
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

    @lru_cache
    def _results(self):
        resp = self._client.get(f"{self._api._base_url}{self.field_uri}/results")
        return resp.json()

    @lru_cache
    def get_files(self) -> List[FileDetailSchema]:
        """Fetch the files details of an Analyses.

        Returns
        -------
        A list of FileDetailSchema
        """
        resp = self._client.get(f"{self._api._base_url}{self.field_uri}/file_details")
        return [FileDetailSchema(**x) for x in resp.json()["files"]]

    # It is almost copy/paste of ResourceDownloadMixin._download
    # I do not want to extract re-usable function just yet I need more than two copies.
    def download_file(
        self,
        filepath: Union[str, FileDetailSchema],
        out_path: Optional[str] = None,
        out_file_obj: Optional[IO] = None,
        progressbar: bool = False,
    ) -> str:
        """Download analysis result file.

        Parameters
        ----------
        filepath: `str` or `FileDetailSchema`
            Must be one of objects or filepathes returned by `get_files`
        out_path: `string`, optional
            Full path to save the file to. If omitted, defaults to the original filename
            in the current working directory.
        out_file_obj: file-like object, optional
            Rather than save the file to a `out_path`, write it to this file-like object.
        progressbar: `bool`, optional
            Display a progress bar using Click for the download?

        Returns
        -------
        `string`
            The path the file was downloaded to, if applicable. Otherwise, None.

        Notes
        -----
        Existing paths will not be overwritten.
        """
        from requests.adapters import HTTPAdapter
        from requests.packages.urllib3.util.retry import Retry

        # You can pass filepath instead of `FileDetailSchema` object.
        if isinstance(filepath, str):
            exists = False
            for x in self.get_files():
                if x.filepath == filepath:
                    filepath = x
                    exists = True
                    break
            if not exists:
                raise OneCodexException(f"Filepath: {filepath} does not exist.")

        if out_path and out_file_obj:
            raise OneCodexException("Please specify only one of: path, file_obj")

        try:
            if out_path is None and out_file_obj is None:
                out_path = os.path.join(os.getcwd(), filepath.filename)

            if out_path and os.path.exists(out_path):
                raise OneCodexException(f"{out_path} already exists. Will not overwrite.")

            session = requests.Session()

            # Retry up to 5 times with backoff timing of 2s, 4s, 8s, 16s, and 32s (applies to all
            # HTTP methods). 404 is included for cases where the file is being asynchronously
            # uploaded to S3 and is expected to be available soon.
            retry_strategy = Retry(
                total=5,
                backoff_factor=2,
                status_forcelist=[404, 429, 500, 502, 503, 504],
                allowed_methods=None,
            )
            adapter = HTTPAdapter(max_retries=retry_strategy)

            session.mount("http://", adapter)
            session.mount("https://", adapter)
            resp = session.get(filepath.url, stream=True)

            if out_path:
                f_out = open(out_path, "wb")
            else:
                f_out = out_file_obj

            if progressbar:
                progress_label = os.path.basename(out_path) if out_path else filepath.filename
                with click.progressbar(length=filepath.size, label=progress_label) as bar:
                    for data in resp.iter_content(chunk_size=1024):
                        bar.update(len(data))
                        f_out.write(data)
            else:
                for data in resp.iter_content(chunk_size=1024):
                    f_out.write(data)

            # do not close the handle if file_obj is used
            if not out_file_obj:
                f_out.close()

        except KeyboardInterrupt:
            if out_path and os.path.exists(out_path):
                os.remove(out_path)
            raise
        except requests.exceptions.HTTPError as exc:
            if exc.response.status_code == 401:
                raise OneCodexException("You must be logged in to download files.")
            elif exc.response.status_code == 402:
                raise OneCodexException(
                    "You must either have a premium platform account or be in "
                    "a notebook environment to download files. Please feel free to contact us "
                    "about your subscription at support@onecodex.com."
                )
            elif exc.response.status_code == 403:
                raise OneCodexException("You are not authorized to download this file.")
            else:
                raise OneCodexException(
                    "Download failed with an HTTP status code {}.".format(exc.response.status_code)
                )

        return out_path


class Analyses(_AnalysesBase, AnalysisSchema):
    _resource_path = "/api/v1/analyses"


class Alignments(_AnalysesBase, AlignmentSchema):
    _resource_path = "/api/v1/alignments"


class Classifications(_AnalysesBase, ClassificationSchema):
    _resource_path = "/api/v1/classifications"
    # root & cellular organisms
    _NONSPECIFIC_TAX_IDS = {"1", "131567"}

    @lru_cache
    def _results(self):
        # results_uri is a pre-signed URL included in the API response; it's None when the
        # analysis isn't available. Fetch directly from the URL when present to avoid the round-trip
        # through the API server. Local paths are supported for testing.
        if self.results_uri is not None:
            try:
                return _load_results_uri(self.results_uri)
            except requests.HTTPError:
                # Pre-signed URL may have expired; fall through to the API.
                pass

        # Fall back to the /results API endpoint (e.g. analysis not yet available, older API,
        # or expired pre-signed URL).
        resp = self._client.get(f"{self._api._base_url}{self.field_uri}/results")
        return resp.json()

    def results(self, json: bool = True) -> dict | pd.DataFrame:
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

    def _table(self) -> pd.DataFrame:
        import pandas as pd

        return pd.DataFrame(self._results()["table"])

    def table(self) -> pd.DataFrame:
        """Return the complete results table for the classification.

        Returns
        -------
        table : `pd.DataFrame`
            A Pandas DataFrame of the classification results.
        """
        return self.results(json=False)

    @property
    def _classification_stats(self) -> dict[str, int]:
        """Return some high-level classification stats needed for normalization."""
        results = self.results()
        host_tax_ids = set(results["host_tax_ids"])

        n_host_reads = 0
        n_nonspecific_reads = 0
        n_mapped_reads = 0
        for row in results["table"]:
            n_mapped_reads += row["readcount"]
            if row["tax_id"] in host_tax_ids:
                n_host_reads += row["readcount"]
            if row["tax_id"] in self._NONSPECIFIC_TAX_IDS:
                n_nonspecific_reads += row["readcount"]

        return {
            "n_reads_total": results["n_reads"],
            "n_host_reads": n_host_reads,
            "n_nonspecific_reads": n_nonspecific_reads,
            "n_mapped_reads": n_mapped_reads,
            "n_mapped_microbial_reads": n_mapped_reads - n_host_reads - n_nonspecific_reads,
        }

    @classmethod
    def where(  # type: ignore[override]
        cls,
        *filters: str | dict,
        sort: str | list[str] | None = None,
        limit: int | None = None,
        public: bool = False,
        filter: Any = None,
        created_at: datetime | DatetimeFilter = UNSET,
        updated_at: datetime | DatetimeFilter = UNSET,
        complete: bool | BoolFilter = UNSET,
        draft: bool | BoolFilter = UNSET,
        success: bool | BoolFilter = UNSET,
        error_msg: str | StrFilter | None = UNSET,
        job: Jobs | str | RefFilter = UNSET,
        sample: Samples | str | RefFilter = UNSET,
    ) -> "SampleCollection":
        """Query classifications and return a :class:`SampleCollection`.

        Classifications are taxonomic results — typically the One Codex
        Database run against each sample. Filters mirror those on
        :meth:`Analyses.where`.

        Examples
        --------
        Find recent successful classifications::

            from datetime import datetime, timedelta, timezone
            since = (datetime.now(timezone.utc) - timedelta(days=1)).isoformat()
            ocx.Classifications.where(
                complete=True,
                success=True,
                created_at={"$gte": since},
            )

        Find the classification for a specific sample::

            cls_run = ocx.Classifications.where(sample=sample)[0]

        See :meth:`OneCodexBase.where` for the full operator reference.
        """
        from onecodex.models.collection import SampleCollection

        classifications = super(Classifications, cls).where(
            *filters,
            sort=sort,
            limit=limit,
            public=public,
            filter=filter,
            created_at=created_at,
            updated_at=updated_at,
            complete=complete,
            draft=draft,
            success=success,
            error_msg=error_msg,
            job=job,
            sample=sample,
        )
        return SampleCollection(classifications, Classifications)


class FunctionalProfiles(_AnalysesBase, FunctionalRunSchema):
    _resource_path = "/api/v1/functional_profiles"

    _FUNCTIONAL_RESULTS_VERSION = 1

    def _condensed_results(self) -> Optional[dict]:
        if self.results_uri is None:
            return None

        try:
            results = _load_results_uri(self.results_uri)
        except requests.HTTPError:
            # presigned URL may have expired
            return None

        # any updates to the condensed results may impact the positional arrays,
        # e.g., [id, name, cpm, rpk] and we should disallow version mismatches to
        # prevent possible issues when rehydrating results
        if results.get("version") != self._FUNCTIONAL_RESULTS_VERSION:
            return None

        return results

    def _results(self) -> dict:
        condensed_results = self._condensed_results()

        if condensed_results is not None:
            return _rehydrate_functional_results(condensed_results)

        return super()._results()

    def _filtered_results(
        self,
        annotation: FunctionalAnnotations,
        metric: FunctionalAnnotationsMetric,
        taxa_stratified: bool,
    ):
        condensed_results = self._condensed_results()

        if condensed_results is not None:
            return _rehydrate_filtered_functional_results(
                condensed_results,
                annotation=annotation,
                metric=metric,
                taxa_stratified=taxa_stratified,
            )

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
            columns = {
                "id": pd.Series(dtype="str"),
                "name": pd.Series(dtype="str"),
                "value": pd.Series(dtype="float"),
            }
            if taxa_stratified:
                # ensure taxon_id and taxon_name are present for downstream processing
                columns["taxon_id"] = pd.Series(dtype="str")
                columns["taxon_name"] = pd.Series(dtype="str")
            return pd.DataFrame(columns)
        return pd.DataFrame(result_json["table"])


class Panels(_AnalysesBase, PanelSchema):
    _resource_path = "/api/v1/panels"

    def results(self, json=True):
        if json is True:
            return self._results()
        else:
            raise NotImplementedError("Panel results only available as JSON at this time.")


class Mlsts(_AnalysesBase, MlstSchema):
    _resource_path = "/api/v1/mlsts"

    def results(self, json=True):
        if json is True:
            return self._results()
        else:
            raise NotImplementedError("MLST results only available as JSON at this time.")


class Workflows(_AnalysesBase, WorkflowSchema):
    _resource_path = "/api/v1/workflows"

    def results(self, json=True):
        if json is True:
            return self._results()
        else:
            raise NotImplementedError("Workflow results only available as JSON at this time.")
