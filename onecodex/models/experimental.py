from onecodex.models import OneCodexBase, ResourceList
from onecodex.models.helpers import ResourceDownloadMixin
from onecodex.models.analysis import Analyses
from onecodex.lib.enums import FunctionalAnnotations


class AnnotationSets(OneCodexBase, ResourceDownloadMixin):
    _resource_path = "/api/v1_experimental/annotation_sets"

    def download(self, path=None, file_obj=None, progressbar=False):
        """Download an AnnotationSet in GenBank format.

        Parameters
        ----------
        path : `string`, optional
            Full path to save the file to. If omitted, defaults to the original filename
            in the current working directory.
        file_obj : file-like object, optional
            Rather than save the file to a path, write it to this file-like object.
        progressbar : `bool`
            Display a progress bar using Click for the download?

        Returns
        -------
        `string`
            The path the file was downloaded to, if applicable. Otherwise, None.

        Notes
        -----
        If no arguments specified, defaults to download the file as the original filename
        in the current working directory. If `file_obj` given, will write data into the
        passed file-like object. If `path` given, will download the file to the path provided,
        but will not overwrite any existing files.
        """
        return self._download(
            "download_uri",
            "annotation_set_" + self.id + ".gbk.gz",
            use_potion_session=False,
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )

    def download_csv(self, path=None, file_obj=None, progressbar=False):
        """Download an AnnotationSet in CSV format.

        Includes Annotation coordinates and sequences in both amino acid and nucleotide space.

        Parameters
        ----------
        path : `string`, optional
            Full path to save the file to. If omitted, defaults to the original filename
            in the current working directory.
        file_obj : file-like object, optional
            Rather than save the file to a path, write it to this file-like object.
        progressbar : `bool`
            Display a progress bar using Click for the download?

        Returns
        -------
        `string`
            The path the file was downloaded to, if applicable. Otherwise, None.

        Notes
        -----
        If no arguments specified, defaults to download the file as the original filename
        in the current working directory. If `file_obj` given, will write data into the
        passed file-like object. If `path` given, will download the file to the path provided,
        but will not overwrite any existing files.
        """
        return self._download(
            "download_csv",
            "annotation_set_" + self.id + ".csv",
            use_potion_session=True,
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )


class Assemblies(OneCodexBase, ResourceDownloadMixin):
    _resource_path = "/api/v1_experimental/assemblies"

    def download(self, path=None, file_obj=None, progressbar=False):
        """Download an Assembly in FASTA format.

        Parameters
        ----------
        path : `string`, optional
            Full path to save the file to. If omitted, the file will be saved in the current working
            directory. The filename will include the taxon name (if available), followed by the
            genome name (if available), followed by the genome UUID. If this assembly is not
            associated with a genome, the filename will be ``assembly_<UUID>.fasta``.
        file_obj : file-like object, optional
            Rather than save the file to a path, write it to this file-like object.
        progressbar : `bool`
            Display a progress bar using Click for the download?

        Returns
        -------
        `string`
            The path the file was downloaded to, if applicable. Otherwise, None.

        Notes
        -----
        Existing paths will not be overwritten.
        """
        return self._download(
            "download_uri",
            # Set `_filename` to `None` to use the filename provided by the server.
            _filename=None,
            use_potion_session=False,
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )


class Genomes(OneCodexBase):
    _resource_path = "/api/v1_experimental/genomes"

    def __repr__(self):
        return "<Genome {} {} ({})>".format(self.id, self.taxon.name, self.name)


class Taxa(OneCodexBase):
    _resource_path = "/api/v1_experimental/taxa"

    def __repr__(self):
        return "<Taxa {} {} ({})>".format(self.taxon_id, self.name, self.rank)

    def genomes(self):
        """Return a list of all Genomes belonging to descendants of this Taxon."""
        return ResourceList(self._resource.genomes(), Genomes)

    def parents(self):
        """Return a list of all parents of this Taxon, at all ranks."""
        return ResourceList(self._resource.parents(), Taxa)


class FunctionalProfiles(Analyses):
    _resource_path = "/api/v1_experimental/functional_profiles"

    def _filtered_results(self, functional_group, metric, taxa_stratified):
        return self._resource.filtered_results(
            functional_group=functional_group, metric=metric, taxa_stratified=taxa_stratified
        )

    def results(self, json=True):
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

    def table(self, annotation="all", taxa_stratified=True):
        """Return a results table for the functional analysis.

        Parameters
        ----------
        annotation : {'all', onecodex.lib.enum.FunctionalAnnotation}, optional
            Either return a table with all annotations or one of `onecodex.lib.enum.FunctionalAnnotation`
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
        if annotation != "all":
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
        return (
            results_df[results_df["taxa_stratified"]]
            if taxa_stratified
            else results_df[~results_df["taxa_stratified"]]
        )

    def filtered_table(self, annotation, metric, taxa_stratified=True):
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
            functional_group=annotation, metric=metric, taxa_stratified=taxa_stratified
        )
        if not result_json["table"]:
            return pd.DataFrame(
                {
                    "id": pd.Series(dtype="str"),
                    "name": pd.Series(dtype="str"),
                    "value": pd.Series(dtype="float"),
                    "taxon_id": pd.Series(dtype="str"),
                    "taxon_name": pd.Series(dtype="str"),
                }
            )
        return pd.DataFrame(result_json["table"])
