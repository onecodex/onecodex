from onecodex.models.base import OneCodexBase
from onecodex.models.helpers import ResourceDownloadMixin
from onecodex.models.schemas.genome import (
    TaxonSchema,
    GenomeSchema,
    AnnotationSetSchema,
    AssemblySchema,
)

# NOTE: these models are only accessible via `X-OneCodex-Api-Experimental` as of 10/2/2025


class AnnotationSets(OneCodexBase, AnnotationSetSchema, ResourceDownloadMixin):
    _resource_path = "/api/v1/annotation_sets"

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
            "annotation_set_" + self.id + ".gbk",
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
            use_client_session=True,
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )


class Assemblies(OneCodexBase, AssemblySchema, ResourceDownloadMixin):
    _resource_path = "/api/v1/assemblies"

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
            path=path,
            file_obj=file_obj,
            progressbar=progressbar,
        )


class Genomes(OneCodexBase, GenomeSchema):
    _resource_path = "/api/v1/genomes"

    def __repr__(self):
        return f"<Genome {self.id} {self.taxon.name} ({self.name})>"


class Taxa(OneCodexBase, TaxonSchema):
    _resource_path = "/api/v1/taxa"

    def __repr__(self):
        return f"<Taxa {self.taxon_id} {self.name} ({self.rank})>"

    def genomes(self):
        """Return a list of all Genomes belonging to descendants of this Taxon."""
        resp = self._client.get(
            f"{self._api._base_url}{self._resource_path}/{self.id}/genomes?expand=all",
        )
        resp.raise_for_status()
        return [Genomes.model_validate(genome) for genome in resp.json()]

    def parents(self):
        """Return a list of all parents of this Taxon, at all ranks."""
        resp = self._client.get(
            f"{self._api._base_url}{self._resource_path}/{self.id}/parents?expand=all",
        )
        resp.raise_for_status()
        return [Taxa.model_validate(taxon) for taxon in resp.json()]
