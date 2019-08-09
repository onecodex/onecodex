from onecodex.models import OneCodexBase, ResourceList
from onecodex.models.helpers import ResourceDownloadMixin


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
