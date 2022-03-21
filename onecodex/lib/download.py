import logging
import os
import os.path
import warnings

import click

from onecodex.utils import FakeProgressBar
from onecodex.exceptions import OneCodexException

log = logging.getLogger("onecodex")


def get_project(ocx, project_name_or_id):
    project = ocx.Projects.get(project_name_or_id)
    if not project:
        try:
            project = ocx.Projects.where(name=project_name_or_id)[0]
        except IndexError:
            raise OneCodexException(
                "No project exists with name or ID '{}'.".format(project_name_or_id)
            )
    return project


def filter_samples_by_tags(ocx, samples, tag_names):
    tags = set()
    for tag_name in tag_names:
        log.info("Filtering samples by tag '{}'...".format(tag_name))
        try:
            tag = ocx.Tags.where(name=tag_name)[0]
        except IndexError:
            warnings.warn("No tag found: '{}'.".format(tag_name))
            continue
        tags.add(tag)

    if tags:
        # Filter for samples with any of the tags above (union)
        return [sample for sample in samples if len(set(sample.tags) & tags)]
    return samples


def get_download_filename(sample):
    pieces = sample.filename.split(".")
    return ".".join(["{}_{}".format(pieces[0], sample.id)] + pieces[1:])


def download_samples(
    ocx,
    outdir,
    project_name_or_id=None,
    tag_names=None,
    prompt=False,
    min_samples_for_prompt=50,
    progressbar=False,
):
    if project_name_or_id:
        log.info("Fetching samples in project '{}'...".format(project_name_or_id))
        project = get_project(ocx, project_name_or_id)
        samples = ocx.Samples.where(project=project)
    else:
        log.info("Fetching all samples in your account...")
        samples = ocx.Samples.all()

    if tag_names:
        samples = filter_samples_by_tags(ocx, samples, tag_names)

    num_samples = len(samples)
    if num_samples == 0:
        log.info("No samples match the filter criteria; nothing to download.")
        return []

    if prompt and num_samples >= min_samples_for_prompt:
        if not click.confirm("Are you sure you'd like to download {} samples?".format(num_samples)):
            return []

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if progressbar:
        label = "Downloading {} sample{}...".format(num_samples, "" if num_samples == 1 else "s")
        progressbar = click.progressbar(length=num_samples, label=label)
    else:
        progressbar = FakeProgressBar()

    filepaths = []
    with progressbar as bar:
        for sample in samples:
            filepath = os.path.join(outdir, get_download_filename(sample))
            try:
                filepaths.append(sample.download(path=filepath, progressbar=False))
            except OneCodexException as e:
                warnings.warn(
                    "Skipping download of sample {} to {}: {}".format(sample.id, filepath, str(e))
                )
            bar.update(1)
    return filepaths
