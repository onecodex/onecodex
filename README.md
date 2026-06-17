# One Codex API - Python Client Library and CLI

[![test](https://github.com/onecodex/onecodex/actions/workflows/test.yml/badge.svg)](https://github.com/onecodex/onecodex/actions/workflows/test.yml)
[![pre-commit](https://github.com/onecodex/onecodex/actions/workflows/pre-commit.yml/badge.svg)](https://github.com/onecodex/onecodex/actions/workflows/pre-commit.yml)

Command line interface (CLI) and Python client library for interacting with the One Codex v1 API.

**Documentation: [API](https://developer.onecodex.com) | [Python](https://onecodex.github.io/onecodex/)**

MAINTAINERS: [@clausmith](https://github.com/clausmith), [@boydgreenfield](https://github.com/boydgreenfield)

## Installation

This package provides 3 major pieces of functionality: (1) a core Python client library; (2) a simple CLI for interacting with the One Codex platform that uses that core library; and (3) optional extensions to the client library, which offers many features aimed at advanced users and provides functionality for use in interactive notebook environments (e.g., IPython notebooks).

Python 3.10 or later is required. **Python 2 is no longer supported.**

## _Basic installation_

The CLI (and core Python library) may be simply installed using `pip`. To download a minimal installation (#1 and #2), simply run:

```bash
pip install onecodex
```

### _Installation with optional extensions_

To also download the optional extensions to the client library, and all of their dependencies, run:

```bash
pip install 'onecodex[all]'
```

## Using the CLI

### Logging in

The CLI supports authentication using either your One Codex API key or your One Codex username and password.
To log in using your username and password:

```bash
onecodex login
```

This command will save a credentials file at `~/.onecodex`, which will then automatically be used for authentication the next time the CLI or Python client library are used (OS X/Linux only). You can clear this file and remove your API key from your machine with `onecodex logout`.

In a shared environment, we recommend directly using your One Codex API key, rather than logging in and storing it in a credentials file. To use API key authentication, simply pass your key as an argument to the `onecodex` command:

```bash
onecodex --api-key=YOUR_API_KEY samples
```

Your API key can be found on the [One Codex settings page](https://app.onecodex.com/settings) and should be 32 character string. You may also generate a new API key on the settings page in the web application. _Note_: Because your API key provides access to all of the samples and metadata in your account, you should immediately reset your key on the website if it is ever accidentally revealed or saved (e.g., checked into a GitHub repository).

### Uploading files

The CLI supports uploading FASTA or FASTQ files (optionally gzip compressed) via the `upload` command.

```bash
onecodex upload bacterial_reads_file.fq.gz
```

Multiple files can be uploaded in a single command as well:

```bash
onecodex upload file1.fq.gz file2.fq.gz ...
```

You can also upload files using the Python client library:

```python
uploaded_sample1 = ocx.Samples.upload("/path/to/file.fastq")

# Or upload a tuple of paired end files
uploaded_sample2 = ocx.Samples.upload(("/path/to/R1.fastq", "/path/to/R2.fastq"))
```

Which returns a `Samples` resource (as of `0.5.0`). Samples can be associated with tags, metadata, and projects at upload timing using those respective keyword arguments:

```python
# Note format must match the schema defined for our API, with arbitrary
# metadata allowed as a single-level dictionary in the `custom` field.
# See https://developer.onecodex.com/api-reference/metadata-resource for details.
metadata = {
    "platform": "Illumina NovaSeq 6000",
    "date_collected": "2019-04-14T00:51:54.832048+00:00",
    "external_sample_id": "my-lims-ID-or-similar",
    "custom": {
        "my-string-field": "A most interesting sample...",
        "my-boolean-field": True,
        "my-number-field-1": 1,
        "my-number-field-2": 2.0,
    }
}
```

Uploads can be made in parallel using Python threads (or multiple processes), e.g.:

```python
import concurrent.futures
uploaded_samples = []

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    futures = {executor.submit(ocx.Samples.upload, file) for file in LIST_OF_FILES}
    for future in concurrent.futures.as_completed(futures):
        try:
            uploaded_samples.append(future.result())
        except Exception as e:
            print("An execption occurred during your upload: {}".format(e))
```

### Resources (CLI)

The CLI supports retrieving One Codex samples, analyses and other resources. For a complete list, see [the documentation](https://onecodex.github.io/onecodex/models.html).

- Your samples (`Samples`)
- Sample metadata (`Metadata`)
- `Analyses`, which include several subtypes with additional functionality and fields:
  - `Classifications`, which are basic metagenomic classification results for your samples
  - `Panels`, which are _in silico_ panels for particular genes or other functional markers ([example on One Codex](https://app.onecodex.com/panel/sample))
- `Jobs`, which provide information on the name, version, and type of analysis which was performed for a given `Analyses`

Simply invoke the `onecodex` command, using one of the above resource names as a subcommand (all lowercase). For example:

```bash
# fetch all your samples
onecodex samples

# fetch a list of panels based on their ids
onecodex panels 0123456789abcdef 0987654321fdecba
```

## Using the Python client library

### Initialization

To load the API, use the following import:

```python
from onecodex.api import Api
```

Instantiate an API client either by passing your API key or automatically fetching your credentials from `~/.onecodex` if you've previously called `onecodex login`.

```python
from onecodex.api import Api

# Instantiate a One Codex API object, will attempt to get credentials from ~/.onecodex
ocx = Api()

# Instantiate an API object, manually specifying an API key
ocx = Api(api_key="YOUR_API_KEY_HERE")
```

### Resources

Resources are exposed as attributes on the API object. You can fetch a resource directly by its ID or you can fetch it using the query interface. Currently you can access resources using either `get()` or `where()`. If you need help finding the ID for a sample, its identifier is part of its url on our webpage: e.g. for an analysis at `https://app.onecodex.com/analysis/public/1d9491c5c31345b6`, the ID is `1d9491c5c31345b6`. IDs are all short unique identifiers, consisting of 16 hexadecimal characters (`0-9a-f`). For a complete list of resource models, see [the documentation](https://onecodex.github.io/onecodex/models.html).

```python
sample_analysis = ocx.Classifications.get("1d9491c5c31345b6")   # Fetch an individual classification
sample_analysis.results()  # Returns classification results as JSON object
sample_analysis.table()    # Returns a pandas dataframe
```

## Custom Workflows

One Codex supports creating, running, and retrieving results from custom bioinformatics workflows, including Nextflow pipelines, from the web app as well as programmatically through the API and client library. See ["Your Authored Workflows on One Codex"](https://docs.onecodex.com/en/collections/9673389-your-authored-workflows-on-one-codex) for more information.

The relevant models are:

- `ocx.Jobs` — a runnable workflow definition (script, image, resource requirements, dependencies).
- `ocx.Assets` — files (e.g. references, databases) that can be attached to a job and made available at run time.
- `ocx.Analyses` — the result of running a job against a sample. Subtypes (`Classifications`, `FunctionalProfiles`, `Panels`, `Mlsts`, `Alignments`, `Workflows`) expose result-specific accessors.

### Running a workflow

Workflows are defined by the `ocx.Jobs` model, and can be run on a sample to generate a new analysis (`ocx.Analyses`). Jobs can also take arguments, passed in via the `job_args` keyword argument.

You can launch a workflow against an uploaded sample directly from the Python client. `Jobs.run()` returns the freshly-created `Analyses` instance, which you can then poll for completion:

```python
job = ocx.Jobs.get("0123456789abcdef")
# or look up by name:
job = ocx.Jobs.where(name="my-custom-job")[0]
sample = ocx.Samples.get("fedcba9876543210")

analysis = job.run(sample, job_args={"min_quality": 30})
analysis.await_completion()                # block until the job finishes

if analysis.success:
    results = analysis.results()
else:
    print(f"Job failed: {analysis.error_msg}")
```

`await_completion()` returns the (refreshed) analysis once it reaches a terminal state — including failure. It does _not_ raise on a failed job; check `analysis.success` and `analysis.error_msg` to distinguish success from failure. The only exception it raises is `TimeoutError`, when a `timeout=` is set and exceeded.

The CLI exposes the same functionality:

```bash
onecodex jobs run <job_id> <sample_id> -a min_quality=30
onecodex analyses await <analysis_id>

# Or block in a single step:
onecodex jobs run <job_id> <sample_id> --arg min_quality=30 --await

# Reuse a prior analysis as a dependency, optionally under a relative path:
onecodex jobs run <job_id> <sample_id> -d <analysis_id>
onecodex jobs run <job_id> <sample_id> -d <analysis_id>=parent_out
```

#### Passing arguments

```shell
onecodex jobs run <job_id> <sample_id> --args-json '{"min_quality": 30, "trim": true}'

# To read arguments from a file, use shell substitution:
onecodex jobs run <job_id> <sample_id> --args-json "$(cat args.json)"
```

`-a/--arg` only supports string values — every `key=value` is sent to the server as a
string. If a job argument expects another type (integer, float, boolean, array, object),
use `--args-json` to pass the full argument set as a JSON object, which preserves types:

`--args-json` is mutually exclusive with `-a/--arg`.

### Creating and updating jobs

You can create and update custom jobs from the client.

```python
asset = ocx.Assets.upload("reference.fa.gz")
parent = ocx.Jobs.get("0123456789abcdef")

job = ocx.Jobs.create(
    name="my-custom-job",
    script=open("run.sh").read(),
    image_uri="docker.io/library/python:3.12",
    job_type="shell_script",  # or "nextflow"
    cpu=1, ram_gb=1, storage_gb=1,
    assets=[asset],
    dependencies=[{"job": parent, "output_dir": "parent_out"}],
)

job.update(name="renamed", description="now with a description")
```

The CLI mirrors this:

```bash
onecodex jobs create \
    --name my-custom-job \
    --script ./run.sh \
    --image-uri docker.io/library/python:3.12 \
    --cpu 1 --ram-gb 1 --storage-gb 1 \
    --asset-id <asset_id> \
    -d <parent_job_id>=parent_out

onecodex jobs update <job_id> --name renamed
```

For long-running analyses, `await_completion()` polls until the analysis reaches a terminal state (`complete=True`). The cadence backs off over time, so failures surface in seconds while longer jobs poll on the order of minutes:

```python
analysis = ocx.Analyses.get("0123456789abcdef")
analysis.await_completion()                # block indefinitely
analysis.await_completion(timeout=600)     # raise TimeoutError after 10 minutes
```

For custom workflow runs, `.logs()` returns the job run logs as a string:

```python
analysis = ocx.Analyses.get("0123456789abcdef")
print(analysis.logs())                     # full log
print(analysis.logs(tail=200))             # last 200 lines
```

The method refreshes `analysis` in place and returns it; check `analysis.success` to see whether it finished cleanly. `analysis.refresh()` is also available if you just need to re-fetch the current state without blocking.

In addition to methods on individual instances of a given resource (e.g., a `Sample` or an `Analysis`), the library also provides methods for aggregating sets of samples or analyses:

```python
all_completed_analyses = ocx.Classifications.where(complete=True)
all_completed_analyses.to_otu()   # Returns a BIOM v1 OTU table as an OrderedDict (JSON-serializable)
all_completed_analyses.to_df()    # Returns a pandas dataframe
```

### Awaiting an analysis

To block until an analysis reaches a terminal state, use the `analyses await` subcommand. Polling starts at a few seconds and backs off, so failures surface quickly while long-running jobs don't get hammered:

```bash
onecodex analyses await 0123456789abcdef
onecodex analyses await 0123456789abcdef --timeout 600
```

The command exits non-zero if the analysis finishes unsuccessfully or times out.

### Dependencies

To re-use the output of a previous run as an input to a new one, pass `dependency_overrides`:

```python
from onecodex.models.misc import DependencyOverride

prior = ocx.Analyses.get("abcdef0123456789")
analysis = job.run(
    sample,
    job_args={"k": 31},
    dependency_overrides=[DependencyOverride(analysis=prior)],
)
```

### Fetching analysis logs

To view the job run logs for a custom workflow analysis, use `analyses logs`:

```bash
onecodex analyses logs 0123456789abcdef
onecodex analyses logs 0123456789abcdef --tail 200
```

`--tail` defaults to the last 1000 lines. Logs are only available for custom
workflow runs.

### Fetching results (files) from an analysis

```python
analysis = ocx.Analyses.get("0123456789abcdef")

output_files = analysis.get_files()

for file in output_files:
    analysis.download_file(file, progressbar=True)
```

## Upgrading from 0.19.x to 1.0

In 1.0, `SampleCollection` no longer takes `metric`, `normalize`, or `rank` at construction time. These are now passed directly to `.to_df()` and the `.plot_*()` functions instead, so you can switch metrics or ranks without rebuilding the collection. `normalize=True` has been removed; use the explicit `normalized_*` metric value instead (e.g. `normalized_readcount_w_children`). In alpha- and beta-diversity functions, the old `metric` argument was also renamed to `diversity_metric` / `distance_metric`, with `metric` now referring to the underlying abundance metric.

### Major changes

- `SampleCollection(...)` no longer accepts `metric`, `normalize`, or `rank`. Pass these to `.to_df()` and `.plot_*()` instead.
- `normalize=True` is gone. Use the matching `normalized_*` metric (e.g. `normalized_readcount_w_children`).
- In alpha-/beta-diversity functions, `metric=` was renamed to `diversity_metric=` / `distance_metric=`; `metric=` now refers to the abundance metric.
- `ClassificationDataframe.ocx` has been removed.

### Examples

```python
# Before:

phylum = SampleCollection(samples, metric='readcount_w_children', normalize=True, rank='phylum')
phylum.plot_bargraph()
phylum_df = phylum.to_df()

species = SampleCollection(samples, metric='readcount_w_children', normalize=True, rank='species')
species.plot_heatmap()

abundance = SampleCollection(samples, metric='abundance_w_children', rank='genus')
abundance.plot_pca()

```

In 1.0.x, build the collection once and pick the metric and rank per call. `normalize=True` becomes the corresponding `normalized_*` metric, and the diversity-vs-abundance ambiguity is gone. Diversity functions now take `diversity_metric` or `distance_metric` alongside an abundance `metric`:

```python
# After:

samples = SampleCollection(samples)

samples.plot_bargraph(metric='normalized_readcount_w_children', rank='phylum')
samples.plot_bargraph(metric='normalized_readcount_w_children', rank='species')
phylum_df = samples.to_df(metric='normalized_readcount_w_children', rank='phylum')

samples.plot_heatmap(metric='normalized_readcount_w_children', rank='species')
samples.plot_pca(metric='abundance_w_children', rank='genus')
```

### Mapping `normalize` to a metric

Prior to 1.0, `metric` and `normalize` were separate arguments. Now, they've been merged into a single
`metric` argument with values corresponding to normalized and un-normalized variants:

| 0.19.x                                            | 1.0.x                                       |
| ------------------------------------------------- | ------------------------------------------- |
| `metric='readcount', normalize=True`              | `metric='normalized_readcount'`             |
| `metric='readcount_w_children', normalize=True`   | `metric='normalized_readcount_w_children'`  |
| `metric='readcount_w_children', normalize=False`  | `metric='readcount_w_children'`             |
| `metric='abundance_w_children'`                   | `metric='abundance_w_children'`             |

For a full list of supported metrics and their definitions, see [the documentation](https://onecodex.github.io/onecodex/enums.html#metric).

## Development

## Environment Setup

Before developing, `git` and `python` version >=3.10 are needed. We recommend using [uv](https://github.com/astral-sh/uv) for Python version management and dependency installation.

To download the client library from GitHub:

```bash
git clone https://github.com/onecodex/onecodex.git
cd onecodex/
```

To set up the project, install dependencies using `uv`:

```bash
# If you are on a M1 Macbook, run the line below, adjusting the version as needed
export HDF5_DIR=/opt/homebrew/Cellar/hdf5/1.12.1_1/

uv sync --all-extras --dev --locked
```

To activate the virtual environment:

```bash
source .venv/bin/activate
```

Tests are run via `pytest` while code formatting and linting is done using [`ruff`](https://github.com/astral-sh/ruff):

```bash
make lint
make test
```

We use [`pre-commit`](https://pre-commit.com) for automated linting using [`ruff`](https://github.com/astral-sh/ruff) and various whitespace and newline formatters during development.

## Writing Unit Tests

We use [pytest](https://docs.pytest.org/) as our unit testing framework. Tests should be able to run without an internet connection, and One Codex API calls must be mocked. We use [responses](https://github.com/getsentry/responses) to mock API responses.

> **Tip:** Any API calls that do not have a matching mock will raise an error. You can figure out which API calls need to be mocked by writing a test, running it, and inspecting the error message to see which route(s) are missing.
>
> **Warning:** Mocked URLs _without_ a query string will ignore query strings in any matching requests. If the mocked URL _includes_ a query string, it will be used when matching requests.

### Fixtures

These pytest fixtures may be helpful when writing unit tests:

- `ocx`: this is a mocked `Api` object that uses the One Codex v1 API schema.
- `api_data`: this mocks some v1 API data.

### Mocking API Data

API data are stored in `tests/data/api/`:

```text
tests/data/api
└── v1  # the API version
    ├── ...
    ├── analyses
    │   └── index.json  # payload for accessing GET::api/v1/analyses. Will also be used to mock each resource instance, e.g. GET::api/v1/analyses/<uuid>
    ├── classifications
    │   ├── 0f4ee4ecb3a3412f
    │   │   └── results
    │   │       └── index.json  # payload for accessing GET::api/v1/classifications/0f4ee4ecb3a3412f/results
    │   └── index.json  # payload for accessing GET::api/v1/classifications. Instance routes are also auto-mocked
    └── ...
```

The directory structure mirrors the One Codex API. For example:

- The payload for API route `api/v1/classifications` is stored at `tests/data/api/v1/classifications/index.json`.
- API route `api/v1/classifications/0f4ee4ecb3a3412f/results` has its payload stored at `tests/data/api/v1/classifications/0f4ee4ecb3a3412f/results/index.json`.

This idea can be extended to arbitrary nesting/depths within the API.

> **Note:** If the payload is large, you can gzip it and name it `index.json.gz`.

A resource's instance list payload (e.g. `api/v1/analyses` gives you a list of analyses) is used to auto-mock each resource instance (e.g. `api/v1/analyses/<uuid>`). You don't need to create an `index.json` for each instance.

### conftest.py

API data is loaded in `tests/conftest.py`. If you need to mock API calls in a way that's not supported by this framework, you can add custom mocked calls in `conftest.py`.

Things that are _not_ supported by mocking in `tests/data/api/`:

- Non-GET requests (e.g. DELETE)
- Query parameters

## Jupyter Notebook Custom Exporters

We also package custom Jupyter notebook [`nbconvert`](https://nbconvert.readthedocs.io/en/latest/index.html) exporters. These can be tested with the following snippets and the provided `example.ipynb` report.

Our `OneCodexHTMLExporter`:

```sh
ONE_CODEX_REPORT_FILENAME=example.html jupyter nbconvert --execute --to onecodex_html --ExecutePreprocessor.timeout=-1 --output="$ONE_CODEX_REPORT_FILENAME" --output-dir="." notebook_examples/example.ipynb && open example.html
```

And using the `OneCodexPDFExporter`:

```sh
ONE_CODEX_REPORT_FILENAME=example.pdf jupyter nbconvert --execute --to onecodex_pdf --ExecutePreprocessor.timeout=-1 --output="$ONE_CODEX_REPORT_FILENAME" --output-dir="." notebook_examples/example.ipynb && open example.pdf
```

Note that `OneCodexPDFExporter` requires the `vl-convert-python` package to be installed.

## Docker

Docker images are built against the `master` branch automatically, and available at:

- `ghcr.io/onecodex/onecodex:latest`

For a complete list of version-tagged images, see [Packages](https://github.com/orgs/onecodex/packages?repo_name=onecodex)

The image contains only the base + all dependencies (equivalent of `pip install onecodex[all]`)
