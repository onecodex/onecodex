# One Codex API - Python Client Library and CLI

[![test](https://github.com/onecodex/onecodex/actions/workflows/test.yml/badge.svg)](https://github.com/onecodex/onecodex/actions/workflows/test.yml)
[![pre-commit](https://github.com/onecodex/onecodex/actions/workflows/pre-commit.yml/badge.svg)](https://github.com/onecodex/onecodex/actions/workflows/pre-commit.yml)

Command line interface (CLI) and Python client library for interacting with the One Codex v1 API ([API docs](https://docs.onecodex.com)).

MAINTAINERS: [@clausmith](https://github.com/clausmith), [@boydgreenfield](https://github.com/boydgreenfield)

# Installation

This package provides 3 major pieces of functionality: (1) a core Python client library; (2) a simple CLI for interacting with the One Codex platform that uses that core library; and (3) optional extensions to the client library, which offers many features aimed at advanced users and provides functionality for use in interactive notebook environments (e.g., IPython notebooks).

Python 3.10 or later is required. **Python 2 is no longer supported.**


### _Basic installation_
The CLI (and core Python library) may be simply installed using `pip`. To download a minimal installation (#1 and #2), simply run:
```shell
pip install onecodex
```

### _Installation with optional extensions_
To also download the optional extensions to the client library, and all of their dependencies, run:

```shell
pip install 'onecodex[all]'
```

# Using the CLI

## Logging in
The CLI supports authentication using either your One Codex API key or your One Codex username and password.
To log in using your username and password:

```shell
onecodex login
```

This command will save a credentials file at `~/.onecodex`, which will then automatically be used for authentication the next time the CLI or Python client library are used (OS X/Linux only). You can clear this file and remove your API key from your machine with `onecodex logout`.

In a shared environment, we recommend directly using your One Codex API key, rather than logging in and storing it in a credentials file. To use API key authentication, simply pass your key as an argument to the `onecodex` command:
```shell
onecodex --api-key=YOUR_API_KEY samples
```

Your API key can be found on the [One Codex settings page](https://app.onecodex.com/settings) and should be 32 character string. You may also generate a new API key on the settings page in the web application. _Note_: Because your API key provides access to all of the samples and metadata in your account, you should immediately reset your key on the website if it is ever accidentally revealed or saved (e.g., checked into a GitHub repository).

## Uploading files
The CLI supports uploading FASTA or FASTQ files (optionally gzip compressed) via the `upload` command.
```shell
onecodex upload bacterial_reads_file.fq.gz
```

Multiple files can be uploaded in a single command as well:
```shell
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


## Resources
The CLI supports retrieving your One Codex samples and analyses. The following resources may be queried:

* Your samples (`Samples`)

* Sample metadata (`Metadata`)

* `Analyses`, which include several subtypes with additional functionality and fields:
    * `Classifications`, which are basic metagenomic classification results for your samples
    * `Panels`, which are _in silico_ panels for particular genes or other functional markers ([example on One Codex](https://app.onecodex.com/panel/sample))

* `Jobs`, which provide information on the name, version, and type of analysis which was performed for a given `Analyses`


Simply invoke the `onecodex` command, using one of the above resource names as a subcommand (all lowercase). For example:
```shell
# fetch all your samples
onecodex samples

# fetch a list of panels based on their ids
onecodex panels 0123456789abcdef 0987654321fdecba
```

# Using the Python client library

## Initialization
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

## Resources

Resources are exposed as attributes on the API object. You can fetch a resource directly by its ID or you can fetch it using the query interface. Currently you can access resources using either `get()` or `where()`. If you need help finding the ID for a sample, its identifier is part of its url on our webpage: e.g. for an analysis at `https://app.onecodex.com/analysis/public/1d9491c5c31345b6`, the ID is `1d9491c5c31345b6`. IDs are all short unique identifiers, consisting of 16 hexadecimal characters (`0-9a-f`).

```python
sample_analysis = ocx.Classifications.get("1d9491c5c31345b6")   # Fetch an individual classification
sample_analysis.results()  # Returns classification results as JSON object
sample_analysis.table()    # Returns a pandas dataframe
```

In addition to methods on individual instances of a given resource (e.g., a `Sample` or an `Analysis`), the library also provides methods for aggregating sets of samples or analyses:


```python
all_completed_analyses = ocx.Classifications.where(complete=True)
all_completed_analyses.to_otu()   # Returns classification results as JSON object
all_completed_analyses.to_df()    # Returns a pandas dataframe
```

# Development

## Environment Setup

Before developing, `git` and `python` version >=3.10 are needed. We recommend using [uv](https://github.com/astral-sh/uv) for Python version management and dependency installation.

To download the client library from GitHub:

```shell
git clone https://github.com/onecodex/onecodex.git
cd onecodex/
```

To set up the project, install dependencies using `uv`:

```shell
# If you are on a M1 Macbook, run the line below, adjusting the version as needed
export HDF5_DIR=/opt/homebrew/Cellar/hdf5/1.12.1_1/

uv sync --all-extras --dev --locked
```

To activate the virtual environment:

```shell
source .venv/bin/activate
```

Tests are run via `pytest` while code formatting and linting is done using [`ruff`](https://github.com/astral-sh/ruff):

```shell
make lint
make test
```

We use [`pre-commit`](https://pre-commit.com) for automated linting using [`ruff`](https://github.com/astral-sh/ruff) and various whitespace and newline formatters during development.

## Writing Unit Tests

We use [pytest](https://docs.pytest.org/) as our unit testing framework. Tests should be able to run without an internet connection, and One Codex API calls must be mocked. We use [responses](https://github.com/getsentry/responses) to mock API responses.

> **Tip:** Any API calls that do not have a matching mock will raise an error. You can figure out which API calls need to be mocked by writing a test, running it, and inspecting the error message to see which route(s) are missing.

> **Warning:** Mocked URLs *without* a query string will ignore query strings in any matching requests. If the mocked URL *includes* a query string, it will be used when matching requests.

### Fixtures

These pytest fixtures may be helpful when writing unit tests:

- `ocx`: this is a mocked `Api` object that uses the One Codex v1 API schema.
- `api_data`: this mocks some v1 API data.

### Mocking API Data

API data are stored in `tests/data/api/`:

```
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

Things that are *not* supported by mocking in `tests/data/api/`:

- Non-GET requests (e.g. DELETE)
- Query parameters

# Jupyter Notebook Custom Exporters

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
