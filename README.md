# One Codex API - Python Client Library and CLI

[![Circle CI](https://circleci.com/gh/onecodex/onecodex.png?style=shield&circle-token=d86a8fc55e54a645ee515387db9acee32068a6ad)](https://circleci.com/gh/onecodex/onecodex) [![Coverage Status](https://coveralls.io/repos/github/onecodex/onecodex/badge.svg)](https://coveralls.io/github/onecodex/onecodex)

Command line interface (CLI) and Python client library for interacting with the One Codex v1 API ([API docs](https://docs.onecodex.com)).

MAINTAINERS: [@bovee](https://github.com/bovee), [@boydgreenfield](https://github.com/boydgreenfield)

# Installation

This package provides 3 major pieces of functionality: (1) a core Python client library; (2) a simple CLI for interacting with the One Codex platform that uses that core library; and (3) optional extensions to the client library, which offers many features aimed at advanced users and provides functionality for use in interactive notebook environments (e.g., IPython notebooks).


### _Basic installation_
The CLI (and core Python library) may be simply installed using `pip`. To download a minimal installation (#1 and #2), simply run:
```shell
pip install onecodex
```


### _Installation with optional extensions_
To also download the optional extensions to the client library - and all of their dependencies - use the command `pip install onecodex[all]`. **Warning:** Because other packages used in the extensions rely upon `numpy` being present during their installation, `numpy` must be installed seperately first. So if you do not have `numpy` installed, and you are going to install `onecodex[all]` please do the following:
```shell
# If numpy is not installed in your environment
pip install numpy

# Once you have numpy installed
pip install onecodex[all]
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

## Initalization
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
all_completed_analyses = ocx.Classifications.where(complete=True)
all_completed_analyses.to_otu()   # Returns classification results as JSON object
all_completed_analyses.to_df()    # Returns a pandas dataframe
```

# Development

Before developing, `git` and `python` (version >=2.7 and/or >=3.4) are needed. We recommend using [pyenv](https://github.com/yyuu/pyenv) for Python version management.

To download the client library from GitHub:

```shell
git clone https://github.com/onecodex/onecodex.git
cd onecodex/
```

To set up the project, first create a virtual environment and then install dependencies:

```shell
virtualenv venv
source venv/bin/activate
pip install numpy  # numpy must be installed before any of its dependencies
pip install -r requirements.txt
```

Test are run through the makefile, and call tox. Note this may take awhile at first because of installing dependencies:

```shell
make lint
make test
```
