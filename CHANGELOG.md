# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added support for pandas 3. pandas 2 remains supported.
- Added filtered readcount metrics (`filtered_readcount`,
  `filtered_readcount_w_children`, and their normalized variants), which are
  computed after [artifact
  filtering](https://docs.onecodex.com/en/articles/3761205-one-codex-database)
  and are consistent across samples with and without abundance estimates.

### Changed

- Renamed the `Metric.is_filtered_readcount_metric` property to
  `Metric.is_abundance_sensitive` to better describe what it flags.

### Removed

- Removed the `unfiltered_readcount` metrics (`Metric.UnfilteredReadcount`,
  `Metric.NormalizedUnfilteredReadcount`, etc.) in favor of the new
  `filtered_readcount` metrics. Note the semantics differ: unfiltered
  readcounts were computed *before* artifact filtering, whereas filtered
  readcounts are computed *after*.

### Fixed

- `plot_functional_heatmap()` no longer raises a `MergeError` when the sample
  collection contains more than one sample without functional profile results.
  Samples lacking functional results are now dropped from the plot (they were
  never plotted anyway).
- Plots no longer raise a `ValueError` when a metadata field name contains a colon character (`:`).

## [v1.1.0] - 2026-06-18

### Added

- Added `updated_at` to `Analyses` schema
- `onecodex.models.Workflows` (a type of `onecodex.models.Analyses`) for [Custom Workflows](https://docs.onecodex.com/en/collections/9673389-your-authored-workflows-on-one-codex)
- Added `Jobs.run()` for running analyses including Workflows
- `onecodex.models.Assets` for querying and creating Workflow assets
- `Analyses.await_completion()` blocks until an analysis reaches a terminal
  state using adaptive polling.
- `onecodex analyses await <analysis_id>` CLI command, with `--timeout`,
  `--initial-interval`, and `--max-interval` options. Also available as
  `onecodex await <analysis_id>` for backward compatibility.
- `Analyses.logs()` fetches the job run logs for an analysis, with an optional
  `tail` parameter to limit the output to the last N lines.
- `onecodex analyses logs <analysis_id>` CLI command, with `--tail` to limit
  the number of log lines returned (defaults to 1000). Only available for
  custom workflow runs.
- `Jobs.create()` and `Jobs.update()` for creating and updating custom jobs
  from the Python client. `assets=[...]` accepts `Assets` instances directly;
  `dependencies=[{"job": parent_job, "output_dir": "..."}]` accepts `Jobs`
  instances directly.
- `onecodex jobs create` and `onecodex jobs update <job_id>` CLI commands.
- `onecodex download samples` accepts `-s/--sample` (repeatable) to download
  specific samples by ID. Mutually exclusive with `--project` and `--tags`.
- `onecodex jobs run` accepts `--args-json` to pass runtime arguments as a JSON
  object, preserving non-string types (integers, floats, booleans, arrays,
  objects). `-a/--arg` only supports string values; use `--args-json` when a job
  argument expects another type. Mutually exclusive with `-a/--arg`.
- `ocx.Jobs.get("...").details()` for retrieving CPU/memory/script/etc for Custom Workflows.

### Changed

- `ApiRef` now accepts any object with a `field_uri` attribute (e.g. resource
  model instances), in addition to the existing `$ref`/`$uri` dict shapes. This
  lets users pass fetched resources directly into create/update calls.
- `OneCodexBase.create()` now validates kwargs against the schema declared in
  `_allowed_methods["create"]` (when one is set), mirroring how `update()`
  already worked. Surfaces server error messages instead of failing on
  response parsing.

### Fixed

- `plot_pcoa` no longer flips axis orientation between runs; PCoA eigenvector
  signs are now pinned so the entry with the largest absolute value on each
  axis is positive.
- CLI commands that take One Codex IDs (`onecodex samples`, `analyses`,
  `classifications`, `panels`, `workflows`, `documents download`,
  `jobs run`/`update`, `upload --sample-id`, `download samples --sample`) now
  validate the format up front and fail with a helpful error instead of an
  opaque 404.

## [v1.0.2] - 2026-05-06

### Added

- Unfiltered readcount metrics added to classification results (#621)
- Creating a results dataframe from a `SampleCollection` with samples with and
  without abundance estimates now raises a warning
- A new model `onecodex.models.Analysis.Mlsts` for Multi-Locus Sequence Type (MLST) analyses

### Fixed

- Fixed `plot_bargraph` "Other" value when using a non-normalized metric (#614)
- Fixed wrong exception type raised when all taxa are filtered out (#594)
- Fixed limit functionality in sample queries (#612)
- Samples without abundances are no longer incorrectly dropped during alpha diversity calculations (#622)

## [v1.0.1] - 2026-03-04

### Fixed

- Samples without abundances were filtered out of heatmaps and stats functions
  when using non-abundance metrics

## [v1.0.0] - 2026-02-09

This release contains breaking changes to plotting and classification results functions such as how `SampleCollection` is instantiated and how `metric`, `diversity_metric`, and `rank` are passed to `SampleCollection.to_df` and plotting functions.

Previously, `metric` and `normalize` were determined when a `SampleCollection` was instantiated. Now, those parameters are passed to `.to_df` and plotting functions separately, allowing users to change metrics at plot time.

**Example:**

```python
##
# in 0.19.x
collection = onecodex.models.SampleCollection(samples, metric='readcount_w_children', normalize=True)
# produce a dataframe using readcount_w_children as the metric
collection.to_df()
# to plot a different metric, create a new SampleCollection
collection = onecodex.models.SampleCollection(samples, metric='abundance_w_children')
collection.plot_bargraph()

##
# in 1.0.0
collection = onecodex.models.SampleCollection(samples)
collection.to_df(metric='normalized_readcount_w_children')
# switch metrics without re-instantiating
collection.plot_bargraph(metric='abundance_w_children')
# for backwards compatibility, the following still works but prints a warning
collection = onecodex.models.SampleCollection(samples, metric='readcount_w_children')
# defaults to the metric provided when the collection was created
collection.to_df()
```

### Added

- `metric` is now an argument for all plot and stats functions. By default, it is `auto`, which is determined by `SampleCollection.automatic_metric` based on existing logic
- `SampleCollection.__repr__` has been updated to differentiate it from lists: `<SampleCollection[Samples] length=60 : [<Samples 0b2d0b5397324841: "SRR4408293.fastq">, ...`
- Added type annotations to `onecodex.Api()` so that models are now known to static analyzers (e.g., `ocx.Samples`)
- `SampleCollection.automatic_rank` now returns the best default `onecodex.enums.Rank` to use based on the contained results
- `SampleCollection.automatic_metric` now returns the best `onecodex.enums.Metric` value to use based on the contained results

### Changed

- You can change plotting metrics without re-instantiating the `SampleCollection`
- In alpha- and beta-diversity related functions (e.g., `plot_mds`), the argument `metric` was changed to either `diversity_metric` or `distance_metric`. The `metric` argument in those functions now corresponds to the abundance metric (see `onecodex.enums.Metric`)
- `onecodex.notebooks.report` functions no longer raise a `OneCodexException` if used outside of IPython.

### Removed

- `ClassificationDataframe.ocx` has been removed
- `SampleCollection()` no longer takes `metric` or `normalize` arguments. Those have been merged into a single `metric=` argument that takes a string or `onecodex.enums.Metric` enum value, which is now passed into dataframe and plotting functions directly
- `normalize` has been removed. Normalized metrics are now available as individual `Metric` enum values
- Removed `test-python3-w-simple-json` from CI. It randomly started breaking. We don't think we need it anymore since we've moved away from Potion

### Fixed

- `readcount` being used instead of the intended metric when clustering samples for PCA/Heatmap when plotting from a `SampleCollection`
- Plot title and table export column names have been updated to be more specific and match custom plots

## [v0.19.5] - 2026-02-06

### Changed

- `ocx.Samples` and `ocx.Metadata` now use a more performant, cursor-based pagination method

## [v0.19.4] - 2026-01-21

### Added

- Adds `onecodex interleave` command for generating interleaved fastq files
  from a pair of R1+R2 files

## [v0.19.3] - 2026-01-09

### Fixed

- Fixes `cost.amount` should be an integer

## [v0.19.2] - 2026-01-08

### Added

- Adds `Analyses.get_files` and `Analyses.download_file` methods to help with analyses' output files.
- Adds `cost` to `Analyses`
- Adds `dependencies` to `Analyses`

## [v0.19.1] - 2025-10-14

### Fixed

- Fixes `TypeError` in `ocx.Api()` constructor

### Changed

- Lowers minimum `pydantic` version to 2.10.6

## [v0.19.0] - 2025-10-06

### Added

- Adds Python 3.13 support
- Adds numpy v2 support

### Changed

- Upgrades Altair to 5.5.0
- Upgrades scikit-bio to 0.7.0. **IMPORTANT:** This changes the log base used in Shannon alpha diversity calculations from base 2 to `e` (natural log)
- Upgrades urllib3 to >=2.5.0
- Modernizes packaging and tooling by migrating to `pyproject.toml`, `uv`, and `ruff`
- Updates `pytest` to fail on warnings
- Replaces Potion-Client resources with Pydantic v2 models
- Changes `metric=auto` behavior based on the number of classification results with abundance estimates

### Fixed

- Fixes `SampleCollection` metadata fetching to only match canonical ranks
- Fixes issue when plotting fields that look like JavaScript paths (e.g. `field.name` or `field[name]`)
- Fixes various documentation issues

### Removed

- Removes `Samples.init_multipart_upload()` in favor of `Samples.upload()`

## [v0.18.0] - 2025-04-18

### Added

- Adds `ocx.FunctionalProfiles` to fetch [functional analyses](https://docs.onecodex.com/en/articles/6293070-functional-analysis) results
- Adds `SampleCollection.to_classification_df()` and
  `SampleCollection.to_functional_df()` methods to return classification or
  functional results as DataFrames
- Adds option to omit taxonomic labels from plot mappings

### Changed

- Ensures singleton FASTQ filenames follow the same format as chunked uploads
  (Oxford Nanopore)

### Fixed

- Fixes issue with `onecodex login` requiring an API key for certain accounts

## [v0.17.0] - 2024-12-03

### Added

- Adds `SampleCollection.beta_diversity_stats()` to test for significant differences between groups of samples based on their beta diversity distances (global and posthoc pairwise PERMANOVA tests are supported)
- Adds categorical data palette to `onecodex` Altair theme
- Adds support for passing both `facet_by` *and* `secondary_haxis` to `plot_metadata()`
- Adds filtering and warning if samples are missing abundances in `plot_metadata()`, `plot_pcoa()`, `plot_mds()`, `plot_pca()`, `alpha_diversity_stats()`, and `beta_diversity_stats()`
- Adds docs about Shannon diversity being calculated using log base 2 instead of base `e`

### Changed

- Upgrades Altair to 5.4.1 and replaces `altair_saver` with `vl-convert`
- Upgrades `notebook` to 6.5.7
- Upgrades `WeasyPrint` to 63.0
- Removes `node` and `vega` npm dependencies
- Removes `jupyter-client` dependency
- Removes `selenium` dependency
- Removes pinned `pydyf` version

### Fixed

- Fixes `Api()` call when Altair is not installed
- Fixes scikit-bio and pandas warnings
- Improves error message and docstring in `SampleCollection.to_df(normalize=False)`

### Deprecated

- Deprecates `AlphaDiversityStatsResults.posthoc_df` in favor of `AlphaDiversityStatsResults.posthoc.adjusted_pvalues`

## [v0.16.0] - 2024-08-15

### Added

- Adds optional `secondary_haxis` parameter to `plot_metadata()`, which can be a field name or tuple of field names, if a second grouping is desired
- Adds alpha diversity stats tests, which are run via `SampleCollection.alpha_diversity_stats()`
- Adds support for Python 3.12
- BIOM export now includes canonical taxonomic lineage
- Adds documentation via GitHub Pages at <https://onecodex.github.io/onecodex/>
- Adds support for concatenating ONT files

### Changed

- Adds `scikit-posthocs` dependency
- Pins `scikit-bio` to 0.6.0
- Updates minimum required version of scipy to 1.11.0
- Updates minimum required version of numpy to 1.21.6
- Pins `pydyf<0.11.0`

### Fixed

- Fixes a bug associated with updating custom metadata on a sample
- Fixes metadata non-unique column name bug
- Provides a more useful error message if `plot_heatmap` is called with any samples without abundances calculated
- Fixes abundance chart for groups of samples where all samples lack abundance estimates

### Removed

- Removes Python 3.8 support

## [v0.15.1] - 2024-01-12

### Changed

- Downgrades `notebook` version pin to 6.4.10

## [v0.15.0] - 2024-01-11

### Added

- Adds `coerce_haxis_dates=True` parameter to `plot_metadata()`

### Fixed

- Fixes date field type coercion in `plot_metadata()`

## [v0.14.0] - 2024-01-05

### Added

- Adds new plotting method `plot_functional_heatmap` to plot a heatmap with specified annotations and metrics of functional profile data

### Changed

- Pins `urllib3<2`
- Updates `jupyter-client` version pin to 8.6.0
- Updates `notebook` version pin to 7.0.6
- Updates `nbconvert` minimum version to 6.4.3
- Updates `WeasyPrint` version pin to 60.1
- Removes `jinja2` version pin
- Replaces `imp` module usage with `importlib` in preparation for Python 3.12 support

### Fixed

- Fixes `FutureWarning`s generated by pandas, scikit-learn, etc.
- `make test` no longer generates warnings
- Fixes notebook exporter commands in `README.md`

## [v0.13.0] - 2023-09-25

### Added

- Adds `name` option to `assets upload` command in experimental API mode

### Fixed

- Fixes heatmap plotting for samples without a specified `haxis` value
- Fixes installation by monkeypatching altair-saver for npm >=9

## [v0.12.0] - 2023-08-30

### Added

- Adds support for Python 3.11
- Adds `assets upload` command in experimental API mode

### Changed

- Handles heatmap display of samples with no abundances calculated
- Excludes samples with no abundances calculated from distance heatmap with abundances metric
- Pins altair to 4.2.2
- Requires jsonschema version >= 3.0
- Improves retry handling
- Pins version of altair-saver that includes a fix for npm >= 9

### Fixed

- Fixes unnecessary global loading of numpy and pandas
- Fixes chart concatenation

### Removed

- Removes support for Python 3.7
- Removes opening links in new tab functionality (can be added manually to a returned chart with `usermeta`)

## [v0.11.0] - 2022-11-29

### Added

- Adds support for testing experimental API
- Adds functionality to `SampleCollection` objects to generate tabular data for `FunctionalRun`s using the `.to_df()` method
- Adds `include_taxa_missing_rank` parameter to `to_df()` method for including taxa that do not have a designated parent at `rank` (will be grouped into a "No <rank>" column) (metrics `readcount_w_children` and `abundance_w_children` are supported)
- Adds support for `readcount_w_children` metric to `plot_bargraph()` when passing `include_taxa_missing_rank` (only `abundance_w_children` was supported previously)
- Adds increased customization of plot legend; `plot_bargraph()`'s `legend` parameter now accepts an `altair.Legend` instance
- Adds `group_by` support to `plot_bargraph`
- Adds support for linking to NCBI taxonomy browser from plots

### Changed

- Pins selenium to <4.3.0
- Pins scipy<1.9
- Default plot legend is bigger (displaying up to 40 items)

### Fixed

- Fixes metadata fetching to support fields that match neither metadata nor taxa instead of raising an error
- Fixes box plot display bug for samples with same value
- Fixes pandas `SettingWithCopyWarning` previously raised by `to_df` and `plot_bargraph`
- Fixes `plot_bargraph` display bug for empty samples (total count or abundance of zero)
- Fixes `None` artifact in PDFs

## [v0.10.0] - 2022-04-05

### Added

- Adds `download_samples()` function and `onecodex download samples` CLI command for downloading batches of samples as FASTA/Q files
- Adds support for Python 3.9 and 3.10
- Adds support for Apple M1 processors

### Changed

- Improves `SampleCollection._collate_results()` runtime by at least 20x
- Improves missing value coercion in `SampleCollection._collate_results()`
- Updates some dependency versions to address security vulnerabilities, performance enhancements, expanded environment support, and easier installation
- `numpy` is no longer required to be installed prior to installing the `onecodex` package
- Relaxes `filelock` dependency version pin
- Changes `help@onecodex.com` to `support@onecodex.com` in user messaging

### Fixed

- Fixes PCoA/MDS plot URL bug
- Fixes `metric` parameter behavior in `SampleCollection` constructor (it is no longer ignored)

### Removed

- Removes Python 2 support
- Removes bash completion
- Removes `Classifications._append_abundance_rollups()` in favor of retrieving `abundance_w_children` metric from One Codex API

## [v0.9.6] - 2021-11-15

### Added

- Adds support for functional analysis profiles as experimental feature
- Adds improved error messaging for `Sample` methods

### Changed

- Pin jupyter-client to version 6.1.12

### Fixed

- Fixes issue with credential file lock

### Removed

- Removes `Classifications.abundances` method

## [v0.9.5] - 2021-03-16

### Added

- Adds support for concatenating sequencing data from multiple lanes
- Adds support for faceting metadata plots

### Changed

- Upgrades minimum sentry-sdk version
- Adds explicit legend title on distance heatmaps
- Adds support for zero-abundance samples in normalization calculations
- Warns when empty boxes are present in faceted boxplots
- Supports more paired file naming schemes for automatic interleaving
- Improves `SampleCollection._collate_results` runtime
- Adds support for an empty domain in `interleave_palette`
- Improves error message handling for samples in an importing state

### Fixed

- Fixed a race condition when running multiple instances of the CLI
- Fixed issues handling symmetry and NaN's in beta diversity calculations
- Fixed beta diversity heatmap layout issues
- Fixed PCA/MDS/PCoA colouring issues
- Fixed broken taxa bargraph links

## [v0.9.4] - 2020-10-27

### Added

- Adds Aitchison distance and supports `manhattan` as a synonym for `cityblock`

### Changed

- Ignores vega-lite warning about boxplots not supporting selection
- Displays facet field name and values below x-axis for taxa barplot and heatmap
- Filters out NaNs from alpha diversity plots

### Fixed

- Improves support for networks with a custom CA
- Fixes support for responsive plots (`width=None` and `width="container"` in `plot_metadata`)

## [v0.9.3] - 2020-08-06

### Changed

- Replaces chao1 with observed_taxa
- Changes the experimental assembly download api to support async retry logic
- Upgrades from deprecated Raven library to sentry_sdk

## [v0.9.2] - 2020-07-23

### Fixed

- Fix bug preventing calculation of weighted Unifrac metric due to taxonomy structure assumptions in `scikit-bio`

## [v0.9.1] - 2020-07-10

### Fixed

- Include fonts for report generation in PyPI distribution

## [v0.9.0] - 2020-07-10

### Changed

- Overhaul styling of Jupyter notebook PDF exports
- Add One Codex theme for Altair plots (enabled by default)
- Refactor how `notebooks.report.references` helper works internally
- Improve upload command to display sample filenames on cancelation
- Bump Altair minimum required library version to 4.1.0

## [v0.8.2] - 2020-06-11

### Fixed

- Fix upload callback retry handling broken in v0.8.1 (c87eb76f0dea544b742c5eacf0a4dfcdcebb87bc)

## [v0.8.1] - 2020-06-08

### Changed

- Improved retry handling for confirm callback POSTs during Sample uploads

### Fixed

- Fix upload and export of PDF reports to One Codex documents portal
- Fix example notebook links to follow updated v0.8.0 plotting conventions

## [v0.8.0] - 2020-06-02

### Added

- Adds more pythonic support for truthy/falsey values to SampleCollection.filter()
- Adds support for passing lists to the various sort_x and sort_y plotting arguments
- Adds a property to SampleCollection to track whether we think the collection is all WGS/metagenomic data
- Adds support for specifying weighted_unifrac and unweighted_unifrac as metrics in the SampleCollection.beta_diversity() method
- Adds support for calculating alpha and beta diversities on normalized data. This is important since we'll now be using abundances by default for most datasets
- Adds support for specifying chart width and height directly in the plot\_\* functions
- Adds option to plot "Other" bars on bargraphs with normalized read counts and abundances
- Adds option to plot "No \<level\>" bars on bargraphs with abundances
- Adds abundance rollups so we can use abundances at all taxonomic ranks
- Adds a project column to the metadata DataFrame of a SampleCollection

### Changed

- Changes default alpha diversity metric from simpson to shannon, since shannon is generally a more appropriate default
- rank="auto" now defaults to species if the field="abundances" or the dataset is metagenomic, instead of genus
- Switches to using a normalized tree for Weighted Unifrac calculations, which gives us Unifrac values in a [0, 1] range
- Defaults to using abundances data for metagenomic datasets instead of readcount_w_children

### Deprecated

- The `field` kwarg on the `SampleCollection` constructor has been renamed to `metric`. We still support passing either, but show a `DeprecationWarning` when passing `field`

### Removed

- Removes support for just specifying `unifrac` as a metric
