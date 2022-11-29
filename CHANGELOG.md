# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.11.0] - 2022-11-28

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

-   Adds more pythonic support for truthy/falsey values to SampleCollection.filter()
-   Adds support for passing lists to the various sort_x and sort_y plotting arguments
-   Adds a property to SampleCollection to track whether we think the collection is all WGS/metagenomic data
-   Adds support for specifying weighted_unifrac and unweighted_unifrac as metrics in the SampleCollection.beta_diversity() method
-   Adds support for calculating alpha and beta diversities on normalized data. This is important since we'll now be using abundances by default for most datasets
-   Adds support for specifying chart width and height directly in the plot\_\* functions
-   Adds option to plot "Other" bars on bargraphs with normalized read counts and abundances
-   Adds option to plot "No \<level\>" bars on bargraphs with abundances
-   Adds abundance rollups so we can use abundances at all taxonomic ranks
-   Adds a project column to the metadata DataFrame of a SampleCollection

### Changed

-   Changes default alpha diversity metric from simpson to shannon, since shannon is generally a more appropriate default
-   rank="auto" now defaults to species if the field="abundances" or the dataset is metagenomic, instead of genus
-   Switches to using a normalized tree for Weighted Unifrac calculations, which gives us Unifrac values in a [0, 1] range
-   Defaults to using abundances data for metagenomic datasets instead of readcount_w_children

### Deprecated

-   The `field` kwarg on the `SampleCollection` constructor has been renamed to `metric`. We still support passing either, but show a `DeprecationWarning` when passing `field`

### Removed

-   Removes support for just specifying `unifrac` as a metric
