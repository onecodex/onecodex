---
name: One Codex
description: The One Codex Python Client provides an object-oriented API for accessing sample data, analyzing and visualizing results, and generating branded reports via Jupyter notebooks rendered to HTML or PDF.
license: MIT
compatibility: Requires Python 3.11+
metadata:
  author: One Codex
  version: "1.0.0"
---

# How to Use This Document

This file is a **SKILL.md** — a structured reference document designed to be loaded into
an LLM-powered coding assistant to give it accurate, up-to-date knowledge of the One Codex
Python client. Without it, assistants may hallucinate incorrect method names, unsupported
filter operators, or non-existent API patterns.

## Loading into your coding tool

**Claude Code** — place this file in your project root. Claude Code will
automatically read files named `SKILL.md` in the working directory.

**Claude.ai / ChatGPT / other chat interfaces** — paste the contents directly into the
conversation, or attach the file at the start of a session before asking questions.

**OpenAI Codex CLI** — pass the file as part of your instructions file or prepend it to
your prompt.

---

# Getting Started

## Installation & Setup

### Basic Installation
```bash
# Core functionality
pip install onecodex==1.0.0

# Full installation with visualization
pip install 'onecodex[all]==1.0.0'
```

### Development Setup
```bash
git clone https://github.com/onecodex/onecodex.git
cd onecodex/
uv sync --all-extras --dev --locked
```

## Authentication

### Interactive Login (Development)
```bash
onecodex login   # Saves credentials to ~/.onecodex
onecodex logout  # Clear credentials (will remove ~/.onecodex)
```

### API Key (Production)
```bash
# Set environment variable
export ONE_CODEX_API_KEY="your_api_key_here"

# Or pass directly
onecodex --api-key your_api_key_here samples
```

### Programmatic Access
```python
import onecodex

# Auto-detect from ~/.onecodex or environment
ocx = onecodex.Api()

# Specify API key directly
ocx = onecodex.Api(api_key="your_api_key_here")

# Custom endpoint (enterprise)
ocx = onecodex.Api(base_url="https://app.onecodex.com")
```
---

# Quick Reference

## Key Resource Types

| Resource | Purpose | Example Usage |
|----------|---------|---------------|
| `Samples` | Sequencing data files | `ocx.Samples.get("sample_id")` |
| `Classifications` | Taxonomic analysis results | `classification.results()` |
| `Analyses` | All analysis types | `ocx.Analyses.where(complete=True)` |
| `Projects` | Sample organization | `sample.project` |
| `Tags` | Sample tags | `ocx.Tags.where(name="trimmed")` |
| `SampleCollection` | Multi-sample operations | `samples.to_classification_df()` |

Experimental models require `ocx = onecodex.Api(experimental=True)`:

| Resource | Purpose | Example Usage |
|----------|---------|---------------|
| `Assets` | File storage and management | `ocx.Assets.upload("file.txt")` |
| `Genomes` | Reference genome browser | `ocx.Genomes.all()` |
| `Assemblies` | Genomic assemblies | `ocx.Assemblies.get("assembly_id")` |
| `AnnotationSets` | Genomic annotations | `ocx.AnnotationSets.all()` |
| `Taxa` | Taxonomic hierarchy browser | `ocx.Taxa.get("taxon_id")` |

## Essential Methods
| Operation | Method | Description |
|-----------|--------|-------------|
| Get DataFrame | `samples.to_classification_df()` | Convert to analysis DataFrame |
| Plot PCA | `samples.plot_pca()` | Sample similarity visualization |
| Plot Composition | `samples.plot_bargraph()` | Taxonomic composition |
| Alpha Diversity | `samples.alpha_diversity()` | Within-sample diversity |
| Beta Diversity | `samples.beta_diversity()` | Between-sample distances |

## Common Enums
```python
from onecodex.lib.enums import Rank, Metric, BetaDiversityMetric, AlphaDiversityMetric

# Taxonomic ranks
Rank.Kingdom, Rank.Phylum, Rank.Class, Rank.Order,
Rank.Family, Rank.Genus, Rank.Species

# Abundance metrics — use as `metric=` argument in plotting/analysis functions
Metric.Readcount, Metric.Abundance, Metric.AbundanceWChildren,
Metric.PropReadcount, Metric.NormalizedReadcount

# Beta-diversity metrics — use as `diversity_metric=` in plot_distance / beta_diversity
BetaDiversityMetric.BrayCurtis, BetaDiversityMetric.Jaccard, BetaDiversityMetric.Aitchison,
BetaDiversityMetric.UnweightedUnifrac, BetaDiversityMetric.WeightedUnifrac

# Alpha-diversity metrics — use as `diversity_metric=` in alpha_diversity
AlphaDiversityMetric.Shannon, AlphaDiversityMetric.Simpson, AlphaDiversityMetric.ObservedTaxa
```

## Query Patterns
```python
# Basic filtering
samples = ocx.Samples.where(limit=50)

# Date filtering — must be a full RFC 3339 timestamp with timezone offset
from datetime import datetime, timedelta, timezone

recent = ocx.Samples.where(updated_at={"$gte": "2024-01-01T00:00:00+00:00"})

# Using Python datetime (recommended)
thirty_days_ago = (datetime.now(timezone.utc) - timedelta(days=30)).isoformat()
recent = ocx.Samples.where(updated_at={"$gte": thirty_days_ago})
older = ocx.Samples.where(created_at={"$lt": "2024-06-01T00:00:00+00:00"})

# Filename partial string match — use $icontains/$contains, NOT $regex (not supported)
samples = ocx.Samples.where(filename={"$icontains": "patient"})   # case-insensitive substring
samples = ocx.Samples.where(filename={"$istartswith": "SRR"})     # case-insensitive prefix
samples = ocx.Samples.where(filename={"$iendswith": ".fastq.gz"}) # case-insensitive suffix

# Filename exclusion
# $ne: exact filename exclusion (server-side)
samples = ocx.Samples.where(filename={"$ne": "excluded_sample.fastq.gz"})

# No server-side "does not contain" operator exists — must be done client-side:
samples = ocx.Samples.where(public=True, limit=100)
results = [s for s in samples if "bad_prefix" not in s.filename.lower()]

# Tag filtering — samples that HAVE a tag (server-side, efficient)
trimmed = ocx.Samples.where(tags=["trimmed"])               # by tag name
trimmed = ocx.Samples.where(tags=["5c1e9e41043e4435"])      # by tag ID
trimmed = ocx.Samples.where(tags=[ocx.Tags.where(name="trimmed")[0]])  # by Tag object

# Tag filtering — samples that DO NOT have a tag (client-side, requires fetching all)
# There is no server-side NOT operator for tags ($containsnone does not exist).
all_samples = ocx.Samples.where()
untrimmed = [s for s in all_samples if not any(t.name == "trimmed" for t in s.tags)]

# Complex queries
filtered = ocx.Samples.where({
    "metadata.sample_type": {"$in": ["clinical", "environmental"]},
    "size": {"$gte": 1000000}
})

# Pagination — Samples uses cursor-based pagination; pages are fetched automatically.
# .where() and .all() retrieve all matching results across all pages by default.
samples = ocx.Samples.where()        # all samples, all pages
samples = ocx.Samples.where(limit=100)  # hard cap at 100 results
```

**Filterable `Samples` fields:** `created_at`, `updated_at`, `filename`, `size`, `status`, `visibility`, `project`, `owner`, `tags`, `primary_classification`, `metadata.*` (own samples only)

**Supported operators:**

| Operator | Accepts | Notes |
|----------|---------|-------|
| `$eq` | str, int, float, bool, datetime, ref, `None` | Implicit when no operator given. `None` → IS NULL |
| `$ne` | str, int, float, bool, datetime, ref, `None` | `None` → IS NOT NULL |
| `$lt` | int, float, datetime | |
| `$lte` | int, float, datetime | |
| `$gt` | int, float, datetime | |
| `$gte` | int, float, datetime | |
| `$between` | `[int, int]` \| `[float, float]` \| `[datetime, datetime]` | Inclusive on both ends |
| `$in` | list of scalars or refs | Scalar and relationship fields |
| `$contains` | str | Case-sensitive substring (`LIKE %val%`) |
| `$icontains` | str | Case-insensitive substring (`ILIKE %val%`) |
| `$startswith` | str | Case-sensitive prefix (`LIKE val%`) |
| `$istartswith` | str | Case-insensitive prefix (`ILIKE val%`) |
| `$endswith` | str | Case-sensitive suffix (`LIKE %val`) |
| `$iendswith` | str | Case-insensitive suffix (`ILIKE %val`) |
| `$containsall` | list of refs | Relationship fields only — all must match (AND) |
| `$containsany` | list of refs | Relationship fields only — any must match (OR) |

### Fetching data in bulk

It is often faster to fetch the data in bulk:

```python
# FAST
# performs filtering on the server
# data is returned in pages (usually ~100 per page)
# related data (such as the primary classification) are often also included by default
samples = ocx.Samples.where(...)

# SLOW
samples = [ ocx.Samples.get(_id) for _id in subset_ids ]

# NOTE: filtering by id with $in is not supported; fetch by ID individually if needed
# samples = ocx.Samples.where({'id': { '$in': subset_ids }})  # does NOT work

classifications = [s.primary_classification for s in samples]

# SLOW
classifications = [ocx.Classifications.get(s.primary_classification.id) for s in subset_ids]

# FAST
samples = ocx.Samples.where(project=project)

# Client-side filtering via SampleCollection.filter() — less efficient, use sparingly
clinical_samples = samples.filter(lambda s: s.metadata.get('sample_type') == 'clinical')
```

## Basic Usage

### Initialize and Test Connection
```python
import onecodex

# Initialize client
ocx = onecodex.Api()

# NOTE: There is no ocx.get_user() method. Do NOT generate code like:
#   ocx.get_user().username  — this does not exist

# List your samples to verify connection
samples = ocx.Samples.where(limit=10)
print(f"You have {len(samples)} recent samples")
```

### Quick Analysis Example
```python
# Load samples for analysis
samples = ocx.Samples.where(project="your_project_id")

# Get classification data as DataFrame
df = samples.to_classification_df()
print(f"Classification data shape: {df.shape}")

# Create visualizations
samples.plot_pca(title="Sample Similarity")
samples.plot_bargraph(rank="phylum", top_n=10)
```

---

# Core API Reference

## Resource Management

### Samples - Sequencing Data
```python
# List and filter samples
samples = ocx.Samples.where(limit=50)

# Get specific sample
sample = ocx.Samples.get("sample_id")

# Access sample properties
print(f"Sample: {sample.filename}")
print(f"Upload date: {sample.created_at}")
print(f"File size: {sample.size}")
```

### Classifications - Analysis Results
```python
# List classifications
classifications = ocx.Classifications.where(complete=True, limit=10)

# Get classification results
classification = ocx.Classifications.get("classification_id")
results = classification.results()

# Access classification data
print(f"Total reads: {results['n_reads']}")
print(f"Taxa found: {len(results['table'])}")
```

### Analyses - All Analysis Types
```python
# List all completed analyses
analyses = ocx.Analyses.where(complete=True)

# Get analysis by ID (auto-detects type)
analysis = ocx.Analyses.get("analysis_id")

print(f"Analysis type: {analysis.analysis_type}")  # e.g. "classification", "functional"
```

### Projects - Sample Organization
```python
# List projects
projects = ocx.Projects.all()

# Get project samples
project_samples = ocx.Samples.where(project="project_id")
```

## File Operations

### Single File Upload
```python
# Basic upload
sample = ocx.Samples.create("path/to/sample.fastq.gz")

# Upload with metadata
sample = ocx.Samples.create(
    "sample.fastq.gz",
    metadata={
        "sample_type": "environmental",
        "date_collected": "2024-01-15",
        "location": "Urban soil"
    }
)
```

### Paired-End Upload
```python
# Upload paired reads
sample = ocx.Samples.create_paired(
    "sample_R1.fastq.gz",
    "sample_R2.fastq.gz",
    metadata={"library_type": "paired-end"}
)
```

### Batch Upload
```python
import os

# Upload multiple files
sample_files = ["sample1.fastq", "sample2.fastq", "sample3.fastq"]

samples = []
for filename in sample_files:
    if os.path.exists(filename):
        sample = ocx.Samples.create(filename)
        samples.append(sample)
        print(f"Uploaded: {sample.filename} (ID: {sample.id})")
```

## Data Access Patterns

### Sample Collections
```python
# Create collection from query
samples = ocx.Samples.where(project="project_id")  # Returns SampleCollection

# Collection operations
print(f"Collection size: {len(samples)}")
df = samples.to_classification_df()  # Combined DataFrame
taxonomy = samples.taxonomy  # Unified taxonomy

# Iterate through samples
for sample in samples:
    print(f"Sample: {sample.filename}")
```

### Resource Relationships

**Standard models:**

| Model | Field | Related model | Notes |
|-------|-------|---------------|-------|
| `Samples` | `.owner` | `Users` | |
| `Samples` | `.metadata` | `Metadata` | |
| `Samples` | `.primary_classification` | `Classifications` | optional |
| `Samples` | `.project` | `Projects` | optional |
| `Samples` | `.tags` | `[Tags]` | list |
| `Metadata` | `.sample` | `Samples` | |
| `Projects` | `.owner` | `Users` | |
| `Documents` | `.uploader` | `Users` | |
| `Documents` | `.downloaders` | `[Users]` | list |
| `Classifications` | `.sample` | `Samples` | |
| `Classifications` | `.job` | `Jobs` | |
| `Classifications` | `.dependencies` | `[Analyses]` | list |
| `Alignments` | `.sample` | `Samples` | |
| `Alignments` | `.job` | `Jobs` | |
| `FunctionalProfiles` | `.sample` | `Samples` | |
| `FunctionalProfiles` | `.job` | `Jobs` | |
| `Panels` | `.sample` | `Samples` | |
| `Panels` | `.job` | `Jobs` | |
| `Analyses` | `.sample` | `Samples` | |
| `Analyses` | `.job` | `Jobs` | |

**Experimental models** (`ocx = Api(experimental=True)`):

| Model | Field | Related model | Notes |
|-------|-------|---------------|-------|
| `Assets` | `.uploaded_by` | `Users` | |
| `Assemblies` | `.owner` | `Users` | |
| `Assemblies` | `.genome` | `Genomes` | optional |
| `Assemblies` | `.input_samples` | `[Samples]` | optional list |
| `Assemblies` | `.job` | `Jobs` | optional |
| `Assemblies` | `.primary_annotation_set` | `AnnotationSets` | optional |
| `AnnotationSets` | `.assembly` | `Assemblies` | |
| `AnnotationSets` | `.job` | `Jobs` | optional |
| `Genomes` | `.assemblies` | `[Assemblies]` | list |
| `Genomes` | `.primary_assembly` | `Assemblies` | optional |
| `Genomes` | `.tags` | `[Tags]` | list |
| `Genomes` | `.taxon` | `Taxa` | |
| `Taxa` | `.parent` | `Taxa` | optional, self-referential |

### Model Fields

**Key:** scalar types (`str`, `int`, `float`, `bool`, `datetime`) support their natural filter operators (see table above). `→ Model` = single-ref field (use `$eq`/`$ne`/`$in`). `→ [Model]` = list-ref field (use `$containsall`/`$containsany`). Fields marked † are not filterable.

**`Samples`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `updated_at` | datetime | optional |
| `filename` | str | optional |
| `size` | int | optional |
| `status` | str | |
| `visibility` | str | |
| `error_msg` | str | optional |
| `metadata` | → `Metadata` | |
| `owner` | → `Users` | |
| `primary_classification` | → `Classifications` | optional |
| `project` | → `Projects` | optional |
| `tags` | → `[Tags]` | |

**`Metadata`** (own samples only for filtering)

| Field | Type | Notes |
|-------|------|-------|
| `sample` | → `Samples` | |
| `starred` | bool | |
| `updated_at` | datetime | optional |
| `date_collected` | datetime | optional |
| `date_sequenced` | datetime | optional |
| `description` | str | optional |
| `external_sample_id` | str | optional |
| `library_type` | str | optional |
| `location_lat` | float | optional |
| `location_lon` | float | optional |
| `location_string` | str | optional |
| `name` | str | optional |
| `platform` | str | optional |
| `sample_type` | str | optional |
| `custom` | dict | † not filterable |

**`Classifications` / `Alignments` / `FunctionalProfiles` / `Panels`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `complete` | bool | |
| `draft` | bool | |
| `success` | bool | optional |
| `error_msg` | str | optional |
| `job` | → `Jobs` | |
| `sample` | → `Samples` | |
| `dependencies` | → `[Analyses]` | |
| `job_args` | dict | † not filterable |
| `cost` | object | † not filterable |

**`Analyses`** (adds one field to the above)

| Field | Type | Notes |
|-------|------|-------|
| `analysis_type` | str | |

**`Projects`**

| Field | Type | Notes |
|-------|------|-------|
| `name` | str | optional |
| `project_name` | str | optional |
| `description` | str | optional |
| `external_id` | str | optional |
| `public` | bool | |
| `owner` | → `Users` | |
| `permissions` | list | † not filterable |

**`Tags`**

| Field | Type |
|-------|------|
| `name` | str |

**`Users`**

| Field | Type |
|-------|------|
| `email` | str |

**`Jobs`**

| Field | Type |
|-------|------|
| `created_at` | datetime |
| `name` | str |
| `analysis_type` | str |
| `public` | bool |

**`Documents`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `filename` | str | |
| `size` | int | optional |
| `uploader` | → `Users` | |
| `downloaders` | → `[Users]` | |

---

Experimental models (`ocx = Api(experimental=True)`):

**`Assets`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `name` | str | |
| `filename` | str | |
| `status` | str | |
| `organization_id` | int | |
| `uuid` | str | |
| `uploaded_by` | → `Users` | |

**`Assemblies`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `filename` | str | optional |
| `size` | int | optional |
| `visibility` | str | |
| `genome` | → `Genomes` | optional |
| `input_samples` | → `[Samples]` | optional |
| `job` | → `Jobs` | optional |
| `owner` | → `Users` | |
| `primary_annotation_set` | → `AnnotationSets` | optional |

**`AnnotationSets`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `assembly` | → `Assemblies` | |
| `job` | → `Jobs` | optional |

**`Genomes`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `name` | str | optional |
| `description` | str | optional |
| `assemblies` | → `[Assemblies]` | |
| `primary_assembly` | → `Assemblies` | optional |
| `tags` | → `[Tags]` | |
| `taxon` | → `Taxa` | |

**`Taxa`**

| Field | Type | Notes |
|-------|------|-------|
| `created_at` | datetime | |
| `name` | str | optional |
| `taxon_id` | str | |
| `rank` | str | optional |
| `parent` | → `Taxa` | optional, self-referential |

---

# Data Analysis

## DataFrames & Data Structures

### Classification DataFrames
```python
from onecodex.lib.enums import Rank, Metric

# Basic DataFrame - species abundance
df = samples.to_classification_df()
print(f"Shape: {df.shape}")  # (n_samples, n_taxa)
print(f"Data type: abundance values")

# Specify rank and metric
df_genus = samples.to_classification_df(
    rank=Rank.Genus,
    metric=Metric.Readcount,
    top_n=20
)

# Long format for statistical analysis
df_long = samples.to_classification_df(table_format="long")
# Columns: classification_id, tax_id, abundance, tax_name
```

### Understanding DataFrame Structure
```python
# Get DataFrame
df = samples.to_classification_df(rank=Rank.Species, top_n=10)

# Structure:
# - Rows: samples (indexed by classification_id)
# - Columns: taxonomic IDs (e.g., '821', '853')
# - Values: abundance/count data
# - ocx_metadata: attached metadata DataFrame

# Access components
abundance_data = df.iloc[:, :-1]  # Just abundance columns
metadata = df.ocx_metadata         # Sample metadata
taxonomy_df = df.ocx_taxonomy     # Taxonomic information
```

### Working with Taxonomic Names
```python
# Add taxonomic names to columns
df = samples.to_classification_df(rank=Rank.Species, top_n=10)

# Method 1: Replace column names with names
tax_names = {}
for tax_id in df.columns:
    if tax_id.isdigit():
        try:
            name = samples.taxonomy.loc[tax_id, 'name']
            tax_names[tax_id] = f"{name} ({tax_id})"
        except KeyError:
            tax_names[tax_id] = f"Unknown ({tax_id})"

df_named = df.rename(columns=tax_names)

# Method 2: Create lookup table
taxonomy_lookup = samples.taxonomy.reset_index()[['tax_id', 'name']]
print(taxonomy_lookup.head())
```

### Sample Metadata Integration
```python
# Add sample metadata to DataFrame
df = samples.to_classification_df()

# Access metadata
metadata = df.ocx_metadata
print(f"Available metadata: {list(metadata.columns)}")

# Add specific metadata columns
if 'sample_type' in metadata.columns:
    df['sample_type'] = metadata['sample_type']
    df['collection_date'] = metadata['date_collected']
```

## Taxonomic Analysis

### Taxonomy Library Integration
```python
import json
import taxonomy

# Create taxonomy instance from One Codex data
def create_taxonomy_from_samples(samples):
    """Create a Taxonomy instance from One Codex sample collection."""

    # Build tax_id to index mapping
    tax_id_to_ix = {}
    taxonomy_df = samples.taxonomy.reset_index()

    for ix, row in enumerate(taxonomy_df.to_dict(orient="records")):
        tax_id_to_ix[row["tax_id"]] = ix

    # Build taxonomy data structure
    tax_data = {"nodes": [], "links": []}

    for ix, row in enumerate(taxonomy_df.to_dict(orient="records")):
        tax_data["nodes"].append({
            "id": int(row["tax_id"]),
            "name": row["name"],
            "rank": row["rank"]
        })

        tax_data["links"].append({
            "source": tax_id_to_ix[row["tax_id"]],
            "target": (
                tax_id_to_ix[row["parent_tax_id"]]
                if row["parent_tax_id"] and row["parent_tax_id"] in tax_id_to_ix
                else tax_id_to_ix["1"]  # Root node
            ),
        })

    # Create Taxonomy instance
    tax = taxonomy.Taxonomy.from_json(json.dumps(tax_data))
    return tax

# Create taxonomy instance for your sample collection
tax = create_taxonomy_from_samples(samples)
print(f"Created taxonomy with {len(tax)} nodes")
```

### Taxonomic Filtering

```python
def filter_taxa_by_parent(samples, tax, parent_name, target_rank="species"):
    """
    Universal function to filter taxa at or below a specific parent.

    Args:
        samples: One Codex SampleCollection
        tax: Taxonomy instance
        parent_tax_id: Tax ID of the parent node (e.g., "10239" for viruses)
        target_rank: Taxonomic rank to filter to (default: "species")

    Returns:
        List of tax_ids for taxa matching criteria
    """

    # Find the parent node
    parent_node = None
    parent_nodes = tax.find_all_by_name(parent_name)
    if not len(parent_nodes) == 1:
        raise Exception(f"Missing or ambiguous name: {parent_name}")

    parent_node = parent_nodes[0]

    print(f"Found parent: {parent_node.name} (ID: {parent_node.tax_id})")

    # Get all descendants at or below this parent
    descendants = tax.descendants(parent_node.tax_id)

    # Filter to target rank (if specified)
    if target_rank:
        target_taxa = [node for node in descendants if node.rank == target_rank]
        print(f"Found {len(target_taxa)} {target_rank} taxa in {parent_node.name}")
    else:
        target_taxa = descendants
        print(f"Found {len(target_taxa)} total taxa in {parent_node.name}")

    return [node.tax_id for node in target_taxa]

# Filter examples
viral_species = filter_taxa_by_parent(samples, tax, 'Viruses')
bacterial_species = filter_taxa_by_parent(samples, tax, 'Bacteria')
fungal_species = filter_taxa_by_parent(samples, tax, 'Fungi')
```

## Statistical Operations

### Alpha Diversity
```python
# Built-in alpha diversity calculation
alpha_diversity = samples.alpha_diversity(rank=Rank.Species)
print(f"Alpha diversity metrics: {alpha_diversity.columns.tolist()}")

# diversity_metric options (AlphaDiversityMetric):
# - AlphaDiversityMetric.Shannon (default)
# - AlphaDiversityMetric.Simpson
# - AlphaDiversityMetric.ObservedTaxa
```

### Beta Diversity
```python
# Calculate beta diversity distance matrix
from onecodex.lib.enums import Metric

distance_matrix = samples.beta_diversity(
    rank=Rank.Species,
    diversity_metric=BetaDiversityMetric.BrayCurtis
)
print(f"Distance matrix shape: {distance_matrix.shape}")
```

---

# Visualization

## Built-in Plotting Functions

### PCA Plots
```python
# Basic PCA
chart = samples.plot_pca(
    rank=Rank.Species,
    title="Sample Similarity Analysis"
)

# Customized PCA
chart = samples.plot_pca(
    rank=Rank.Genus,
    metric=Metric.BrayCurtis,
    title="Genus-level PCA",
    return_chart=True  # Returns Altair chart object
)
chart.show()
```

### Taxonomic Composition
```python
# Bar chart of taxonomic composition
chart = samples.plot_bargraph(
    rank=Rank.Phylum,
    top_n=10,
    title="Phylum-level Composition"
)

# Customized bargraph
chart = samples.plot_bargraph(
    rank=Rank.Family,
    top_n=15,
    metric=Metric.ReadCount,
    title="Family Abundance",
    return_chart=True
)
```

### Abundance Heatmaps
```python
# Heatmap of abundant taxa
chart = samples.plot_heatmap(
    rank=Rank.Genus,
    top_n=20,
    title="Genus Abundance Heatmap"
)

# Heatmap with more taxa
chart = samples.plot_heatmap(
    rank=Rank.Species,
    top_n=25,
    title="Species Abundance Heatmap"
)
```

### Distance Matrices
```python
# Beta diversity visualization
chart = samples.plot_distance(
    rank=Rank.Species,
    diversity_metric=BetaDiversityMetric.UnweightedUnifrac,
    title="UniFrac Distance Matrix"
)

# Bray-Curtis distances
chart = samples.plot_distance(
    rank=Rank.Genus,
    diversity_metric=BetaDiversityMetric.BrayCurtis,
    title="Bray-Curtis Distances"
)
```

## Chart Customization

### Working with Altair Charts
```python
import altair as alt

# Get chart object
chart = samples.plot_pca(return_chart=True)

# Customize chart
customized_chart = chart.properties(
    width=600,
    height=400
).configure_axis(
    labelFontSize=12,
    titleFontSize=14
).configure_title(
    fontSize=16
)

customized_chart.show()
```

### Color Schemes and Themes
```python
# Custom color palette
chart = samples.plot_bargraph(
    rank=Rank.Phylum,
    top_n=8,
    return_chart=True
)

# Apply custom colors
chart = chart.configure_range(
    category=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
)
```

## Export Options

### Saving Charts
```python
# Save as PNG
chart = samples.plot_pca(return_chart=True)
chart.save('pca_analysis.png', scale_factor=2.0)

# Save as SVG
chart.save('pca_analysis.svg')

# Save as HTML
chart.save('pca_analysis.html')
```

---

# Advanced Features

## Experimental Models

The One Codex API provides access to experimental models and features through the `X-OneCodex-Api-Experimental` header. These features are subject to change without notice and should not be relied upon in a production environment.

### Enabling Experimental API Access

```python
import warnings
from onecodex import Api

# Enable experimental API access
ocx = Api(experimental=True)

# This will show a warning about experimental features
```

**Warning:** Experimental API mode enables access to features that are subject to change without notice and should not be relied upon in production environments.

### Available Experimental Models

#### Assets Model
The Assets model provides file storage and management capabilities:

```python
# Upload a file as an asset
asset = ocx.Assets.upload("path/to/file.txt", name="My Data File")
print(f"Uploaded asset: {asset.id}")

# List assets
assets = ocx.Assets.all()

# Update asset metadata
asset.name = "Updated File Name"
asset.save()

# Download an asset
asset.download("downloaded_file.txt")

# Delete an asset
asset.delete()
```

**Asset Properties:**
- `id`: Unique asset identifier
- `name`: User-defined asset name
- `filename`: Original filename
- `s3_uri`: Internal storage URI
- `status`: Asset processing status
- `organization_id`: Associated organization
- `uploaded_by`: User who uploaded the asset
- `uuid`: Asset UUID

#### Genome Models
Access genome-related data including assemblies, annotations, and taxonomic information:

```python
# Browse available genomes
genomes = ocx.Genomes.all()

# Get a specific genome
genome = ocx.Genomes.get("genome_id")
print(f"Genome: {genome.name} ({genome.taxon.name})")

# Access genome assemblies
assemblies = genome.assemblies
primary_assembly = genome.primary_assembly

# Download assembly in FASTA format
if assemblies:
    assembly = assemblies[0]
    assembly.download("assembly.fasta")
```

#### Assembly Models
Work with genomic assemblies:

```python
# List assemblies
assemblies = ocx.Assemblies.all()

# Get assembly details
assembly = ocx.Assemblies.get("assembly_id")
print(f"Assembly size: {assembly.size} bytes")

# Download assembly
assembly.download("my_assembly.fasta")

# Access associated genome and annotations
if assembly.genome:
    genome = assembly.genome
    annotation_set = assembly.primary_annotation_set
```

#### Annotation Set Models
Access genomic annotations:

```python
# List annotation sets
annotation_sets = ocx.AnnotationSets.all()

# Download annotations in GenBank format
annotation_set = ocx.AnnotationSets.get("annotation_set_id")
annotation_set.download("annotations.gbk")

# Download annotations as CSV
annotation_set.download_csv("annotations.csv")
```

#### Taxa Models
Browse taxonomic information:

```python
# Search for taxa
taxa = ocx.Taxa.all()

# Get specific taxon (by Tax ID)
taxon = ocx.Taxa.where(taxon_id='821')[0]
print(f"Taxon: {taxon.name} ({taxon.rank})")

# Get all genomes for a taxon and its descendants
genomes = taxon.genomes()

# Get taxonomic hierarchy
parents = taxon.parents()
```

---

## CLI Tools

### Basic Commands
```bash
# List samples
onecodex samples

# Upload files
onecodex upload sample.fastq.gz
onecodex upload --metadata '{"sample_type": "clinical"}' *.fastq

# Check analysis status
onecodex analyses --complete
```

### Advanced CLI Usage
```bash
# Query with filters
onecodex samples --filter '{"metadata.sample_type": "environmental"}'

# Download results
onecodex download --sample-id SAMPLE_ID --format csv

# Batch operations
onecodex upload --directory /path/to/samples/ --recursive
```

---

*This documentation provides comprehensive coverage of the One Codex Python
client library. For the most current information, refer to the [official
documentation](https://github.com/onecodex/onecodex) and repository.*
