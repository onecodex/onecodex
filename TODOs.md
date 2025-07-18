# Type Checking TODOs for One Codex

## Overview
This document tracks the progress of adding type annotations to the One Codex library using the `ty` type checker.

## Initial Status
- Total diagnostics found: 101
- Python version target: 3.10+
- Excluded directories: tests, notebook_examples, docs, onecodex/vendored

## Type Error Categories

### 1. Import Issues (High Priority)
- [ ] Fix `requests.packages.urllib3.util.retry` import (onecodex/api.py:9)
- [ ] Fix `sentry_sdk._compat.string_types` import (onecodex/utils.py:278)

### 2. Enum Type Issues (High Priority)
- [ ] Fix `FunctionalAnnotations` enum usage (onecodex/analyses.py:346)
- [ ] Fix `FunctionalAnnotationsMetric` enum usage (onecodex/analyses.py:348)
- [ ] Fix `Rank` enum usage (onecodex/analyses.py:389)
- [ ] Fix `Metric` enum attribute access (onecodex/analyses.py:501-502)

### 3. Exception Handling Issues
- [ ] Fix `Exception.errno` attribute access (onecodex/auth.py:84, 132, 151, 153)
- [ ] Fix exception type annotations throughout

### 4. Missing Type Annotations by Module

#### onecodex/analyses.py
- [ ] Add type annotations for class methods and attributes
- [ ] Fix enum parameter defaults

#### onecodex/api.py
- [ ] Add type annotations for API client methods
- [ ] Fix import issues

#### onecodex/auth.py
- [ ] Add proper exception type annotations
- [ ] Add type annotations for auth functions

#### onecodex/cli.py
- [ ] Add type annotations for CLI commands
- [ ] Fix click decorator type issues

#### onecodex/dataframes.py
- [ ] Add pandas DataFrame type annotations
- [ ] Add type annotations for data manipulation functions

#### onecodex/models/
- [ ] Add type annotations for all model classes
- [ ] Ensure Pydantic v2 models are properly typed

#### onecodex/viz/
- [ ] Add type annotations for visualization functions
- [ ] Fix altair chart type issues

### 5. Possibly Unresolved References
- [ ] Fix `chart` variable definition flow (onecodex/viz/_metadata.py:282, 284)
- [ ] Fix `idx` variable definition flow (onecodex/viz/_primitives.py:138)

### 6. Library-Specific Issues
- [ ] Fix concurrent.futures.thread access (onecodex/utils.py:400)
- [ ] Fix ProgressBar._update attribute (onecodex/utils.py:406)
- [ ] Fix requests.exceptions access patterns

## Progress Tracking
- [ ] Phase 1: Fix all import errors
- [ ] Phase 2: Fix enum type issues
- [ ] Phase 3: Add basic type annotations to all public APIs
- [ ] Phase 4: Fix all possibly unresolved references
- [ ] Phase 5: Add comprehensive internal type annotations
- [ ] Phase 6: Achieve 0 type errors

## Notes
- Do not change code behavior, only add type annotations
- Do not add broad type ignores without discussion
- Focus on public API typing first
- Maintain compatibility with Python 3.10+
