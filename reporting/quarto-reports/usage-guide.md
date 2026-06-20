# Quarto Reports - Usage Guide

## Overview

Quarto is a language-agnostic scientific publishing system (R, Python, Julia, Observable) and the successor to R Markdown. The decisions that actually matter when using it are which computational engine runs the document, the difference between caching and freezing computations, why relative paths break on render, how native cross-references must be labeled, and that the document captures code but not the environment that produced it.

## Prerequisites

```bash
# Install Quarto - download from https://quarto.org/docs/download/
quarto check

# For Python documents (jupyter engine)
pip install jupyter matplotlib

# For R documents (knitr engine)
install.packages(c('knitr', 'rmarkdown'))

# PDF output
quarto install tinytex
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a Quarto report of my scRNA-seq analysis in Python"
- "Set up freeze so CI can render my Quarto website without R installed"
- "My @fig-x cross-reference renders as ?@fig-x - fix the labels"
- "My report serves stale results after I updated the data"
- "Render to both HTML and PDF from one source"

## Example Prompts

### Basic Documents

> "Create a Quarto document for my Python analysis with code folding and a table of contents"

> "Convert my Jupyter notebook into a Quarto report"

### Engine, Cache, and Freeze

> "Set up freeze: auto and commit _freeze/ so the GitHub Action renders without a Python kernel"

> "My cached chunk doesn't update when the input CSV changes - bind the cache to the data"

### Cross-References and Parameters

> "Add a cross-referenced figure and table to my report"

> "Parameterize my Python report by sample using a parameters cell"

### Reproducibility

> "Pin the environment with renv so this Quarto report reproduces next year"

## What the Agent Will Do

1. Choose or confirm the engine (knitr if any R chunk, jupyter for Python-only, knitr+reticulate for mixed)
2. Set `execute-dir: project` (or use `here::here()`) so relative paths resolve from the project root
3. Label figures/tables with the required type prefix (`fig-`/`tbl-`) and captions for cross-references
4. Parameterize correctly for the engine (`params:` for knitr, a `parameters`-tagged cell for jupyter)
5. Configure `freeze: auto` for CI and pair the report with a lockfile/container for real reproducibility

## Quarto vs R Markdown

| Feature | R Markdown | Quarto |
|---------|------------|--------|
| Languages | R-centric (others as knitr engines) | R, Python, Julia, Observable |
| Cross-references | needs bookdown `*_document2` | native (`@fig-`/`@tbl-`) |
| Caching | knitr `cache` | knitr `cache` + project `freeze` |
| Websites / books / slides | blogdown / bookdown / xaringan | built-in |

## Tips

- The engine is chosen by the document's languages; mixed R+Python runs on knitr+reticulate in one process, not jupyter.
- Cache speeds up the local render; freeze lets someone else's CI render with no runtime. Neither invalidates on a raw-data change under identical code.
- `cache.extra = tools::md5sum('data.csv')` is the fix for stale cache after a data edit.
- Freeze only acts on full-project renders; a single-file `quarto render file.qmd` always executes.
- A label is a cross-reference only with its type prefix (`fig-`/`tbl-`/`sec-`), all lower-case, plus a caption.
- The document does not pin packages - add renv.lock / conda / Docker for actual reproducibility.

## Related Skills

- reporting/rmarkdown-reports - R Markdown for R-only workflows
- reporting/jupyter-reports - Parameterized notebooks Quarto can consume
- data-visualization/ggplot2-fundamentals - Creating visualizations for reports
