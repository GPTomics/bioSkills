# R Markdown Reports - Usage Guide

## Overview

R Markdown combines R code, results, and narrative into reproducible documents (HTML, PDF, Word, slides) via knitr then pandoc. The decisions that actually matter are making the document self-sufficient so it knits the same way a colleague's fresh session does, invalidating the cache when the data changes, resolving relative paths from the right directory, and reaching for bookdown when cross-references are needed.

## Prerequisites

```r
install.packages(c('rmarkdown', 'knitr', 'bookdown', 'DT', 'kableExtra', 'here'))
# For PDF output:
tinytex::install_tinytex()
```

## Quick Start

Tell your AI agent what you want to do:
- "Create an R Markdown template for my DESeq2 RNA-seq analysis"
- "My report knits differently than it runs interactively - make it self-sufficient"
- "My cached chunk won't update after I changed the counts file"
- "Add cross-referenced figures and tables to my report"
- "Generate a parameterized report for each sample"

## Example Prompts

### Basic Reports

> "Create an R Markdown document for my DESeq2 results with code folding and a TOC"

> "Generate a PDF report summarizing my variant-calling pipeline"

### Reproducibility and Caching

> "Make this report render in a fresh session the way the Knit button does"

> "Cache the long DESeq2 step but invalidate it when the count matrix changes"

### Cross-References and Parameters

> "Add numbered, cross-referenced figures using bookdown"

> "Set up a parameterized template that accepts a sample ID and FDR threshold"

### Tables

> "Add an interactive DT table for exploration and a static kableExtra table for the PDF"

## What the Agent Will Do

1. Structure the `.Rmd` with a setup chunk, YAML output format, and self-sufficient chunks
2. Bind expensive cached chunks to their input files with `cache.extra`
3. Resolve paths with `here::here()` so the report knits from any directory
4. Switch to a bookdown `*_document2` format when cross-references are requested
5. Parameterize with `params:` and pair the report with renv for environment pinning

## Document Types

| Output | Best for |
|--------|----------|
| `html_document` / `bookdown::html_document2` | interactive reports; the `*2` form adds cross-refs |
| `pdf_document` / `bookdown::pdf_document2` | publication, archiving; needs LaTeX |
| `word_document` | collaborator editing |
| `ioslides_presentation` / `xaringan` | slides |

## Tips

- `render()` from the console sees the caller's global objects; the Knit button uses a fresh session. Make every object chunk-created, and test with `envir=new.env()`.
- `cache=TRUE` keys on chunk code, not data; add `cache.extra=tools::md5sum('file')` so a data edit invalidates it.
- knit working directory is the `.Rmd` folder, not the project root; use `here::here()` for relative paths.
- Base rmarkdown cannot cross-reference; use bookdown `*_document2` with a labeled, captioned chunk.
- DT tables are interactive HTML widgets, not for print; use kableExtra for static publication tables.
- The document does not pin packages - add renv.lock (and a container for full reproducibility); end with `sessionInfo()`.

## Related Skills

- reporting/quarto-reports - Successor with native cross-references and multi-language support
- reporting/publication-tables - Formatted static tables for reports
- data-visualization/ggplot2-fundamentals - Creating publication-quality figures
