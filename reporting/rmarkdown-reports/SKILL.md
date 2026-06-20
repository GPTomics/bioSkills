---
name: bio-reporting-rmarkdown-reports
description: Creates reproducible R Markdown analysis reports (HTML, PDF, Word) with knitr, covering the render pipeline, the interactive-vs-knit session trap, cache invalidation, bookdown cross-references, parameterization, and environment pinning. Use when generating an R-based analysis report, debugging a report that knits differently than it runs interactively, or fixing caching or cross-references.
tool_type: r
primary_tool: rmarkdown
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: rmarkdown 2.25+, knitr 1.45+, bookdown 0.37+, DESeq2 1.42+, ggplot2 3.5+, DT 0.31+, kableExtra 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package (`?rmarkdown::render`, `?knitr::opts_chunk`) and adapt the example to the actual API rather than retrying.

# R Markdown Reports

**"Create an R Markdown report"** -> Write an R-centric document combining code chunks, results, and narrative that knits to HTML/PDF/Word.
- R: `rmarkdown::render('report.Rmd')`, or the Knit button in RStudio

## The Pipeline: knitr, Then Pandoc

An `.Rmd` always renders in two stages: **knitr** executes the chunks and weaves the results into an intermediate `.md`, then **pandoc** converts that `.md` into the target format (LaTeX via a TeX engine for PDF). knitr is the only execution engine for `.Rmd` (other languages run only as knitr engines); it is the successor to Sweave, adding caching, hooks, and markdown hosting. `rmarkdown::render()` orchestrates both stages. Knowing the split explains most failures: a chunk error is knitr; a formatting or cross-reference problem is usually pandoc/bookdown.

## The "Works Interactively, Fails on Knit" Trap

This is the single most common reproducibility surprise. `rmarkdown::render()` defaults to `envir = parent.frame()`, so calling it from the console evaluates chunks in the caller's environment - it can SEE objects sitting in the interactive global env. The RStudio **Knit button does NOT**: it spawns a fresh, clean R session. So a report that relies on a `df` created interactively renders fine via `render()` from the console, then fails when a colleague clicks Knit or CI runs it, because the fresh session has no `df`.

Guards:
- Treat the document as self-sufficient - every object must be CREATED in a chunk, never assumed present.
- To mimic the button before trusting a report, render in isolation: `rmarkdown::render('r.Rmd', envir = new.env())`, or in a fresh process via `callr::r(...)` / `xfun::Rscript_call(rmarkdown::render, ...)`.

## Caching: cache Keys on Code, Not Data

`cache=TRUE` stores a chunk's result in a `*_cache/` dir and reloads it on re-knit if the chunk is "unchanged" - where the cache key is an MD5 of the chunk CODE plus evaluating options. The footgun: if a chunk reads `data.csv` and the FILE changes but the chunk code is byte-identical, the hash is unchanged and knitr serves the STALE cached result. Bind the data into the key:

````r
```{r de-analysis, cache=TRUE, cache.extra=tools::md5sum('counts.csv')}
dds <- DESeq(DESeqDataSetFromMatrix(counts, metadata, ~ condition))
```
````

Cross-chunk dependencies are not tracked automatically either: if chunk B uses an object from chunk A, editing A does not invalidate B's cache by default - declare `dependson='de-analysis'` (or `autodep=TRUE`, best-effort).

## The Working-Directory Trap

knitr evaluates chunks with the working directory set to the directory of the `.Rmd`, NOT the project root. So `read.csv('data/x.csv')` works when run interactively from the project root but breaks on knit if the `.Rmd` lives in `reports/`. Fixes, in order of preference: `here::here('data/x.csv')` (anchors to the project root, most robust); `knitr::opts_knit$set(root.dir = '...')` in the setup chunk (note `opts_knit`, not `opts_chunk`); or `rmarkdown::render('r.Rmd', knit_root_dir = '...')`. Never `setwd()` in a chunk - it desyncs figure/cache file placement.

## Cross-References Require bookdown

Base rmarkdown CANNOT cross-reference figures, tables, sections, or equations. Use a bookdown output format - `bookdown::html_document2`, `bookdown::pdf_document2`, `bookdown::word_document2` - which add numbering and `\@ref(type:label)`. Two hard requirements: the figure/table chunk must be LABELED, and it must have a CAPTION (`fig.cap=`); a captionless figure is emitted unnumbered and cannot be referenced.

````yaml
output:
  bookdown::html_document2:
    toc: true
````
````r
```{r volcano, fig.cap="Volcano plot of differential expression"}
plot(res$log2FoldChange, -log10(res$pvalue))
```
See Figure \@ref(fig:volcano).
````

(Quarto has native cross-references without bookdown - see reporting/quarto-reports.)

## Parameterized Reports

Declare defaults in YAML and read them as a read-only list:

````yaml
params:
  count_file: "counts.csv"
  fdr_threshold: 0.05
````
````r
counts <- read.csv(params$count_file)
```
````

Override per render and loop over samples:

```r
rmarkdown::render('report.Rmd', params = list(count_file = 'sampleB.csv'),
                  output_file = 'sampleB_report.html')
```

`rmarkdown::render(..., params = 'ask')` launches the "Knit with Parameters" UI.

## Document Basics, Tables, and Output

```yaml
---
title: "RNA-seq Report"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    code_folding: hide
    self_contained: true   # base64-embed assets into one portable HTML
---
```

A `setup` chunk with `knitr::opts_chunk$set(echo=TRUE, message=FALSE, warning=FALSE, fig.width=10)` sets document-wide defaults. Section tabs use `## Results {.tabset}`. Inline results splice with `` `r ...` ``. For tables: `knitr::kable()` + `kableExtra` for STATIC publication tables; `DT::datatable()` for INTERACTIVE HTML exploration - DT is a JavaScript widget, not for print/PDF, and it inflates the HTML (see reporting/publication-tables for the formatted-table decision). `self_contained: true` (default for `html_document`) embeds all assets into one portable file at a size cost; htmlwidgets get inlined too.

## The Document Captures Code, Not the Environment

rmarkdown does not pin package versions or R itself. A report that knits perfectly today can change output next year when a dependency updates. The document gives byte-reproducible output only if code, data, AND versions are unchanged - and versions are not in the repo unless pinned. Add `renv::snapshot()` (`renv.lock`, commit it) for package pinning, and a container (Docker/Apptainer) when the OS, TeX, and pandoc must also be fixed. End the report with `sessionInfo()` / `sessioninfo::session_info()` - provenance for the reader, not a restore mechanism. Seed any stochastic step (`set.seed`).

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Renders from console, fails on Knit | `render` sees globals (`parent.frame`); Knit uses a fresh session | make every object chunk-created; test with `envir=new.env()` |
| Stale results after editing data | `cache` keys on code, not data | `cache.extra=tools::md5sum('data.csv')` |
| `read.csv('data/..')` fails on knit | working dir = `.Rmd` folder, not project root | `here::here()` or `knit_root_dir=` |
| `\@ref(fig:x)` shows as `??` | base rmarkdown can't cross-ref, or no caption/label | bookdown `*_document2` + chunk label + `fig.cap` |
| Edited upstream chunk, downstream cache stale | dependencies not tracked | `dependson=` or `autodep=TRUE` |
| Report changes output months later | environment not pinned | renv.lock + container; seed RNGs |
| PDF knit fails | no LaTeX | `tinytex::install_tinytex()` |

## Related Skills

- reporting/quarto-reports - Successor with native cross-references and multi-language support
- reporting/publication-tables - Formatted static tables (gt/gtsummary/flextable) for reports
- reporting/figure-export - Exporting the report's figures for publication
- differential-expression/de-results - The analysis these reports typically present

## References

- Xie Y. Dynamic Documents with R and knitr. 2nd ed. Chapman & Hall/CRC; 2015
- Xie Y. bookdown: Authoring Books and Technical Documents with R Markdown. Chapman & Hall/CRC; 2016
- Xie Y, Allaire JJ, Grolemund G. R Markdown: The Definitive Guide. Chapman & Hall/CRC; 2018
- Xie Y, Dervieux C, Riederer E. R Markdown Cookbook. Chapman & Hall/CRC; 2020
