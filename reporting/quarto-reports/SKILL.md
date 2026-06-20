---
name: bio-reporting-quarto-reports
description: Builds reproducible Quarto reports, presentations, and websites across R, Python, and Julia, with correct engine selection, cache-vs-freeze semantics, native cross-references, parameters, and environment pinning. Use when creating a Quarto report of an analysis, setting up freeze for CI, or debugging cross-references, caching, or working-directory issues.
tool_type: mixed
primary_tool: Quarto
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: Quarto 1.4+, knitr 1.45+, pandoc 3.1+ (bundled), scanpy 1.10+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `quarto --version`, `quarto check`, `quarto render --help`

Some flags and project keys move between Quarto releases (e.g. the file-based `--execute-params`); confirm against `quarto render --help`. If a render fails, run `quarto check` and adapt to the installed version rather than retrying.

# Quarto Reports

**"Create a Quarto analysis report"** -> Write a document mixing code (R/Python/Julia), narrative, and figures that executes through a computational engine and renders to HTML/PDF/Word.
- CLI: `quarto render report.qmd --to html`

## The Pipeline: Engine, Then Pandoc

Both Quarto and R Markdown end at pandoc; what differs is what runs before it. Quarto first picks a computational ENGINE, then pandoc converts to the target format. The engine is a property of the document's languages, and it determines what runtime the rendering machine needs:

- Any `{r}` chunk present -> **knitr** engine (same knit -> md -> pandoc path as R Markdown).
- Only `{python}`/`{julia}` chunks -> **jupyter** engine (executes via a Jupyter kernel, then pandoc).
- **Both R and Python -> knitr + reticulate in ONE process**, so R and Python share a session and can pass objects back and forth. (This is why mixed-language docs "just work" through knitr, not jupyter.)
- Override in YAML: `engine: knitr` / `engine: jupyter`, or pin a kernel with `jupyter: python3`.

The consequence: a Python-only `.qmd` on the jupyter engine needs a registered Jupyter kernel; switched to knitr+reticulate it needs R+reticulate instead. Freeze (below) lets CI skip needing either.

## Cache vs Freeze (the load-bearing distinction)

These solve DIFFERENT problems and are constantly conflated:

> knitr **cache** makes a SINGLE render faster by skipping unchanged chunks. Quarto **freeze** lets a DIFFERENT machine (CI / a website build) render with NO language runtime installed, by reusing stored results.

| | cache (`execute: cache`) | freeze (`execute: freeze`) |
|---|---|---|
| Granularity | per-chunk (knitr) / per-notebook (jupyter-cache) | per-document |
| Problem solved | skip unchanged chunks during a render | skip ALL execution on publish/CI |
| Key | MD5(code + evaluating options); data only via `cache.extra` | source-file hash (`auto`) or never re-run (`true`) |
| Lives in | `*_cache/` (per-doc) | `_freeze/` (project - **commit it**) |
| Runtime needed to render? | yes (still renders, skips some chunks) | **no** - CI renders with no R/Python |
| Invalidates on upstream DATA change? | NO unless `cache.extra` | only via source change (`auto`); data not auto-tracked |
| Scope | within one render | only FULL project renders |

Two edges that trip everyone:
- **The stale-cache footgun:** `cache=TRUE` keys on chunk CODE, not the data it reads. If `data.csv` changes but the chunk code is byte-identical, the cached (stale) result is served. Bind the data into the key: `cache.extra = tools::md5sum('data.csv')`. Cross-chunk dependencies need `dependson='chunkA'` (or `autodep=TRUE`, best-effort).
- **Freeze only acts on FULL project renders.** `quarto render onefile.qmd` and `quarto render subdir/` always execute, ignoring `freeze:`. Arrange CI to do a whole-project `quarto render` so frozen results are honored. Commit `_freeze/` so others render without reproducing the environment.
- **They compose, not conflict.** Freeze decides whether the project re-executes at all; when it does (source changed), knitr cache still skips unchanged chunks within that run.

## The Working-Directory Trap

Chunks execute with the working directory set to the document's folder, NOT the project root (default `execute-dir: file`). So `pd.read_csv('data/x.csv')` works interactively from the project root but breaks on render when the `.qmd` lives in `reports/`. Set `project: execute-dir: project` in `_quarto.yml` to run all chunks from the project root, or use root-anchored paths (`here::here(...)` in R). Never `setwd()` in a chunk - it desyncs figure/cache file placement.

## Cross-References Need the Type Prefix

A Quarto label is a cross-reference ONLY if it starts with a reserved lower-case type prefix: `fig-`, `tbl-`, `sec-`, `eq-`, `lst-`, theorem/callout families. `#| label: scatter` is a dead anchor; `#| label: fig-scatter` is referenceable as `@fig-scatter`. This is the #1 cause of a reference rendering as `?@fig-x`.

````markdown
```{python}
#| label: fig-umap
#| fig-cap: "UMAP embedding colored by cluster"
sc.pl.umap(adata, color='leiden')
```
See @fig-umap. Methods are in @sec-methods.
````

A figure/table from a code cell needs both the prefixed `label` and a `fig-cap`/`tbl-cap`. Section refs need `{#sec-methods}` on the heading AND `number-sections: true`. (Base R Markdown cannot cross-reference at all - that requires bookdown; see reporting/rmarkdown-reports.)

## Parameters: knitr vs jupyter Differ

- **knitr engine:** YAML `params:` block, accessed read-only as `params$x`. Override: `quarto render doc.qmd -P alpha:0.2`.
- **jupyter engine:** there is NO `params:` block. Designate a cell tagged `parameters` (papermill convention) with default assignments; variables are then top-level names. A `params:` YAML block on a jupyter-engine document is silently ignored - a common bug.

````markdown
```{python}
#| tags: [parameters]
input_file = "adata.h5ad"
n_top_genes = 2000
```
````

`-P key:val` overrides on the CLI for both engines.

## Document Basics and Layout

```yaml
---
title: "Analysis Report"
date: today
format:
  html:
    toc: true
    code-fold: true
    embed-resources: true   # one portable self-contained HTML
execute:
  warning: false
  freeze: auto
---
```

Per-cell options use the `#|` hash-pipe (`#| echo: false`, `#| fig-width: 8`, `#| cache: true`). Tabsets group alternative views under `::: {.panel-tabset}`; callouts (`::: {.callout-note}`) flag notes/warnings/tips. Render multiple formats by listing them under `format:` and `quarto render` (or `--to pdf`); PDF needs a TeX engine (`quarto install tinytex`).

## Self-Contained Output

`embed-resources: true` base64-inlines images, CSS, and JS into one portable HTML (maps to pandoc `--embed-resources --standalone`; the older `--self-contained` is deprecated since pandoc 2.19). htmlwidgets (plotly, DT) get inlined too, so an interactive report is one openable file - but each widget library inflates the size.

## The Document Captures Code, Not the Environment

Quarto does not pin package versions or the interpreter. A `.qmd` that renders perfectly today can silently change output next year when a dependency updates. The document gives byte-reproducible output only if code, data, AND versions are unchanged - and versions are not in the repo unless pinned. For real reproducibility add a lockfile/container: `renv::snapshot()` (`renv.lock`) for R, `environment.yml`/`requirements.txt` for Python, Docker/Apptainer when the OS, TeX, and pandoc must also be pinned. **Freeze is not reproducibility** - `_freeze/` lets CI skip execution, but the frozen results came from an uncaptured environment. Record provenance with `sessionInfo()` / `sessioninfo::session_info()` (provenance, not a restore mechanism). For journal submission, Quarto manuscript/journal templates (`quarto-journals/...`) produce article-formatted output from the same source.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `@fig-x` renders as `?@fig-x` | label missing the type prefix | name it `fig-x`/`tbl-x` and give it a caption |
| `params:` ignored on a Python doc | jupyter engine uses a `parameters`-tagged cell, not `params:` | tag a cell `parameters`, or use the knitr engine |
| Stale results after editing data | `cache` keys on code, not data | `cache.extra = tools::md5sum('data.csv')` |
| CI re-runs everything despite freeze | single-file/subdir render ignores freeze | do a full-project `quarto render`; commit `_freeze/` |
| `read_csv('data/..')` fails on render | working dir = doc folder, not project root | `execute-dir: project` or `here::here()` |
| Report reproduces differently months later | environment not pinned | renv.lock / conda env / container |
| PDF render fails | no TeX engine | `quarto install tinytex` |

## Related Skills

- reporting/rmarkdown-reports - R-focused alternative; needs bookdown for cross-references
- reporting/jupyter-reports - Parameterized notebook execution Quarto can consume
- reporting/publication-tables - Formatted tables to embed in the report
- data-visualization/ggplot2-fundamentals - Figures for R-engine reports
- data-visualization/interactive-visualization - Interactive dashboards (Quarto `format: dashboard` for static/self-contained, Shiny when a running server is acceptable)

## References

- Knuth DE. Literate Programming. Comput J. 1984;27(2):97-111. doi:10.1093/comjnl/27.2.97
- Xie Y, Allaire JJ, Grolemund G. R Markdown: The Definitive Guide. Chapman & Hall/CRC; 2018
- Xie Y, Dervieux C, Riederer E. R Markdown Cookbook. Chapman & Hall/CRC; 2020
- Quarto documentation: quarto.org (cache/freeze, cross-references, parameters, execution engine)
