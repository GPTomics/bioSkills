# reporting

## Overview

Reproducible report generation for bioinformatics analyses: literate-programming documents, QC aggregation, and publication-ready figures and tables. The cross-cutting principle is that a report, figure, or table is only as reproducible as the environment that produced it - the document captures code, not the package versions, so pin the environment (lockfile/container) and seed stochastic steps alongside.

**Tool type:** mixed | **Primary tools:** RMarkdown, Quarto, Jupyter/papermill, MultiQC, matplotlib, gtsummary

## Skills

| Skill | Description |
|-------|-------------|
| rmarkdown-reports | Reproducible R Markdown reports: render pipeline, the interactive-vs-knit trap, cache invalidation, bookdown cross-references |
| quarto-reports | Multi-language Quarto reports: engine selection, cache vs freeze, native cross-references, parameters |
| jupyter-reports | Parameterized notebook execution with papermill: clean-kernel reproducibility, batch reports, aggregation |
| automated-qc-reports | Aggregate QC across samples with MultiQC: scrapes-not-computes, sample-name resolution, report-to-gate |
| figure-export | Publication figures: hybrid rasterization, editable fonts, RGB/CMYK, perceptual colormaps, journal specs |
| publication-tables | Table 1 and result tables: descriptive-vs-inferential statistics, the baseline p-value fallacy, Word/Excel export, gene-symbol safety |

## Method Selection

Report framework:

| Situation | Use |
|-----------|-----|
| R-only analysis | rmarkdown-reports (or quarto-reports with the knitr engine) |
| Python or mixed R/Python, multi-format, websites, native cross-refs | quarto-reports |
| Batch per-sample reports from a Python notebook template | jupyter-reports (papermill) |
| Heavy multi-hour compute | a workflow manager (Snakemake/Nextflow), with a notebook/Quarto report on top |

Deliverable:

| Deliverable | Use |
|-------------|-----|
| Cross-tool QC summary across samples | automated-qc-reports (MultiQC) |
| Publication figure (PDF/TIFF) | figure-export |
| Table 1 / results / supplementary table | publication-tables |

## Example Prompts

- "Create an RMarkdown report for my RNA-seq analysis"
- "Set up a Quarto document with freeze so CI renders without a kernel"
- "Run my analysis notebook on all samples with papermill"
- "Generate a MultiQC report from my pipeline outputs and gate on mapping rate"
- "Export my UMAP figure for Nature without making an unopenable vector file"
- "Make a Table 1 by treatment arm with SMD instead of p-values"
- "Write a gene-symbol-safe supplementary table to Excel"

## Requirements

```bash
# R packages
install.packages(c('rmarkdown', 'knitr', 'bookdown', 'gtsummary', 'flextable', 'gt', 'kableExtra', 'ggrastr', 'ragg', 'showtext'))

# Quarto
# Download from https://quarto.org/docs/download/

# Python
pip install jupyter papermill nbconvert multiqc matplotlib great_tables openpyxl
```

## Related Skills

- **differential-expression** - Analysis to report
- **pathway-analysis** - Enrichment results
- **data-visualization** - Plots feeding figure-export
- **read-qc** - QC outputs MultiQC aggregates
- **clinical-biostatistics** - CONSORT trial reporting context for Table 1
