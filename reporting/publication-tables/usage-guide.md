# Publication Tables Usage Guide

## Overview

This guide covers generating publication-ready tables programmatically - descriptive Table 1, regression and differential-expression result tables, and supplementary tables. The decisions that matter are choosing the right statistics (descriptive vs inferential, and never a baseline p-value on a randomized Table 1), exporting to the journal's actual target format (Word fidelity is the chronic pain point), and protecting gene symbols from Excel auto-conversion. A table, like a figure, is a deterministic function of data and code - generate it, never hand-edit it.

## Prerequisites

```r
# R
install.packages(c('gtsummary', 'gt', 'flextable', 'kableExtra', 'tableone'))
```

```bash
# Python
pip install great_tables pandas tableone openpyxl
```

## Quick Start

Tell your AI agent what you want to do:
- "Make a Table 1 of baseline characteristics by treatment arm"
- "Export my DESeq2 results as a formatted Word table"
- "Build a regression table with odds ratios and 95% CIs"
- "Write a supplementary gene table to Excel without SEPT9 turning into a date"
- "Use SMD instead of p-values to show baseline balance"

## Example Prompts

### Descriptive Tables

> "Create a Table 1 by arm with median (IQR) for continuous variables and show missing counts"

> "Add standardized mean differences instead of a p-value column to my baseline table"

### Result Tables

> "Turn my logistic regression into a table of odds ratios with confidence intervals"

> "Merge the univariable and multivariable regression results side by side"

### Export

> "Export this gtsummary table to Word with the formatting intact"

> "Write my full DE results as a gene-symbol-safe CSV and Excel supplement"

## What the Agent Will Do

1. Decide descriptive vs inferential and pick the right summary statistics (median/IQR for skewed data, n(%) for categorical)
2. For a randomized Table 1, use SMD rather than p-values, and show missingness explicitly
3. Generate the table with gtsummary/gt/flextable (R) or great_tables/tableone (Python)
4. Export to the requested target format (flextable for Word, gt/kableExtra for LaTeX/HTML)
5. Protect gene-symbol columns when writing Excel/CSV supplements

## Tips

- A baseline p-value in a randomized trial tests a null known to be true; report SMD for balance, not p-values.
- gtsummary defaults continuous variables to median (IQR) because biomedical variables are usually skewed - override to mean(SD) only after checking normality.
- Show missingness (`missing="ifany"`); never compute percentages that silently hide dropped rows.
- flextable `save_as_docx()` is the most reliable Word path; great_tables cannot export Word or LaTeX.
- Report to significant figures, not the float default; p-values to 2-3 sig figs or "<0.001".
- For Excel gene tables, prefer CSV (import as text) or force the gene column to text format `'@'`; verify by reopening.
- DT tables are for interactive exploration only, never the published static table.

## Related Skills

- reporting/figure-export - The figure counterpart to table export
- reporting/rmarkdown-reports - Embedding tables in R reports
- clinical-biostatistics/trial-reporting - CONSORT context for Table 1
- differential-expression/de-results - Result tables these formatters present
