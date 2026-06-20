---
name: bio-reporting-publication-tables
description: Builds publication-ready tables - descriptive Table 1, regression and differential-expression result tables, and supplementary tables - with gtsummary, gt, flextable, and kableExtra (R) or great_tables, pandas, and tableone (Python), choosing the right statistics and the right export format. Use when making a Table 1, exporting a formatted results table for a paper, or writing a gene-symbol-safe supplementary table.
tool_type: mixed
primary_tool: gtsummary
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: gtsummary 2.0+, gt 0.10+, flextable 0.9+, kableExtra 1.4+, great_tables 0.13+, pandas 2.2+, tableone 0.9+, openpyxl 3.1+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('gtsummary')` then `?tbl_summary` (gtsummary had a major API refresh at v2.0)
- Python: `pip show great_tables` then `help(great_tables.GT.save)`

If code throws an error, introspect the installed package and adapt the example to the actual API rather than retrying.

# Publication-Ready Tables

**"Make my Table 1"** / **"export this results table for the paper"** -> Generate the table programmatically with the right statistics and export it to the journal's target format.
- R: `gtsummary::tbl_summary(data, by=arm)` then `as_flex_table()` -> Word
- Python: `great_tables.GT(df)` -> HTML/PNG; `tableone.TableOne(...)` for a descriptive table

## The Load-Bearing Idea: A Table Is Structure + Precision + the Right Statistics

Three orthogonal concerns, and conflating them is where tables go wrong:

- **Structure** - rows are units (subjects, genes, models), columns are variables/groups, spanners group columns, footnotes/source-notes carry the apparatus. This is the grammar of tables that gt, great_tables, and flextable all encode.
- **Precision** - report to MEANINGFUL precision, not the float default. P-values to 2-3 significant figures or "<0.001"; estimates to the precision the CI supports; percentages to 0-1 decimal. A `3.14159265` mean is noise.
- **The right statistics** - a table is DESCRIPTIVE (summarize the sample: n(%), mean(SD) or median(IQR)) or INFERENTIAL (estimate + CI + test statistic, with the effect size primary and the p-value never alone). Declare which.

The deepest framing, shared with figures: a table is a deterministic function of data + code. Same input + same code -> the same table, byte for byte. Manual edits in Word break this; every number must trace to a line of code. That principle dictates the tooling - generate programmatically and never hand-edit the output.

## Tool Decision (by language and target format)

Target format dominates: Word for most biomedical journals, LaTeX for some genomics/physics venues, HTML for web/preprints/Quarto, CSV/Excel for machine-readable supplements.

| Need | Tool | Path |
|------|------|------|
| R, descriptive Table 1 or regression results | **gtsummary** | `tbl_summary` / `tbl_regression` -> `as_flex_table()` / `as_gt()` |
| R, going to Word/PowerPoint | **flextable** | `save_as_docx()` (most reliable Word fidelity; pairs with officer) |
| R, going to LaTeX/PDF | **gt** or **kableExtra** | `gtsave('x.tex')` / `kbl(format='latex')` |
| R, going to HTML/Quarto | **gt** or **kableExtra** | `as_raw_html()` / `save_kable()` |
| Python, display HTML/image | **great_tables** | `GT(df)` -> `save('x.png')` / `as_raw_html()` |
| Python, going to LaTeX | **pandas Styler** | `df.style.format(...).to_latex()` (pandas 1.3+) |
| Either, classic Table 1 with SMD | **tableone** | `CreateTableOne` (R) / `TableOne` (Python) |
| Machine-readable supplement | **CSV** (preferred) or Excel | gene-symbol-safe export (below) |

`gt` has the broadest R export (HTML/PNG/PDF/RTF/LaTeX/Word). `great_tables.save()` is image+PDF only (HTML via `as_raw_html()`/`write_raw_html()`) - NO native Word or LaTeX, a real limitation vs the R stack. **DT is for interactive exploration and online-only/interactive HTML supplements** - never for a static print/PDF table.

## Table 1 Is Descriptive, Not Inferential

Table 1 reports baseline characteristics so the reader can judge who was studied and how comparable the groups are. It describes the sample; it is not a place to test hypotheses.

**The p-value fallacy (randomized trials):** adding a p-value column comparing arms in a randomized trial is discouraged by CONSORT and statisticians. In a properly randomized trial any baseline imbalance is by definition due to chance, so the test asks whether a difference could have arisen by chance when the assignment WAS by chance - it tests a null already known true. A "significant" baseline p-value is a Type I error by construction; a non-significant one tells nothing new. The right response to a worrying imbalance on a prognostic covariate is to adjust for it (pre-specified ANCOVA covariate), not test it (Senn 1994; CONSORT 2010 item 15). gtsummary's documentation cautions against `add_p()` on a randomized Table 1 for this reason.

- **Randomized Table 1:** no p-value column. To convey balance, use **standardized mean differences (SMD)** - they describe the magnitude of imbalance (|SMD| > 0.1 is a common "notable" rule of thumb) without the inferential fallacy. If a journal or regulator nonetheless requires a baseline comparison column, report it but interpret per Senn: a "significant" baseline difference in a properly randomized trial is a Type I error, not evidence of confounding. Note `add_difference()` compares exactly two groups; for >2 arms, SMD is defined pairwise, so report reference-group or all-pairs SMDs rather than one omnibus value.
- **Observational studies:** a comparison column can be defensible (the groups genuinely may differ), but multiplicity (many rows -> many tests) and "significant does not mean important" still bite, and SMD is the standard balance diagnostic in propensity-score/causal contexts. SMD is preferred over p-values for balance in essentially all cases.

## Continuous Summaries and Missingness

- **mean(SD) vs median(IQR) is distribution-driven.** Symmetric -> mean(SD); skewed/heavy-tailed (most biomarkers, counts, lab values, length-of-stay) -> median(IQR), because the mean is pulled by the tail. gtsummary defaults continuous variables to `median (p25, p75)` - defensible because biomedical variables are usually skewed and normality cannot be assumed column by column. Override per-variable (`statistic = list(age ~ "{mean} ({sd})")`) only after checking normality. Categorical: n(%), and state row% vs column% (baseline tables want column%).
- **Missingness must be SHOWN, not silently dropped.** The cardinal sin is computing percentages on complete cases with no indication rows were dropped - a reader cannot tell 90% from 90%-of-the-60%-with-data. gtsummary's `missing = "ifany"` (default) shows a missing row when any value is absent; `missing_text = "Unknown"` labels it. Never compute denominators that hide missingness; if rows are dropped, report the analyzed N in the table or a footnote.

## Getting Formatting Into the Target Format

- **Word** (the chronic pain point): flextable `save_as_docx()` is the most reliable path; gtsummary `as_flex_table()` then `save_as_docx()` gives the best Word fidelity; gt `gtsave('x.docx')` works but routes through rmarkdown and supports fewer Word styles. great_tables has NO native Word export.
- **LaTeX:** `gtsave('x.tex')`, `kableExtra::kbl(format='latex', booktabs=TRUE)`, or pandas `Styler.to_latex()`.
- **HTML:** gt `as_raw_html()`, great_tables `as_raw_html()`/`write_raw_html()`, kableExtra `save_kable()`.

## The Excel Gene-Symbol Hazard

Excel, with default settings, auto-converts gene symbols and IDs when it PARSES them - opening a CSV, double-clicking, or typing: `SEPT2` -> `2-Sep`, `MARCH1` -> `1-Mar`; RIKEN IDs like `2310009E13` -> `2.31E+13` (precision lost irreversibly); long numeric accessions lose trailing digits to float rounding. (The corruption is a parse behavior, not a write behavior - a string written by openpyxl stays intact until Excel re-interprets it, which is why forcing text format matters.) Ziemann et al. 2016 found ~19.6% of papers with supplementary Excel gene lists affected; Abeysooriya et al. 2021 showed it persisted at 30.9% and drove HGNC to rename the families (SEPT->SEPTIN, MARCH->MARCHF) in 2020 - biology changed its nomenclature to defend against a spreadsheet bug.

When a gene table must reach Excel:
- **Prefer CSV** and tell the consumer to import the gene column as Text (not double-click).
- **If writing .xlsx**, set the gene column to Excel text format `'@'` (openpyxl `cell.number_format = '@'`; XlsxWriter `add_format({'num_format': '@'})` or `write_string()`). Note pandas `Styler.format` is IGNORED by `to_excel` - set the number format via the writer, not the styler.
- **Verify by reopening** - the only sure check.

## Precision and Locale

Report to significant figures, not the float default (`fmt_number(decimals=)`, `pvalue_fun`, Styler `.format('{:.2f}')`). Watch the decimal-comma locale trap: a CSV written in a `,`-decimal locale becomes unparseable elsewhere and Excel may re-misinterpret columns. Write numeric supplements with `.`-decimal and document the locale rather than relying on the system default.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| p-value column on a randomized Table 1 | testing a null known to be true | drop it; use SMD for balance |
| Percentages do not add up / hide dropped rows | missingness silently excluded | `missing="ifany"`; report analyzed N |
| `SEPT2` became a date in the supplement | Excel auto-conversion | CSV + import-as-text, or `'@'` text format in .xlsx |
| mean(SD) misleads on a skewed variable | wrong summary statistic | median(IQR) for skewed data |
| Word table lost its formatting | exported HTML/LaTeX into Word | flextable `save_as_docx()` |
| great_tables won't save to Word/LaTeX | not supported (image/HTML/PDF only) | use the R stack, or export PNG/HTML |
| Excel export ignored my number format | `Styler.format` is dropped by `to_excel` | set `number_format` via the ExcelWriter |

## Related Skills

- reporting/figure-export - The figure counterpart to table export
- reporting/rmarkdown-reports - Embedding kable/gt tables in R reports
- reporting/quarto-reports - Embedding tables in Quarto reports
- clinical-biostatistics/trial-reporting - CONSORT trial reporting context for Table 1
- differential-expression/de-results - Result tables these formatters present

## References

- Senn S. Testing for baseline balance in clinical trials. Stat Med. 1994;13(17):1715-1726. doi:10.1002/sim.4780131703
- Schulz KF, Altman DG, Moher D; CONSORT Group. CONSORT 2010 Statement: updated guidelines for reporting parallel group randomised trials. BMC Med. 2010;8:18 (item 15, baseline table). doi:10.1186/1741-7015-8-18
- Ziemann M, Eren Y, El-Osta A. Gene name errors are widespread in the scientific literature. Genome Biol. 2016;17:177. doi:10.1186/s13059-016-1044-7
- Abeysooriya M, Soria M, Kasu MS, Ziemann M. Gene name errors: Lessons not learned. PLoS Comput Biol. 2021;17(7):e1008984. doi:10.1371/journal.pcbi.1008984
