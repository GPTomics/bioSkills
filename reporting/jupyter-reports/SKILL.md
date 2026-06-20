---
name: bio-reporting-jupyter-reports
description: Runs parameterized Jupyter notebooks as reproducible batch report generators with papermill, renders them to HTML/PDF with nbconvert, aggregates results across samples, and makes notebook outputs trustworthy. Use when generating per-sample analysis reports, executing a notebook template across many datasets, or fixing notebooks that do not reproduce.
tool_type: python
primary_tool: papermill
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: papermill 2.6+, nbconvert 7.16+, nbclient 0.10+, jupyter-client 8+, scrapbook 0.5+, jupytext 1.16+, nbstripout 0.7+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show papermill` then `help(papermill.execute_notebook)`
- CLI: `papermill --help`, `jupyter nbconvert --help`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Jupyter Reports with Papermill

**"Generate a reproducible analysis report"** -> Execute a parameterized notebook in a clean kernel, producing an executed-notebook artifact, then render it to a code-hidden HTML/PDF report.
- Python: `papermill.execute_notebook(input, output, parameters={...})`
- CLI: `jupyter nbconvert --execute --to html notebook.ipynb`

## The Load-Bearing Idea: A Notebook Has Hidden State; the Report Is the Executed Artifact

A `.ipynb` is JSON holding cell source, **saved outputs**, and per-cell **execution_count** integers - and those three are decoupled. The kernel is a long-lived process with a mutable namespace, so a user can run cell 5, edit cell 2, rerun it, delete cell 3, run cell 7, and save. The stored outputs then reflect a kernel state that NO top-to-bottom rerun reproduces. Non-monotonic execution counts are the forensic tell.

This is intrinsic to the REPL-on-a-document model, not a bug. The discipline answer is mechanical: **Restart Kernel and Run All** before sharing, so saved outputs equal one clean linear run. papermill and `nbconvert --execute` are that discipline automated - they always spin up a FRESH kernel and run cells strictly in document order, so the output is by construction the record of one clean run. The value is not "trust the saved outputs," it is "regenerate them from a known-empty state."

The empirical case (Pimentel et al.): of ~1.4M GitHub notebooks, only ~24% of those executed finished without error and only ~4% reproduced their stored outputs; ~36% had out-of-order cells, and the dominant failure was `ImportError` - environment, not logic. (It is a public-GitHub corpus, so the absolute rates carry selection bias, but the failure mode is the lesson.) Two lessons: saved outputs are not reproducible (always re-execute), and the #1 fix is pinning the environment (see below).

## Two Artifacts People Conflate

- **Executed notebook** - papermill's output, or `nbconvert --execute --to notebook`. A record of the run: same cells, freshly computed outputs, an injected-parameters cell. For audit/debug/aggregation, not for humans to read as prose.
- **Human report** - `nbconvert --to html/pdf`, optionally `--no-input` to hide code. The rendered deliverable.

A papermill pipeline does both: execute the parameterized notebook -> output `.ipynb` (the evidence) -> convert to HTML/PDF (the report).

## Parameterizing with papermill

Tag ONE cell `parameters` holding defaults. At execution papermill inserts a NEW cell tagged `injected-parameters` immediately AFTER it, containing only the overrides; because Python runs top-to-bottom, the injected cell shadows the defaults.

```python
import papermill as pm
pm.execute_notebook('template.ipynb', 'out/sampleA.ipynb',
                    parameters={'sample_id': 'sampleA', 'fdr_threshold': 0.05})
```

- If NO cell is tagged `parameters`, the injected cell goes at the TOP and downstream references to a param raise NameError. Keep all parameters in one tagged cell.
- **Type coercion is the #1 CLI gotcha.** `-p name value` YAML-parses the value (`-p n 5` -> int, `-p flag true` -> bool); `-r name value` keeps it a raw string (`-r chrom 1` -> `"1"`, needed for sample IDs like `007` or chromosome `"1"`). `-y`/`--parameters_yaml` and `-f`/`--parameters_file` pass lists/dicts. The Python API passes native objects directly with no YAML round-trip - prefer it in pipelines to avoid coercion surprises.

## Execution Controls That Matter in Pipelines

- `execution_timeout` (CLI `--execution-timeout`) - seconds per cell; **default is forever (None)**. A long bioinformatics cell hangs a pipeline silently without this. Set it.
- `kernel_name` / `--kernel` - must be a REGISTERED kernelspec; a mismatch is a top failure cause. The kernel must point at the pinned environment.
- `log_output=True` / `--log-output` - stream each cell's stdout/stderr for CI visibility.
- On a cell exception papermill raises `PapermillExecutionError`, STILL writes the output notebook with the traceback captured (evidence preserved), then exits non-zero. Tag a cell `raises-exception` to allow it to fail without aborting. Fail loud, keep the artifact.
- papermill reads/writes notebooks from `s3://`, `gs://`, `adl://`/`abs://`, and `http(s)://` out of the box (cloud connectors are extras: `pip install papermill[s3]`), so serverless per-sample execution works.
- The fresh-kernel guarantee holds when papermill owns the kernel lifecycle (the default `execute_notebook` path). Advanced patterns that reuse a kernel across notebooks forgo the from-empty-state guarantee - let papermill create the kernel.

## Rendering to a Report

`nbconvert --execute` runs the notebook top-to-bottom in a fresh kernel and writes regenerated outputs - that re-execution, not the saved outputs, is what makes a converted report reproducible.

```bash
jupyter nbconvert --to html --no-input out/sampleA.ipynb     # code-hidden stakeholder report
jupyter nbconvert --execute --to html template.ipynb         # execute then render in one step
jupyter nbconvert --to webpdf out/sampleA.ipynb              # PDF without LaTeX
```

- `--to pdf` needs a LaTeX toolchain (xelatex + pandoc) - the classic CI pain. `--to webpdf` renders HTML then prints via headless Chromium (`pip install nbconvert[webpdf]`, `--allow-chromium-download`) and avoids TeX entirely. The tradeoff: webpdf uses Chromium page geometry, so wide tables and long code lines can clip at the page edge, whereas LaTeX `--to pdf` paginates and wraps better for table-heavy reports.
- `--no-input` hides all code; `--no-prompt` drops the `In[ ]:`/`Out[ ]:` prompts; `TagRemovePreprocessor.remove_cell_tags` strips cells by tag (tag setup cells to remove them). Custom branded layouts use the nbconvert 6+ directory-template system (`--template <name>`).

## Aggregating Results Across Samples

`papermill.record` is deprecated (since papermill 1.0). Use **scrapbook**: in the template, `import scrapbook as sb; sb.glue('auc', 0.91)` records named scraps (and `sb.glue('fig', obj, display=True)` for figures). Downstream, `sb.read_notebooks('reports/').papermill_dataframe` aggregates every executed notebook's scraps into one tidy table - the canonical pattern for looping papermill over a sample sheet then collecting per-sample QC into a cohort summary.

## Version Control: Never Commit Outputs

`.ipynb` is JSON with embedded base64 outputs and `execution_count`, so committing it raw gives giant unreviewable diffs, brutal merge conflicts, and leaked data. Three complementary fixes:

- **nbstripout** - a git clean filter (`*.ipynb filter=nbstripout` in `.gitattributes`) that strips outputs and execution counts on `git add`; the working copy keeps its outputs. Highest-leverage single fix.
- **jupytext** - pairs `.ipynb` with a text twin (`py:percent` is runnable and diff-friendly, or Markdown/`.qmd`); commit the text file as source of truth, regenerate outputs. Lets a notebook be code-reviewed as a normal PR.
- **jupyter-cache** - caches executed outputs keyed by a code-cell hash; the engine Jupyter Book/MyST-NB and Quarto use to skip re-running expensive cells. Right tool when notebooks are doc sources with long compute.

## Reproducible Execution Is Not Reproducible Science

State this plainly: papermill and `nbconvert --execute` guarantee EXECUTION ORDER and a clean kernel - they do NOT guarantee the ENVIRONMENT or numerical determinism. The same parameterized notebook rerun against newer numpy/scanpy, a different BLAS, or an unseeded RNG runs flawlessly and produces DIFFERENT numbers. Order-reproducibility is necessary, not sufficient. Pair papermill with: a pinned environment (conda lockfile / `requirements.txt` with exact versions, ideally a container; the registered kernel must point at it), seeded RNGs for any stochastic step (clustering, UMAP, splits, bootstraps), and pinned reference/DB versions. To catch silent drift, **nbval** (`pytest --nbval`) re-executes and compares new outputs against stored ones in CI.

## When a Notebook Report Is the Wrong Tool

The notebook is the REPORT, not the PIPELINE. Heavy compute (alignment, variant calling, large scanpy integration, anything multi-hour or needing a scheduler, retries, or parallel fan-out) belongs in a workflow manager (Snakemake/Nextflow/WDL). papermill notebooks shine as the final per-sample or per-cohort summary plus figures over already-computed results. Rule of thumb: if `--retry`, `--cluster`, or a DAG is wanted, it is a pipeline; if a parameterized HTML/PDF of results is wanted, it is a papermill report. Note Quarto can consume a `.ipynb` directly, so papermill (parameterize/execute) and Quarto (render) interoperate.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Saved outputs do not match a rerun | Out-of-order interactive execution / hidden state | Restart-and-run-all; let papermill/`nbconvert --execute` regenerate |
| Parameter is a bool when a string was wanted | `-p` YAML-parses values | use `-r name value` for raw strings, or the Python API |
| NameError on a parameter | No `parameters`-tagged cell | tag one cell `parameters`; keep all params there |
| Pipeline hangs on one cell | `execution_timeout` defaults to forever | set `execution_timeout` |
| `--to pdf` fails in CI | no LaTeX toolchain | use `--to webpdf` (headless Chromium) |
| Giant notebook diffs / leaked data in git | committing outputs | nbstripout filter or jupytext pairing |
| Reruns months later give different numbers | environment/seed not pinned | pin env + container, seed RNGs, nbval in CI |

## Related Skills

- reporting/quarto-reports - Document-first reporting that can render the same notebooks
- reporting/rmarkdown-reports - R-based literate reports
- workflows/scrnaseq-pipeline - Heavy compute that belongs in a pipeline, with the notebook as the summary report
- single-cell/preprocessing - Analysis embedded inside per-sample notebook templates

## References

- Pimentel JF, Murta L, Braganholo V, Freire J. A large-scale study about quality and reproducibility of Jupyter notebooks. MSR '19, IEEE/ACM. 2019:507-517. doi:10.1109/MSR.2019.00077
- Pimentel JF, Murta L, Braganholo V, Freire J. Understanding and improving the quality and reproducibility of Jupyter notebooks. Empir Softw Eng. 2021;26:65. doi:10.1007/s10664-021-09961-9
- Rule A, Birmingham A, Zuniga C, et al. Ten simple rules for writing and sharing computational analyses in Jupyter Notebooks. PLoS Comput Biol. 2019;15(7):e1007007. doi:10.1371/journal.pcbi.1007007
- Kluyver T, Ragan-Kelley B, Pérez F, et al. Jupyter Notebooks - a publishing format for reproducible computational workflows. ELPUB 2016, IOS Press. 2016:87-90. doi:10.3233/978-1-61499-649-1-87
