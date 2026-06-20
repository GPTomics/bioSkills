---
name: bio-reporting-automated-qc-reports
description: Aggregates per-tool QC metrics (FastQC, fastp, alignment, quantification, variant calling, single-cell) into one interactive MultiQC report, and guides module scoping, sample-name resolution, large-cohort behavior, and turning the report into an actual QC gate. Use when summarizing QC across many samples, building a shareable quality report, or wiring automated QC into a pipeline.
tool_type: cli
primary_tool: multiqc
---

## Version Compatibility

Reference examples tested with: MultiQC 1.21+ (Plotly era), FastQC 0.12+, STAR 2.7.11+, Subread 2.0+, salmon 1.10+, samtools 1.19+, Picard 3.1+, fastp 0.23+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

Config keys and defaults move between MultiQC releases (the plotting backend changed from HighCharts to Plotly at 1.20; flat-plot and AI thresholds shifted). When a default matters, confirm it against the installed version: `python -c "import multiqc; print(multiqc.__version__)"` and check that version's `config_defaults.yaml`.

If code throws an error, run `multiqc --help` and adapt flags to the installed version rather than retrying.

# Automated QC Reports with MultiQC

**"Aggregate QC results into one report"** -> Walk a directory of tool outputs, parse the metrics those tools already wrote, and render one interactive HTML report plus a parseable `multiqc_data/` directory.
- CLI: `multiqc <dir>` (scans for recognized tool outputs)

## The Load-Bearing Idea: MultiQC Aggregates, It Does Not Measure

MultiQC computes nothing. It SCRAPES the log/metrics files that FastQC, STAR, Picard, salmon, bcftools, etc. already wrote, re-tabulates those numbers, and renders them. Every value in a report traces back to an upstream tool's output file. Four consequences that drive every real decision below:

- **The report is a triage SNAPSHOT, not a pass/fail gate.** There is no "MultiQC said FAIL -> stop the pipeline." MultiQC has no fail-on-threshold and exits 0 on bad QC. Green/amber/red is either the upstream tool's own status (FastQC writes PASS/WARN/FAIL; MultiQC just displays it) or a threshold a human configured. Gating is a SEPARATE step (see From Report to Gate).
- **Garbage upstream becomes a clean-looking report.** Run a tool with the wrong strandedness, wrong reference, or dedup on amplicon data, and MultiQC faithfully aggregates the wrong numbers into a polished HTML. Polish is not evidence of correctness.
- **An empty report means "nothing matched," not "QC passed."** Point MultiQC at the wrong directory or over-filter modules and it emits a near-empty report with a log warning and exit 0. Always check the sample count in the header against the roster expected.
- **It is only as current as its parsers.** Each module is a hand-written parser keyed to a specific log format. An upstream version bump that changes a header can silently drop a file or mis-map a column.

## Basic Usage

```bash
multiqc results/ -o qc_report/            # scan results/, write qc_report/multiqc_report.html
multiqc results/ -n project_qc -o qc/     # custom report name
multiqc results/ -m fastqc -m star        # ONLY these modules (see scoping below)
multiqc results/ -c multiqc_config.yaml   # reproducible config-driven report
```

## Supported Tools

MultiQC ships parsers for 100+ tools. Common assay groupings:

| Stage | Tools with modules |
|-------|--------------------|
| Read QC | FastQC, fastp, Cutadapt, falco |
| Alignment | STAR, HISAT2, BWA, Bowtie2, samtools, Qualimap, Picard |
| Quantification | featureCounts, Salmon, kallisto, RSeQC |
| Variant calling | bcftools, GATK, Picard, SnpEff, VEP |
| Single-cell | Cell Ranger, STARsolo |

## Module Detection Is Regex - Scope It

Detection runs off `search_patterns.yaml`: each module declares a filename glob/regex (`fn`/`fn_re`) and/or a file-content match (`contents`/`contents_re`, bounded by `num_lines`). Loose patterns (`*.txt`, `*.log`, `*.json`) in a messy directory cause FALSE module matches and PHANTOM samples - a file that is not really that tool's output gets parsed as one. A single file can also satisfy two modules.

Scope explicitly rather than trusting auto-detection across thousands of samples:

```bash
multiqc results/ --ignore "*_tmp/" --ignore "work/"   # drop paths from the search
multiqc results/ -m fastqc -m star -m salmon          # run ONLY named modules
multiqc results/ -e snippy -e custom_content          # run all EXCEPT named modules
```

Tighten an over-loose pattern by overriding `sp:` in the config (`sp: {mytool: {fn: 'real_name_*.txt'}}`). Production configs pin `sp:` and `module_order` instead of relying on detection.

## Sample Names Are Derived, Not Declared

Sample names are NOT read from a manifest. MultiQC derives each name from the matched filename (or a sample column inside the file), then "cleans" it by trimming a ~100-entry default list of extensions (`fn_clean_exts`: `.gz`, `.fastq`, `.bam`, `_fastqc`, ...). This is how `sampleA_R1.fastq.gz`, `sampleA.sorted.bam`, and `sampleA.salmon/` all collapse to one `sampleA` row gathering read, alignment, and quant metrics.

The same mechanism is the #1 large-cohort bug:
- **Merge** - two genuinely different inputs clean to the same name (two lanes both reduce to `sampleA`) and silently overwrite each other's metrics.
- **Split** - one sample appears as several rows because different tools cleaned its name differently (one kept `_L001`, another stripped it).

`multiqc_data/multiqc_sources.txt` maps every parsed file to the sample name it produced - read it first when diagnosing duplicate/missing rows. Controls:

| Need | Control |
|------|---------|
| Add suffixes to strip (keep defaults) | `extra_fn_clean_exts:` in config (do NOT override `fn_clean_exts`, which replaces the defaults) |
| Use the log filename as the name | `--fn_as_s_name` (config `use_filename_as_sample_name`) |
| Disambiguate by directory | `--dirs` / `-d`, `--dirs-depth N` |
| Keep full names, no cleaning | `--fullnames` / `-s` |
| Rename at report time | `--replace-names map.tsv` (pattern -> replacement, two columns) |
| Offer toggleable name sets | `--sample-names headered.tsv` (relabel buttons, does not merge rows) |

## General Statistics and Conditional Formatting Are Configured, Not Authoritative

The General Statistics table is one row per sample with columns each module contributes. Cell colors come from `table_cond_formatting_rules` (numeric `gt`/`lt`/`eq`/`ge`/`le`, string `s_eq`/`s_contains`/`s_ne`). A red ">10% duplication" cell is red because someone wrote that rule (or because a module ships a built-in default rule), not because biology says 10% is bad. Treat formatting as a configured convenience; absence of red is not a pass, and presence of red is not a biological verdict. Column visibility/order/naming are config too (`table_columns_visible`, `table_columns_placement`, `table_columns_name`).

## Large Cohorts: MultiQC Downgrades Automatically

To keep the single HTML openable, MultiQC silently changes rendering as series counts grow. The exact thresholds have moved across versions - verify against the installed `config_defaults.yaml` - but the behaviors are:

| Behavior | Config key | Effect |
|----------|-----------|--------|
| Table -> violin/beeswarm plot | `max_table_rows` (~500) | above the limit the General Stats "table" becomes a distribution plot; per-cell view is lost |
| Interactive plot deferred | `plots_defer_loading_numseries` (~100) | viewer must click to render |
| Interactive -> flat image | `plots_flat_numseries` (moved across versions; HighCharts-era 100, current default much higher) | plots render as static PNG/SVG |

Force a mode for reproducible visuals across cohort sizes: `--flat` / `--interactive` (config `plots_force_flat` / `plots_force_interactive`). At tens of thousands of samples, also scope with `-m`/`--ignore` or split into per-batch reports - MultiQC holds all parsed data in memory before rendering.

## Custom Content (Injecting Custom Metrics)

Two mechanisms; `--custom-data-file` does NOT exist.

- **`_mqc` suffix** - any file named `*_mqc.{tsv,csv,txt,yaml,json,png,...}` is auto-discovered and rendered with no config. The suffix is what makes it findable.
- **`custom_data` in the config** - define a section with `plot_type` (`bargraph`, `linegraph`, `table`, `generalstats`, `image`, ...) and inline data or a search pattern. `plot_type: generalstats` injects columns straight into General Statistics.

## From Report to Gate (the decision MultiQC does not make)

MultiQC is a viewer; QC GATING is separate. The machine-readable truth lives in `multiqc_data/`: `multiqc_data.json` (all parsed values), per-module `multiqc_*.txt` tables, and `multiqc_general_stats.txt`. Build a gate ON TOP of that file, not by scraping the HTML:

```bash
multiqc results/ -o qc/ --data-format json     # write multiqc_data.json
# a downstream script parses qc/multiqc_data/multiqc_data.json,
# applies thresholds, and exits non-zero / quarantines failing samples.
```

This is the correct division of labor: MultiQC presents; the pipeline (nf-core modules, a purpose-built gater like CheckQC, a Nextflow/Snakemake check, or a parse-and-exit script) decides. Building fail-on-threshold logic inside MultiQC is a category error.

## AI Summaries Send Data Off-Network

MultiQC (1.27+) can prepend an LLM-written natural-language summary (`--ai` / `--ai-summary`, `--ai-summary-full`; providers via `ai_provider`: `seqera`, `openai`, `anthropic`, `aws_bedrock`, `custom`; keys via `OPENAI_API_KEY` / `ANTHROPIC_API_KEY` / `SEQERA_ACCESS_TOKEN`). It is OFF by default. When enabled it transmits the aggregated QC metrics - and, unless `ai_anonymize_samples` is set, the SAMPLE NAMES - to an external API over the internet. For clinical, patient, or embargoed data this can be a data-governance violation; use `--no-ai` to strip AI controls from a shared report, or the in-browser on-demand mode (summary stays in browser local storage, not baked into the distributed HTML). Confirm the exact key spelling against the installed version.

## Reproducibility

- **Pin both the MultiQC version and the upstream tool versions.** Parsers, default thresholds, and column sets shift between releases; a report from 1.30 is not byte-stable against one from 1.14.
- **Pin a config, treat the report as a deliverable.** The nf-core pattern: ship `multiqc_config.yml` locking `title`, `module_order`, sample-name cleaning, and `sp:` patterns; expose `--multiqc_config` to layer a user config on top (both apply, user wins). nf-core also emits a `methods_description_template.yml` so the report carries auto-generated methods text and citations for only the tools that ran, plus a consolidated software-versions table.
- **Cross-check the roster.** If 3 of 100 BAMs failed earlier and produced no metrics, MultiQC reports a clean 97-sample report with no indication 3 are missing. It aggregates what exists and has no notion of an expected sample set.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Near-empty report, exit 0 | No files matched a search pattern (wrong dir, over-filtered) | Check header sample count; read `multiqc_sources.txt`; relax `-m`/`--ignore` |
| Two samples merged into one row | Names collide after `fn_clean_exts` cleaning | `extra_fn_clean_exts`, `--dirs`, or `--replace-names`; verify in `multiqc_sources.txt` |
| One sample split across rows | Tools cleaned the name differently | `--fn_as_s_name` or `extra_fn_clean_exts` to normalize |
| Phantom sample / wrong module | Loose pattern matched an unrelated file | `--ignore` the path or tighten `sp:`; restrict with `-m` |
| Metric missing after a tool upgrade | Upstream log-format drift broke the parser | Pin tool + MultiQC versions; check the module changelog |
| "table" rendered as a violin plot | Rows exceeded `max_table_rows` | Raise the limit or split the cohort |
| Sensitive sample names left the network | AI summary enabled | `--no-ai`, or `ai_anonymize_samples`; default is off |

## Related Skills

- read-qc/quality-reports - Generate the FastQC/falco inputs MultiQC aggregates
- read-qc/fastp-workflow - Preprocessing QC that feeds the report
- read-qc/rnaseq-qc - Post-alignment RNA-seq metrics surfaced in MultiQC
- workflows/rnaseq-to-de - Full pipeline that emits a MultiQC report as a deliverable

## References

- Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. doi:10.1093/bioinformatics/btw354
- MultiQC documentation: docs.seqera.io/multiqc (post-2024 features incl. Plotly plots and AI summaries are documented here and in the GitHub CHANGELOG, not a separate paper)
