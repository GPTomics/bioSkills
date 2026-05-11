# nf-core/rnaseq on LANTA and RAPTOR

Source of truth: nf-core/rnaseq usage docs, https://nf-co.re/rnaseq/usage

Use this reference with either `lanta.md` or `raptor.md`. Do not use it alone for cluster-specific paths.

## Pipeline Scope

`nf-core/rnaseq` is an RNA-seq analysis pipeline supporting STAR, RSEM, HISAT2, and Salmon, with gene/isoform counts and extensive QC. Current public docs showed version `3.25.0` when this reference was written; verify the desired release before production runs.

Prefer a pinned release:

```bash
-r 3.25.0
```

Update the version only after checking nf-core release notes and cluster compatibility.

## Samplesheet Contract

The required columns are:

```csv
sample,fastq_1,fastq_2,strandedness
```

Optional columns such as `seq_platform` and `seq_center` can be used when read-group metadata matters. For paired-end data, provide both FASTQ paths. For single-end data, leave `fastq_2` empty if supported by the pipeline version in use.

Use `strandedness` values accepted by the pipeline, commonly:

- `auto`
- `unstranded`
- `forward`
- `reverse`

Use `auto` only when the sample type and pipeline behavior make automatic inference acceptable.

## Reference Strategy

Prefer explicit reference files over relying on implicit defaults:

```bash
--fasta <REFERENCE_DIR>/genome.fa
--gtf <REFERENCE_DIR>/genes.gtf
```

Keep genome FASTA, GTF/GFF, transcriptome, indexes, annotation, and downstream DE assumptions on the same genome build. Do not invent reference paths.

## Cluster Launch Pattern

Use the matching site config:

```bash
# LANTA
-c nextflow.lanta.config
-work-dir /scratch/<PROJECT_ID-SHORTNAME>/<USER>/nextflow-work/rnaseq
--outdir /project/<PROJECT_ID-SHORTNAME>/results/rnaseq

# RAPTOR
-c nextflow.raptor.config
-work-dir <SCRATCH_DIR>/nextflow-work/rnaseq
--outdir <PROJECT_DIR>/results/rnaseq
```

Include run metadata:

```bash
-resume \
-with-report <OUTDIR>/report.html \
-with-trace <OUTDIR>/trace.txt \
-with-timeline <OUTDIR>/timeline.html \
-with-dag <OUTDIR>/flow.svg
```

## Common Parameters to Decide

Ask or verify before final production commands:

- paired-end vs single-end FASTQ layout
- strandedness or whether `auto` is acceptable
- organism and genome build
- FASTA/GTF/reference index location
- whether reads need trimming, contaminant removal, rRNA removal, or UMI handling
- desired aligner/quantifier strategy if the default is not acceptable
- expected sample count and biological replicate structure
- output location and retention policy

## Expected Outputs

Expect a structured `--outdir` with QC, alignment/quantification outputs, MultiQC reports, and pipeline execution metadata. Point users to the pipeline output docs for exact paths because paths can change between releases.

## Troubleshooting

Start with cluster-specific checks from `lanta.md` or `raptor.md`, then check these workflow items:

1. Validate the samplesheet paths and required columns.
2. Confirm all FASTQ files are visible from compute nodes.
3. Confirm reference FASTA/GTF paths are visible from compute nodes.
4. Confirm container engine and cache settings match the cluster.
5. Inspect `.nextflow.log` and failed work directories.
6. Use `trace.txt` to identify processes needing more memory, CPUs, or time.
7. Rerun with `-resume` after fixing the cause.
