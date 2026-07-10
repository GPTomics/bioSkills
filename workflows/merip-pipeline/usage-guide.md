# MeRIP-seq Pipeline - Usage Guide

## Overview

Complete workflow from MeRIP-seq raw FASTQ to differential m6A peaks, DRACH motif sanity check, and Guitar transcript-feature metagene with stop-codon enrichment as the biological QC anchor. Chains fastp trimming, STAR splice-aware genome alignment (no deduplication for non-UMI MeRIP), deepTools replicate-concordance and IP-enrichment QC, PreSeq saturation curves, exomePeak2 peak calling (transcript-aware GC-corrected GLM), the four-BAM-vector differential interface (`bam_ip` + `bam_input` + `bam_treated_ip` + `bam_treated_input`), ChIPseeker feature annotation, and Guitar visualisation. Defers per-skill deep treatment to the four `epitranscriptomics/` skills.

## Prerequisites

```bash
conda install -c bioconda star samtools deeptools preseq fastp multiqc macs3 homer bedtools
```

```r
# NOTE: exomePeak2 exists only in Bioconductor 3.18 (deprecated 3.19, removed 3.20). On current
# Bioconductor install it from the Bioc 3.18 archive (or substitute a successor m6A caller).
BiocManager::install(c('exomePeak2', 'GenomicFeatures', 'BSgenome.Hsapiens.UCSC.hg38',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene', 'rtracklayer',
                       'ChIPseeker', 'Guitar'))
```

## Quick Start

- "Run the complete MeRIP-seq pipeline from FASTQ to differential m6A peaks"
- "Analyse my MeRIP-seq data for m6A peaks with QC, differential, and metagene plot"
- "Process my epitranscriptomics data end-to-end with exomePeak2 + Guitar"
- "Build a Snakemake / Nextflow wrapper around the four epitranscriptomics skills"

## Example Prompts

### Full Pipeline

> "Run the complete MeRIP-seq pipeline on 3 control + 3 treatment paired IP/input libraries; trim with fastp, align with STAR to GRCh38, do NOT deduplicate; QC with deepTools plotFingerprint + replicate Spearman + PreSeq saturation; call peaks with exomePeak2; cross-check with MACS3 broad; confirm DRACH with HOMER; run exomePeak2 differential via the four-BAM-vector interface; annotate with ChIPseeker; render Guitar stop-codon metagene as biological QC anchor."

> "Build a Snakemake workflow wrapping fastp -> STAR -> deepTools QC -> PreSeq -> exomePeak2 -> exomePeak2 differential -> ChIPseeker -> Guitar."

### Specific Steps

> "Just call exomePeak2 peaks from my aligned IP/Input BAMs."

> "Run exomePeak2 differential between control and treatment using bam_ip + bam_treated_ip; export the differential peaks with FDR < 0.05 and |log2FC| > 0.5."

> "Render the canonical Guitar metagene from my exomePeak2 peaks; confirm stop-codon enrichment."

## What the Agent Will Do

1. Adapter-trim IP and Input FASTQ with fastp (minimum read length 25; no UMI handling)
2. Align IP and Input separately with STAR splice-aware to the GENOME (not transcriptome)
3. Sort, index, and run flagstat / idxstats per BAM
4. Skip deduplication (standard non-UMI MeRIP convention)
5. Run deepTools multiBamSummary + plotCorrelation for replicate Spearman
6. Run deepTools plotFingerprint for IP enrichment QC; flag failed IPs (JS distance < 0.5)
7. Run PreSeq lc_extrap for library complexity; recommend rarefaction depth for cross-condition comparison
8. Call peaks with exomePeak2 against the matched TxDb + BSgenome
9. Cross-check with MACS3 `--broad --keep-dup all` (optional)
10. Confirm DRACH motif enrichment via HOMER `findMotifsGenome.pl -rna` (sanity check on peak set, not a per-peak filter)
11. Run exomePeak2 differential by populating `bam_treated_ip` + `bam_treated_input` (no `mode=` or `experiment_design=` argument; those don't exist)
12. Apply effect-size filter (|log2FC| >= 0.5) AND FDR < 0.05 AND replicate-direction concordance
13. Annotate peaks against transcript features with ChIPseeker
14. Render canonical Guitar transcript-feature metagene; expect stop-codon enrichment as biological QC anchor
15. Flag peaks within 50 nt of TSS as m6A-or-m6Am ambiguous
16. Generate IP-over-Input log2 bigWig and per-locus pyGenomeTracks browser figures for figure assembly

## Tips

- Standard MeRIP has NO UMI: do NOT deduplicate. Picard MarkDuplicates collapses real coverage at high-expression transcripts.
- Align to the GENOME for MeRIP peak calling. Transcriptome alignment is for m6anet-analysis only (ONT DRS).
- Peak counts are library-size-dependent. Rarefy to common depth or report saturation curves alongside.
- exomePeak2 has NO `mode=` or `experiment_design=` argument. Populate `bam_treated_ip` + `bam_treated_input` for differential.
- For batch / antibody-lot covariate adjustment, fall through to featureCounts-on-peaks then DESeq2 (see `epitranscriptomics/m6a-differential`).
- DRACH motif is a sanity check on the peak set (HOMER P-value < 1e-50). NEVER post-hoc filter individual peaks by DRACH content.
- Stop-codon enrichment in the Guitar metagene is the biological QC anchor (Dominissini 2012 *Nature* 485:201; Meyer 2012 *Cell* 149:1635). Absence means failed IP or non-m6A modification.
- 5'UTR peaks (within ~50 nt of TSS) are m6A-or-m6Am ambiguous because anti-m6A antibodies cross-react with PCIF1-deposited cap m6Am.
- N >= 3 biological replicates per condition; N=2 is under-powered for differential (McIntyre 2020 *Sci Rep* 10:6590).
- MACS3 `--keep-dup all` is mandatory; default `--keep-dup 1` collapses MeRIP signal.

## Related Skills

- epitranscriptomics/merip-preprocessing - Per-step preprocessing detail
- epitranscriptomics/m6a-peak-calling - exomePeak2 / MeTPeak / MACS3 deep treatment
- epitranscriptomics/m6a-differential - Differential methods and batch / lot handling
- epitranscriptomics/modification-visualization - Guitar metagene, heatmaps, browser tracks
- epitranscriptomics/m6anet-analysis - ONT DRS alternative for orthogonal stoichiometry
- chip-seq/peak-calling - Sibling IP-vs-input peak-calling framework
- chip-seq/chipseq-qc - IP enrichment QC concepts
- read-alignment/star-alignment - Splice-aware alignment
- workflow-management/snakemake-workflows - Snakemake orchestration
- workflow-management/nextflow-pipelines - Nextflow orchestration
- workflows/rnaseq-to-de - General RNA-seq pipeline patterns
