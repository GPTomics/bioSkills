# MeRIP-seq Preprocessing - Usage Guide

## Overview

End-to-end preprocessing of methylated-RNA-immunoprecipitation (MeRIP / m6A-seq) IP and matched input libraries. Covers adapter trimming with MeRIP-appropriate defaults, splice-aware genome alignment (STAR or HISAT2), the explicit do-NOT-deduplicate convention for non-UMI MeRIP, replicate-concordance QC via deepTools, IP enrichment QC via plotFingerprint, library-complexity saturation curves via PreSeq, antibody-lot metadata recording for cross-batch reconciliation, and IP-over-Input log2 bigWig generation for downstream visualisation. Produces the paired BAM files required by exomePeak2, MeTPeak, and MACS3.

## Prerequisites

```bash
conda install -c bioconda star hisat2 samtools deeptools preseq fastp trim-galore picard multiqc
```

For HISAT2-only environments:

```bash
conda install -c bioconda hisat2 samtools deeptools preseq fastp multiqc
```

Reference inputs:

- Genome FASTA (e.g., GRCh38 primary assembly)
- Matching GENCODE / Ensembl GTF (intronic / exonic annotation)
- Paired-end FASTQ files for each IP and matched input library
- Sample sheet recording: sample ID, condition, biological replicate, antibody clone, antibody lot

## Quick Start

- "Trim and align my MeRIP IP and input FASTQ pairs with STAR"
- "Align MeRIP data with HISAT2 on a memory-constrained node"
- "Check replicate concordance and IP enrichment across my MeRIP libraries"
- "Build a saturation curve so I can compare peak counts across conditions of different depth"
- "Generate IP-over-Input log2 bigWig tracks for downstream metagene plots"
- "Tell me which replicate has a failed IP before I run peak calling"

## Example Prompts

### Alignment

> "Align paired-end MeRIP-seq IP and input libraries to GRCh38 with STAR using splice-aware defaults; retain multi-mappers up to 20 per read; produce coordinate-sorted indexed BAMs ready for exomePeak2."

> "Use HISAT2 instead of STAR because the compute node has only 16 GB of memory; align IP and input libraries with --dta for downstream-compatibility; pipe through samtools sort."

### Trimming

> "Adapter-trim my MeRIP FASTQ pairs with fastp using --detect_adapter_for_pe; require minimum read length 25 nt; do NOT pass any UMI flags because the library is standard non-UMI MeRIP."

### Replicate Concordance

> "Compute Spearman correlation across all IP and input BAMs at 10 kb bins; render as a clustered heatmap so I can flag a divergent replicate before peak calling."

### IP Enrichment QC

> "Run deepTools plotFingerprint on all IP and input libraries; report JS distance and synthetic JS distance per sample; flag any IP whose Lorenz curve sits on the diagonal as a failed IP."

### Library Complexity

> "Run PreSeq lc_extrap on each sorted BAM; subsample all BAMs to the lowest common unique-read depth so peak counts are comparable across conditions."

### Downstream Tracks

> "Build per-replicate log2 (IP / Input) bigWig tracks at 25 bp bin resolution with pseudocount 1; CPM-normalise; output ready for downstream Guitar metagene and pyGenomeTracks browser plots."

### Failed IP Triage

> "Tell me which of my 6 IP replicates is most likely a failed IP given the fingerprint metrics and per-transcript IP/input ratio distributions."

### Dedup Decision

> "Confirm whether to deduplicate this MeRIP library; the protocol used standard non-UMI poly(A) selection followed by Synaptic Systems anti-m6A IP."

## What the Agent Will Do

1. Inspect FASTQ headers and metadata to confirm UMI status (almost always: no UMI for standard MeRIP)
2. Adapter-trim with fastp or Trim Galore; minimum length 25 nt
3. Build STAR or HISAT2 index against the matched genome and GTF (one-time)
4. Loop IP and input libraries through alignment with identical parameters
5. Sort coordinate, index, and run flagstat + idxstats per BAM
6. Run deepTools multiBamSummary + plotCorrelation for replicate concordance
7. Run deepTools plotFingerprint for IP enrichment QC; report JS distance per sample
8. Run PreSeq c_curve and lc_extrap for library complexity / saturation
9. Skip dedup (unless library is explicitly UMI-MeRIP, which is rare)
10. Generate per-replicate log2 (IP/Input) bigWig tracks for downstream visualisation
11. Aggregate all QC outputs into a MultiQC HTML report
12. Record antibody clone + lot in sample metadata for downstream design matrix
13. Recommend rarefaction depth for cross-condition peak-count comparison based on saturation curves

## Tips

- Standard MeRIP has NO UMI: do NOT deduplicate. The single most common preprocessing mistake in the field. dedup collapses real coverage at high-expression transcripts.
- IP and input MUST come from the same biological replicate; pair propagates through the downstream differential design matrix.
- STAR needs ~30 GB for the human index; HISAT2 needs ~12 GB. Choose by memory budget; both give comparable accuracy for MeRIP.
- Align to the GENOME for downstream MeRIP peak calling. Align to the TRANSCRIPTOME only for downstream m6Anet (ONT DRS).
- Peak counts are library-size-dependent. Never compare peak counts across libraries of different depth without rarefaction or saturation correction.
- Anti-m6A antibodies are polyclonal and lot-variable. Use one lot for an entire study. Record clone + lot in metadata.
- 3 biological replicates per condition is the practical minimum (McIntyre 2020 *Sci Rep* 10:6590). 4-5 preferred for differential downstream analysis.
- The NEB EpiMark kit includes Gluc (m6A-modified) and Cluc (unmodified) spike-in control RNAs. Use them for per-sample IP-efficiency QC; most users discard the controls.
- A failed IP can produce peaks. Always inspect plotFingerprint AND the per-transcript IP/input ratio distribution BEFORE peak calling.

## Related Skills

- m6a-peak-calling - Immediate downstream consumer of the BAM pairs produced here
- m6a-differential - Downstream differential analysis; uses IP/input pairing recorded here
- m6anet-analysis - ONT direct-RNA alternative; uses TRANSCRIPTOME alignment, NOT the genome BAMs from here
- modification-visualization - Uses the IP-over-Input bigWig tracks produced here
- read-qc/quality-reports - Upstream FastQC / MultiQC before trimming
- read-alignment/star-alignment - General STAR splice-aware alignment patterns
- read-alignment/hisat2-alignment - HISAT2 graph-based alternative
- alignment-files/sam-bam-basics - General BAM mechanics, samtools fundamentals
- alignment-files/duplicate-handling - General dedup philosophy (NOT applicable to non-UMI MeRIP)
- chip-seq/chipseq-qc - ChIP-seq IP QC concepts transfer directly to MeRIP
- workflows/rnaseq-to-de - End-to-end pipeline orchestration
