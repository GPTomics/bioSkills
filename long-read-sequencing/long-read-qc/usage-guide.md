# Long-Read Quality Control - Usage Guide

## Overview
Long-read QC assesses Oxford Nanopore and PacBio run quality and filters reads for the downstream goal. It differs fundamentally from short-read QC: there is no per-cycle quality plot, the headline metric is read N50 (length-weighted), and the reported per-read Qscore is an uncalibrated basecaller posterior - real accuracy (percent identity) requires aligning to a reference. Run-health metrics (pore activity, yield-over-time, translocation speed) come from the basecaller's sequencing_summary.txt, not the FASTQ. The right filter is intent-dependent: assembly preserves the long, low-quality reads and small replicons; variant calling filters almost nothing; HiFi is already Q20+ and is filtered on its rq tag. This skill covers NanoPlot, cramino, NanoComp, pycoQC/toulligQC, seqkit, chopper, and Filtlong, plus the chimera trap that fabricates SVs and the run-health red flags an expert reads.

## Prerequisites
```bash
conda install -c bioconda nanoplot nanocomp cramino chopper filtlong seqkit pycoqc
# toulligQC, fastcat, Porechop_ABI optionally for ONT-native QC and adapter discovery
# Keep the basecaller's sequencing_summary.txt for run-health QC
```

## Quick Start
Tell your AI agent what you want to do:
- "Run NanoPlot QC on my Nanopore reads and report the read N50"
- "Get the real percent identity of my reads against the reference"
- "Check my run health from the sequencing summary"
- "Filter my reads appropriately for assembly"

## Example Prompts

### Run overview and real accuracy
> "Give me an overview QC of my ONT FASTQ (length N50, yield), then align to the reference and tell me the real gap-compressed identity - I know the FASTQ Qscore overstates accuracy."

### Run health
> "From my sequencing_summary.txt, tell me whether the run had pore death, translocation-speed drift, or a high unclassified-barcode fraction."

### Intent-conditioned filtering
> "Filter my deep Nanopore data for a bacterial assembly. Subsample by quality to about 100x and make sure you do not erase small plasmids with a length cut."

### Comparing barcodes
> "Compare read length, quality, and identity across my four barcodes and flag any outlier."

### Chimera check
> "My SV calls have suspicious translocations. Check whether chimeric reads (internal adapters) could be the cause."

## What the Agent Will Do
1. Summarize length distribution, read N50, and yield from the FASTQ (NanoPlot/seqkit).
2. Align (or use an existing BAM) and report gap-compressed identity (cramino / NanoPlot --bam) as the real accuracy.
3. Read run-health from the sequencing_summary.txt (pycoQC/toulligQC): pore activity, yield-over-time, translocation speed, barcode breakdown.
4. Choose a filter conditioned on the downstream goal (assembly vs variant calling vs HiFi vs cDNA).
5. For assembly, subsample by quality with Filtlong rather than a hard length cut.
6. Flag run-health red flags and possible chimeras.

## Tips
- The FASTQ Qscore is an uncalibrated posterior; align to a reference for real percent identity.
- Always keep the sequencing_summary.txt; pycoQC/toulligQC need it and FASTQ-only hand-off loses run-health forever.
- For assembly, subsample by quality (Filtlong --target_bases) and never apply a hard length floor above your smallest replicon - the longest reads are the lowest-Q and small plasmids vanish.
- For variant calling, filter almost nothing; the caller models per-base Q and wants depth.
- HiFi is Q20+ already; filter on rq >= 0.99, not on Phred quality.
- Chimeras (internal adapters) masquerade as SVs; check whether Dorado already trimmed/split, and use Porechop_ABI for unknown adapters.
- NanoFilt and the rrwick Porechop are deprecated; use chopper and Porechop_ABI.

## Related Skills

- basecalling - Produces reads and the sequencing_summary.txt
- long-read-alignment - Produces the BAM for real percent identity
- structural-variants - Chimeras flagged here fabricate SVs
- genome-assembly/long-read-assembly - Subsample by quality before assembling
- read-qc/quality-reports - General short-read-oriented QC
- sequence-io/sequence-statistics - FASTA/FASTQ summary statistics

## Resources
- [NanoPack / NanoPlot](https://github.com/wdecoster/nanopack)
- [cramino](https://github.com/wdecoster/cramino)
- [pycoQC](https://github.com/a-slide/pycoQC)
- [Filtlong](https://github.com/rrwick/Filtlong)
