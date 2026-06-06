# Structural Variant Detection - Usage Guide

## Overview
Long-read structural variant calling detects deletions, insertions, inversions, duplications, and translocations from Oxford Nanopore and PacBio alignments - events that short reads largely miss because a 150 bp read cannot span an SV breakpoint or its repeat context. The central idea this skill teaches is that an SV call is a representation artifact as much as a biological fact: in tandem repeats the same event has many valid encodings, so the tandem-repeat BED supplied to the caller, the aligner, and the Truvari benchmark parameters determine precision and recall as much as the caller does. It covers Sniffles2 (the germline workhorse and its .snf cohort workflow), cuteSV's mandatory per-platform parameters, assembly-based calling, the somatic boundary to Severus/nanomonsv, and Truvari benchmarking against GIAB.

## Prerequisites
```bash
conda install -c bioconda sniffles cutesv svim minimap2 samtools truvari
# Reference-matched tandem-repeat BED (ships with Sniffles annotations/ for GRCh38/hs37d5)
# GIAB SV Tier1 / CMRG truth sets for benchmarking
```

## Quick Start
Tell your AI agent what you want to do:
- "Call structural variants from my Nanopore BAM with Sniffles2 and a tandem-repeat BED"
- "Joint-genotype SVs across my cohort"
- "Call SVs with cuteSV using the correct parameters for HiFi"
- "Benchmark my SV calls against GIAB HG002"

## Example Prompts

### Germline single sample
> "Call SVs from my ONT R10 BAM with Sniffles2. Supply the GRCh38 tandem-repeat BED and the reference so insertions carry sequence, and tell me why the TR BED matters."

### Cohort joint genotyping
> "I have 20 Nanopore samples. Build per-sample .snf signatures and merge them into one jointly-genotyped cohort VCF so absent samples still get real genotypes."

### Platform-correct cuteSV
> "Call SVs with cuteSV from my PacBio HiFi BAM using the HiFi-recommended cluster-bias and merge-ratio parameters, with genotyping enabled."

### Somatic
> "I have tumor-normal Nanopore data and want somatic SVs. Is Sniffles --mosaic appropriate, or should I use a paired caller?"

### Benchmarking
> "Benchmark my Sniffles VCF against GIAB Tier1 with Truvari, run truvari refine, and explain why the region set changes the F1."

## What the Agent Will Do
1. Map (or check the mapping) with the platform preset and `-Y` so supplementary alignments keep breakpoint sequence.
2. Call with a TR-aware caller, supplying `--tandem-repeats` and `--reference`.
3. For cohorts, build per-sample `.snf` and merge; for known panels, force-call with `--genotype-vcf`.
4. For cuteSV, apply the platform-matched parameter set and enable `--genotype`.
5. Route somatic tumor-normal work to Severus or nanomonsv.
6. Benchmark with Truvari (plus `refine`), reporting the region set, TR BED, and params.

## Tips
- Always supply a reference-matched tandem-repeat BED to the caller; it is the single biggest FP-reduction lever in repeats.
- Always pass `--reference` to Sniffles so insertions carry their sequence (needed for MEI classification and sequence-aware benchmarking).
- Map with minimap2 `-Y` (soft-clip supplementaries) - split-read callers reconstruct breakpoints from that sequence. NGMLR is a higher-precision/slower legacy niche.
- cuteSV parameters differ by platform (ONT vs HiFi vs CLR); the defaults are not platform-appropriate, and `--genotype` is off by default.
- The Sniffles `.snf` is a binary signature index, not a VCF; use it as Sniffles input.
- Sniffles `--mosaic` is single-sample low-VAF, not a tumor-normal caller; use Severus/nanomonsv for paired somatic.
- When benchmarking, state the region set, TR BED, and Truvari params, and run `truvari refine`; Tier1 F1 overstates whole-genome performance.

## Related Skills

- long-read-alignment - SV-ready mapping with `-Y`
- basecalling - Read accuracy/length affecting breakpoint precision
- clair3-variants - Small variants (<50 bp) belong to Clair3
- haplotype-phasing - Haplotag for phased / haplotype-specific SVs
- genome-assembly/hifi-assembly - Phased assembly for assembly-based SV calling
- variant-calling/structural-variant-calling - The variant-calling-side SV view
- variant-calling/vcf-manipulation - Filter and merge SV VCFs

## Resources
- [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)
- [cuteSV](https://github.com/tjiangHIT/cuteSV)
- [Truvari](https://github.com/ACEnglish/truvari)
- [GIAB benchmarks](https://www.nist.gov/programs-projects/genome-bottle)
