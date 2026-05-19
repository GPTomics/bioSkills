# CLIP-seq Alignment - Usage Guide

## Overview

Align preprocessed CLIP-seq reads (eCLIP, iCLIP, iCLIP2, PAR-CLIP) to genome with crosslink-preserving parameters. The defining constraints: end-to-end alignment (no 5' soft-clip), strict mismatch ceiling (except PAR-CLIP T->C), and unique-mapper-only by default - with multi-mapper EM rescue (CLAM) reserved for RBPs that bind repetitive elements (Alu, LINE-1, LTR). STAR is the ENCODE standard; HISAT2 is the low-memory alternative; bowtie2 is acceptable only for non-spliced targets.

## Prerequisites

```bash
conda install -c bioconda star hisat2 bowtie2 samtools umi_tools
# CLAM for multi-mapper EM rescue (repeat-binding RBPs)
pip install CLAM
```

## Quick Start

Tell your AI agent what you want to do:
- "Align my eCLIP with the ENCODE STAR parameter block"
- "Raise mismatch ceiling for PAR-CLIP to keep T-to-C reads"
- "Use HISAT2 because my server only has 16 GB RAM"
- "My RBP binds LINE-1 repeats - run STAR with --outFilterMultimapNmax 100 and rescue with CLAM"
- "Filter MAPQ >= 10 and UMI-dedupe with method=unique"
- "Allele-specific binding study - add WASP filter to STAR"

## Example Prompts

### ENCODE-style eCLIP alignment

> "Align my eCLIP paired-end reads with STAR using --alignEndsType EndToEnd, --outFilterMultimapNmax 1, --outFilterMismatchNoverReadLmax 0.04"

> "Sort, index, and UMI-dedupe with `umi_tools dedup --method=unique --paired`"

### PAR-CLIP exception

> "Raise STAR mismatch ceiling to 0.07 so the T-to-C signal does not get filtered as sequencing error"

### Repeat-binding RBP rescue

> "MATR3 binds LINE-1 - emit up to 100 multi-mappers from STAR and feed to CLAM for EM-based reassignment"

> "How much does multi-mapper rescue increase my peak count for FUS? Compare unique-only vs CLAM"

### Low-memory alternative

> "STAR needs 30 GB; my machine has 16 - switch to HISAT2 with --no-softclip"

### Allele-specific

> "BEAPR-compatible alignment for allele-specific binding - add --waspOutputMode SAMtag with the heterozygous SNP VCF"

### Validation

> "Confirm my BAM did not soft-clip the truncation site - count reads with S in CIGAR"

> "Check that R2 5' positions cluster at expected truncation sites, not at uniform offsets"

## What the Agent Will Do

1. Pick STAR (default), HISAT2 (low memory), or bowtie2 (non-splice target) based on the scenario
2. Apply the ENCODE eCLIP parameter block: `--alignEndsType EndToEnd`, `--outFilterMultimapNmax 1`, `--outFilterMismatchNoverReadLmax 0.04`, `--outFilterType BySJout`, `--outFilterScoreMinOverLread 0.66`, `--outFilterMatchNminOverLread 0.66`
3. For PAR-CLIP: override mismatch ceiling to 0.07
4. For repeat-binding RBPs: raise multi-mapper ceiling to 100 and add CLAM EM rescue
5. For allele-specific: add WASP filter against het VCF
6. Index BAM, MAPQ filter (>=10 STAR, >=30 bowtie2/HISAT2), UMI-dedupe with `--method=unique`
7. Validate post-alignment: soft-clip rate < 1%, R2 5' positions at expected truncation sites

## Tips

- **End-to-end is sacred.** Soft-clip destroys the 5' truncation = CL site -1 base.
- **STAR default is `Local`.** Always explicitly set `--alignEndsType EndToEnd` for CLIP.
- **0.04 mismatch ceiling is for iCLIP/eCLIP only.** PAR-CLIP needs 0.07 to retain T->C signal.
- **Unique-mappers only by default.** Multi-mapper rescue is opt-in for repeat-binding RBPs.
- **CLAM is the EM solution.** Recovers 10-30% additional peaks in repeat-rich regions.
- **STAR memory hurts.** 30 GB for human. HISAT2 (8 GB) is the practical fallback.
- **bowtie2 misses splice junctions.** Use it for snoRNA/7SK targets only.
- **WASP filter is mandatory for ASB.** Reference-allele bias inflates REF allele frequency 1-5%.
- **MAPQ thresholds differ by aligner.** STAR 255 = unique; bowtie2 42; HISAT2 60. Use `-q 10` (STAR), `-q 30` (others) conservatively.
- **`umi_tools dedup --method=unique` is the ENCODE convention.** Directional is more conservative but slower.

## Related Skills

- clip-seq/clip-preprocessing - UMI extraction and adapter trimming
- clip-seq/clip-qc - Library complexity and read distribution checks
- clip-seq/crosslink-site-detection - Why 5' preservation matters
- clip-seq/clip-peak-calling - Peak calling on the dedup BAM
- read-alignment/star-alignment - General STAR usage
- read-alignment/bowtie2-alignment - General bowtie2 usage
- alignment-files/sam-bam-basics - BAM manipulation
