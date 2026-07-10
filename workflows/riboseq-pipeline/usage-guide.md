# Ribo-seq Pipeline - Usage Guide

## Overview

Complete workflow from Ribo-seq FASTQ through periodicity QC, P-site calibration, ORF detection, and differential translation efficiency. The pipeline gates each downstream analysis on library quality and records the two upstream decisions (harvest/drug and UMIs) that determine which conclusions are valid.

## Prerequisites

```bash
conda install -c bioconda cutadapt umi_tools bowtie2 star plastid samtools ribocode
```

```r
BiocManager::install(c("riboWaltz", "riborex", "DESeq2"))
```

## Quick Start

- "Analyze my Ribo-seq data from FASTQ to translation efficiency"
- "Run the ribosome profiling pipeline with periodicity QC"
- "Detect translated ORFs and differential TE from my Ribo-seq"

## Example Prompts

### Full Pipeline

> "Run the complete Ribo-seq pipeline from FASTQ"

> "Process my footprints and paired RNA-seq through to differential TE"

### Specific Steps

> "Just calibrate P-site offsets with riboWaltz"

> "Gate ORF detection on the periodicity QC result"

> "My library has UMIs - include extraction and deduplication"

## What the Agent Will Do

1. Extract UMIs, trim, deplete rRNA, align end-to-end, and deduplicate (UMIs only)
2. Run periodicity QC and calibrate per-length P-site offsets
3. Gate downstream analyses on the frame-0 fraction
4. Detect ORFs with RiboCode
5. Test differential TE with a count-based GLM
6. Optionally run stalling or initiation-site mapping

## Tips

- **Harvest gates dwell-time** - CHX pre-treatment invalidates stalling analysis
- **UMIs gate dedup** - deduplicate only with UMIs
- **End-to-end alignment** - soft-clipping corrupts P-site offsets
- **Periodicity is a hard gate** - no periodicity means gene-level counts only
- **RiboCode read lengths** come from metaplots, not -l
- **Differential TE** uses count GLMs (riborex/Xtail/anota2seq), not ratio tests

## Related Skills

- ribo-seq/riboseq-preprocessing - Preprocessing decisions
- ribo-seq/ribosome-periodicity - Periodicity QC and offsets
- ribo-seq/orf-detection - ORF calling
- ribo-seq/translation-efficiency - Differential TE
- ribo-seq/initiation-site-mapping - Start-codon mapping
- differential-expression/deseq2-basics - Count-based testing
