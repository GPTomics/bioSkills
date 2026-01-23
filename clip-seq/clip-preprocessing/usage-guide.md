# CLIP-seq Preprocessing - Usage Guide

## Overview

Preprocess CLIP-seq data with UMI extraction, adapter trimming, and PCR duplicate removal specific to CLIP protocols.

## Prerequisites

```bash
pip install umi_tools
conda install -c bioconda cutadapt
```

## Quick Start

Tell your AI agent:
- "Extract UMIs from my eCLIP data"
- "Trim adapters from CLIP reads"
- "Deduplicate using UMIs"
- "Process my iCLIP FASTQ files"

## Example Prompts

### UMI Handling

> "Extract 10-nt UMIs from read 1"

> "Deduplicate my aligned BAM using UMIs"

> "What is my library complexity?"

### Adapter Trimming

> "Trim Illumina adapters from my CLIP reads"

> "Run two-pass trimming for eCLIP protocol"

> "Remove inline barcodes"

### Quality Control

> "Check duplication rate in my CLIP library"

> "What percentage of reads have UMIs?"

> "Assess library complexity"

## What the Agent Will Do

1. Extract UMIs from read names (umi_tools extract)
2. Trim 3' and 5' adapters
3. Filter by minimum length
4. After alignment, deduplicate with UMIs
5. Report QC metrics

## Tips

- **UMI pattern** varies by protocol - check your library prep
- **eCLIP** typically has 10-nt UMIs in read 1
- **Deduplication** is critical - CLIP has high PCR duplication
- **Library complexity** should be > 0.3 for good libraries
- **Two-pass trimming** needed for read-through adapters
