# Methylation Calling - Usage Guide

## Overview
bismark_methylation_extractor processes Bismark BAM files to extract per-cytosine methylation information, producing various output formats for downstream analysis and visualization.

## Prerequisites
```bash
conda install -c bioconda bismark

bismark_methylation_extractor --version
```

## Quick Start
Tell your AI agent what you want to do:
- "Extract methylation calls from my Bismark BAM files"
- "Generate bedGraph files for IGV visualization"
- "Create coverage files for methylKit analysis"

## Example Prompts
### Basic Extraction
> "Run methylation extraction on my deduplicated Bismark BAM file"

> "Extract CpG methylation with bedGraph output for visualization"

### M-Bias Analysis
> "Check M-bias plots for my samples and recommend ignore parameters"

> "My M-bias shows end bias, help me rerun extraction with trimming"

### Output Formats
> "Generate a cytosine report for bsseq analysis"

> "Create coverage files compatible with methylKit"

## What the Agent Will Do
1. Run bismark_methylation_extractor with appropriate flags
2. Generate M-bias plots to check for positional bias
3. Review M-bias and recommend --ignore parameters if needed
4. Produce coverage files, bedGraph, and/or cytosine reports
5. Explain output files and their downstream uses

## Output Files

| File | Content | Downstream Use |
|------|---------|----------------|
| CpG_context_*.txt | Per-read CpG calls | Custom analysis |
| *.bismark.cov | Per-CpG summary | methylKit input |
| *.bedGraph | Methylation track | IGV/UCSC |
| *.CpG_report | All genome CpGs | bsseq input |

## Tips
- Always use --no_overlap for paired-end data to avoid double-counting
- Review M-bias plots before downstream analysis; ideal is a flat line around 70-80%
- Use --ignore and --ignore_3prime if M-bias shows sharp changes at read ends
- Include --cytosine_report and --genome_folder if you need all CpG positions (for bsseq)
- For memory issues, increase --buffer_size (e.g., --buffer_size 20G)
- Use --gzip to compress output files and save disk space
- Coverage files (*.bismark.cov) are the most common input for differential analysis
