# Methylation Calling - Usage Guide

## Overview

bismark_methylation_extractor processes Bismark BAM files to extract per-cytosine methylation information. It reads the XM tag containing methylation calls and produces various output formats for downstream analysis.

## When to Use This Skill

- You have aligned Bismark BAM files
- You need per-CpG methylation levels
- You want to generate bedGraph files for visualization
- You need coverage files for methylKit/bsseq

## Workflow

### 1. Run Methylation Extraction

```bash
bismark_methylation_extractor \
    --paired-end \
    --no_overlap \
    --gzip \
    --bedGraph \
    --cytosine_report \
    --genome_folder genome/ \
    -o methylation_calls/ \
    sample.deduplicated.bam
```

### 2. Check M-bias Plot

Review the M-bias plots to check for systematic bias at read ends:
```bash
# If bias found at positions 1-3:
bismark_methylation_extractor --ignore 3 --ignore_r2 3 ...
```

### 3. Use Output for Analysis

```bash
# Coverage file for methylKit
methylation_calls/sample.bismark.cov.gz

# bedGraph for visualization
methylation_calls/sample.bedGraph.gz

# Full report for bsseq
methylation_calls/sample.CpG_report.txt.gz
```

## Understanding M-Bias

The M-bias plot shows methylation level by read position. Ideal: flat line around 70-80%. Problems:
- Sharp increase/decrease at ends = adapter contamination or end-repair bias
- Solution: Use --ignore and --ignore_3prime parameters

## Output Files Explained

| File | Content | Downstream Use |
|------|---------|----------------|
| CpG_context_*.txt | Per-read CpG calls | Custom analysis |
| *.bismark.cov | Per-CpG summary | methylKit input |
| *.bedGraph | Methylation track | IGV/UCSC |
| *.CpG_report | All genome CpGs | bsseq input |

## Common Issues

### Memory Errors

```bash
# Increase buffer size
bismark_methylation_extractor --buffer_size 20G sample.bam
```

### Missing cytosine_report

Requires --genome_folder to know all CpG positions in genome.

## Resources

- [Bismark User Guide](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html)
