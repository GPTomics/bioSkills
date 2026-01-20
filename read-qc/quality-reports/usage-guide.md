# Quality Reports - Usage Guide

## Overview
Quality reports are the first step in any NGS analysis. FastQC generates per-sample reports showing quality scores, adapter content, GC bias, and duplication levels. MultiQC aggregates multiple FastQC reports into a single interactive summary.

## Prerequisites
```bash
# Conda (recommended)
conda install -c bioconda fastqc multiqc

# pip (MultiQC only)
pip install multiqc
```

## Quick Start
Tell your AI agent what you want to do:
- "Run FastQC on all my FASTQ files"
- "Generate a MultiQC report from my QC results"
- "Check the quality of my sequencing data"

## Example Prompts

### Initial QC
> "Run FastQC on all FASTQ files in my data directory and create a summary report"

> "Check the quality of my raw reads before trimming"

### Aggregating Reports
> "Combine all my FastQC reports into a single MultiQC summary"

> "Generate a project-wide QC report with custom sample names"

### Interpreting Results
> "My FastQC shows adapter contamination, what should I do?"

> "Explain the quality metrics in my MultiQC report"

## What the Agent Will Do
1. Create an output directory for QC reports
2. Run FastQC on all FASTQ files with appropriate threading
3. Generate MultiQC summary from individual reports
4. Identify any quality issues (adapter contamination, low quality, duplication)

## Key Quality Metrics

| Metric | Good | Warning | Action |
|--------|------|---------|--------|
| Per base quality | >Q30 throughout | <Q20 at ends | Trim |
| Adapter content | <1% | >5% | Trim adapters |
| Duplication | <20% | >50% | Dedup |
| GC content | Normal curve | Secondary peak | Investigate |

## Tips
- Run FastQC both before and after trimming to verify improvement
- Use MultiQC config files to customize report titles and sample name cleaning
- Check "Overrepresented sequences" to identify unknown adapters
- High duplication may be normal for RNA-seq or amplicon data
- Secondary peaks in GC distribution often indicate contamination

## Resources
- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC Documentation](https://multiqc.info/)
- [FastQC Example Reports](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
