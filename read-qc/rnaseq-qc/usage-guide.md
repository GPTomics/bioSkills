# RNA-seq QC - Usage Guide

## Overview
RNA-seq QC goes beyond general read quality to assess RNA-specific metrics including rRNA contamination, library strandedness, gene body coverage, and transcript integrity.

## Prerequisites
```bash
# RSeQC
pip install RSeQC

# SortMeRNA
conda install -c bioconda sortmerna

# Picard
conda install -c bioconda picard

# MultiQC
pip install multiqc
```

## Quick Start
Tell your AI agent what you want to do:
- "Check the RNA-seq specific quality metrics for my samples"
- "Verify the strandedness of my RNA-seq library"
- "Assess gene body coverage and RNA integrity"

## Example Prompts

### Strandedness
> "Determine the strandedness of my RNA-seq library using RSeQC"

> "Check if my library is stranded or unstranded before running alignment"

### rRNA Contamination
> "Measure the rRNA contamination level in my samples"

> "Filter rRNA reads from my RNA-seq data using SortMeRNA"

### RNA Integrity
> "Calculate gene body coverage to assess RNA degradation"

> "Compute TIN scores for all my samples to check RNA integrity"

### Comprehensive QC
> "Run full RNA-seq QC including strandedness, coverage, and rRNA metrics"

## What the Agent Will Do
1. Run RSeQC tools on aligned BAM files
2. Determine library strandedness for correct quantification settings
3. Calculate gene body coverage to detect 3'/5' bias
4. Compute TIN scores per transcript
5. Measure rRNA contamination levels
6. Aggregate results with MultiQC

## Key Metrics

| Metric | Good | Warning | Action |
|--------|------|---------|--------|
| rRNA contamination | <10% | >20% | Filter rRNA or re-prep |
| Gene body coverage | Even | 3' bias | Check RNA quality |
| Mean TIN | >70 | <50 | Exclude degraded samples |

## Tips
- Run strandedness check early; wrong settings cause major quantification errors
- 3' bias in gene body coverage indicates RNA degradation
- High rRNA despite depletion may require computational filtering with SortMeRNA
- Use salmon `-l A` to auto-detect strandedness if unsure
- TIN scores below 50 suggest significant degradation; consider excluding those samples

## Resources
- [RSeQC Documentation](http://rseqc.sourceforge.net/)
- [SortMeRNA GitHub](https://github.com/biocore/sortmerna)
- [Picard RNA Metrics](https://broadinstitute.github.io/picard/)
