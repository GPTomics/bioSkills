# Quality Reports - Usage Guide

## Overview
Quality reports are the first step in any NGS analysis. FastQC (or the faster drop-in falco) generates per-file reports of Phred quality, per-base composition, GC, duplication, overrepresented sequences, and adapter content; MultiQC aggregates them into one cross-sample summary. The skill teaches reading the PLOTS against the assay rather than trusting FastQC's pass/warn/fail, which are calibrated to whole-genome DNA and false-fail on RNA-seq, amplicon, bisulfite, small-RNA, single-cell, and long-read data.

## Prerequisites
```bash
# Conda (recommended)
conda install -c bioconda fastqc multiqc falco seqkit nanoplot

# pip (MultiQC, NanoPlot)
pip install multiqc nanoplot
```

## Quick Start
Tell your AI agent what you want to do:
- "Run FastQC on all my FASTQ files and aggregate with MultiQC"
- "Check the quality of my sequencing data before trimming"
- "My FastQC shows adapter contamination, what should I do?"
- "Is this duplication level a problem for RNA-seq?"
- "QC my Nanopore reads"

## Example Prompts

### Initial QC
> "Run FastQC on all FASTQ files in my data directory and build a MultiQC summary"

> "Check the quality of my raw reads before trimming and flag any outlier samples"

### Aggregating and comparing
> "Combine all my FastQC reports into one MultiQC report with custom sample names"

> "Compare QC before and after trimming in a single report"

### Interpreting results
> "My FastQC fails per-base sequence content in the first 12 bases, is that a problem?"

> "I see a 3' rise in G content on NovaSeq data, what is that?"

> "Explain the duplication plot in my MultiQC report for an RNA-seq library"

## What the Agent Will Do
1. Run FastQC (or falco) per file with appropriate threading, or NanoPlot for long reads
2. Aggregate the per-file reports with MultiQC and read the General Statistics table for outliers
3. Interpret each flagged module against the assay (distinguish expected false-fails from real defects)
4. Route remediation to the right skill (adapter-trimming, quality-filtering, fastp-workflow, contamination-screening)

## Tips
- Treat a FastQC FAIL as a hypothesis about a WGS library; for any other assay, first ask whether the module is expected to deviate for that chemistry.
- On 2-color instruments (NextSeq/NovaSeq), high-quality poly-G tails appear as a 3' G-content rise, NOT as low quality; remove them with a chemistry-aware poly-G trim, not a quality cutoff.
- Blocky/quantized quality boxes on NovaSeq are binned quality scores (expected), not bad data.
- High duplication is normal for RNA-seq, amplicon, and ChIP/ATAC; it is a prompt to assess library complexity (preseq), not an automatic reason to remove duplicates.
- The "Overrepresented sequences" module prints the offending sequence so it can be BLASTed to identify the source.
- MultiQC scrapes logs rather than re-analyzing data; check multiqc_sources.txt if two samples merge into one row.
- Run FastQC before and after trimming and overlay them in MultiQC to verify improvement.

## Resources
- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC Documentation](https://multiqc.info/)
- [falco](https://github.com/smithlabcode/falco)
- [NanoPlot](https://github.com/wdecoster/NanoPlot)

## Related Skills
read-qc/adapter-trimming - Remove read-through adapter flagged by the adapter-content panel
read-qc/quality-filtering - Drop low-quality reads and trim ends
read-qc/fastp-workflow - All-in-one QC + trim, including 2-color poly-G
read-qc/contamination-screening - Resolve bimodal-GC or unexpected overrepresented sequences
read-qc/rnaseq-qc - Transcriptome QC on the aligned BAM
sequence-io/sequence-statistics - Programmatic per-file sequence summaries
