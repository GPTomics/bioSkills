# alignment-files

## Overview

Working with SAM/BAM/CRAM alignment files using samtools and pysam. Covers the standard NGS workflow: viewing, sorting, indexing, filtering, marking duplicates, and preparing data for variant calling.

**Tool type:** cli | **Primary tools:** samtools, pysam

## Skills

| Skill | Description |
|-------|-------------|
| sam-bam-basics | View, convert SAM/BAM/CRAM, understand format structure |
| alignment-indexing | Create BAI/CSI indices, enable random region access |
| alignment-sorting | Sort by coordinate or name, merge BAM files, collate pairs |
| duplicate-handling | Mark and remove PCR/optical duplicates |
| bam-statistics | Flagstat, depth, coverage, QC metrics |
| alignment-validation | Insert size, proper pairing, strand balance, MAPQ distribution |
| alignment-filtering | Filter by flags, quality, regions, subsample reads |
| reference-operations | Index FASTA, create dictionaries, generate consensus |
| pileup-generation | Generate pileup for variant calling |

## Example Prompts

- "View the first 100 alignments in my BAM file"
- "Convert this BAM file to CRAM format"
- "Show me the header of this BAM file"
- "How do I decode SAM FLAG 147?"
- "Sort this BAM file by coordinate"
- "Merge these BAM files into one"
- "Create an index for my BAM file"
- "Get reads from chromosome 1, positions 1000000-2000000"
- "Get alignment statistics for this BAM file"
- "What is the mapping rate?"
- "Check insert size distribution"
- "Validate my alignment quality"
- "Check proper pairing rate"
- "Calculate coverage across my target regions"
- "What is the duplicate rate?"
- "Keep only properly paired reads"
- "Remove duplicates and low-quality alignments"
- "Extract reads from regions in my BED file"
- "Subsample to 10% of reads"
- "Mark duplicates in my BAM file"
- "Generate pileup for variant calling"

## Requirements

```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Related Skills

- **read-qc** - Quality control before alignment
- **variant-calling** - bcftools for VCF/BCF operations
- **genome-intervals** - BED file operations for region filtering
- **database-access** - Download reference sequences from NCBI
