# Alignment Files

Skills for working with SAM/BAM/CRAM alignment files using samtools and pysam.

## Overview

This category covers the standard NGS alignment file workflow: viewing, sorting, indexing, filtering, and preparing data for variant calling.

**Tool type:** `cli` (with Python alternatives via pysam)
**Primary tools:** samtools, pysam

## Skills

| Skill | Description |
|-------|-------------|
| [sam-bam-basics](sam-bam-basics/) | View, convert SAM/BAM/CRAM; understand format structure |
| [alignment-indexing](alignment-indexing/) | Create BAI/CSI indices, enable random region access |
| [alignment-sorting](alignment-sorting/) | Sort by coordinate or name, merge BAM files, collate pairs |
| [duplicate-handling](duplicate-handling/) | Mark and remove PCR/optical duplicates |
| [alignment-statistics](alignment-statistics/) | Flagstat, depth, coverage, QC metrics |
| [alignment-filtering](alignment-filtering/) | Filter by flags, quality, regions |
| [reference-operations](reference-operations/) | Index FASTA, create dictionaries, generate consensus |
| [pileup-generation](pileup-generation/) | Generate pileup for variant calling |

## Workflow

```
FASTQ reads
    |
    v
[Aligner: bwa/bowtie2/STAR]
    |
    v
SAM/BAM (unsorted)
    |
    +---> [sam-bam-basics] - View, inspect, convert formats
    |
    v
[alignment-sorting] - Sort by coordinate
    |
    v
[alignment-indexing] - Create BAI index
    |
    +---> [alignment-statistics] - QC metrics (flagstat, stats)
    |
    v
[duplicate-handling] - Mark/remove PCR duplicates
    |
    +---> [alignment-filtering] - Filter by quality, flags, regions
    |
    v
[pileup-generation] - Generate pileup
    |
    v
[variant-calling] - bcftools call (next category)
```

## Example Prompts

### Viewing and Converting
- "View the first 100 alignments in my BAM file"
- "Convert this BAM file to CRAM format"
- "Show me the header of this BAM file"
- "How do I decode SAM FLAG 147?"

### Sorting, Merging, and Indexing
- "Sort this BAM file by coordinate"
- "Merge these BAM files into one"
- "Create an index for my BAM file"
- "Get reads from chromosome 1, positions 1000000-2000000"

### Quality Control
- "Get alignment statistics for this BAM file"
- "What is the mapping rate?"
- "Calculate coverage across my target regions"
- "What is the duplicate rate?"

### Filtering
- "Keep only properly paired reads"
- "Remove duplicates and low-quality alignments"
- "Extract reads from regions in my BED file"
- "Subsample to 10% of reads"

### Variant Calling Prep
- "Mark duplicates in my BAM file"
- "Generate pileup for variant calling"
- "Prepare this BAM for GATK"

## Installation

### samtools
```bash
# Conda
conda install -c bioconda samtools

# macOS
brew install samtools

# Ubuntu
sudo apt install samtools
```

### pysam
```bash
pip install pysam
```

## Related Skills

- **sequence-io** - FASTQ input files
- **variant-calling** - bcftools for VCF/BCF operations
- **database-access** - Download reference sequences from NCBI

## References

- [samtools documentation](http://www.htslib.org/doc/samtools.html)
- [pysam documentation](https://pysam.readthedocs.io/)
- [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
