# alignment-files

## Overview

Working with SAM/BAM/CRAM alignment files using samtools and pysam. Covers the standard NGS workflow: viewing, sorting, indexing, filtering, marking duplicates, and preparing data for variant calling.

**Tool type:** cli | **Primary tools:** samtools, pysam

## Skills

| Skill | Description |
|-------|-------------|
| sam-bam-basics | View, convert SAM/BAM/CRAM, FLAG semantics, MAPQ-by-aligner, CRAM REF_PATH |
| alignment-indexing | Create BAI/CSI indices, enable random region access, idxstats caveats |
| alignment-sorting | Sort by coordinate or name, collate vs sort -n, merge with header dedup |
| duplicate-handling | Assay-aware markdup, optical-distance per platform, UMI-aware alternatives |
| bam-statistics | flagstat / stats / mosdepth, assay-specific QC thresholds |
| alignment-validation | File integrity (quickcheck, ValidateSamFile, M5), insert-size by library, contamination |
| alignment-filtering | Aligner-aware MAPQ, expression filter, assay-aware variant-prep recipes |
| alignment-amplicon-clipping | Trim PCR primers from amplicon BAMs (samtools ampliconclip) |
| reference-operations | Index FASTA, sequence dict, consensus generation, contig-naming, GRCh38 flavors |
| pileup-generation | bcftools mpileup, BAQ semantics, library-typed flag recipes |

## Example Prompts

- "View the first 100 alignments in my BAM file"
- "Convert this BAM file to CRAM format with the right reference"
- "Show me the header and PG-chain of this BAM file"
- "Decode SAM FLAG 147"
- "Why is `-q 30` filtering out all my STAR-aligned reads?"
- "Sort this BAM file by coordinate"
- "Merge these BAM files; check sort-order consistency first"
- "Create a CSI index for my wheat genome BAM"
- "Get reads from chr1:1000000-2000000"
- "What is the mapping rate, and what does flagstat NOT tell me?"
- "Check insert size distribution for this ATAC-seq library"
- "Run a file-integrity check (quickcheck + ValidateSamFile)"
- "Cross-validate my BAM's @SQ M5 tags against the reference dict"
- "Estimate cross-sample contamination with verifybamid2 / somalier"
- "Calculate coverage across my target regions with mosdepth"
- "What duplicate rate is normal for an exome panel?"
- "Mark optical duplicates on a NovaSeq run"
- "Skip dedup -- this is bulk RNA-seq without UMIs"
- "Run UMI-aware dedup on a 10x Cell Ranger BAM"
- "Trim ARTIC SARS-CoV-2 primers from my amplicon BAM"
- "Filter for SV calling -- keep supplementary alignments"
- "Filter for HaplotypeCaller germline calling"
- "Subsample to 10% of reads, deterministically with a seed"
- "Generate a per-position pileup with BAQ off for long reads"
- "Modern variant calling pipeline (bcftools mpileup, not the deprecated samtools -g)"
- "Build a viral consensus FASTA with IUPAC codes"

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
