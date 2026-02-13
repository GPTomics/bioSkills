---
name: bio-clip-seq-clip-alignment
description: Align CLIP-seq reads to the genome with crosslink site awareness. Use when mapping preprocessed CLIP reads for peak calling.
tool_type: cli
primary_tool: STAR
---

## Version Compatibility

Reference examples tested with: Bowtie2 2.5.3+, STAR 2.7.11+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# CLIP-seq Alignment

**"Align my CLIP-seq reads to the genome"** → Map preprocessed CLIP reads with splice-aware alignment and crosslink site extraction for downstream peak calling.
- CLI: `STAR` with CLIP-optimized parameters (no multi-mappers, short reads)
- CLI: `bowtie2` for unspliced protocols

## STAR Alignment

**Goal:** Align CLIP-seq reads to the genome with splice awareness and strict uniqueness filtering.

**Approach:** Run STAR with single-mapping only, low mismatch tolerance, and end-to-end alignment to maximize crosslink site precision.

```bash
STAR --runMode alignReads \
    --genomeDir STAR_index \
    --readFilesIn trimmed.fq.gz \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 1 \
    --alignEndsType EndToEnd \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix clip_
```

## Bowtie2 Alternative

```bash
bowtie2 -x genome_index \
    -U trimmed.fq.gz \
    --very-sensitive \
    -p 8 \
    | samtools view -bS - \
    | samtools sort -o aligned.bam
```

## Post-Alignment Processing

**Goal:** Index aligned reads and remove PCR duplicates using UMI information.

**Approach:** Index the BAM with samtools and run umi_tools dedup to collapse UMI-duplicate reads.

```bash
# Index
samtools index aligned.bam

# Deduplicate with UMIs
umi_tools dedup \
    --stdin=aligned.bam \
    --stdout=deduped.bam
```

## Related Skills

- clip-preprocessing - Prepare reads
- clip-peak-calling - Call peaks
