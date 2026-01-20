# Alignment Indexing - Usage Guide

## Overview
Create and use indices for random access to BAM and CRAM files, enabling fast region queries without reading entire files.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "Index my BAM file for random access"
- "Get per-chromosome read counts from the index"
- "Check if my BAM file has an index"
- "Create a CSI index for large chromosomes"

## Example Prompts

### Creating Indices
> "Create an index for sample.bam"

> "Index my BAM file using 8 threads"

> "Create CSI index for my genome with large chromosomes"

### Using Indices
> "Count reads per chromosome using idxstats"

> "Extract reads from chr1:1000000-2000000"

> "Get reads overlapping regions in my BED file"

### FASTA Indexing
> "Index my reference FASTA for random access"

> "Extract sequence for chr1:1000-2000 from reference"

### Index Statistics
> "Show per-chromosome read distribution"

> "Calculate mitochondrial contamination percentage"

> "Check X/Y ratio for sex determination"

## What the Agent Will Do

1. Verify the BAM file is coordinate-sorted (required for indexing)
2. Create the appropriate index type (BAI, CSI, or CRAI)
3. Place the index file alongside the BAM file
4. Verify the index was created successfully
5. Use the index for efficient region queries or statistics

## Index Types

### BAI (BAM Index)
Standard index for BAM files. Chromosomes must be < 512 Mbp.
```bash
samtools index input.bam  # Creates input.bam.bai
```

### CSI (Coordinate-Sorted Index)
For genomes with large chromosomes (> 512 Mbp).
```bash
samtools index -c input.bam  # Creates input.bam.csi
```

### CRAI (CRAM Index)
Index for CRAM files.
```bash
samtools index input.cram  # Creates input.cram.crai
```

## Common Commands

### Creating Indices
```bash
samtools index sample.bam           # Basic indexing
samtools index -@ 8 large.bam       # Multi-threaded
samtools index -c input.bam         # CSI for large chromosomes
```

### Batch Indexing
```bash
for bam in *.bam; do
    if [ ! -f "${bam}.bai" ]; then
        samtools index "$bam"
    fi
done
```

### Region Queries
```bash
samtools view input.bam chr1:1000000-2000000       # Single region
samtools view input.bam chr1:1000-2000 chr2:3000-4000  # Multiple regions
samtools view -L targets.bed input.bam             # From BED file
samtools view -c input.bam chr1:1000000-2000000    # Count in region
```

### Index Statistics
```bash
samtools idxstats input.bam  # Per-chromosome counts

# Total mapped reads
samtools idxstats input.bam | awk '{sum += $3} END {print sum}'

# Mitochondrial fraction
samtools idxstats input.bam | awk '/^chrM/ {mt=$3} {total+=$3} END {printf "MT: %.2f%%\n", mt/total*100}'
```

### FASTA Indexing
```bash
samtools faidx reference.fa                    # Create index
samtools faidx reference.fa chr1:1000-2000     # Extract region
```

## Python with pysam

### Region Queries
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000000, 2000000):
        print(f'{read.query_name}: {read.reference_start}-{read.reference_end}')
    count = bam.count('chr1', 1000000, 2000000)
    print(f'Total reads: {count}')
```

### Index Statistics
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    stats = bam.get_index_statistics()
    total_mapped = sum(s.mapped for s in stats)
    for stat in stats:
        pct = stat.mapped / total_mapped * 100 if total_mapped > 0 else 0
        print(f'{stat.contig}: {stat.mapped:,} ({pct:.1f}%)')
```

### FASTA Access
```python
import pysam

with pysam.FastaFile('reference.fa') as ref:
    seq = ref.fetch('chr1', 999, 2000)  # 0-based coordinates
    length = ref.get_reference_length('chr1')
```

## Troubleshooting

### "random alignment retrieval only works for indexed BAM"
BAM file lacks an index:
```bash
samtools index input.bam
```

### "file is not coordinate sorted"
BAM must be sorted before indexing:
```bash
samtools sort -o sorted.bam unsorted.bam
samtools index sorted.bam
```

### "could not retrieve index file"
Index file is missing or in wrong location:
```bash
ls -la input.bam*  # Should see input.bam and input.bam.bai
```

### Chromosome Name Mismatch
Names must match exactly between query and BAM:
```bash
samtools view -H input.bam | grep "^@SQ"  # Check names
# If BAM uses "chr1", query must use "chr1", not "1"
```

## Tips
- Always index after sorting - they go together as a pair
- Keep the index file alongside the BAM file with the same prefix
- Rebuild the index after any BAM modification
- Use CSI format for genomes with very large chromosomes
- Check that index exists before using fetch() in pysam
- Multi-threading (`-@ 8`) significantly speeds up large file indexing
- Use idxstats for quick per-chromosome summaries without reading alignments
