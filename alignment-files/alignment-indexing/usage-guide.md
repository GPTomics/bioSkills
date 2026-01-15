# Alignment Indexing Usage Guide

This guide covers creating and using indices for random access to BAM/CRAM files.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools`)
- pysam installed (`pip install pysam`)
- Coordinate-sorted BAM/CRAM file

## Why Index?

Without an index, accessing a specific genomic region requires reading the entire file. An index allows jumping directly to any region, making queries O(log n) instead of O(n).

## Index Types

### BAI (BAM Index)

The standard index format for BAM files.

```bash
samtools index input.bam
```

Creates `input.bam.bai` alongside the BAM file.

**Limitations**: Reference sequences must be < 512 Mbp (2^29 bp).

### CSI (Coordinate-Sorted Index)

For genomes with large chromosomes or when using custom bin sizes.

```bash
samtools index -c input.bam
```

Creates `input.bam.csi`.

**When to use**: Chromosomes > 512 Mbp, or when you need different binning.

### CRAI (CRAM Index)

Index for CRAM files (automatically created with correct extension).

```bash
samtools index input.cram
```

Creates `input.cram.crai`.

## Creating Indices

### Basic Indexing

```bash
# Index a BAM file
samtools index sample.bam

# Verify index was created
ls -la sample.bam*
```

### Multi-threaded Indexing

For large files, use multiple threads:

```bash
samtools index -@ 8 large.bam
```

### Batch Indexing

Index multiple files:

```bash
for bam in *.bam; do
    if [ ! -f "${bam}.bai" ]; then
        samtools index "$bam"
    fi
done
```

## Using Indices

### Region Queries with samtools

```bash
# Single region
samtools view input.bam chr1:1000000-2000000

# Multiple regions
samtools view input.bam chr1:1000-2000 chr2:3000-4000

# Count reads in region
samtools view -c input.bam chr1:1000000-2000000

# Regions from BED file
samtools view -L targets.bed input.bam
```

### Region Queries with pysam

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    # Fetch all reads overlapping region
    for read in bam.fetch('chr1', 1000000, 2000000):
        print(f'{read.query_name}: {read.reference_start}-{read.reference_end}')

    # Count reads in region
    count = bam.count('chr1', 1000000, 2000000)
    print(f'Total reads: {count}')
```

### Iterate by Chromosome

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for chrom in bam.references:
        count = bam.count(chrom)
        print(f'{chrom}: {count} reads')
```

## Index Statistics

### samtools idxstats

Quick per-chromosome read counts without reading alignment data:

```bash
samtools idxstats input.bam
```

Output:
```
chr1    248956422    5000000    1000
chr2    242193529    4500000    800
*       0            0          50000
```

Columns:
1. Reference name
2. Reference length
3. Mapped read count
4. Unmapped read count (placed on this reference)

The `*` row shows reads not mapped to any reference.

### Parse idxstats

```bash
# Total mapped reads
samtools idxstats input.bam | awk '{sum += $3} END {print sum}'

# Reads per chromosome as percentage
total=$(samtools idxstats input.bam | awk '{sum += $3} END {print sum}')
samtools idxstats input.bam | awk -v total="$total" '{if ($3 > 0) printf "%s\t%.2f%%\n", $1, $3/total*100}'
```

### pysam Index Statistics

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    stats = bam.get_index_statistics()
    total_mapped = sum(s.mapped for s in stats)

    for stat in stats:
        pct = stat.mapped / total_mapped * 100 if total_mapped > 0 else 0
        print(f'{stat.contig}: {stat.mapped:,} ({pct:.1f}%)')
```

## FASTA Indexing

Index reference FASTA for random access (separate from BAM indexing):

### Create FASTA Index

```bash
samtools faidx reference.fa
```

Creates `reference.fa.fai` with format:
```
chr1    248956422    6    60    61
chr2    242193529    253404903    60    61
```

Columns: name, length, offset, line bases, line width

### Fetch FASTA Region

```bash
samtools faidx reference.fa chr1:1000-2000
```

### pysam FastaFile

```python
import pysam

with pysam.FastaFile('reference.fa') as ref:
    # Fetch sequence (0-based coordinates)
    seq = ref.fetch('chr1', 999, 2000)  # Gets bases 1000-2000
    print(seq)

    # Get chromosome length
    length = ref.get_reference_length('chr1')
    print(f'chr1 length: {length}')
```

## Troubleshooting

### "random alignment retrieval only works for indexed BAM"

The BAM file lacks an index:

```bash
samtools index input.bam
```

### "file is not coordinate sorted"

BAM must be sorted by coordinate before indexing:

```bash
samtools sort -o sorted.bam unsorted.bam
samtools index sorted.bam
```

### "could not retrieve index file"

Index file is missing or in wrong location. Check:

```bash
ls -la input.bam*
# Should see input.bam and input.bam.bai (or input.bai)
```

### Chromosome name mismatch

Chromosome names must match exactly between query and BAM:

```bash
# Check what chromosomes are in the file
samtools view -H input.bam | grep "^@SQ"

# Common issue: "chr1" vs "1"
# If BAM uses "chr1", query must use "chr1"
samtools view input.bam chr1:1000-2000  # Not "1:1000-2000"
```

## Best Practices

1. **Always index after sorting**: Sort and index are typically done together
2. **Keep index with BAM**: Index file should stay alongside BAM
3. **Rebuild after modification**: Any BAM modification requires re-indexing
4. **Use CSI for large genomes**: If working with very large chromosomes
5. **Check before fetching**: Verify index exists before using fetch()

## See Also

- [samtools index documentation](http://www.htslib.org/doc/samtools-index.html)
- [samtools idxstats documentation](http://www.htslib.org/doc/samtools-idxstats.html)
- [pysam documentation](https://pysam.readthedocs.io/)
