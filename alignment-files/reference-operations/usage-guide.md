# Reference Operations Usage Guide

This guide covers working with reference genomes and generating consensus sequences.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools`)
- pysam installed (`pip install pysam`)

## Reference File Types

| File | Extension | Purpose |
|------|-----------|---------|
| FASTA | `.fa`, `.fasta` | Reference sequences |
| FAI Index | `.fa.fai` | Random access index |
| Dictionary | `.dict` | SAM header format (GATK/Picard) |

## Indexing Reference FASTA

Before most operations, the reference needs to be indexed.

### Create FAI Index

```bash
samtools faidx reference.fa
```

This creates `reference.fa.fai` which enables:
- Random access to any region
- Quick chromosome length lookup
- Efficient memory usage (doesn't load entire file)

### FAI File Format

```
chr1    248956422    6    60    61
chr2    242193529    253404903    60    61
```

| Column | Description |
|--------|-------------|
| 1 | Sequence name |
| 2 | Sequence length |
| 3 | Byte offset of sequence start |
| 4 | Bases per line |
| 5 | Bytes per line (including newline) |

## Extracting Sequences

### Single Region

```bash
# 1-based coordinates
samtools faidx reference.fa chr1:1000-2000
```

Output:
```
>chr1:1000-2000
ACGTACGTACGT...
```

### Multiple Regions

```bash
samtools faidx reference.fa chr1:1000-2000 chr2:3000-4000 > regions.fa
```

### Entire Chromosome

```bash
samtools faidx reference.fa chr1 > chr1.fa

# Don't forget to re-index
samtools faidx chr1.fa
```

### Reverse Complement

```bash
samtools faidx -i reference.fa chr1:1000-2000
```

### Regions from File

```bash
# regions.txt: one region per line (chr:start-end)
while read region; do
    samtools faidx reference.fa "$region"
done < regions.txt > regions.fa
```

## Sequence Dictionary

A sequence dictionary is a SAM-format header containing reference information. Required by GATK and Picard tools.

### Create Dictionary

```bash
samtools dict reference.fa -o reference.dict
```

### With Metadata

```bash
samtools dict \
    -a GRCh38 \
    -s "Homo sapiens" \
    reference.fa -o reference.dict
```

### Dictionary Format

```
@HD VN:1.6 SO:unsorted
@SQ SN:chr1 LN:248956422 M5:6aef897c3d6ff0c78aff06ac189178dd UR:file:reference.fa
@SQ SN:chr2 LN:242193529 M5:f98db672eb0993dcfdabafe2a882905c UR:file:reference.fa
```

| Tag | Description |
|-----|-------------|
| SN | Sequence name |
| LN | Sequence length |
| M5 | MD5 checksum of sequence |
| UR | URI of sequence file |

## Generating Consensus

Create a consensus sequence from aligned reads.

### Basic Consensus

```bash
samtools consensus input.bam -o consensus.fa
```

### From Specific Region

```bash
samtools consensus -r chr1:1000000-2000000 input.bam -o region_consensus.fa
```

### Output as FASTQ

Include quality scores:

```bash
samtools consensus -f fastq input.bam -o consensus.fq
```

### Minimum Depth

Only call bases with sufficient coverage:

```bash
# Require at least 10 reads to call a base
samtools consensus -d 10 input.bam -o consensus.fa
```

### Call All Positions

Include positions with no coverage (as N):

```bash
samtools consensus -a input.bam -o consensus.fa
```

## Working with pysam

### Fetch Sequences

```python
import pysam

with pysam.FastaFile('reference.fa') as ref:
    # 0-based, half-open coordinates
    seq = ref.fetch('chr1', 999, 2000)  # Gets bases 1000-2000
    print(seq)
```

### Get Reference Info

```python
import pysam

with pysam.FastaFile('reference.fa') as ref:
    print(f'Chromosomes: {ref.nreferences}')

    for name in ref.references:
        length = ref.get_reference_length(name)
        print(f'{name}: {length:,} bp')
```

### Extract to FASTA

```python
import pysam

regions = [('chr1', 0, 10000), ('chr2', 5000, 15000)]

with pysam.FastaFile('reference.fa') as ref:
    for chrom, start, end in regions:
        seq = ref.fetch(chrom, start, end)
        print(f'>{chrom}:{start+1}-{end}')
        for i in range(0, len(seq), 60):
            print(seq[i:i+60])
```

### Build Simple Consensus

```python
import pysam
from collections import Counter

def simple_consensus(bam_path, chrom, start, end, min_depth=5):
    consensus = []

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup in bam.pileup(chrom, start, end, truncate=True):
            bases = Counter()
            for read in pileup.pileups:
                if read.is_del or read.is_refskip:
                    continue
                qpos = read.query_position
                if qpos is not None:
                    base = read.alignment.query_sequence[qpos]
                    bases[base.upper()] += 1

            total = sum(bases.values())
            if total >= min_depth:
                consensus.append(bases.most_common(1)[0][0])
            else:
                consensus.append('N')

    return ''.join(consensus)

seq = simple_consensus('input.bam', 'chr1', 1000000, 1001000)
print(f'>consensus')
print(seq)
```

### Compare to Reference

```python
import pysam

def compare_to_ref(bam_path, ref_path, chrom, start, end):
    consensus = simple_consensus(bam_path, chrom, start, end)

    with pysam.FastaFile(ref_path) as ref:
        reference = ref.fetch(chrom, start, end)

    differences = []
    for i, (c, r) in enumerate(zip(consensus, reference)):
        if c != r and c != 'N':
            differences.append((start + i, r, c))

    return differences

diffs = compare_to_ref('input.bam', 'reference.fa', 'chr1', 1000000, 1001000)
for pos, ref_base, cons_base in diffs:
    print(f'{pos}: {ref_base} -> {cons_base}')
```

## Common Workflows

### Prepare Reference for Analysis

```bash
#!/bin/bash
REF=$1

if [ ! -f "$REF" ]; then
    echo "Reference file not found: $REF"
    exit 1
fi

echo "Indexing for samtools..."
samtools faidx "$REF"

echo "Creating sequence dictionary..."
samtools dict "$REF" -o "${REF%.fa}.dict"

echo "Indexing for BWA..."
bwa index "$REF"

echo "Done. Files created:"
ls -la "${REF}"*
```

### Get Chromosome Sizes

```bash
# For bedtools, UCSC tools, etc.
cut -f1,2 reference.fa.fai > chrom.sizes
```

### Create Subset Reference

```bash
# Extract only main chromosomes
samtools faidx reference.fa \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
    chr21 chr22 chrX chrY chrM > main_chroms.fa

samtools faidx main_chroms.fa
samtools dict main_chroms.fa -o main_chroms.dict
```

### Validate Reference Setup

```bash
# Check all required files exist
REF=reference.fa

echo "Checking $REF..."
[ -f "$REF" ] && echo "  FASTA: OK" || echo "  FASTA: MISSING"
[ -f "${REF}.fai" ] && echo "  FAI: OK" || echo "  FAI: MISSING"
[ -f "${REF%.fa}.dict" ] && echo "  DICT: OK" || echo "  DICT: MISSING"

# Test fetch
echo "Testing fetch..."
samtools faidx "$REF" chr1:1-100 > /dev/null && echo "  Fetch: OK" || echo "  Fetch: FAILED"
```

## Troubleshooting

### "faidx reference file not found"

Index the reference first:

```bash
samtools faidx reference.fa
```

### "invalid region" error

Check chromosome names match exactly:

```bash
# See what chromosomes are in reference
grep "^>" reference.fa | head

# Or from index
cut -f1 reference.fa.fai
```

### CRAM requires reference

When working with CRAM files, always specify the reference:

```bash
samtools view -T reference.fa input.cram
```

## See Also

- [samtools faidx documentation](http://www.htslib.org/doc/samtools-faidx.html)
- [samtools dict documentation](http://www.htslib.org/doc/samtools-dict.html)
- [samtools consensus documentation](http://www.htslib.org/doc/samtools-consensus.html)
