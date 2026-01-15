# Pileup Generation Usage Guide

This guide covers generating pileup data for variant calling and position-level analysis.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools`)
- bcftools installed (`conda install -c bioconda bcftools`)
- pysam installed (`pip install pysam`)
- Indexed reference FASTA (`samtools faidx reference.fa`)

## What is Pileup?

A pileup represents all reads covering each position in the reference genome. At each position, you can see:

- The reference base
- Read depth (coverage)
- All bases present in reads
- Base qualities
- Strand information

This is the foundation for variant calling.

## Text Pileup Format

### Generating Text Pileup

```bash
samtools mpileup -f reference.fa input.bam > pileup.txt
```

### Understanding the Output

```
chr1    1000    A    15    ....,,..,..G...    FFFFFFFGFFFFFFFF
```

| Column | Value | Description |
|--------|-------|-------------|
| 1 | chr1 | Chromosome |
| 2 | 1000 | Position (1-based) |
| 3 | A | Reference base |
| 4 | 15 | Depth |
| 5 | ....,,..,..G... | Read bases |
| 6 | FFF... | Base qualities (ASCII - 33) |

### Read Bases Encoding

The pileup string shows what each read has at this position:

| Symbol | Meaning |
|--------|---------|
| `.` | Matches reference (forward strand) |
| `,` | Matches reference (reverse strand) |
| `A/C/G/T` | Mismatch (uppercase = forward strand) |
| `a/c/g/t` | Mismatch (lowercase = reverse strand) |
| `*` | Deletion |
| `+NNNbases` | Insertion of NNN bases |
| `-NNNbases` | Deletion of NNN bases |
| `^Q` | Start of read (Q is MAPQ as ASCII) |
| `$` | End of read |

### Example Interpretation

```
chr1    1000    A    20    ..,,.+2AT.,.G,..*..,
```

- 20 reads cover this position
- Most match reference A (`.` and `,`)
- One read has a 2bp insertion (`+2AT`)
- One read has G instead of A
- One read has a deletion (`*`)

## Quality Filtering

### Minimum Mapping Quality

Only include reads with MAPQ >= threshold:

```bash
samtools mpileup -f reference.fa -q 20 input.bam
```

### Minimum Base Quality

Only count bases with quality >= threshold:

```bash
samtools mpileup -f reference.fa -Q 20 input.bam
```

### Combined (Recommended)

```bash
samtools mpileup -f reference.fa -q 20 -Q 20 input.bam
```

## Region-Specific Pileup

### Single Region

```bash
samtools mpileup -f reference.fa -r chr1:1000000-2000000 input.bam
```

### From BED File

```bash
samtools mpileup -f reference.fa -l targets.bed input.bam
```

## Variant Calling Pipeline

The standard variant calling workflow pipes mpileup to bcftools:

### Simple Pipeline

```bash
samtools mpileup -f reference.fa input.bam | bcftools call -mv -o variants.vcf
```

### With Quality Filtering

```bash
samtools mpileup -f reference.fa -q 20 -Q 20 input.bam | \
    bcftools call -mv -Oz -o variants.vcf.gz

bcftools index variants.vcf.gz
```

### BCF Intermediate

For large datasets, BCF (binary VCF) is more efficient:

```bash
# Generate BCF
samtools mpileup -f reference.fa -g -o raw.bcf input.bam

# Call variants
bcftools call -mv raw.bcf -o variants.vcf
```

### Multi-Sample Calling

```bash
samtools mpileup -f reference.fa sample1.bam sample2.bam sample3.bam | \
    bcftools call -mv -o variants.vcf
```

## Working with pysam

### Basic Pileup Iteration

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000):
        pos = pileup_column.pos
        depth = pileup_column.n
        print(f'chr1:{pos+1} depth={depth}')
```

### Truncate to Exact Region

By default, pysam pileup shows all positions covered by reads that overlap the region. Use `truncate=True` for exact boundaries:

```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000, truncate=True):
        # Only positions 1000000-1000999
        pass
```

### Access Individual Reads

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1000001, truncate=True):
        print(f'Position {pileup_column.pos}, depth {pileup_column.n}')

        for pileup_read in pileup_column.pileups:
            # Check for deletion/skip
            if pileup_read.is_del:
                print(f'  {pileup_read.alignment.query_name}: DELETION')
                continue
            if pileup_read.is_refskip:
                print(f'  {pileup_read.alignment.query_name}: REFSKIP')
                continue

            # Get the base
            qpos = pileup_read.query_position
            base = pileup_read.alignment.query_sequence[qpos]
            qual = pileup_read.alignment.query_qualities[qpos]
            strand = '-' if pileup_read.alignment.is_reverse else '+'

            print(f'  {pileup_read.alignment.query_name}: {base} Q{qual} {strand}')
```

### Count Alleles

```python
import pysam
from collections import Counter

def count_alleles(bam_path, chrom, pos, min_qual=20):
    alleles = Counter()

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, pos, pos + 1,
                                         truncate=True,
                                         min_base_quality=min_qual):
            if pileup_column.pos != pos:
                continue

            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del:
                    alleles['DEL'] += 1
                elif pileup_read.is_refskip:
                    pass  # Ignore reference skips
                else:
                    qpos = pileup_read.query_position
                    base = pileup_read.alignment.query_sequence[qpos].upper()
                    alleles[base] += 1

    return dict(alleles)

# Example: Check for SNP at position
counts = count_alleles('input.bam', 'chr1', 1000000)
print(f'Allele counts: {counts}')

total = sum(counts.values())
for base, count in sorted(counts.items(), key=lambda x: -x[1]):
    freq = count / total * 100
    print(f'  {base}: {count} ({freq:.1f}%)')
```

### Find Variants

```python
import pysam
from collections import Counter

def find_variants(bam_path, ref_path, chrom, start, end, min_depth=10, min_alt_freq=0.1):
    variants = []

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        with pysam.FastaFile(ref_path) as ref:
            for pileup_column in bam.pileup(chrom, start, end,
                                             truncate=True,
                                             min_base_quality=20):
                pos = pileup_column.pos
                ref_base = ref.fetch(chrom, pos, pos + 1).upper()

                alleles = Counter()
                for pileup_read in pileup_column.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    qpos = pileup_read.query_position
                    base = pileup_read.alignment.query_sequence[qpos].upper()
                    alleles[base] += 1

                total = sum(alleles.values())
                if total < min_depth:
                    continue

                for base, count in alleles.items():
                    if base != ref_base:
                        freq = count / total
                        if freq >= min_alt_freq:
                            variants.append({
                                'chrom': chrom,
                                'pos': pos + 1,  # 1-based
                                'ref': ref_base,
                                'alt': base,
                                'depth': total,
                                'alt_count': count,
                                'freq': freq
                            })

    return variants

variants = find_variants('input.bam', 'reference.fa', 'chr1', 1000000, 1001000)
for v in variants:
    print(f'{v["chrom"]}:{v["pos"]} {v["ref"]}>{v["alt"]} '
          f'({v["alt_count"]}/{v["depth"]} = {v["freq"]:.1%})')
```

## Performance Considerations

### Max Depth

High-coverage regions can use excessive memory. Limit with `-d`:

```bash
samtools mpileup -f reference.fa -d 1000 input.bam
```

### Output Format

BCF is faster to process than text pileup:

```bash
# Faster: output BCF
samtools mpileup -f reference.fa -g input.bam -o output.bcf

# Slower: text pileup
samtools mpileup -f reference.fa input.bam > output.txt
```

### Parallel Processing

Process chromosomes in parallel:

```bash
for chr in chr1 chr2 chr3; do
    samtools mpileup -f reference.fa -r "$chr" input.bam | \
        bcftools call -mv -o "${chr}.vcf" &
done
wait
```

## Troubleshooting

### "No sequences in common"

Reference doesn't match the BAM:

```bash
# Check BAM chromosomes
samtools view -H input.bam | grep "^@SQ"

# Check reference chromosomes
grep "^>" reference.fa
```

### Empty Output

- Check region exists: `samtools view input.bam chr1:1000-2000 | head`
- Check reference is indexed: `ls reference.fa.fai`
- Check for chromosome name mismatch (chr1 vs 1)

### Memory Issues

Reduce max depth:

```bash
samtools mpileup -f reference.fa -d 500 input.bam
```

## See Also

- [samtools mpileup documentation](http://www.htslib.org/doc/samtools-mpileup.html)
- [bcftools call documentation](http://www.htslib.org/doc/bcftools.html#call)
- [pysam pileup documentation](https://pysam.readthedocs.io/)
