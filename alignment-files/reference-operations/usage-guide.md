# Reference Operations - Usage Guide

## Overview
Work with reference genomes including indexing FASTA files, extracting sequences, creating sequence dictionaries, and generating consensus sequences from alignments.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "Index my reference FASTA file"
- "Extract sequence for chr1:1000-2000"
- "Create a sequence dictionary for GATK"
- "Generate consensus sequence from my alignments"

## Example Prompts

### FASTA Indexing
> "Index my reference genome for random access"

> "Prepare the reference for variant calling with GATK"

> "Create all required index files for my reference"

### Sequence Extraction
> "Extract the sequence for chr1:1000000-2000000"

> "Get multiple genomic regions and save to FASTA"

> "Extract the reverse complement of a region"

### Sequence Dictionary
> "Create a sequence dictionary for my reference"

> "Prepare reference files for Picard tools"

### Consensus Generation
> "Generate a consensus sequence from my BAM file"

> "Create consensus for a specific region with minimum depth"

> "Compare consensus to reference to find differences"

## What the Agent Will Do

1. Check if the reference file exists and is readable
2. Create the appropriate index files (FAI for samtools, dict for GATK/Picard)
3. Enable random access to any genomic region
4. Extract sequences or generate consensus as requested
5. Verify output files were created successfully

## Reference File Types

| File | Extension | Purpose |
|------|-----------|---------|
| FASTA | `.fa`, `.fasta` | Reference sequences |
| FAI Index | `.fa.fai` | Random access index |
| Dictionary | `.dict` | SAM header format (GATK/Picard) |

## Common Commands

### Create FASTA Index
```bash
samtools faidx reference.fa
```

Creates `reference.fa.fai` enabling:
- Random access to any region
- Quick chromosome length lookup
- Efficient memory usage

### Create Sequence Dictionary
```bash
samtools dict reference.fa -o reference.dict

# With metadata
samtools dict -a GRCh38 -s "Homo sapiens" reference.fa -o reference.dict
```

### Extract Sequences
```bash
samtools faidx reference.fa chr1:1000-2000           # Single region
samtools faidx reference.fa chr1:1000-2000 chr2:3000-4000 > regions.fa  # Multiple
samtools faidx reference.fa chr1 > chr1.fa            # Entire chromosome
samtools faidx -i reference.fa chr1:1000-2000        # Reverse complement
```

### Generate Consensus
```bash
samtools consensus input.bam -o consensus.fa                    # Basic
samtools consensus -r chr1:1000000-2000000 input.bam -o region.fa  # Specific region
samtools consensus -f fastq input.bam -o consensus.fq           # With quality
samtools consensus -d 10 input.bam -o consensus.fa              # Min depth 10
samtools consensus -a input.bam -o consensus.fa                 # Include N for no coverage
```

### Prepare Reference for Analysis
```bash
# Complete setup script
samtools faidx reference.fa
samtools dict reference.fa -o reference.dict
bwa index reference.fa  # If using BWA

# Get chromosome sizes for bedtools/UCSC
cut -f1,2 reference.fa.fai > chrom.sizes
```

### Create Subset Reference
```bash
samtools faidx reference.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
    chrX chrY chrM > main_chroms.fa
samtools faidx main_chroms.fa
samtools dict main_chroms.fa -o main_chroms.dict
```

## Python with pysam

### Fetch Sequences
```python
import pysam

with pysam.FastaFile('reference.fa') as ref:
    seq = ref.fetch('chr1', 999, 2000)  # 0-based, half-open
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

### Extract Multiple Regions
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
```

## FAI File Format

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

## Troubleshooting

### "faidx reference file not found"
Index the reference first:
```bash
samtools faidx reference.fa
```

### "invalid region" error
Check chromosome names match exactly:
```bash
grep "^>" reference.fa | head      # Check names in FASTA
cut -f1 reference.fa.fai           # Check names in index
# Use matching names: "chr1" not "1" if file uses "chr1"
```

### CRAM requires reference
Always specify reference for CRAM operations:
```bash
samtools view -T reference.fa input.cram
```

### Validate Reference Setup
```bash
REF=reference.fa
[ -f "$REF" ] && echo "FASTA: OK" || echo "FASTA: MISSING"
[ -f "${REF}.fai" ] && echo "FAI: OK" || echo "FAI: MISSING"
[ -f "${REF%.fa}.dict" ] && echo "DICT: OK" || echo "DICT: MISSING"
samtools faidx "$REF" chr1:1-100 > /dev/null && echo "Fetch: OK" || echo "Fetch: FAILED"
```

## Tips
- Always create both FAI index and dict file when setting up a new reference
- pysam uses 0-based coordinates while samtools uses 1-based
- Keep index files alongside the FASTA with the same prefix
- Use consensus with minimum depth to avoid calling bases from noise
- Chromosome names must match exactly between reference and BAM files
- For CRAM files, always provide the reference with `-T` flag
- Store chromosome sizes file for bedtools and UCSC tools
