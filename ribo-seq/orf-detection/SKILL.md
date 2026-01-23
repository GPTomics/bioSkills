---
name: bio-ribo-seq-orf-detection
description: Detect actively translated open reading frames from Ribo-seq data including canonical ORFs, uORFs, and novel ORFs. Use when identifying translated regions beyond annotated coding sequences.
tool_type: python
primary_tool: RiboCode
---

# ORF Detection

## RiboCode Workflow

```bash
# Step 1: Prepare annotation
prepare_transcripts \
    -g annotation.gtf \
    -f genome.fa \
    -o ribocode_annot

# Step 2: Run RiboCode
RiboCode \
    -a ribocode_annot \
    -c config.txt \
    -l 27,28,29,30 \
    -o output_prefix

# config.txt format:
# SampleName  AlignmentFile  Stranded
# sample1     sample1.bam    yes
```

## One-Step RiboCode

```bash
# All-in-one command
RiboCode_onestep \
    -g annotation.gtf \
    -r riboseq.bam \
    -f genome.fa \
    -l 27,28,29,30 \
    -o output_dir
```

## RiboCode Output

| File | Description |
|------|-------------|
| *_ORF_result.txt | Detected ORFs with coordinates |
| *_ORF_result.html | Interactive visualization |
| *_binomial_test.txt | Statistical test results |

## Parse RiboCode Results

```python
import pandas as pd

def load_ribocode_orfs(filepath):
    '''Load RiboCode ORF predictions'''
    df = pd.read_csv(filepath, sep='\t')

    # ORF categories
    categories = {
        'annotated': df[df['ORF_type'] == 'annotated'],
        'uORF': df[df['ORF_type'] == 'uORF'],
        'dORF': df[df['ORF_type'] == 'dORF'],
        'novel': df[df['ORF_type'].isin(['novel', 'noncoding'])]
    }

    return df, categories
```

## Alternative: RibORF

```bash
# RibORF uses random forest classifier
RibORF.py \
    -f genome.fa \
    -r riboseq.bam \
    -g annotation.gtf \
    -o output_dir
```

## Manual ORF Detection

```python
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence, min_length=30):
    '''Find all ORFs in a sequence'''
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']

    orfs = []
    seq = str(sequence).upper()

    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon == start_codon:
                # Find next stop codon
                for j in range(i + 3, len(seq) - 2, 3):
                    if seq[j:j+3] in stop_codons:
                        orf_length = j - i + 3
                        if orf_length >= min_length:
                            orfs.append({
                                'start': i,
                                'end': j + 3,
                                'frame': frame,
                                'length': orf_length,
                                'sequence': seq[i:j+3]
                            })
                        break

    return orfs

def detect_translated_orfs(orfs, coverage_data, min_coverage=10):
    '''Filter ORFs by Ribo-seq coverage'''
    translated = []
    for orf in orfs:
        cov = coverage_data[orf['start']:orf['end']]
        if sum(cov) >= min_coverage:
            translated.append(orf)
    return translated
```

## uORF Analysis

```python
def find_uorfs(transcript, cds_start):
    '''Find upstream ORFs before main CDS'''
    utr5 = transcript[:cds_start]
    uorfs = find_orfs(utr5)

    # Classify uORFs
    for uorf in uorfs:
        if uorf['end'] <= cds_start:
            uorf['type'] = 'contained'  # Fully in 5' UTR
        else:
            uorf['type'] = 'overlapping'  # Overlaps main CDS

    return uorfs
```

## ORF Categories

| Type | Description |
|------|-------------|
| annotated | Known CDS in annotation |
| uORF | Upstream of main CDS |
| dORF | Downstream of main CDS |
| internal | Within CDS, different frame |
| noncoding | In annotated non-coding RNA |
| novel | Unannotated region |

## Related Skills

- **ribosome-periodicity** - Validate ORF calling
- **translation-efficiency** - Quantify ORF translation
- **differential-expression** - Compare ORF expression
