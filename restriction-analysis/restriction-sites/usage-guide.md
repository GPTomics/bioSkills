# Finding Restriction Sites Usage Guide

Bio.Restriction allows you to search DNA sequences for restriction enzyme recognition sites.

## Prerequisites

```bash
pip install biopython
```

## Basic Workflow

### 1. Load Your Sequence

```python
from Bio import SeqIO

record = SeqIO.read('plasmid.fasta', 'fasta')
seq = record.seq
```

### 2. Search for Sites

```python
from Bio.Restriction import EcoRI

sites = EcoRI.search(seq)
print(f'EcoRI cuts at: {sites}')
```

### 3. Search Multiple Enzymes

```python
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII

batch = RestrictionBatch([EcoRI, BamHI, HindIII])
analysis = Analysis(batch, seq)
results = analysis.full()

for enzyme, sites in results.items():
    if sites:
        print(f'{enzyme}: {sites}')
```

## Understanding Cut Positions

Positions returned are **1-based** and indicate where the enzyme cuts:

```
EcoRI: G^AATTC
       |
       Cut position = 1 (after G)
```

For a sequence:
```
Position: 1234567890...
Sequence: ATGAATTCGC...
EcoRI cuts at position 4 (between G and A)
```

## Linear vs Circular

For circular DNA (plasmids), sites near the origin matter:

```python
# Linear DNA (default)
sites = EcoRI.search(seq, linear=True)

# Circular DNA
sites = EcoRI.search(seq, linear=False)
```

## Enzyme Properties

```python
from Bio.Restriction import EcoRI

# Recognition site
EcoRI.site        # 'GAATTC'

# Cut characteristics
EcoRI.is_blunt()      # False (makes sticky ends)
EcoRI.is_5overhang()  # True (5' overhang)
EcoRI.is_3overhang()  # False

# Overhang details
EcoRI.ovhg        # 4 (overhang length)
EcoRI.ovhgseq     # 'AATT' (overhang sequence)
```

## Common Enzyme Collections

| Collection | Description |
|------------|-------------|
| AllEnzymes | All ~800 enzymes in database |
| CommOnly | Commercially available only |
| Custom RestrictionBatch | Your selected enzymes |

## Filtering Results

```python
from Bio.Restriction import Analysis, CommOnly

analysis = Analysis(CommOnly, seq)

# Get specific categories
cutters = analysis.only_cut()         # All that cut
once = analysis.once_cutters()        # Cut exactly once
twice = analysis.twice_cutters()      # Cut exactly twice
none = analysis.only_dont_cut()       # Don't cut at all
```

## Common Issues

**Empty results:**
- Check sequence is DNA (not protein)
- Verify sequence has the recognition site
- Check linear/circular setting

**Module not found:**
- Install biopython: `pip install biopython`
- Import from Bio.Restriction, not Bio.restriction
