# Finding Restriction Sites - Usage Guide

## Overview
Search DNA sequences for restriction enzyme recognition sites using Bio.Restriction.

## Prerequisites
```bash
pip install biopython
```

## Quick Start
Tell your AI agent what you want to do:
- "Find all EcoRI sites in my plasmid sequence"
- "Search for multiple restriction enzyme sites in this DNA sequence"

## Example Prompts

### Single Enzyme Search
> "Find all EcoRI cut sites in plasmid.fasta"

> "Where does BamHI cut in my sequence?"

### Multiple Enzyme Search
> "Search for EcoRI, BamHI, and HindIII sites in my plasmid"

> "Find all commercially available enzymes that cut my sequence exactly once"

### Filtering Results
> "Which enzymes cut my sequence twice?"

> "Find all enzymes that don't cut my insert sequence"

## What the Agent Will Do
1. Load your DNA sequence from file
2. Search for specified enzyme recognition sites
3. Report cut positions (1-based)
4. Optionally filter results by cut frequency

## Code Patterns

### Basic Search
```python
from Bio import SeqIO
from Bio.Restriction import EcoRI

record = SeqIO.read('plasmid.fasta', 'fasta')
sites = EcoRI.search(record.seq)
print(f'EcoRI cuts at: {sites}')
```

### Multiple Enzymes
```python
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII

batch = RestrictionBatch([EcoRI, BamHI, HindIII])
analysis = Analysis(batch, seq)
results = analysis.full()

for enzyme, sites in results.items():
    if sites:
        print(f'{enzyme}: {sites}')
```

### Filtering Results
```python
from Bio.Restriction import Analysis, CommOnly

analysis = Analysis(CommOnly, seq)
once = analysis.once_cutters()       # Cut exactly once
twice = analysis.twice_cutters()     # Cut exactly twice
none = analysis.only_dont_cut()      # Don't cut at all
```

## Understanding Cut Positions

Positions returned are 1-based and indicate where the enzyme cuts:
```
EcoRI: G^AATTC
       |
       Cut position = 1 (after G)
```

## Linear vs Circular DNA
```python
# Linear DNA (default)
sites = EcoRI.search(seq, linear=True)

# Circular DNA (plasmids)
sites = EcoRI.search(seq, linear=False)
```

## Enzyme Properties
```python
EcoRI.site            # 'GAATTC' - recognition site
EcoRI.is_blunt()      # False (makes sticky ends)
EcoRI.is_5overhang()  # True (5' overhang)
EcoRI.ovhgseq         # 'AATT' (overhang sequence)
```

## Common Enzyme Collections

| Collection | Description |
|------------|-------------|
| AllEnzymes | All ~800 enzymes in database |
| CommOnly | Commercially available only |
| Custom RestrictionBatch | Your selected enzymes |

## Tips
- Use `linear=False` for plasmid sequences
- Use CommOnly for practical cloning applications
- Check sequence is DNA (not protein) if getting empty results
- Import from Bio.Restriction (capital R)
