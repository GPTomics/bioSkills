# Fragment Analysis - Usage Guide

## Overview
Predict and analyze DNA fragments produced by restriction enzyme digestion.

## Prerequisites
```bash
pip install biopython
```

## Quick Start
Tell your AI agent what you want to do:
- "Predict the fragment sizes from an EcoRI digest"
- "What bands will I see on a gel after cutting with BamHI?"

## Example Prompts

### Single Digest
> "What fragment sizes will I get from an EcoRI digest of my plasmid?"

> "Predict the gel pattern for HindIII digestion of sequence.fasta"

### Double Digest
> "Calculate fragment sizes for an EcoRI + BamHI double digest"

> "What bands will I see from digesting with both PstI and SalI?"

### Gel Comparison
> "Compare my predicted fragments to a 1kb ladder"

> "My gel shows bands at 3000, 2000, and 1000 bp - does this match EcoRI digestion?"

### Verification
> "Verify my digest worked by comparing observed vs expected fragments"

## What the Agent Will Do
1. Load your DNA sequence
2. Find enzyme cut positions
3. Calculate fragment sizes
4. Sort fragments for gel comparison
5. Optionally compare to observed gel results

## Code Patterns

### Basic Fragment Prediction
```python
from Bio import SeqIO
from Bio.Restriction import EcoRI

record = SeqIO.read('plasmid.fasta', 'fasta')
fragments = EcoRI.catalyze(record.seq)[0]
sizes = sorted([len(f) for f in fragments], reverse=True)
print(f'Fragment sizes: {sizes}')
```

### Understanding catalyze()
```python
# catalyze() returns a tuple
five_prime_frags, three_prime_frags = EcoRI.catalyze(seq)
# [0]: 5' fragments (most common use)
# [1]: 3' fragments (for asymmetric cuts)
```

### Linear vs Circular DNA

| DNA Type | n cuts | Fragments |
|----------|--------|-----------|
| Linear | n | n + 1 |
| Circular | n | n |

```python
# Plasmid (circular)
fragments = EcoRI.catalyze(seq, linear=False)[0]
```

### Double Digest
```python
from Bio.Restriction import EcoRI, BamHI

ecori_sites = EcoRI.search(seq)
bamhi_sites = BamHI.search(seq)
all_sites = sorted(set(ecori_sites + bamhi_sites))

def calc_fragments(seq_len, positions, linear=True):
    if not positions:
        return [seq_len]
    positions = sorted(positions)
    frags = []
    if linear:
        frags.append(positions[0])
        for i in range(len(positions) - 1):
            frags.append(positions[i + 1] - positions[i])
        frags.append(seq_len - positions[-1])
    else:
        for i in range(len(positions) - 1):
            frags.append(positions[i + 1] - positions[i])
        frags.append((seq_len - positions[-1]) + positions[0])
    return frags

sizes = calc_fragments(len(seq), all_sites, linear=True)
```

### Gel Simulation
```python
def gel_pattern(sizes, ladder=[10000, 5000, 3000, 2000, 1500, 1000, 500]):
    all_bands = sorted(set(sizes + ladder), reverse=True)
    for band in all_bands:
        marker = 'L' if band in ladder else ' '
        sample = '=' * (sizes.count(band) * 4) if band in sizes else ''
        print(f'{band:>6} {marker} | {sample}')
```

### Comparing Predicted vs Observed
```python
predicted = [3000, 2000, 1000]
observed = [3050, 1980, 1020]  # From gel image
tolerance = 100  # bp

for pred in predicted:
    matches = [obs for obs in observed if abs(pred - obs) <= tolerance]
    if matches:
        print(f'{pred} bp matches {matches[0]} bp')
```

## Tips
- Check linear vs circular setting if fragment count is wrong
- For linear DNA: sum of fragments should equal sequence length
- Small fragments (<100 bp) may run off the gel
- Allow ~5-10% tolerance when comparing to gel measurements
