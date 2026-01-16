# Fragment Analysis Usage Guide

Predict and analyze DNA fragments produced by restriction enzyme digestion.

## Prerequisites

```bash
pip install biopython
```

## Basic Fragment Prediction

```python
from Bio import SeqIO
from Bio.Restriction import EcoRI

record = SeqIO.read('plasmid.fasta', 'fasta')
seq = record.seq

# Get fragments
fragments = EcoRI.catalyze(seq)[0]

# Get sizes
sizes = sorted([len(f) for f in fragments], reverse=True)
print(f'Fragment sizes: {sizes}')
```

## Understanding catalyze()

The `catalyze()` method returns a tuple:
- `[0]`: 5' fragments (most common use)
- `[1]`: 3' fragments (for asymmetric cuts)

```python
five_prime_frags, three_prime_frags = EcoRI.catalyze(seq)
```

## Linear vs Circular DNA

| DNA Type | n cuts | Fragments |
|----------|--------|-----------|
| Linear | n | n + 1 |
| Circular | n | n |

```python
# Plasmid (circular)
fragments = EcoRI.catalyze(seq, linear=False)[0]
```

## Double Digest

For two enzymes, combine their cut positions:

```python
from Bio.Restriction import EcoRI, BamHI

# Get all cut positions
ecori_sites = EcoRI.search(seq)
bamhi_sites = BamHI.search(seq)
all_sites = sorted(set(ecori_sites + bamhi_sites))

# Calculate fragments from positions
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

## Gel Simulation

Match predicted fragments to a DNA ladder:

```python
def gel_pattern(sizes, ladder=[10000, 5000, 3000, 2000, 1500, 1000, 500]):
    all_bands = sorted(set(sizes + ladder), reverse=True)
    for band in all_bands:
        marker = 'L' if band in ladder else ' '
        sample = '=' * (sizes.count(band) * 4) if band in sizes else ''
        print(f'{band:>6} {marker} | {sample}')
```

## Comparing Results

When comparing predicted vs gel-observed sizes:

```python
predicted = [3000, 2000, 1000]
observed = [3050, 1980, 1020]  # From gel image

# Allow tolerance for gel measurement error
tolerance = 100  # bp

for pred in predicted:
    matches = [obs for obs in observed if abs(pred - obs) <= tolerance]
    if matches:
        print(f'{pred} bp matches {matches[0]} bp')
```

## Common Issues

**Wrong number of fragments:**
- Check linear vs circular setting
- Verify enzyme recognition site exists

**Fragments don't add up:**
- For linear DNA: sum should equal sequence length
- For circular: same, but n cuts give n fragments

**Missing small fragments:**
- Small fragments (<100 bp) may run off gel
- Include them in calculations
