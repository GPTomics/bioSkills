# Fragment Analysis - Usage Guide

## Overview
Predict the DNA fragments a restriction digest produces -- sizes and sequences, for single or double digests on linear or circular DNA -- and interpret them against an agarose gel. Useful for planning a diagnostic digest to verify a clone and for matching observed bands to an expected pattern.

## Prerequisites
```bash
pip install biopython
```
A DNA sequence in FASTA or GenBank. For plasmids, prefer GenBank so topology is known.

## Quick Start
Tell your AI agent what you want to do:
- "Predict the fragment sizes from an EcoRI digest of my plasmid"
- "What bands will I see on a gel after cutting with BamHI?"
- "Calculate fragments for an EcoRI + BamHI double digest"
- "My gel shows 3000, 2000, 1000 bp -- does this match an EcoRI digest?"

## Example Prompts

### Single Digest
> "What fragment sizes will I get from an EcoRI digest of my plasmid?"

> "Predict the gel pattern for HindIII digestion of sequence.fasta"

### Double Digest
> "Calculate fragment sizes for an EcoRI + BamHI double digest"

> "What bands will I see from digesting with both PstI and SalI?"

### Gel Comparison
> "Compare my predicted fragments to a 1 kb ladder and tell me which bands co-migrate"

> "My gel shows bands at 3000, 2000, and 1000 bp - does this match EcoRI digestion?"

### Verification
> "Verify my digest worked by comparing observed vs expected fragments"

## What the Agent Will Do
1. Load the DNA sequence and determine topology (linear or circular).
2. Find cut positions and compute fragment sizes with the correct topology.
3. Run the sum check (fragments must total the molecule length).
4. Sort fragments and flag any too close to resolve on a gel.
5. Optionally compare predicted vs observed bands.

## Code Patterns

### Basic Fragment Prediction
```python
from Bio import SeqIO
from Bio.Restriction import EcoRI

record = SeqIO.read('plasmid.fasta', 'fasta')
fragments = EcoRI.catalyze(record.seq, linear=False)   # tuple of fragment Seqs; no [0]
sizes = sorted([len(f) for f in fragments], reverse=True)
print(f'Fragment sizes: {sizes}')
```

### Understanding catalyze()
```python
# catalyze() returns the tuple of fragments DIRECTLY:
#   EcoRI.catalyze(seq) -> (Seq(...), Seq(...), Seq(...))
# Do NOT write catalyze(seq)[0] expecting "the fragments" -- that is only the FIRST fragment,
# and iterating it then counts its bases and gives nonsense sizes.
fragments = EcoRI.catalyze(seq, linear=True)
sizes = [len(f) for f in fragments]
```

### Linear vs Circular DNA

| DNA Type | n cuts | Fragments |
|----------|--------|-----------|
| Linear | n | n + 1 |
| Circular | n | n |

```python
fragments = EcoRI.catalyze(seq, linear=False)   # plasmid -> n fragments from n cuts
```

### Double Digest
```python
from Bio.Restriction import EcoRI, BamHI, RestrictionBatch, Analysis

def calc_fragments(seq_len, positions, linear=True):
    cuts = sorted(set(positions))
    if not cuts:
        return [seq_len]
    spans = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
    if linear:
        return [cuts[0]] + spans + [seq_len - cuts[-1]]
    return spans + [(seq_len - cuts[-1]) + cuts[0]]

positions = [p for s in Analysis(RestrictionBatch([EcoRI, BamHI]), seq).with_sites().values() for p in s]
sizes = calc_fragments(len(seq), positions, linear=True)
```

### Gel Simulation
```python
def gel_pattern(sizes, ladder=(10000, 5000, 3000, 2000, 1500, 1000, 500)):
    for band in sorted(set(sizes) | set(ladder), reverse=True):
        marker = 'L' if band in ladder else ' '
        sample = '=' * (sizes.count(band) * 4)   # doubled width = co-migrating doublet
        print(f'{band:>6} {marker} | {sample}')
```

### Comparing Predicted vs Observed
```python
predicted = [3000, 2000, 1000]
observed = [3050, 1980, 1020]   # measured off the gel
tolerance = 100                 # bp; ~5-10% of band size is realistic

for pred in predicted:
    match = next((o for o in observed if abs(pred - o) <= tolerance), None)
    print(f'{pred} bp -> {"matches " + str(match) if match else "NOT seen (co-migration? partial digest?)"}')
```

## Reading The Gel (Physics, Not BioPython)
- Migration is ~linear in -log10(size) only within a gel's resolving window; size only against a ladder on the same gel.
- Agarose percent sets the window: ~0.7% for ~0.8-12 kb, ~1% for ~0.5-10 kb, ~2% for ~0.1-2 kb.
- Above ~20-50 kb, fragments co-migrate (reptation); resolving them needs pulsed-field (PFGE).
- Two similar sizes co-migrate as one brighter band; a "missing" band is often a doublet -- check the size sum.
- Uncut plasmid runs as supercoiled (fast) + nicked (slow): the same molecule, two bands. Only a linearized plasmid sizes correctly off a linear ladder.

## Tips
- Always run the sum check: for linear DNA the fragment sizes must total the sequence length.
- Set `linear=False` for plasmids or the fragment count is off by one.
- Small fragments (<100 bp) run off a low-percent gel; very large ones need PFGE.
- Extra, evenly spaced bands usually mean a partial digest -- drive the reaction to completion before blaming the map.

## Related Skills

- restriction-sites - Find the cut positions
- restriction-mapping - Order sites and compute distances
- enzyme-selection - Choose enzymes for a resolvable, diagnostic pattern
- sequence-manipulation/seq-objects - Work with the fragment Seq objects
