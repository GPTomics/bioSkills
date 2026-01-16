# Restriction Mapping Usage Guide

Create visual representations of restriction enzyme cut sites on DNA sequences.

## Prerequisites

```bash
pip install biopython
```

## Basic Map Generation

```python
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq

batch = RestrictionBatch([EcoRI, BamHI, HindIII])
analysis = Analysis(batch, seq)

# Print map
analysis.print_as('map')
```

## Output Formats

### Map Format
Visual representation with positions:
```
analysis.print_as('map')
```

### Linear Format
Simple list of sites:
```
analysis.print_as('linear')
```

### Tabular Format
Structured data:
```
analysis.print_as('tabulate')
```

## Calculating Fragment Distances

For cloning, you often need distances between sites:

```python
results = analysis.full()

# Collect all cut positions
all_positions = []
for enzyme, sites in results.items():
    for site in sites:
        all_positions.append((site, str(enzyme)))

# Sort by position
all_positions.sort()

# Calculate distances
for i in range(len(all_positions) - 1):
    pos1, enz1 = all_positions[i]
    pos2, enz2 = all_positions[i + 1]
    distance = pos2 - pos1
    print(f'{enz1}({pos1}) to {enz2}({pos2}): {distance} bp')
```

## Linear vs Circular DNA

For circular plasmids, fragment calculations must account for the wrap-around:

```python
# For circular DNA
analysis = Analysis(batch, seq, linear=False)

# Manual circular distance calculation
def circular_fragments(sites, seq_len):
    sites = sorted(sites)
    fragments = []
    for i in range(len(sites) - 1):
        fragments.append(sites[i + 1] - sites[i])
    # Wrap-around
    fragments.append((seq_len - sites[-1]) + sites[0])
    return fragments
```

## Integrating with GenBank Features

Find which features overlap with cut sites:

```python
from Bio import SeqIO

record = SeqIO.read('plasmid.gb', 'genbank')

for enzyme, sites in results.items():
    for site in sites:
        for feature in record.features:
            if feature.location.start <= site <= feature.location.end:
                print(f'{enzyme} at {site} is within {feature.type}')
```

## Exporting Maps

### To Text File

```python
with open('map.txt', 'w') as f:
    f.write(analysis.format_as('map'))
```

### To CSV

```python
import csv

with open('sites.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Enzyme', 'Position'])
    for enzyme, sites in results.items():
        for site in sites:
            writer.writerow([str(enzyme), site])
```

## Tips

- Use `linear=False` for plasmids
- Sort all sites for proper distance calculations
- Check overlaps with important features before cloning
- Consider using GenomeDiagram for publication-quality figures
