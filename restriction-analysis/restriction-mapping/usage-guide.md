# Restriction Mapping - Usage Guide

## Overview
Create visual representations of restriction enzyme cut sites on DNA sequences.

## Prerequisites
```bash
pip install biopython
```

## Quick Start
Tell your AI agent what you want to do:
- "Generate a restriction map of my plasmid"
- "Show me where enzymes cut relative to each other"

## Example Prompts

### Basic Mapping
> "Create a restriction map for my plasmid with EcoRI, BamHI, and HindIII"

> "Print a visual map of cut sites in sequence.fasta"

### Distance Calculations
> "Calculate the distances between restriction sites in my plasmid"

> "How far apart are the EcoRI and BamHI sites?"

### Feature Integration
> "Show which restriction sites overlap with annotated features in my GenBank file"

### Export
> "Export the restriction map to a CSV file"

> "Save the restriction site positions to a text file"

## What the Agent Will Do
1. Load your DNA sequence
2. Search for specified enzyme cut sites
3. Generate visual map or calculate distances
4. Optionally check overlaps with sequence features
5. Export results if requested

## Code Patterns

### Basic Map Generation
```python
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII

record = SeqIO.read('sequence.fasta', 'fasta')
batch = RestrictionBatch([EcoRI, BamHI, HindIII])
analysis = Analysis(batch, record.seq)

analysis.print_as('map')      # Visual representation
analysis.print_as('linear')   # Simple list
analysis.print_as('tabulate') # Structured data
```

### Calculating Fragment Distances
```python
results = analysis.full()

all_positions = []
for enzyme, sites in results.items():
    for site in sites:
        all_positions.append((site, str(enzyme)))

all_positions.sort()

for i in range(len(all_positions) - 1):
    pos1, enz1 = all_positions[i]
    pos2, enz2 = all_positions[i + 1]
    print(f'{enz1}({pos1}) to {enz2}({pos2}): {pos2 - pos1} bp')
```

### Circular DNA Distance
```python
def circular_fragments(sites, seq_len):
    sites = sorted(sites)
    fragments = []
    for i in range(len(sites) - 1):
        fragments.append(sites[i + 1] - sites[i])
    fragments.append((seq_len - sites[-1]) + sites[0])
    return fragments
```

### Check Feature Overlaps
```python
record = SeqIO.read('plasmid.gb', 'genbank')

for enzyme, sites in results.items():
    for site in sites:
        for feature in record.features:
            if feature.location.start <= site <= feature.location.end:
                print(f'{enzyme} at {site} is within {feature.type}')
```

### Export to CSV
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
