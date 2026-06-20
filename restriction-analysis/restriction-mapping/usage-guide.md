# Restriction Mapping - Usage Guide

## Overview
Build restriction maps that place enzyme cut sites in order along a molecule, with the distances between them and the wrap-around fragment for plasmids. Produce text or graphical maps, overlay GenBank features, and order sites from single and double digests when the sequence is unknown.

## Prerequisites
```bash
pip install biopython matplotlib    # matplotlib only needed for graphical maps
```
A DNA sequence in FASTA or GenBank. Use GenBank for plasmids so topology and features are available.

## Quick Start
Tell your AI agent what you want to do:
- "Generate a restriction map of my plasmid with EcoRI, BamHI, and HindIII"
- "List cut sites in order with the distance between each"
- "Show which restriction sites fall inside annotated features"
- "Draw a graphical restriction map I can put in a figure"

## Example Prompts

### Basic Mapping
> "Create a restriction map for my plasmid with EcoRI, BamHI, and HindIII"

> "Print a visual map of cut sites in sequence.fasta"

### Distance Calculations
> "Calculate the distances between restriction sites in my plasmid, treating it as circular"

> "How far apart are the EcoRI and BamHI sites?"

### Feature Integration
> "Show which restriction sites overlap annotated features in my GenBank file"

### Export
> "Save the ordered restriction sites to a CSV file I name"

> "Draw a graphical map and save it to map.png"

## What the Agent Will Do
1. Load the DNA sequence and topology (linear or circular).
2. Find cut sites for the requested enzymes.
3. Order the sites and compute inter-site distances (closing the circle for plasmids).
4. Render a text or graphical map, optionally overlaying features.
5. Export to a caller-named file if asked.

## Code Patterns

### Basic Map Generation
```python
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII

record = SeqIO.read('sequence.fasta', 'fasta')
analysis = Analysis(RestrictionBatch([EcoRI, BamHI, HindIII]), record.seq)

analysis.print_as('map'); analysis.print_that()       # print_as SETS format; print_that() renders
analysis.print_as('linear'); analysis.print_that()
report = analysis.format_output()                     # capture as a string (not format_as)
```

### Calculating Distances (and closing a circle)
```python
cuts = sorted((pos, str(enz)) for enz, sites in analysis.full().items() for pos in sites)
seq_len = len(record.seq)
for i, (pos, enz) in enumerate(cuts):
    nxt = cuts[(i + 1) % len(cuts)][0]
    span = (nxt - pos) if nxt > pos else (seq_len - pos) + nxt   # wrap on circular DNA
    print(f'{enz}({pos}) -> next in {span} bp')
```

### Circular DNA Distance
```python
def circular_fragments(sites, seq_len):
    s = sorted(sites)
    spans = [s[i + 1] - s[i] for i in range(len(s) - 1)]
    return spans + [(seq_len - s[-1]) + s[0]]   # wrap-around fragment closes the circle
```

### Check Feature Overlaps
```python
record = SeqIO.read('plasmid.gb', 'genbank')
analysis = Analysis(RestrictionBatch([EcoRI, BamHI]), record.seq, linear=False)
for enzyme, sites in analysis.with_sites().items():
    for pos in sites:
        for feature in record.features:
            if int(feature.location.start) <= pos <= int(feature.location.end):
                print(f'{enzyme} at {pos} is within {feature.type}')
```

### Export to a Caller-Named CSV
```python
import csv
def export_sites(analysis, out_path):           # caller supplies out_path; do not hard-code
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['Enzyme', 'Position'])
        for enzyme, sites in analysis.full().items():
            for pos in sites:
                w.writerow([str(enzyme), pos])
```

### Ordering Sites From a Gel (Classical Mapping)
When the sequence is unknown, deduce order from fragment sizes: single digests give sizes only; comparing single vs double digests locates one enzyme's sites inside the other's fragments; partial-digest end labeling (Smith-Birnstiel) measures each site's distance from a labeled end directly. Use the sum-of-fragments check and multiple enzymes to resolve ambiguity.

## Tips
- Use `linear=False` for plasmids, and close the circle with `(seq_len - last) + first`.
- `print_as(...)` only sets the format; call `print_that()` to actually render.
- Capture text with `format_output()`; `format_as()` does not exist.
- Write graphical/CSV outputs only to a path the caller names, so scripts do not litter files.
- For richer annotated figures, see data-visualization/genome-tracks.

## Related Skills

- restriction-sites - Find the cut positions a map is built from
- fragment-analysis - Fragment sizes and gel interpretation behind classical mapping
- enzyme-selection - Choose informative enzymes for a map
- data-visualization/genome-tracks - Publication-quality annotated map figures
- sequence-io/read-sequences - Load FASTA or GenBank input
