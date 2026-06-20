# Finding Restriction Sites - Usage Guide

## Overview
Search DNA sequences for restriction enzyme recognition sites with Bio.Restriction, for a single enzyme, a curated panel, or a whole commercial set, on linear or circular DNA. The result is a list of cut positions per enzyme, plus the ability to filter enzymes by how many times they cut.

## Prerequisites
```bash
pip install biopython
```
A DNA sequence in FASTA or GenBank format. Use GenBank when the molecule is a plasmid so the agent can read its circular topology and features.

## Quick Start
Tell your AI agent what you want to do:
- "Find all EcoRI sites in my plasmid sequence"
- "Search for EcoRI, BamHI, and HindIII sites in this DNA"
- "Which commercial enzymes cut my insert exactly once?"
- "Does NotI cut this sequence at all?"

## Example Prompts

### Single Enzyme Search
> "Find all EcoRI cut sites in plasmid.fasta"

> "Where does BamHI cut in my sequence, treating it as circular?"

### Multiple Enzyme Search
> "Search for EcoRI, BamHI, and HindIII sites in my plasmid"

> "Find all commercially available enzymes that cut my sequence exactly once"

### Filtering Results
> "Which enzymes cut my sequence twice?"

> "Find all enzymes that don't cut my insert sequence"

## What the Agent Will Do
1. Load the DNA sequence from file (FASTA or GenBank).
2. Choose enzyme scope: a named enzyme, a `RestrictionBatch`, `CommOnly`, or `AllEnzymes`.
3. Search with the correct topology (`linear=False` for plasmids).
4. Report cut positions (1-based, the base just 3' of the cut), and optionally filter enzymes by cut count.

## Code Patterns

### Basic Search
```python
from Bio import SeqIO
from Bio.Restriction import EcoRI

record = SeqIO.read('plasmid.fasta', 'fasta')
sites = EcoRI.search(record.seq, linear=False)   # circular plasmid
print(f'EcoRI cuts at: {sites}')
```

### Multiple Enzymes
```python
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII

batch = RestrictionBatch([EcoRI, BamHI, HindIII])
analysis = Analysis(batch, seq)

for enzyme, sites in analysis.with_sites().items():   # only enzymes that cut
    print(f'{enzyme}: {sites}')
```

### Filtering Results
```python
from Bio.Restriction import Analysis, CommOnly

analysis = Analysis(CommOnly, seq)
once = analysis.with_N_sites(1)   # cut exactly once (linearization)
twice = analysis.with_N_sites(2)  # cut exactly twice (excision)
none = analysis.without_site()    # do not cut at all (safe in a digest)
```
Note: the `with_N_sites(n)` / `with_sites()` / `without_site()` names are the current BioPython API. Older guides used `once_cutters()` / `twice_cutters()` / `only_dont_cut()`, which raise `AttributeError` in current versions.

## Understanding Cut Positions

A position is 1-based and equals the first base of the downstream fragment, i.e. the base immediately 3' of the cut on the top strand. It is NOT the start of the recognition site:
```
seq:    ...G A A T T C...      EcoRI site GAATTC starts at position p
cut:       ^                   top-strand cut is between G and A
search():  returns p+1         (the A after the cut), not p
```
Read the cut geometry directly with `elucidate()`:
```python
EcoRI.elucidate()   # 'G^AATT_C'  -> ^ top-strand cut, _ bottom-strand cut, 5' overhang AATT
```

## Linear vs Circular DNA
```python
sites = EcoRI.search(seq, linear=True)    # linear molecule (PCR product, genomic fragment)
sites = EcoRI.search(seq, linear=False)   # plasmid; finds sites spanning the origin
```

## Enzyme Properties
```python
EcoRI.site            # 'GAATTC' recognition site
EcoRI.is_blunt()      # False
EcoRI.is_5overhang()  # True
EcoRI.ovhg            # -4  (negative = 5' overhang, positive = 3', zero = blunt)
EcoRI.ovhgseq         # 'AATT'
EcoRI.is_ambiguous()  # False (True for N-spacer sites like BstXI; NOT for IUPAC sites like HincII)
```

## Common Enzyme Collections

| Collection | Description |
|------------|-------------|
| AllEnzymes | All known enzymes (1088 in current REBASE build), including non-commercial |
| CommOnly | Commercially available only (623); the practical default |
| Custom RestrictionBatch | A panel of enzymes chosen for the experiment |

## Tips
- Use `linear=False` for plasmid sequences, or sites near the origin are missed.
- Default to `CommOnly` so the answer is an enzyme that can actually be purchased.
- A position is the cut site, not the recognition-site start; do not subtract or add to "fix" it.
- Negative `ovhg` means a 5' overhang; confirm any end with `elucidate()`.
- Import from `Bio.Restriction` (capital R). Empty results usually mean the sequence is protein or the wrong topology.

## Related Skills

- restriction-mapping - Order cut sites and draw a map
- enzyme-selection - Choose enzymes by cut frequency, overhang, methylation, or compatible ends
- fragment-analysis - Turn cut positions into fragment sizes and gel patterns
- golden-gate-assembly - Screen a part for internal Type IIS sites
- sequence-io/read-sequences - Load the sequence to search
