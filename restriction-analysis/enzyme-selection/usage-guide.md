# Enzyme Selection - Usage Guide

## Overview
Choose restriction enzymes based on various criteria for cloning experiments.

## Prerequisites
```bash
pip install biopython
```

## Quick Start
Tell your AI agent what you want to do:
- "Find enzymes that cut my plasmid exactly once"
- "Which enzymes don't cut my insert sequence?"

## Example Prompts

### Single-Cutters
> "Find all enzymes that cut my plasmid exactly once for linearization"

> "Which commercial enzymes are single-cutters in pUC19?"

### Non-Cutters
> "Find enzymes that don't cut my insert but do cut the vector"

> "Which enzymes should I avoid for my insert sequence?"

### Compatible Ends
> "Find enzymes with sticky ends compatible with BamHI"

> "Which enzymes produce the same overhang as EcoRI?"

### Cloning Strategy
> "Find enzyme pairs for directional cloning of my insert into pET28a"

> "Suggest enzymes for Golden Gate assembly of my construct"

### Methylation
> "Which enzymes are blocked by Dam methylation?"

> "Find methylation-insensitive alternatives to MboI"

## What the Agent Will Do
1. Analyze your vector and/or insert sequences
2. Search for enzymes matching your criteria
3. Filter by commercial availability
4. Check compatibility and methylation sensitivity
5. Recommend enzyme choices for your cloning strategy

## Code Patterns

### Single-Cutters (Linearization)
```python
from Bio import SeqIO
from Bio.Restriction import Analysis, CommOnly

record = SeqIO.read('plasmid.fasta', 'fasta')
analysis = Analysis(CommOnly, record.seq, linear=False)
once_cutters = analysis.once_cutters()
print(f'Found {len(once_cutters)} single-cutters')
```

### Non-Cutters (Insert Protection)
```python
analysis = Analysis(CommOnly, insert_seq)
non_cutters = analysis.only_dont_cut()
```

### Compatible Enzymes
```python
from Bio.Restriction import BamHI

compatible = BamHI.compatible_end()
```

### Directional Cloning Selection
```python
vec_once = set(Analysis(CommOnly, vector_seq).once_cutters().keys())
ins_none = set(Analysis(CommOnly, insert_seq).only_dont_cut())
candidates = vec_once & ins_none

five_prime = [e for e in candidates if e.is_5overhang()]
three_prime = [e for e in candidates if e.is_3overhang()]
blunt = [e for e in candidates if e.is_blunt()]
```

### Rare Cutters (8-base)
```python
eight_cutters = [e for e in CommOnly if len(e.site) == 8]
```

### Methylation Sensitivity
```python
enzyme.is_dam_methylable()  # True if blocked by Dam
enzyme.is_dcm_methylable()  # True if blocked by Dcm
```

### Golden Gate Compatibility
```python
from Bio.Restriction import BsaI

sites = BsaI.search(insert_seq)
if not sites:
    print('Insert is Golden Gate compatible')
```

## Overhang Types

| Type | Example | Use Case |
|------|---------|----------|
| 5' overhang | EcoRI, BamHI | Most common cloning |
| 3' overhang | PstI, KpnI | Specific strategies |
| Blunt | EcoRV, SmaI | When no compatible sites |

## Recognition Site Length

| Length | Frequency | Use |
|--------|-----------|-----|
| 4 bp | ~256 bp | Frequent cutting |
| 6 bp | ~4096 bp | Standard cloning |
| 8 bp | ~65536 bp | Rare cutting |

## Methylation Sensitivity

| Methylation | Enzymes Blocked | Alternative |
|-------------|-----------------|-------------|
| Dam (GATC) | MboI, DpnII, Sau3AI | Sau3AI (partially resistant) |
| Dcm (CCWGG) | BstNI, EcoRII | BsaWI |

Note: DpnI requires Dam methylation to cut.

## Type IIS Enzymes (Golden Gate)

| Enzyme | Recognition | Overhang |
|--------|-------------|----------|
| BsaI | GGTCTC | 4 bp |
| BsmBI | CGTCTC | 4 bp |
| BbsI | GAAGAC | 4 bp |
| SapI | GCTCTTC | 3 bp |

## Tips
- Use 6-cutters for routine cloning
- Use 8-cutters for large constructs
- Check both vector AND insert for cut sites
- Consider methylation sensitivity for genomic DNA
- For Golden Gate, verify insert lacks Type IIS sites
- Always use CommOnly for practical applications
