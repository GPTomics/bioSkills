# Enzyme Selection Usage Guide

Choose restriction enzymes based on various criteria for cloning experiments.

## Prerequisites

```bash
pip install biopython
```

## Common Selection Criteria

### 1. Single-Cutters (Linearization)

Find enzymes that cut a plasmid exactly once:

```python
from Bio import SeqIO
from Bio.Restriction import Analysis, CommOnly

record = SeqIO.read('plasmid.fasta', 'fasta')
analysis = Analysis(CommOnly, record.seq, linear=False)

once_cutters = analysis.once_cutters()
print(f'Found {len(once_cutters)} single-cutters')
```

### 2. Non-Cutters (Insert Protection)

Find enzymes that don't cut your insert:

```python
analysis = Analysis(CommOnly, insert_seq)
non_cutters = analysis.only_dont_cut()
```

### 3. Compatible Enzymes

Find enzymes with compatible sticky ends:

```python
from Bio.Restriction import BamHI

compatible = BamHI.compatible_end()
# Returns enzymes that produce compatible overhangs
```

### 4. Rare Cutters

8-base recognition (cut rarely):

```python
eight_cutters = [e for e in CommOnly if len(e.site) == 8]
```

## Cloning Strategy Selection

### Directional Cloning

Need two different enzymes that:
- Cut vector once each
- Don't cut insert
- Produce incompatible ends (for directionality)

```python
# Find candidates
vec_once = set(Analysis(CommOnly, vector_seq).once_cutters().keys())
ins_none = set(Analysis(CommOnly, insert_seq).only_dont_cut())
candidates = vec_once & ins_none

# Pick one 5' and one 3' overhang enzyme
five_prime = [e for e in candidates if e.is_5overhang()]
three_prime = [e for e in candidates if e.is_3overhang()]
```

### Blunt-End Cloning

For blunt-end ligation:

```python
blunt = [e for e in candidates if e.is_blunt()]
```

## Overhang Types

| Type | Example | Use Case |
|------|---------|----------|
| 5' overhang | EcoRI, BamHI | Most common cloning |
| 3' overhang | PstI, KpnI | Specific strategies |
| Blunt | EcoRV, SmaI | When no compatible sites |

## Commercial Availability

Always check if enzyme is commercially available:

```python
from Bio.Restriction import CommOnly

# CommOnly contains only commercial enzymes
analysis = Analysis(CommOnly, seq)  # Use CommOnly, not AllEnzymes
```

## Recognition Site Length

| Length | Frequency | Use |
|--------|-----------|-----|
| 4 bp | ~256 bp | Frequent cutting |
| 6 bp | ~4096 bp | Standard cloning |
| 8 bp | ~65536 bp | Rare cutting |

## Methylation Sensitivity

When working with genomic DNA from bacteria (E. coli), methylation can block certain enzymes:

| Methylation | Enzymes Blocked | Alternative |
|-------------|-----------------|-------------|
| Dam (GATC) | MboI, DpnII, Sau3AI | Sau3AI (partially resistant) |
| Dcm (CCWGG) | BstNI, EcoRII | BsaWI |

```python
# Check methylation sensitivity
enzyme.is_dam_methylable()  # True if blocked by Dam
enzyme.is_dcm_methylable()  # True if blocked by Dcm
```

**Note:** DpnI is the opposite - it *requires* Dam methylation to cut.

## Golden Gate / Type IIS Cloning

Type IIS enzymes cut outside their recognition site, enabling scarless assembly:

| Enzyme | Recognition | Overhang |
|--------|-------------|----------|
| BsaI | GGTCTC | 4 bp |
| BsmBI | CGTCTC | 4 bp |
| BbsI | GAAGAC | 4 bp |
| SapI | GCTCTTC | 3 bp |

For Golden Gate cloning, your insert must be **free of the Type IIS site**:

```python
from Bio.Restriction import BsaI

sites = BsaI.search(insert_seq)
if not sites:
    print('Insert is Golden Gate compatible')
```

## Tips

- Use 6-cutters for routine cloning
- Use 8-cutters for large constructs
- Check both vector AND insert for cut sites
- Consider methylation sensitivity for genomic DNA
- For Golden Gate, verify insert lacks Type IIS sites
- Verify commercial availability before planning
