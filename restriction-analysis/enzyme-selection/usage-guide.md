# Enzyme Selection - Usage Guide

## Overview
Choose restriction enzymes for cloning or diagnostics by composing constraints: cut frequency, overhang type, recognition-site length, commercial availability, compatible ends, and methylation sensitivity. The skill helps pick enzymes to linearize a vector, drop in an insert directionally, set up a diagnostic digest, or avoid an enzyme silenced by Dam/Dcm methylation.

## Prerequisites
```bash
pip install biopython
```
Vector and/or insert sequences in FASTA or GenBank. Know which strain the DNA was prepped in (dam+ dcm+ is standard) when methylation matters.

## Quick Start
Tell your AI agent what you want to do:
- "Find enzymes that cut my plasmid exactly once for linearization"
- "Which commercial enzymes cut the vector once but not my insert?"
- "Find enzymes with sticky ends compatible with BamHI"
- "Which of these enzymes are blocked by Dam methylation?"

## Example Prompts

### Single-Cutters
> "Find all commercial enzymes that cut my plasmid exactly once for linearization"

> "Which single-cutters in pUC19 leave a 5' overhang?"

### Non-Cutters
> "Find enzymes that don't cut my insert but do cut the vector"

> "Which enzymes should I avoid because they cut inside my insert?"

### Compatible Ends
> "Find enzymes with sticky ends compatible with BamHI"

> "Which enzymes leave the same overhang as EcoRI?"

### Cloning Strategy
> "Find enzyme pairs for directional cloning of my insert into pET28a"

> "Suggest a single-cutter and a non-cutting partner for my construct"

### Methylation
> "Which of these enzymes are blocked by Dam methylation in DH5-alpha DNA?"

> "Find a methylation-insensitive alternative to MboI for GATC"

## What the Agent Will Do
1. Analyze the vector and/or insert sequences.
2. Search enzymes against the requested constraints (cut count, overhang, length, availability).
3. Restrict to `CommOnly` so recommendations can be purchased.
4. Check compatibility and flag methylation-blocked enzymes from the curated Dam/Dcm cases.
5. Recommend a concrete enzyme or pair, noting any star-activity / HF consideration.

## Code Patterns

### Single-Cutters (Linearization)
```python
from Bio import SeqIO
from Bio.Restriction import Analysis, CommOnly

record = SeqIO.read('plasmid.fasta', 'fasta')
analysis = Analysis(CommOnly, record.seq, linear=False)
single_cutters = analysis.with_N_sites(1)
print(f'Found {len(single_cutters)} single-cutters')
```

### Non-Cutters (Insert Protection)
```python
non_cutters = Analysis(CommOnly, insert_seq).without_site()
```

### Compatible Enzymes
```python
from Bio.Restriction import BamHI
partners = BamHI.compatible_end()                     # ANY enzyme that can leave a 5'-GATC end
fixed = [e for e in partners if e.is_palindromic()]   # keep Type IIP isocaudomers (BglII, BclI, MboI, Sau3AI)
```
`compatible_end()` also lists Type IIS enzymes (user-defined overhang); filter to fixed-overhang Type IIP enzymes before treating one as a drop-in ligation partner.

### Directional Cloning Selection
```python
vec_once = set(Analysis(CommOnly, vector_seq, linear=False).with_N_sites(1))
ins_none = set(Analysis(CommOnly, insert_seq).without_site())
candidates = vec_once & ins_none

five_prime = [e for e in candidates if e.is_5overhang()]
three_prime = [e for e in candidates if e.is_3overhang()]
blunt = [e for e in candidates if e.is_blunt()]
```
Note: `with_N_sites(n)` / `without_site()` / `with_sites()` are the current API. Older `once_cutters()` / `only_dont_cut()` raise `AttributeError`.

### Rare Cutters (8-base)
```python
eight_cutters = [e for e in CommOnly if len(e.site) == 8]
```

### Methylation Sensitivity
BioPython's `enzyme.is_methylable()` is a coarse REBASE flag (it returns True for the Dam-INSENSITIVE Sau3AI and for EcoRI), does not distinguish Dam/Dcm/CpG, and does not give the direction of the effect. Use the curated tables below and consult REBASE for the specific methyltransferase rather than relying on the flag.

### Golden Gate Compatibility (quick screen)
```python
from Bio.Restriction import BsaI
if not BsaI.search(insert_seq):
    print('No internal BsaI sites: part is Golden Gate ready for BsaI')
# Full Type IIS selection, domestication, and overhang design -> golden-gate-assembly
```

## Overhang Types

| Type | Example | Use Case |
|------|---------|----------|
| 5' overhang | EcoRI, BamHI | Most common cloning |
| 3' overhang | PstI, KpnI | Specific strategies (note PstI is a 3' overhang) |
| Blunt | EcoRV, SmaI | When no compatible sticky sites; inefficient, non-directional |

## Recognition Site Length

| Length | Naive frequency | Use |
|--------|-----------------|-----|
| 4 bp | ~256 bp | Frequent cutting (fingerprinting, RFLP) |
| 6 bp | ~4096 bp | Standard cloning |
| 8 bp | ~65536 bp | Rare cutting; large constructs and mapping |

Real genomes deviate: vertebrate CpG suppression makes CpG-containing sites (NotI especially) far rarer than the naive figure.

## Methylation Sensitivity (Curated)

Standard cloning strains are dam+ dcm+, so plasmid/insert DNA is methylated at GATC (Dam) and CCWGG (Dcm). Re-prep in a dam- dcm- strain (GM2163, JM110, INV110) to use a blocked enzyme.

| Methyltransferase | Site | Blocked | Insensitive | Requires methylation |
|-------------------|------|---------|-------------|----------------------|
| Dam (GATC, m6A) | GATC | MboI, DpnII | Sau3AI | DpnI (cuts only methylated) |
| Dcm (CCWGG, 5mC) | CCWGG | EcoRII | BstNI | - |
| CpG (mammalian, 5mC) | CCGG | HpaII | MspI | - |

## Type IIS Enzymes (Golden Gate)

| Enzyme | Recognition | Overhang |
|--------|-------------|----------|
| BsaI | GGTCTC | 4 bp 5' |
| BsmBI | CGTCTC | 4 bp 5' |
| BbsI | GAAGAC | 4 bp 5' |
| SapI | GCTCTTC | 3 bp 5' |

Selecting, domesticating, and designing fusion overhangs for these is its own analysis -> golden-gate-assembly.

## Tips
- Use 6-cutters for routine cloning, 8-cutters for large constructs, 4-cutters when frequent cuts are wanted.
- Check BOTH vector and insert: cut vector once, leave insert intact.
- Two different sticky ends give directional cloning and block vector self-ligation; dephosphorylate the vector (rSAP / Antarctic phosphatase, heat-inactivatable) when using a single compatible end.
- For DNA from standard strains, confirm the chosen enzyme is not Dam/Dcm-blocked.
- Prefer High-Fidelity (HF) enzymes when a clean overnight digest matters (no star activity).
- When picking a pair for a double digest, choose enzymes active together in one universal buffer (rCutSmart / FastDigest); if none exists, plan a sequential digest (lower-salt enzyme first) rather than forcing a single tube.
- Always restrict to `CommOnly` for practical, buyable recommendations.

## Related Skills

- restriction-sites - Find where the selected enzymes cut
- restriction-mapping - Map the selected enzyme sites
- fragment-analysis - Predict the fragments a chosen digest produces
- golden-gate-assembly - Select and validate Type IIS enzymes for scarless assembly
- primer-design/primer-basics - Add chosen restriction sites to primer tails
