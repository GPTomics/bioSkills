# Golden Gate / Type IIS Assembly - Usage Guide

## Overview
Design and validate Type IIS scarless assembly (Golden Gate, MoClo, Golden Braid). The skill screens parts for internal BsaI/BsmBI/BbsI/SapI sites (domestication), previews the fusion overhang a cut exposes, and validates a fusion-overhang set so parts assemble in exactly one order. The principle is that a Type IIS enzyme cuts outside its recognition site, so the overhang is user-defined and the site is removed from the product, letting one tube digest and ligate at once.

## Prerequisites
```bash
pip install biopython
```
Part sequences in FASTA (or GenBank). For coding parts, know the reading frame so domestication can use synonymous codons.

## Quick Start
Tell your AI agent what you want to do:
- "Check my parts for internal BsaI sites before Golden Gate assembly"
- "Domesticate this CDS by removing its BsmBI site without changing the protein"
- "Preview the fusion overhangs my BsaI digest will leave"
- "Validate that my four junction overhangs assemble uniquely"

## Example Prompts

### Domestication
> "Scan part.fasta for internal BsaI, BsmBI, BbsI, and SapI sites"

> "Remove the internal BsaI site from my coding sequence using a silent mutation"

### Overhang Design and Validation
> "Are these four fusion overhangs distinct and non-palindromic: AATG, GCTT, TACT, AGGA?"

> "Preview the 4-bp overhangs a BsaI digest of my construct will expose"

### Planning
> "Which Type IIS enzyme should I use for MoClo Level 1 vs Level 2?"

> "Compare Golden Gate, classic restriction-ligation, and Gibson for assembling six parts"

## What the Agent Will Do
1. Pick the level's Type IIS enzyme (BsaI/BsmBI/BbsI/SapI).
2. Scan every part for internal sites on both strands and report any to domesticate.
3. Optionally rewrite a coding part with synonymous codons to remove a site.
4. Preview the fusion overhangs and validate the overhang set (distinct, non-palindromic, no cross-ligation).
5. Recommend the one-pot reaction setup (Type IIS enzyme + T4 ligase, 37/16 C cycling).

## Code Patterns

### Domestication Scan
```python
from Bio import SeqIO
from Bio.Restriction import BsaI, BsmBI, BbsI, SapI

record = SeqIO.read('part.fasta', 'fasta')
for enzyme in (BsaI, BsmBI, BbsI, SapI):
    hits = enzyme.search(record.seq)   # both strands
    print(enzyme, 'clean' if not hits else f'internal sites at {hits}')
```

### Preview Forward-Site Overhangs
```python
from Bio.Restriction import BsaI
overhangs = [str(record.seq[p - 1:p + 3]) for p in BsaI.search(record.seq)]   # forward sites
```

### Validate an Overhang Set
```python
from Bio.Seq import Seq

def validate(overhangs):
    bad = []
    if len(set(overhangs)) != len(overhangs):
        bad.append('duplicates')
    for o in overhangs:
        if o == str(Seq(o).reverse_complement()):
            bad.append(f'{o} palindromic')
    return bad or ['OK']
```

## Design Rules (Why They Matter)
- One distinct overhang per junction, or parts assemble in multiple orders.
- No palindromic overhang (it self-ligates); no overhang equal to the reverse complement of another (it cross-ligates).
- Avoid homopolymer overhangs (low ligation fidelity); prefer curated high-fidelity sets (Potapov 2018).
- Domesticate every part against the level's enzyme; orient the enzyme sites so cleavage removes them from the product.
- Alternate enzymes between MoClo levels (e.g. BsaI then BsmBI) so each round removes the previous level's sites.

## Tips
- Run the domestication scan on the full part in context; a site near a trimmed end can be missed because its cut falls off the sequence.
- SapI leaves a 3-nt (codon-length) overhang, useful for in-frame protein fusions.
- Prefer High-Fidelity variants (BsaI-HFv2, BsmBI-v2) for clean one-pot reactions.
- If parts cannot be domesticated, consider Gibson assembly instead (sequence-independent overlaps).

## Related Skills

- enzyme-selection - Choose a classic enzyme when scarless assembly is not needed
- restriction-sites - Find any enzyme's sites in a part
- fragment-analysis - Verify an assembly by predicted digest
- genome-engineering/grna-design - Design the constructs this assembly builds
- sequence-manipulation/transcription-translation - Confirm domestication preserved the protein
