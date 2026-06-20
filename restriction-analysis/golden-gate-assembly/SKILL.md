---
name: bio-restriction-golden-gate-assembly
description: Design and validate Type IIS scarless DNA assembly (Golden Gate, MoClo) using Biopython Bio.Restriction. Screens parts for internal BsaI/BsmBI/BbsI/SapI sites (domestication), previews the fusion overhangs a digest exposes, and validates a fusion-overhang set for distinctness and fidelity. Use when designing a Golden Gate or MoClo assembly, domesticating a part by removing internal Type IIS sites, or choosing and checking fusion overhangs for one-pot assembly.
tool_type: python
primary_tool: Bio.Restriction
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+ (API verified on 1.86)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Restriction.BsaI.search)` to confirm the search API

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Golden Gate / Type IIS Assembly

**"Design (or check) my Golden Gate assembly"** -> Make each part free of the assembly enzyme's internal sites, and give every junction a distinct, well-behaved fusion overhang, so one tube of enzyme plus ligase builds the construct directionally and scarlessly.
- Python: `Bio.Restriction` to find internal Type IIS sites and read the overhang a cut exposes; the overhang-set rules are sequence logic, not a database call.

The whole method rests on one property: a Type IIS enzyme **cuts outside its recognition sequence**, so the 4-nt overhang it leaves is set by the user's flanking DNA, and the recognition site is placed to be **removed from the final product**. Because the ligated junction no longer contains the site, the enzyme cannot re-cut it, so digestion and ligation run together in one pot. Two design obligations follow, and both are what this skill checks: (1) **domestication** -- no part may contain an internal copy of the assembly enzyme's site, or it will be fragmented during assembly; (2) a set of **distinct, non-palindromic fusion overhangs** -- one per junction -- so parts assemble in exactly one order.

## Choosing The Assembly Enzyme

| Enzyme | Recognition | Overhang | Typical role |
|--------|-------------|----------|--------------|
| BsaI (Eco31I) | GGTCTC(1/5) | 4 nt 5' | The default Golden Gate / MoClo Level 1 enzyme; BsaI-HFv2 for fidelity |
| BsmBI (Esp3I) | CGTCTC(1/5) | 4 nt 5' | MoClo Level 0 / Level 2 (alternates with BsaI between levels) |
| BbsI (BpiI) | GAAGAC(2/6) | 4 nt 5' | Alternative when BsaI/BsmBI sites cannot be domesticated out |
| SapI (LguI) | GCTCTTC(1/4) | 3 nt 5' | 3-nt (codon-length) overhangs for reading-frame-preserving fusions |

Hierarchical systems (MoClo, Golden Braid) **alternate enzymes between levels** so each assembly round removes the previous level's sites: assemble Level 0 -> 1 with one enzyme, 1 -> 2 with the other. Pick the level's enzyme first, then domesticate every part against it.

## Method Context

| | Golden Gate (Type IIS) | Classic restriction-ligation (Type IIP) | Gibson assembly |
|--|------------------------|------------------------------------------|-----------------|
| Junction defined by | user-designed 4-nt overhang (3 for SapI) | the enzyme's fixed overhang | ~20-40 bp designed homology |
| Scar | none (site removed from product) | a restriction-site scar at each junction | none |
| Reaction | one-pot, one-step (37/16 C cycling) | sequential digest -> purify -> ligate | isothermal 50 C |
| Fragments per reaction | many (20-30+, more with optimized sets) | few | several |
| Main constraint | domestication; distinct overhang set | needs available compatible sites | terminal homology only |

Choose Golden Gate when assembling several parts repeatedly from a standardized library; classic digestion for a one-off two-piece clone with convenient sites; Gibson when parts cannot be domesticated or no site layout works.

## Domesticate: Find Internal Sites

**Goal:** Confirm a part carries no internal copy of the assembly enzyme's site (on either strand), and locate any that must be removed.

**Approach:** `enzyme.search(seq)` finds Type IIS sites on both strands (the recognition sequence is asymmetric, so its reverse-complement is detected too). Any hit inside a part is a defect to silently mutate away.

```python
from Bio import SeqIO
from Bio.Restriction import BsaI, BsmBI, BbsI, SapI

record = SeqIO.read('part.fasta', 'fasta')

for enzyme in (BsaI, BsmBI, BbsI, SapI):
    hits = enzyme.search(record.seq)            # both strands; recognition site is asymmetric
    status = 'clean' if not hits else f'{len(hits)} internal site(s) at {hits} -> domesticate'
    print(f'{enzyme} ({enzyme.site}): {status}')
```

## Domesticate: Break A Site By A Silent Mutation

**Goal:** Remove an internal site from a coding part without changing the protein.

**Approach:** Anchor on the recognition sequence itself (on both strands), not the cut position -- a Type IIS enzyme cuts *outside* its site, so mutating the codon at the cut would not touch the site. Walk the codons overlapping each recognition-site occurrence, swap one for a synonymous codon that breaks the site, and assert the protein is unchanged. This needs the reading frame.

```python
from Bio.Data import CodonTable
from Bio.Seq import Seq

def domesticate_cds(cds, enzyme, frame=0):
    '''Remove an enzyme's internal sites from a CDS by synonymous codon swaps (reading frame `frame`).'''
    table = CodonTable.unambiguous_dna_by_id[1]
    syn = {}
    for codon, aa in table.forward_table.items():
        syn.setdefault(aa, []).append(codon)
    s = list(str(cds).upper())
    end = frame + 3 * ((len(s) - frame) // 3)
    protein = str(Seq(''.join(s[frame:end])).translate())
    site = str(enzyme.site)
    motifs = (site, str(Seq(site).reverse_complement()))    # both strands; Type IIS sites are unambiguous
    for _ in range(len(s)):
        seqstr = ''.join(s)
        if not enzyme.search(Seq(seqstr)):
            break
        hit = max(seqstr.find(m) for m in motifs)           # a recognition-site start (either strand)
        for ci in range(((hit - frame) // 3) * 3 + frame, hit + len(site), 3):
            if ci < frame or ci + 3 > len(s):
                continue
            codon = ''.join(s[ci:ci + 3])
            alt = next((a for a in syn.get(table.forward_table.get(codon), ())
                        if a != codon and not enzyme.search(Seq(''.join(s[:ci] + list(a) + s[ci + 3:])))), None)
            if alt:
                s = s[:ci] + list(alt) + s[ci + 3:]
                break
    assert str(Seq(''.join(s[frame:end])).translate()) == protein   # protein unchanged
    return Seq(''.join(s))

# A site that overlaps only Met/Trp codons (no synonyms) is rare but cannot be silently broken;
# always confirm the result is clean: assert not enzyme.search(domesticate_cds(cds, enzyme))
```

## Preview The Fusion Overhangs A Digest Exposes

**Goal:** Read the actual 4-nt overhang each Type IIS cut would leave, to confirm junctions match as designed.

**Approach:** Anchor on the literal forward recognition sequence so only forward-oriented sites are read (a reverse-oriented site cuts on the other side and would otherwise return a misleading overhang). The 5' overhang starts `fst5` bases after the recognition-site start.

```python
from Bio.Restriction import BsaI

def forward_overhangs(seq, enzyme=BsaI, width=4):
    '''4-nt overhangs at FORWARD-oriented Type IIS sites, anchored on the recognition sequence.'''
    s, site, out = str(seq).upper(), str(enzyme.site), []
    i = s.find(site)
    while i >= 0:
        cut = i + enzyme.fst5                  # top-strand cut offset from the site start
        out.append(s[cut:cut + width])
        i = s.find(site, i + 1)
    return out

# Reverse-oriented sites cut on the other side; for a full construct, design the overhangs
# explicitly (below) rather than inferring every one from sequence.
```

## Validate A Fusion-Overhang Set

**Goal:** Check that the overhangs chosen for all junctions assemble uniquely and ligate efficiently.

**Approach:** Apply the design rules: every overhang distinct; none palindromic (self-ligates); no overhang equal to the reverse complement of another (cross-ligates); avoid all-identical bases. High-throughput ligation-fidelity data (Potapov 2018) underlies curated high-fidelity sets used for large assemblies.

```python
from Bio.Seq import Seq

def validate_overhang_set(overhangs):
    issues = []
    if len(set(overhangs)) != len(overhangs):
        issues.append('duplicate overhangs (parts assemble ambiguously)')
    for o in overhangs:
        if o == str(Seq(o).reverse_complement()):
            issues.append(f'{o} is palindromic (self-ligates)')
        if len(set(o)) == 1:
            issues.append(f'{o} is a homopolymer (low ligation fidelity)')
    rc = {o: str(Seq(o).reverse_complement()) for o in overhangs}
    for a in overhangs:
        for b in overhangs:
            if a != b and rc[a] == b:
                issues.append(f'{a} is the reverse complement of {b} (cross-ligates)')
    return issues or ['overhang set OK']

print(validate_overhang_set(['AATG', 'GCTT', 'TACT', 'GGGG']))
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Assembly drops or scrambles a part | An internal Type IIS site fragmented it | Domesticate every part against the level's enzyme before assembly |
| Junctions ligate in the wrong order or orientation | Two junctions share an overhang, or one is the reverse complement of another | Use a distinct, non-self-complementary overhang per junction; check with `validate_overhang_set` |
| Empty or low-efficiency assembly | Palindromic or homopolymer overhang, or wrong enzyme/buffer cycling | Avoid palindromic/homopolymer overhangs; cycle 37/16 C with a Type IIS enzyme + T4 ligase |
| `search()` misses a reverse-oriented site | Site too close to the sequence end so the cut falls off it | Domesticate on the full part in context, not a trimmed fragment |
| Recognition site still present in the product | Site placed so the cut does not remove it | Orient Type IIS sites so cleavage excises them from the assembled junction |

## Related Skills

- enzyme-selection - Choose a classic restriction enzyme when scarless assembly is not needed
- restriction-sites - Find any enzyme's sites in a part
- fragment-analysis - Predict fragments to verify an assembly digest
- genome-engineering/grna-design - Design constructs that this assembly will build
- sequence-manipulation/transcription-translation - Confirm domestication kept the reading frame

## References

- Engler C, Kandzia R, Marillonnet S. A one pot, one step, precision cloning method with high throughput capability. PLoS One. 2008;3(11):e3647. doi:10.1371/journal.pone.0003647
- Engler C, Gruetzner R, Kandzia R, Marillonnet S. Golden gate shuffling: a one-pot DNA shuffling method based on type IIs restriction enzymes. PLoS One. 2009;4(5):e5553. doi:10.1371/journal.pone.0005553
- Weber E, Engler C, Gruetzner R, Werner S, Marillonnet S. A modular cloning system for standardized assembly of multigene constructs. PLoS One. 2011;6(2):e16765. doi:10.1371/journal.pone.0016765
- Potapov V, Ong JL, Kucera RB, et al. Comprehensive profiling of four base overhang ligation fidelity by T4 DNA ligase and application to DNA assembly. ACS Synth Biol. 2018;7(11):2665-2674. doi:10.1021/acssynbio.8b00333
