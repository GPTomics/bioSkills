---
name: bio-restriction-mapping
description: Build restriction maps showing enzyme cut positions and inter-site distances along DNA using Biopython Bio.Restriction. Produces text or graphical maps for linear and circular molecules, orders sites from single and double digests, and overlays GenBank features. Use when creating a restriction map of a sequence, ordering cut sites along a plasmid, or relating sites to annotated features.
tool_type: python
primary_tool: Bio.Restriction
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+ (API verified on 1.86), matplotlib 3.7+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Restriction.Analysis.print_as)` to confirm format names

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Restriction Mapping

**"Make a restriction map of my sequence"** -> Place each enzyme's cut sites along the molecule, in order, with the distances between them and (for plasmids) the wrap-around fragment.
- Python: `Bio.Restriction.Analysis(...).print_as('map')` for a quick text map; `search()` positions + `matplotlib` for a graphical one.

A map is more than a list of positions: it is the ordering and spacing of sites, and on a plasmid the ordering is circular. Two things separate a correct map from a wrong one. First, **a circular molecule wraps**: the fragment between the last site and the first crosses the origin, so its length is `(seq_len - last) + first`, not `seq_len - last`. Second, when sites come from a gel rather than a known sequence, **order is deduced, not given** -- single digests give sizes, and only comparing single vs double digests (or partial digests) orders them.

## Choosing The Map Representation

| Need | Representation | How |
|------|----------------|-----|
| Quick look while exploring | Text map | `Analysis.print_as('map')` / `'linear'` |
| Capture to a string/report | Formatted text | `Analysis.format_output()` |
| Publication / slide figure | Graphical map | `search()` positions drawn with matplotlib |
| Sites vs annotated features | Feature overlay | iterate `record.features` against each cut position |
| Ordering sites from a gel | Digest comparison | single vs double (or partial) digest fragment patterns |

## Quick Text Map

```python
from Bio import SeqIO
from Bio.Restriction import EcoRI, BamHI, HindIII, RestrictionBatch, Analysis

record = SeqIO.read('sequence.fasta', 'fasta')
analysis = Analysis(RestrictionBatch([EcoRI, BamHI, HindIII]), record.seq)

analysis.print_as('map'); analysis.print_that()      # visual map to stdout
analysis.print_as('linear'); analysis.print_that()   # linear list
report = analysis.format_output()                    # capture as a string (not format_as)
```

## Ordered Site List With Distances

**Goal:** A single ordered table of every cut, which enzyme made it, and the distance to the next.

**Approach:** Collect `(position, enzyme)` from `Analysis.full()`, sort by position, and walk the list. For circular DNA, close the loop with the wrap-around span.

```python
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII, XhoI, NotI

seq = record.seq
seq_len = len(seq)
circular = False    # set True for a plasmid (and use linear=not circular below)

analysis = Analysis(RestrictionBatch([EcoRI, BamHI, HindIII, XhoI, NotI]), seq, linear=not circular)
cuts = sorted((pos, str(enz)) for enz, sites in analysis.full().items() for pos in sites)

for i, (pos, enz) in enumerate(cuts):
    nxt = cuts[(i + 1) % len(cuts)][0]
    span = (nxt - pos) if nxt > pos else (seq_len - pos) + nxt   # wrap on circular
    last = (i == len(cuts) - 1)
    dist = span if (circular or not last) else seq_len - pos
    print(f'{pos:6d} bp ({pos / seq_len:5.1%}) {enz:8s} -> next in {dist} bp')
```

## Graphical Map (matplotlib)

**Goal:** A figure with the molecule as an axis and a labeled tick per cut site.

**Approach:** Draw the backbone, place a vertical tick at each `search()` position, and stack enzymes on separate rows. Write the figure only to a path the caller names (so running this does not litter the working directory).

```python
import matplotlib
matplotlib.use('Agg')                  # headless; no display needed
import matplotlib.pyplot as plt
from Bio.Restriction import EcoRI, BamHI, HindIII

def draw_map(seq, enzymes, out_path):
    seq_len = len(seq)
    fig, ax = plt.subplots(figsize=(10, 2 + 0.4 * len(enzymes)))
    ax.hlines(0, 0, seq_len, color='black')
    for row, enz in enumerate(enzymes, start=1):
        for pos in enz.search(seq):
            ax.vlines(pos, row - 0.3, row + 0.3, color='C0')
            ax.text(pos, row + 0.35, str(pos), ha='center', va='bottom', fontsize=7)
        ax.text(-0.02 * seq_len, row, str(enz), ha='right', va='center')
    ax.set_xlim(0, seq_len); ax.set_yticks([]); ax.set_xlabel('position (bp)')
    fig.savefig(out_path, dpi=200, bbox_inches='tight'); plt.close(fig)

# draw_map(record.seq, [EcoRI, BamHI, HindIII], 'my_map.png')   # caller supplies the path
```

## Map Against GenBank Features

```python
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI

record = SeqIO.read('plasmid.gb', 'genbank')
analysis = Analysis(RestrictionBatch([EcoRI, BamHI]), record.seq, linear=False)

for enzyme, sites in analysis.with_sites().items():
    for pos in sites:
        hits = [f.qualifiers.get('label', f.qualifiers.get('gene', [f.type]))[0]
                for f in record.features
                if int(f.location.start) <= pos <= int(f.location.end)]
        print(f'{enzyme} at {pos}: {", ".join(hits) or "intergenic"}')
```

## Ordering Sites From A Gel (Classical Mapping)

When the sequence is unknown, a map is reconstructed from fragment sizes, not read off positions. The logic, in order of power:

- **Single digest** gives fragment sizes but not their order. A circular molecule cut n times gives n fragments; a linear one gives n+1.
- **Double digest** orders sites: run enzyme A alone, B alone, and A+B together. A single-digest fragment that disappears and is replaced by two smaller ones in the double digest contains a B site, which is thereby located inside that fragment. Iterate, enforcing that all fragment sizes sum to the molecule length. In practice one enzyme pair often admits several orderings consistent with the same band sizes (co-migrating or symmetric fragments); resolving a unique map needs additional enzymes or the partial-digest method below.
- **Partial digest with end labeling** (Smith-Birnstiel) measures positions directly: label one end, cut only a random subset of sites, and each labeled partial runs from the labeled end to one site, so its size is that site's distance from the labeled end. This sidesteps the combinatorial ambiguity of double digests.

Maps from one enzyme pair are often ambiguous (co-migrating or symmetric fragments fit multiple orderings); resolving a unique map needs several enzymes and the sum-of-fragments constraint.

## Circular Fragment Distances

```python
def circular_distances(sites, seq_len):
    '''Fragment sizes around a circle from sorted cut positions.'''
    s = sorted(sites)
    spans = [s[i + 1] - s[i] for i in range(len(s) - 1)]
    return spans + [(seq_len - s[-1]) + s[0]]    # the wrap-around fragment closes the circle

frags = circular_distances(EcoRI.search(record.seq, linear=False), len(record.seq))
assert sum(frags) == len(record.seq)             # the circle must be fully accounted for
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `AttributeError: ... 'format_as'` | Method is `format_output` | Use `Analysis.format_output()` to get the text as a string |
| Wrap-around fragment is too short on a plasmid | Used `seq_len - last_site` instead of `(seq_len - last) + first` | Close the circle across the origin |
| Site near the origin missing on a plasmid map | Built the map with `linear=True` | Pass `linear=False` for circular DNA |
| Running a mapping script litters PNG/TXT files | Wrote outputs to a hard-coded filename | Write only to a path the caller supplies (or a temp dir) |
| Two enzymes' sites cannot be ordered from one gel | Single digest gives sizes, not order | Add a double digest (or partial-digest end-labeling) and use the sum check |

## Related Skills

- restriction-sites - Find the cut positions a map is built from
- fragment-analysis - Fragment sizes and gel interpretation behind classical mapping
- enzyme-selection - Choose informative enzymes for a map
- data-visualization/genome-tracks - Richer graphical layouts for annotated maps
- sequence-io/read-sequences - Load FASTA or GenBank input

## References

- Smith HO, Birnstiel ML. A simple method for DNA restriction site mapping. Nucleic Acids Res. 1976;3(9):2387-2398. doi:10.1093/nar/3.9.2387
- Roberts RJ, Vincze T, Posfai J, Macelis D. REBASE: a database for DNA restriction and modification: enzymes, genes and genomes. Nucleic Acids Res. 2023;51(D1):D629-D630. doi:10.1093/nar/gkac975
