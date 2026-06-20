---
name: bio-restriction-fragment-analysis
description: Predict restriction digest fragment sizes and gel patterns using Biopython Bio.Restriction. Computes fragment lengths and sequences for single and double digests on linear or circular DNA, and interprets them against an agarose gel. Use when predicting the fragments from a digest, planning a diagnostic digest to verify a clone, or matching observed gel bands to an expected pattern.
tool_type: python
primary_tool: Bio.Restriction
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+ (API verified on 1.86)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython`, then check the return shape with `from Bio.Restriction import EcoRI; from Bio.Seq import Seq; EcoRI.catalyze(Seq('GAATTCGAATTC'))` (a tuple of fragment Seqs)

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Restriction Fragment Analysis

**"What fragments will this digest produce?"** -> Simulate the cut and return fragment lengths (and sequences), for one enzyme or a double digest, on linear or circular DNA.
- Python: `enzyme.catalyze(seq, linear=...)` returns a tuple of fragment `Seq` objects directly.

Two facts decide whether the prediction is right. First, **`catalyze()` returns the tuple of fragments itself** -- `EcoRI.catalyze(seq)` is `(Seq(...), Seq(...), Seq(...))`. Do NOT write `catalyze(seq)[0]` to "get the fragments": that returns only the first fragment, and iterating it then counts its bases, producing nonsense sizes like `[1, 1, 1, ...]`. Second, **topology sets the fragment count**: a linear molecule with n cuts yields n+1 fragments; a circular molecule with n cuts yields n. Passing the default `linear=True` to a plasmid invents one extra fragment that does not exist on the bench.

## Fragment Count By Topology

| Molecule | n cut sites yields | Why |
|----------|--------------------|-----|
| Linear (PCR product, genomic fragment, lambda) | n + 1 fragments | Cut once -> two pieces |
| Circular (plasmid, many viral genomes) | n fragments | Cut once -> one linearized piece; the origin-spanning fragment wraps |

A plasmid showing three bands on a single-enzyme digest has three sites, not two. Counting plasmid bands as if the molecule were linear is the most common digest-interpretation error.

## Predict Fragment Sizes

**Goal:** Turn a digest into the band sizes a gel would show.

**Approach:** Call `catalyze()` with the correct topology, take the lengths of the returned fragments, sort descending. The summed sizes must equal the molecule length -- a built-in correctness check.

```python
from Bio import SeqIO
from Bio.Restriction import EcoRI

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq

fragments = EcoRI.catalyze(seq, linear=True)   # tuple of Seq; NO [0]
sizes = sorted((len(f) for f in fragments), reverse=True)
assert sum(sizes) == len(seq)                  # fragments must account for the whole molecule
print(f'{len(sizes)} fragments: {sizes}')
```

## Linear vs Circular Digestion

```python
from Bio.Restriction import EcoRI

linear_frags   = EcoRI.catalyze(seq, linear=True)
circular_frags = EcoRI.catalyze(seq, linear=False)   # plasmid: one fewer fragment
print(f'Linear: {len(linear_frags)} fragments; Circular: {len(circular_frags)} fragments')
```

## Double Digest

**Goal:** Predict fragments when two enzymes cut the same molecule.

**Approach:** A double digest cuts at the union of both site sets. The robust way is to build a `RestrictionBatch` and let `Analysis` collect every position, then compute fragments from the sorted positions (this generalizes to any number of enzymes and to circular DNA). Sequential `catalyze` calls also work but are easy to get wrong on circular DNA.

```python
from Bio.Restriction import EcoRI, BamHI, RestrictionBatch, Analysis

def fragments_from_positions(seq_len, positions, linear=True):
    '''Fragment sizes (bp) from a set of 1-based cut positions.'''
    cuts = sorted(set(positions))
    if not cuts:
        return [seq_len]
    spans = [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]
    if linear:
        return [cuts[0]] + spans + [seq_len - cuts[-1]]
    return spans + [(seq_len - cuts[-1]) + cuts[0]]    # circular: wrap-around fragment

batch = RestrictionBatch([EcoRI, BamHI])
positions = [p for sites in Analysis(batch, seq).with_sites().values() for p in sites]
sizes = sorted(fragments_from_positions(len(seq), positions, linear=True), reverse=True)
print(f'Double digest: {len(sizes)} fragments: {sizes}')
assert sum(sizes) == len(seq)
```

## Interpreting The Gel: Size, Resolution, And Topology

Fragment sizes are not what a gel directly reports; migration distance is. The decisions below come from gel physics, not from BioPython, and they determine whether predicted bands are actually resolvable.

| Reality | Consequence for interpretation |
|---------|--------------------------------|
| Migration distance is ~linear in -log10(size) only within a gel's resolving window | Sizing is reliable only against a ladder run on the same gel; never reuse a standard curve across gels |
| Agarose percent sets the window (approximate, buffer- and voltage-dependent: ~0.7% resolves ~0.8-12 kb; ~1% ~0.5-10 kb; ~2% ~0.1-2 kb) | Choose the percent for the bands of interest; small fragments run off a low-percent gel, large fragments pile up on a high-percent one |
| Above ~20-50 kb fragments co-migrate (reptation); resolving them needs pulsed-field (PFGE) | Two large predicted bands may appear as one on a conventional gel |
| Two fragments of similar size co-migrate as one band of doubled intensity (intensity ~ mass) | A "missing" predicted band is often a co-migrating doublet, not an error -- check the sum |
| Uncut plasmid runs as supercoiled (fast, anomalous) + nicked (slow), same molecule | Do not size a supercoiled band off a linear ladder; only a linearized (single-cut) plasmid sizes correctly |

## Simulate A Gel Pattern (Text)

**Goal:** Lay predicted digests next to a ladder to see which bands resolve and which co-migrate.

**Approach:** Pool all band sizes and the ladder, sort descending, and mark each lane. Co-migration shows up as multiple marks at one size.

```python
def simulate_gel(digests, ladder=None):
    '''digests: {lane_name: [sizes]} -> print a text gel against a ladder.'''
    ladder = ladder or [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250]
    bands = sorted(set(ladder).union(*[set(s) for s in digests.values()]), reverse=True)
    header = f'{"size":>6} | {"ladder":^6} | ' + ' | '.join(f'{n:^8}' for n in digests)
    print(header); print('-' * len(header))
    for b in bands:
        row = f'{b:>6} | {"---" if b in ladder else "":^6} | '
        row += ' | '.join(f'{"=" * 4 * lane.count(b):^8}' for lane in digests.values())
        print(row)

simulate_gel({'EcoRI': sizes})
```

## Diagnostic Digest Report

**Goal:** Document an expected digest to plan a clone-verification gel.

**Approach:** Report site count, positions, fragment sizes, and the total; flag any pair of fragments too close to resolve.

```python
def fragment_report(seq, enzyme, linear=True, resolution=0.10):
    sites = enzyme.search(seq, linear=linear)
    sizes = sorted((len(f) for f in enzyme.catalyze(seq, linear=linear)), reverse=True)
    print(f'{enzyme} ({enzyme.site}): {len(sites)} site(s) at {sites}')
    print(f'  {len(sizes)} fragments, total {sum(sizes)} bp: {sizes}')
    close = [(a, b) for a, b in zip(sizes, sizes[1:]) if a and (a - b) / a < resolution]
    if close:
        print(f'  WARNING likely co-migrating (<{resolution:.0%} apart): {close}')
    return sizes

fragment_report(seq, EcoRI)
```

## Compare Expected vs Observed Bands

```python
def compare_fragments(expected, observed, tolerance=50):
    '''Match predicted sizes to gel-measured sizes within a tolerance (bp).'''
    obs = list(observed)
    matched, missing = [], []
    for exp in expected:
        hit = next((o for o in obs if abs(exp - o) <= tolerance), None)
        if hit is None:
            missing.append(exp)
        else:
            matched.append((exp, hit)); obs.remove(hit)
    print('matched:', matched)
    if missing: print('missing (predicted, not seen -- co-migration or partial digest?):', missing)
    if obs:     print('extra (seen, not predicted -- star activity, contaminant, or wrong map?):', obs)

compare_fragments(expected=[3000, 2000, 1500, 500], observed=[3050, 2000, 1480, 510, 200])
```

Extra bands and a smeary ladder of intermediate sizes often signal a **partial digest** (not every site cut), which produces sums of adjacent complete-digest fragments. Drive the reaction to completion (more enzyme or time) before concluding the map is wrong.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Fragment sizes are all 1 (or wildly wrong) | `catalyze(seq)[0]` returns the first fragment, then `len(f) for f in it` counts bases | Use `fragments = enzyme.catalyze(seq, linear=...)` with no `[0]` |
| Predicted one fragment too many for a plasmid | Digested a circular molecule with `linear=True` | Pass `linear=False`; circular gives n fragments from n cuts |
| Fragment sizes do not sum to the molecule length | A band was missed or two co-migrate | The `sum(sizes) == len(seq)` check is mandatory; a shortfall means a hidden doublet or run-off fragment |
| Two predicted bands never separate on the gel | Sizes within the gel's resolving limit, or both >20-50 kb | Adjust agarose percent, or use PFGE for very large fragments |

## Related Skills

- restriction-sites - Find the cut positions that drive fragment calculation
- restriction-mapping - Order sites and compute inter-site distances
- enzyme-selection - Choose enzymes that give a resolvable, diagnostic band pattern
- sequence-manipulation/seq-objects - Work with the fragment Seq objects

## References

- Smith HO, Birnstiel ML. A simple method for DNA restriction site mapping. Nucleic Acids Res. 1976;3(9):2387-2398. doi:10.1093/nar/3.9.2387
- Schwartz DC, Cantor CR. Separation of yeast chromosome-sized DNAs by pulsed field gradient gel electrophoresis. Cell. 1984;37(1):67-75. doi:10.1016/0092-8674(84)90301-5
