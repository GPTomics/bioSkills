---
name: bio-restriction-sites
description: Find restriction enzyme cut sites in DNA sequences using Biopython Bio.Restriction. Searches single enzymes, batches, or commercial enzyme sets and returns cut positions for linear or circular DNA. Use when locating where one or more restriction enzymes cut a sequence, screening a sequence for the presence or absence of a site, or counting how often an enzyme cuts.
tool_type: python
primary_tool: Bio.Restriction
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+ (API verified on 1.86)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Restriction.Analysis)` to check method names

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying. The
Analysis "cutters" methods in particular were renamed across versions (see Common Errors).

# Finding Restriction Sites

**"Find where this enzyme cuts my DNA"** -> Return the cut positions for one or more restriction enzymes along a linear or circular sequence.
- Python: `enzyme.search(seq, linear=...)` for one enzyme; `Bio.Restriction.Analysis(batch, seq, linear=...).full()` for many.

The one fact that governs every result: `search()` returns a **1-based position equal to the first base of the downstream fragment** (the base immediately 3' of the cut on the top strand), **not** the start of the recognition site. For EcoRI `G^AATTC` whose site starts at position 4, `search` reports 5 (the base after the cut). Confusing the cut position with the recognition-site start is the single most common bug in restriction code, and it propagates silently into fragment sizes and map coordinates.

## The Decisions That Shape A Site Search

| Decision | Options | When to pick which |
|----------|---------|--------------------|
| Enzyme scope | one enzyme / a curated `RestrictionBatch` / `CommOnly` / `AllEnzymes` | A named enzyme when the assay dictates it; a small batch for a cloning panel; `CommOnly` (623 buyable enzymes) when the answer must be an enzyme one can purchase; `AllEnzymes` (1088, includes non-commercial) only for exhaustive in-silico surveys |
| Topology | `linear=True` (default) / `linear=False` | `linear=False` for any plasmid, viral circle, or BAC. A circular molecule lets a site span the origin and changes fragment counts (see fragment-analysis) |
| Question | does it cut? / how often? / where? | `search()` for positions; `Analysis.with_N_sites(n)` for exact cut counts; `bool(search())` for a yes/no screen |

Use `CommOnly` not `AllEnzymes` by default: proposing an enzyme nobody sells wastes a wet-lab cycle. The README's legacy "800+" figure is stale; the installed database holds 1088 enzymes total, 623 commercially available.

## Search With One Enzyme

```python
from Bio import SeqIO
from Bio.Restriction import EcoRI

record = SeqIO.read('sequence.fasta', 'fasta')
seq = record.seq

sites = EcoRI.search(seq)              # list of 1-based cut positions, e.g. [5, 14]
print(f'EcoRI cuts {len(sites)} time(s) at {sites}')
if not sites:
    print('EcoRI does not cut this sequence')
```

## Search With A Batch Of Enzymes

**Goal:** Screen a sequence against several enzymes at once and keep only those that cut.

**Approach:** Build a `RestrictionBatch`, run `Analysis.full()` to get every enzyme's positions, then filter to cutters with `with_sites()`.

```python
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI, HindIII, XhoI

batch = RestrictionBatch([EcoRI, BamHI, HindIII, XhoI])
analysis = Analysis(batch, seq, linear=True)

cutters = analysis.with_sites()        # {enzyme: [positions]} only enzymes that cut
for enzyme, positions in cutters.items():
    print(f'{enzyme}: {positions}')
```

## Filter By Cut Count

**Goal:** Separate single-cutters (linearize a plasmid), double-cutters (excise an insert), and non-cutters (safe to carry through a digest).

**Approach:** `Analysis` exposes `with_N_sites(n)` for an exact count and `without_site()` for enzymes with no site. (The older `once_cutters()`/`twice_cutters()`/`only_dont_cut()` names do not exist in current BioPython.)

```python
from Bio.Restriction import Analysis, CommOnly

analysis = Analysis(CommOnly, seq, linear=False)   # circular plasmid

single_cutters = analysis.with_N_sites(1)          # {enzyme: [pos]} good for linearization
double_cutters = analysis.with_N_sites(2)          # {enzyme: [pos, pos]} good for excision
non_cutters    = analysis.without_site()           # {enzyme: []} safe in a multi-step digest
all_cutters    = analysis.with_sites()             # any number of sites

print(f'{len(single_cutters)} single-cutters, {len(non_cutters)} non-cutters')

# Pretty-print a chosen subset
analysis.print_as('map')
analysis.print_that(single_cutters)                # formats the dict you pass it
```

## Built-In Enzyme Collections

```python
from Bio.Restriction import AllEnzymes, CommOnly, Analysis

print(f'{len(AllEnzymes)} known enzymes, {len(CommOnly)} commercially available')

analysis = Analysis(CommOnly, seq)                 # default: only buyable enzymes
for enzyme, positions in analysis.with_sites().items():
    print(f'{enzyme}: {positions}')
```

## Linear vs Circular DNA

```python
from Bio.Restriction import EcoRI

sites_linear   = EcoRI.search(seq, linear=True)    # ends are free; no wrap-around
sites_circular = EcoRI.search(seq, linear=False)   # a site may span the origin
```

A circular search can find a site that straddles position 1, which a linear search misses. Always pass `linear=False` for plasmids; the fragment count and map differ (see restriction-analysis/fragment-analysis).

## Read An Enzyme's Cut Geometry

**Goal:** Know what ends an enzyme leaves before designing a ligation.

**Approach:** `elucidate()` draws the cut unambiguously; the boolean predicates and the signed `ovhg` summarize it. The sign convention is the trap: **negative `ovhg` is a 5' overhang, positive is a 3' overhang, zero is blunt.**

```python
from Bio.Restriction import EcoRI, KpnI, EcoRV

for enz in (EcoRI, KpnI, EcoRV):
    print(enz, enz.elucidate())        # EcoRI G^AATT_C ; KpnI G_GTAC^C ; EcoRV GAT^_ATC
    print(f'  site={enz.site} ovhg={enz.ovhg} ovhgseq={enz.ovhgseq!r}'
          f' 5prime={enz.is_5overhang()} 3prime={enz.is_3overhang()} blunt={enz.is_blunt()}')
# EcoRI.ovhg == -4  -> a 5' overhang (NOT +4). In elucidate, ^ = top-strand cut, _ = bottom-strand cut.
```

## Access Enzymes By Name

```python
from Bio.Restriction import AllEnzymes

if 'EcoRI' in AllEnzymes:
    ecori = AllEnzymes.get('EcoRI')
    sites = ecori.search(seq)
```

## Ambiguous And Interrupted Recognition Sites

Not every enzyme has a fixed 6-bp palindrome. Degenerate sites use IUPAC codes (HincII `GTYRAC`), and interrupted palindromes carry an unspecified N spacer (BstXI `CCANNNNNNTGG`, DraIII `CACNNNGTG`). Note what `is_ambiguous()` actually means in BioPython: it is True when the site or cut is ambiguous -- N-spacer / interrupted sites (BstXI, DraIII) and enzymes that cut outside their site -- but it is **False** for a fully IUPAC-degenerate site whose cut is fixed, such as HincII `GTYRAC` (BioPython reports that as `is_defined()`). So `is_ambiguous()` does not detect IUPAC degeneracy; read `enzyme.site` for the actual letters. Either way, the expected cut frequency for a degenerate or N-containing site is not a clean `1/4^n`, so do not estimate cutter rarity from site length alone for these enzymes.

```python
from Bio.Restriction import HincII, BstXI

for enz in (HincII, BstXI):
    print(enz, enz.site, 'ambiguous=', enz.is_ambiguous())
```

## Search Many Sequences

```python
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis, EcoRI, BamHI

batch = RestrictionBatch([EcoRI, BamHI])
for record in SeqIO.parse('sequences.fasta', 'fasta'):
    cutters = Analysis(batch, record.seq).with_sites()
    print(record.id, {str(e): p for e, p in cutters.items()})
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `AttributeError: 'Analysis' object has no attribute 'once_cutters'` | Method renamed across BioPython versions | Use `with_N_sites(1)` / `with_N_sites(2)`; `without_site()` for non-cutters; `with_sites()` for any cutter |
| `AttributeError: ... 'print_that_cut'` / `'esite'` | These names do not exist | Use `print_as(...)` + `print_that(dct)`; read the cut with `elucidate()` |
| Fragment sizes or map coordinates off by a few bases | Treated `search()` output as the recognition-site start | The integer is the cut position = first base of the downstream fragment (1-based) |
| Site near the origin missed on a plasmid | Searched with `linear=True` | Pass `linear=False` for circular DNA |
| Reported a 5' overhang as 3' (or vice versa) | Misread the `ovhg` sign | Negative `ovhg` = 5' overhang, positive = 3' overhang, zero = blunt; confirm with `elucidate()` |
| Proposed enzyme cannot be purchased | Searched `AllEnzymes` | Search `CommOnly` when the answer must be a buyable enzyme |

## Related Skills

- restriction-mapping - Order cut sites and draw a map with inter-site distances
- enzyme-selection - Choose enzymes by cut frequency, overhang, methylation sensitivity, or compatible ends
- fragment-analysis - Turn cut positions into fragment sizes and gel patterns
- golden-gate-assembly - Screen a part for internal Type IIS sites before scarless assembly
- sequence-io/read-sequences - Load the FASTA or GenBank sequence to search

## References

- Roberts RJ, Vincze T, Posfai J, Macelis D. REBASE: a database for DNA restriction and modification: enzymes, genes and genomes. Nucleic Acids Res. 2023;51(D1):D629-D630. doi:10.1093/nar/gkac975
- Roberts RJ, Belfort M, Bestor T, et al. A nomenclature for restriction enzymes, DNA methyltransferases, homing endonucleases and their genes. Nucleic Acids Res. 2003;31(7):1805-1812. doi:10.1093/nar/gkg274
