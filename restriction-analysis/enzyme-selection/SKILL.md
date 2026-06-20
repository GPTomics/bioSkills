---
name: bio-restriction-enzyme-selection
description: Select restriction enzymes for cloning or diagnostics using Biopython Bio.Restriction. Finds enzymes by cut frequency, overhang type, recognition-site length, commercial availability, compatible ends, and methylation sensitivity, and identifies isoschizomers and compatible pairs. Use when choosing which enzymes to use to linearize a vector, drop in an insert, set up a diagnostic digest, or pick a methylation-insensitive enzyme.
tool_type: python
primary_tool: Bio.Restriction
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+ (API verified on 1.86)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Restriction.Analysis)` to confirm method names

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying. The
Analysis cut-count methods were renamed across versions (see Common Errors).

# Restriction Enzyme Selection

**"Pick enzymes to clone my insert into this vector"** -> Search the enzyme database under the constraints that actually matter: cuts the vector once, leaves the insert intact, makes a usable end, is buyable, and is not silenced by methylation.
- Python: `Bio.Restriction.Analysis(CommOnly, seq)` with `with_N_sites`/`without_site`, plus enzyme predicates for overhang and compatibility.

The canonical selection is an intersection, not a single query: an enzyme that **cuts the vector exactly once at the cloning site** AND **does not cut the insert** AND **leaves the intended overhang** AND **is commercially available** AND **is not blocked by the methylation on the source DNA**. Each constraint below is one filter in that intersection; the skill's value is composing them, not running any one alone.

## The Selection Decision Table

| Constraint | What to ask | API / source |
|------------|-------------|--------------|
| Cut frequency | Cut vector once? Leave insert uncut? | `Analysis.with_N_sites(1)`, `Analysis.without_site()` |
| Recognition length | naive 1/4^n spacing: 4-cutter (~256 bp, frequent), 6-cutter (~4 kb, routine cloning), 8-cutter (~65 kb, rare; large constructs, mapping) -- real genomes deviate (see below) | `len(enzyme.site)` |
| Overhang | 5' overhang (most common), 3' overhang, or blunt (non-directional, inefficient ligation) | `is_5overhang()`, `is_3overhang()`, `is_blunt()` |
| Directionality | Two different ends so the insert goes one way and the vector cannot self-ligate | two single-cutters with non-compatible ends |
| Availability | Can it be purchased? | membership in `CommOnly` |
| Methylation | Is the site blocked by Dam/Dcm/CpG on this DNA? | curated table below + REBASE (not the coarse `is_methylable()`) |
| Fidelity | Avoid star activity under forcing conditions | prefer High-Fidelity (HF) enzymes; benchtop, not in BioPython |

Naive cut-frequency intuition (a 6-cutter every 4^6 = 4096 bp) fails on real genomes: vertebrate CpG suppression makes any CpG-containing site -- NotI `GCGGCCGC` above all -- far rarer than 1/4^n, which is exactly why NotI and other 8-cutters are the rare-cutters of choice for large mammalian fragments.

## Find Enzymes By Cut Frequency

**Goal:** Sort candidates into single-cutters (linearize), double-cutters (excise), and non-cutters (safe through the digest).

**Approach:** `Analysis` exposes `with_N_sites(n)` for an exact cut count and `without_site()` for non-cutters. (The older `once_cutters()` / `twice_cutters()` / `only_dont_cut()` / `only_cut()` names do not exist in current BioPython and raise `AttributeError`.)

```python
from Bio import SeqIO
from Bio.Restriction import Analysis, CommOnly

record = SeqIO.read('sequence.fasta', 'fasta')
analysis = Analysis(CommOnly, record.seq)

single_cutters = analysis.with_N_sites(1)   # linearization candidates
double_cutters = analysis.with_N_sites(2)   # excise-an-insert candidates
non_cutters    = analysis.without_site()    # safe to keep in a multi-enzyme digest
all_cutters    = analysis.with_sites()      # any number of sites
print(f'{len(single_cutters)} single-cutters, {len(non_cutters)} non-cutters')
```

## Select A Pair For Directional Cloning

**Goal:** Two enzymes that each cut the vector once, neither cuts the insert, and their ends differ so the ligation is directional and the vector cannot recircularize.

**Approach:** Intersect "cuts vector once" with "does not cut insert", then pair candidates whose ends are mutually INCOMPATIBLE -- that is what makes the cloning directional and stops the vector self-ligating. Two enzymes are an incompatible (directional) pair when neither appears in the other's `compatible_end()`.

```python
from itertools import combinations
from Bio.Restriction import Analysis, CommOnly

def directional_pairs(vector_seq, insert_seq):
    vec_once  = set(Analysis(CommOnly, vector_seq, linear=False).with_N_sites(1))
    ins_clear = set(Analysis(CommOnly, insert_seq).without_site())
    candidates = sorted(vec_once & ins_clear, key=str)
    pairs = []
    for a, b in combinations(candidates, 2):
        if b not in a.compatible_end():          # incompatible ends -> directional, no self-ligation
            pairs.append((a, b))
    return pairs                                  # each pair cuts vector once, leaves insert intact
```

## Find Compatible And Isocaudomer Ends

**Goal:** Identify enzymes whose overhangs ligate together, including enzymes with different recognition sites that leave the same overhang (isocaudomers).

**Approach:** `compatible_end()` returns every enzyme that *can* leave a compatible overhang -- including Type IIS enzymes whose overhang is user-defined, not fixed. For a real ligation partner, filter the result to fixed-overhang Type IIP enzymes (the true isocaudomers, e.g. BamHI/BglII/BclI/Sau3AI all leave 5'-GATC); a Type IIS enzyme listed here is not a drop-in cloning partner.

```python
from Bio.Restriction import BamHI

partners = BamHI.compatible_end()                  # any enzyme that can leave a 5'-GATC end
fixed = [e for e in partners if e.is_palindromic()] # keep Type IIP isocaudomers (BglII, BclI, MboI, Sau3AI...)
print(f'BamHI isocaudomers (fixed overhang): {sorted(str(e) for e in fixed)}')
```

Ligating two different-but-compatible sites usually creates a hybrid junction that **neither enzyme re-cleaves** (BamHI `G^GATCC` + BglII `A^GATCT` -> `GGATCT`/`AGATCC`, which is neither site). This makes the join directional and is used deliberately to destroy one site -- but whether the junction is recut is pair-dependent, so verify the specific pair rather than assuming.

## Filter By Overhang And Recognition Length

```python
from Bio.Restriction import CommOnly, Analysis

cutters = Analysis(CommOnly, record.seq).with_sites()

blunt   = [e for e in cutters if e.is_blunt()]
five_p  = [e for e in cutters if e.is_5overhang()]
three_p = [e for e in cutters if e.is_3overhang()]

six_cutters   = [e for e in CommOnly if len(e.site) == 6]   # routine cloning
eight_cutters = [e for e in CommOnly if len(e.site) == 8]   # rare cutters
```

## Methylation Sensitivity (The Silent-Failure Trap)

Standard E. coli cloning strains (DH5-alpha, JM109, TOP10) are **dam+ dcm+**, so plasmid and insert DNA prepped from them is methylated at GATC (Dam, N6-methyladenine) and CCWGG (Dcm, 5-methylcytosine). An enzyme blocked by that mark will fail or partially cut even though the recognition site is present -- a silent failure. Mammalian genomic DNA additionally carries CpG (5mC) methylation. The fix is to re-propagate the DNA in a **dam- dcm- strain** (GM2163, JM110, INV110) before cutting.

| Site context | Enzyme | Behavior on the methylated site |
|--------------|--------|----------------------------------|
| Dam GATC | DpnI | Cuts ONLY when fully Dam-methylated (methylation-dependent) |
| Dam GATC | DpnII, MboI | Blocked by Dam methylation (cut only unmethylated GATC) |
| Dam GATC | Sau3AI | Insensitive to Dam (cuts methylated or not) |
| CpG CCGG | HpaII | Blocked by CpG methylation of the internal C |
| CpG CCGG | MspI | Cuts regardless of CpG methylation (isoschizomer of HpaII) |
| Dam-overlapping | ClaI `ATCGAT`, XbaI `TCTAGA` | Blocked when flanking bases create an overlapping Dam GATC |

Do NOT rely on BioPython's `enzyme.is_methylable()` to make this decision: it is a coarse REBASE flag (it returns True for Sau3AI, which is actually Dam-insensitive, and for EcoRI), does not distinguish Dam vs Dcm vs CpG, and does not indicate the direction of the effect. Use the curated cases above and consult REBASE for the specific methyltransferase that blocks a given enzyme.

```python
from Bio.Restriction import DpnI, DpnII, Sau3AI, MboI

# Curated, not from is_methylable(): the GATC quartet a cloner must know.
dam_behavior = {
    'DpnI': 'requires Dam methylation to cut',
    'DpnII': 'blocked by Dam methylation',
    'MboI': 'blocked by Dam methylation',
    'Sau3AI': 'insensitive to Dam methylation',
}
for enz in (DpnI, DpnII, MboI, Sau3AI):
    print(f'{enz} ({enz.site}): {dam_behavior[str(enz)]}')
```

## Isoschizomers, Neoschizomers, And Why The Choice Matters

```python
from Bio.Restriction import SmaI, XmaI

print('SmaI isoschizomers:', SmaI.isoschizomers())     # all same-site enzymes (Cfr9I, TspMI, XmaI)
print('SmaI elucidate:', SmaI.elucidate())             # CCC^_GGG  -> blunt
print('XmaI elucidate:', XmaI.elucidate())             # C^CCGG_G  -> 5' overhang
```

- Isoschizomers recognize the same site; pick among them for a different buffer, supplier, or methylation sensitivity (HpaII vs MspI differ only in CpG sensitivity; MboI vs Sau3AI in Dam sensitivity).
- A neoschizomer recognizes the same site but cuts at a different position -- the lever for choosing blunt vs sticky ends from one sequence (SmaI `CCC^GGG` blunt vs XmaI `C^CCGGG` 5' overhang). Note: in current BioPython `neoschizomers()` and `isoschizomers()` overlap (both list all same-site enzymes), so confirm the actual cut difference with `elucidate()` rather than trusting the method name to filter.

## Star Activity And High-Fidelity Enzymes

Under forcing conditions -- >5% glycerol, low ionic strength, high pH (>8), large enzyme excess or over-long incubation, or Mn2+ replacing Mg2+ -- many enzymes relax specificity and cut near-cognate sites ("star activity"; EcoRI* is the classic case). When a clean digest matters, prefer an engineered High-Fidelity (HF) enzyme (e.g. EcoRI-HF), which is selected to show no star activity even in overnight, high-unit digests. This is a benchtop property, not encoded in BioPython; surface it when recommending an enzyme.

## Type IIS / Golden Gate

Type IIS enzymes (BsaI, BsmBI, BbsI, SapI) cut outside their recognition site and enable scarless, directional, one-pot assembly. Selecting and validating them -- including domestication of internal sites and fusion-overhang design -- is its own analysis; route to restriction-analysis/golden-gate-assembly. For plain selection, a part is "Golden Gate ready" for an enzyme when `enzyme.search(seq)` returns no internal sites.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `AttributeError: ... 'once_cutters'` / `'only_dont_cut'` / `'only_cut'` | Methods renamed across BioPython versions | `with_N_sites(1)`/`with_N_sites(2)` for exact counts; `without_site()` for non-cutters; `with_sites()` for any cutter |
| `AttributeError: ... 'is_dam_methylable'` / `'is_dcm_methylable'` | These methods do not exist | Use the curated Dam/Dcm table above; consult REBASE for specifics |
| `AttributeError: ... 'fst3cut'` / `'fst5cut'` | Attribute names are `fst3` / `fst5` | Use `enzyme.fst5` / `enzyme.fst3` (Type IIS cut offsets) |
| Chosen enzyme fails to cut a real prep | Site blocked by Dam/Dcm methylation | Re-prep DNA in a dam- dcm- strain, or pick a methylation-insensitive enzyme |
| Recommended enzyme cannot be bought | Searched `AllEnzymes` | Restrict to `CommOnly` |
| Blunt clone has high background / wrong orientation | Blunt ends ligate inefficiently and non-directionally | Prefer two different sticky ends; dephosphorylate the vector (see usage guide) |
| Double digest only partially cuts | The two chosen enzymes share no buffer where both are fully active | Pick a pair compatible in one universal buffer (rCutSmart / FastDigest); otherwise digest sequentially, lower-salt enzyme first (this is a selection criterion when choosing the pair) |

## Related Skills

- restriction-sites - Find where the selected enzymes cut
- restriction-mapping - Map the selected enzyme sites
- fragment-analysis - Predict the fragments a chosen digest produces
- golden-gate-assembly - Select and validate Type IIS enzymes for scarless assembly
- primer-design/primer-basics - Add chosen restriction sites to PCR primer tails

## References

- Roberts RJ, Vincze T, Posfai J, Macelis D. REBASE: a database for DNA restriction and modification: enzymes, genes and genomes. Nucleic Acids Res. 2023;51(D1):D629-D630. doi:10.1093/nar/gkac975
- Waalwijk C, Flavell RA. MspI, an isoschizomer of HpaII which cleaves both unmethylated and methylated HpaII sites. Nucleic Acids Res. 1978;5(9):3231-3236. doi:10.1093/nar/5.9.3231
- Geier GE, Modrich P. Recognition sequence of the dam methylase of Escherichia coli K12 and mode of cleavage of DpnI endonuclease. J Biol Chem. 1979;254(4):1408-1413.
- Wei H, Therrien C, Blanchard A, Guan S, Zhu Z. The Fidelity Index provides a systematic quantitation of star activity of DNA restriction endonucleases. Nucleic Acids Res. 2008;36(9):e50. doi:10.1093/nar/gkn182
