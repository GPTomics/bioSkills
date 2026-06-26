---
name: bio-sequence-properties
description: Calculate nucleotide and protein sequence properties (GC content, GC skew, molecular weight, melting temperature, isoelectric point, instability, hydropathy) with Biopython. Use when analyzing sequence composition, computing primer Tm, estimating DNA or protein mass, or profiling protein biophysical properties.
tool_type: python
primary_tool: Bio.SeqUtils
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Sequence Properties

Calculate physical and chemical properties of nucleotide and protein sequences using Biopython.

**"Calculate GC content"** -> Compute the fraction of G+C bases in a nucleotide sequence.
- Python: `gc_fraction(seq)` (Bio.SeqUtils) - returns a FRACTION 0-1, multiply by 100 for percent.

**"Compute a primer melting temperature"** -> Estimate Tm for hybridization or PCR.
- Python: `MeltingTemp.Tm_NN(seq)` (Bio.SeqUtils) - nearest-neighbor, the accurate method for primers.

**"Analyze protein properties"** -> Compute MW, pI, stability, hydrophobicity from an amino-acid sequence.
- Python: `ProteinAnalysis(str_seq)` (Bio.SeqUtils.ProtParam).

## The governing principle

Most of these functions return a plausible number for any input, so the danger is silent wrongness, not crashes. Three defaults bite hardest: `gc_fraction` returns a FRACTION (not the percent the legacy `GC()` returned), `molecular_weight` defaults to a SINGLE strand (~half a duplex), and `Tm_Wallace`/`Tm_GC` are composition-only methods that are wrong for real primers. Pick the function to match the question and verify its units, not just that it ran.

## Required Imports

```python
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight, GC123, GC_skew, MeltingTemp, nt_search, seq1, seq3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
```

## DNA/RNA Properties

### GC Content

`gc_fraction()` returns a fraction in [0, 1]. The legacy `GC()` returned a percent in [0, 100] and was REMOVED in BioPython 1.82 (`from Bio.SeqUtils import GC` now raises ImportError).

```python
from Bio.SeqUtils import gc_fraction

seq = Seq('ATGCGATCGATCGATCGATCG')
gc = gc_fraction(seq)        # 0.476... (FRACTION, not percent)
gc_percent = gc * 100        # 47.6 - multiply for percent
```

Factor-of-100 trap: porting `GC(seq)` to `gc_fraction(seq)` without `* 100` silently underreports 100x, so downstream filters like "GC > 40" reject everything.

Ambiguity-default trap: the new default `ambiguous='remove'` strips ambiguity codes before computing, but legacy `GC()` counted them in the length only, which equals the new `ambiguous='ignore'`. The faithful drop-in replacement is `gc_fraction(seq, ambiguous='ignore') * 100`. The modes only diverge on sequences that actually contain ambiguity codes (so clean test fixtures hide the difference, real data exposes it).

```python
gc_fraction(seq, ambiguous='remove')    # default: ambiguity codes stripped (neither numerator nor denominator)
gc_fraction(seq, ambiguous='ignore')    # counts ambiguous in denominator only - matches legacy GC()
gc_fraction(seq, ambiguous='weighted')  # each code contributes its mean GC probability (N/X = 0.5)
```

### GC at Codon Positions (GC123)

`GC123()` returns FOUR PERCENTAGES (0-100) - total GC plus GC at codon positions 1, 2, 3 (position 3 is the wobble base, most free to vary under codon bias). Note the unit inconsistency with `gc_fraction`: these are percentages, not fractions. `GC123` does not handle ambiguity codes.

```python
from Bio.SeqUtils import GC123

gc_total, gc_pos1, gc_pos2, gc_pos3 = GC123(Seq('ATGCGATCGATCGATCGATCG'))  # all 0-100
```

### GC Skew

`GC_skew(seq, window=100)` returns `(G-C)/(G+C)` for each non-overlapping window. A window with no G or C returns 0 (the zero-division is guarded), which can be misread as "no skew" rather than "no data".

```python
from Bio.SeqUtils import GC_skew

skew_values = GC_skew(seq, window=1000)  # list of per-window skew values
```

Biology: cumulative GC skew has a global MINIMUM at the replication origin (oriC) and a MAXIMUM at the terminus on circular bacterial chromosomes - the leading strand is G-enriched from strand-asymmetric mutation/repair. This is the basis of in-silico origin prediction. Compute the cumulative skew by taking the running sum of `GC_skew()`.

`xGC_skew()` is a GRAPHICS routine (its docstring literally says "GRAPHICS !!!") that draws on a Tkinter canvas. It raises a loud TclError/ImportError in headless environments. For headless cumulative-skew analysis, sum `GC_skew()` yourself instead.

### Molecular Weight

`molecular_weight(seq, seq_type='DNA', double_stranded=False, circular=False, monoisotopic=False)`.

```python
from Bio.SeqUtils import molecular_weight

dna = Seq('ATGCGATCG')
mw_ss = molecular_weight(dna)                          # single-stranded (DEFAULT)
mw_ds = molecular_weight(dna, double_stranded=True)    # full duplex mass
mw_circ = molecular_weight(dna, circular=True)         # no terminal phosphate adjustment
mw_rna = molecular_weight(Seq('AUGCGAUCG'), seq_type='RNA')
mw_prot = molecular_weight(Seq('MRCRS'), seq_type='protein')
mw_mono = molecular_weight(dna, monoisotopic=True)     # most-abundant isotope, for high-res MS
```

Double-stranded trap (~2x, silent): the default `double_stranded=False` returns a single-strand mass - roughly half a duplex. For genomic dsDNA pass `double_stranded=True`. It is NOT exactly half, because the complementary strand has a different base composition, so a non-self-complementary single-strand number cannot be "fixed later" by doubling. ng-to-molecule, copy-number, and molarity conversions all come out ~2x wrong.

Monoisotopic vs average: the default is average mass (bulk, spectrophotometric). Pass `monoisotopic=True` to match high-resolution mass spec (ESI/MALDI); the wrong choice is a silent systematic offset that grows with mass. Ambiguous letters raise ValueError (loud - good).

### Melting Temperature

**Goal:** Estimate the Tm of an oligo, choosing a method that matches its length and use.

**Approach:** Use `Tm_NN` for any PCR primer (nearest-neighbor, sequence-order aware). Reserve `Tm_Wallace` for very short probes and `Tm_GC` only when a composition-only estimate is acceptable.

```python
from Bio.SeqUtils import MeltingTemp as mt

primer = Seq('ACGGTCAGGTCAGGTACGGT')

tm = mt.Tm_NN(primer, strict=True)                       # accurate primer Tm
tm_salt = mt.Tm_NN(primer, Na=50, dnac1=250, dnac2=250)  # 50 mM Na+, 250 nM each strand
tm_mg = mt.Tm_NN(primer, Mg=1.5, dNTPs=0.2, saltcorr=7)  # Mg/dNTPs ONLY honored at saltcorr 6 or 7
```

| Method | Model | Use when | Caveat |
|--------|-------|----------|--------|
| `Tm_Wallace` | 4(G+C) + 2(A+T) "2+4 rule" | Oligos <=14 nt only | Ignores order/salt; WRONG for primers (off 5-10 C+) |
| `Tm_GC` | GC-content empirical equation | Longer sequences, rough estimate | Composition-only, no nearest-neighbor info |
| `Tm_NN` | Nearest-neighbor thermodynamics | PCR primers, probes, accurate work | Needs realistic salt/strand conc for absolute values |

`Tm_NN` defaults: `nn_table=None` which selects `DNA_NN3` (Allawi & SantaLucia 1997), `saltcorr=5`, strand concentrations `dnac1=dnac2=25` nM. `saltcorr` ranges 1-7, but only 6 (Owczarzy 2004) and 7 (Owczarzy 2008) actually use Mg2+/dNTPs - setting `Mg`/`dNTPs` with `saltcorr<=5` silently ignores them. Keep `strict=True` for primer work: it raises on ambiguous or unsupported nearest-neighbor pairs, whereas `strict=False` silently skips them and underestimates Tm.

### IUPAC-Aware Search (nt_search)

```python
from Bio.SeqUtils import nt_search

result = nt_search('ATGCGATCGATCGATNGATC', 'GATNGATC')  # ['GAT[GATC]GATC', 4] - result[0] is the EXPANDED regex, result[1:] are 0-based starts
```

## Protein Properties

**Goal:** Compute biophysical properties of a protein from its amino-acid sequence.

**Approach:** Create one `ProteinAnalysis` object and call its methods. Non-standard residues (B, Z, X, U, `*`, `-`) are absent from the parameter tables and raise KeyError, so sanitize first.

```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis

clean = 'MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNG'.replace('*', '').replace('X', '')
protein = ProteinAnalysis(clean)

mw = protein.molecular_weight()              # protein average MW (Daltons)
pi = protein.isoelectric_point()             # pI from linear-sequence pKa tables
charge = protein.charge_at_pH(7.0)           # net charge at a given pH
ii = protein.instability_index()             # Guruprasad: > 40 => predicted unstable
gravy = protein.gravy()                       # mean Kyte-Doolittle hydropathy (neg = hydrophilic)
arom = protein.aromaticity()                 # relative frequency of F + W + Y
helix, turn, sheet = protein.secondary_structure_fraction()
eps_reduced, eps_oxidized = protein.molar_extinction_coefficient()  # 280 nm, (reduced, cystine)
flex = protein.flexibility()                 # per-residue, fixed window of 9
```

Interpretation caveats:
- `isoelectric_point()` and `charge_at_pH()` use fixed pKa tables on the LINEAR sequence - they ignore 3D environment and post-translational modifications, so a measured pI can differ by a full pH unit or more.
- `gravy(scale='KyteDoolitle')` - the default scale literal is MISSPELLED `'KyteDoolitle'` (one 't'). Passing the correctly-spelled `'KyteDoolittle'` raises KeyError. A single whole-protein average also collapses local topology (a TM helix plus hydrophilic loops can average near 0); use a windowed hydropathy profile for membrane topology.
- `secondary_structure_fraction()` is a composition propensity estimate, not a structure prediction.
- `instability_index()` uses Guruprasad's dipeptide method; > 40 predicts an unstable protein.

### Amino-Acid Code Conversion

```python
from Bio.SeqUtils import seq1, seq3

seq1('MetAlaGlyTrp')           # 'MAGW'  (3-letter -> 1-letter)
seq3('MAGW')                   # 'MetAlaGlyTrp'  (1-letter -> 3-letter, no separator)
```

## Code Patterns

### Per-Record GC Across a FASTA

**Goal:** Report length and GC percent for every record in a file.

**Approach:** Stream records with `SeqIO.parse`, compute GC per record, multiply the fraction by 100.

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def analyze_fasta(filename):
    return [{'id': r.id, 'length': len(r.seq), 'gc': gc_fraction(r.seq) * 100} for r in SeqIO.parse(filename, 'fasta')]
```

### Cumulative GC Skew (Headless)

**Goal:** Locate a candidate replication origin without the Tkinter graphics routine.

**Approach:** Take per-window skew from `GC_skew`, accumulate it, and read off the minimum (oriC) and maximum (terminus).

```python
from Bio.SeqUtils import GC_skew

def cumulative_skew(seq, window=10000):
    skew = GC_skew(seq, window=window)
    positions, cumulative, total = [], [], 0
    for i, s in enumerate(skew):
        total += s
        positions.append(i * window)
        cumulative.append(total)
    ori = positions[cumulative.index(min(cumulative))]
    return positions, cumulative, ori
```

### Full Protein Report

**Goal:** Summarize the key biophysical metrics of a protein in one pass.

**Approach:** Sanitize non-standard residues, build one `ProteinAnalysis` object, and collect each metric into a dict.

```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def protein_report(sequence):
    clean = str(sequence).upper().replace('*', '').replace('X', '')
    protein = ProteinAnalysis(clean)
    helix, turn, sheet = protein.secondary_structure_fraction()
    return {
        'length': len(clean),
        'molecular_weight': protein.molecular_weight(),
        'isoelectric_point': protein.isoelectric_point(),
        'charge_at_pH7': protein.charge_at_pH(7.0),
        'instability_index': protein.instability_index(),
        'gravy': protein.gravy(),
        'aromaticity': protein.aromaticity(),
        'helix_fraction': helix, 'turn_fraction': turn, 'sheet_fraction': sheet,
    }
```

### CpG Observed/Expected Ratio

```python
def cpg_ratio(seq):
    s = str(seq).upper()
    expected = (s.count('C') * s.count('G')) / len(s) if s else 0
    return s.count('CG') / expected if expected > 0 else 0
```

## Property Reference

| Property | Function | Units / Notes |
|----------|----------|---------------|
| GC content | `gc_fraction()` | FRACTION 0-1 (multiply by 100 for percent) |
| GC by codon position | `GC123()` | FOUR PERCENTAGES 0-100 (total + pos 1/2/3) |
| GC skew | `GC_skew()` | (G-C)/(G+C) per window; 0 = no G/C in window |
| Molecular weight | `molecular_weight()` | Daltons; single-strand by DEFAULT |
| Melting temp | `MeltingTemp.Tm_NN()` | Celsius; accurate for primers |
| pI / charge | `isoelectric_point()` / `charge_at_pH()` | Linear pKa only, ignores 3D/PTMs |
| Instability | `instability_index()` | > 40 => predicted unstable |
| Hydropathy | `gravy()` | neg = hydrophilic; default scale misspelled `'KyteDoolitle'` |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| GC values look 100x too small; filters reject all | Used `gc_fraction()` (fraction) where percent expected | Multiply by 100 |
| GC differs from legacy `GC()` on real data | Default `ambiguous='remove'` vs legacy `'ignore'` | Use `gc_fraction(seq, ambiguous='ignore') * 100` for a faithful drop-in |
| `GC_skew` returns 0 across a region | Window had no G or C (guarded division), not true zero skew | Treat 0 as "no data"; widen the window |
| `xGC_skew` raises TclError/ImportError | It is a Tkinter graphics routine, fails headless | Sum `GC_skew()` yourself for cumulative skew |
| Primer Tm off by 5-10 C | Used `Tm_Wallace`/`Tm_GC` (composition-only) | Use `Tm_NN` with realistic salt and strand concentration |
| Mg/dNTPs change nothing in `Tm_NN` | `Mg`/`dNTPs` ignored unless `saltcorr` is 6 or 7 | Set `saltcorr=7` (Owczarzy 2008) for divalent correction |
| MW ~half of expected for genomic dsDNA | `molecular_weight` default `double_stranded=False` | Pass `double_stranded=True` |
| `KeyError` from `ProteinAnalysis` | Non-standard residue (B, Z, X, U, `*`, `-`) | Strip or replace before analysis |
| `KeyError` from `gravy('KyteDoolittle')` | Default scale literal is misspelled `'KyteDoolitle'` (one 't') | Omit the argument or pass `'KyteDoolitle'` |

## References

Lobry JR (1996) Asymmetric substitution patterns in the two DNA strands of bacteria. Mol Biol Evol 13(5):660-665.

Lobry JR, Gautier C (1994) Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in 999 Escherichia coli chromosome-encoded genes. Nucleic Acids Res 22(15):3174-3180.

SantaLucia J Jr (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS 95(4):1460-1465.

Kyte J, Doolittle RF (1982) A simple method for displaying the hydropathic character of a protein. J Mol Biol 157(1):105-132.

Guruprasad K, Reddy BVB, Pandit MW (1990) Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. Protein Eng 4(2):155-161.

## Related Skills

- seq-objects - Create and modify Seq objects before property calculation
- codon-usage - GC123 and codon-bias indices for coding-sequence analysis
- transcription-translation - Translate a CDS before protein property analysis
- sequence-io/sequence-statistics - File-level statistics (N50, totals, dataset GC)
- primer-design/primer-basics - Design primers where Tm_NN and GC content drive the choices
- restriction-analysis/restriction-sites - Locate enzyme recognition sites in the same sequence
