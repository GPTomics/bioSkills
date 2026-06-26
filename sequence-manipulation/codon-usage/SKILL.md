---
name: bio-codon-usage
description: Analyze codon usage and calculate CAI (Codon Adaptation Index), RSCU, and Nc with Biopython, and produce naive max-CAI codon-optimized sequences. Use when scoring a gene's codon bias against a host, optimizing a CDS for heterologous expression, or studying synonymous codon selection.
tool_type: python
primary_tool: Bio.SeqUtils.CodonAdaptationIndex
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Codon Usage

Analyze codon usage patterns, score adaptation to a host, and optimize coding sequences for expression.

**"Analyze codon usage"** -> Count codons in a coding sequence, compute frequencies and bias metrics.
- Python: `Counter` on in-frame triplets + RSCU/Nc helpers (BioPython + standard library)

**"Score this gene against a host"** -> Compute the Codon Adaptation Index from a reference set of highly expressed genes.
- Python: `CodonAdaptationIndex(reference_seqs).calculate(query)` (BioPython)

**"Optimize codons for expression"** -> Replace each codon with the host's single most-preferred synonymous codon.
- Python: `CodonAdaptationIndex(reference_seqs).optimize(seq)` (BioPython)

## The governing principle

CAI is meaningless without an expression-biased reference. The relative-adaptiveness weights (w) must be built from the **highly expressed genes of the TARGET organism** (ribosomal proteins, elongation factors). A CAI computed against a whole-genome average, or against the wrong organism, is a number with no biological meaning. There is no bundled reference index in modern Biopython, so the reference set is always the caller's responsibility.

Two silent traps dominate this skill:
- **Out-of-frame input is silently corrupted.** `calculate` blindly steps `range(0, len, 3)` from position 0; it never checks the reading frame. A frame-shifted CDS returns a plausible CAI computed from garbage codons. Frame correctness is the caller's job (length divisible by 3, starts at the first base of codon 1).
- **Naive max-CAI optimization can REDUCE expression.** `optimize()` is the textbook max-CAI output (single most-frequent codon per amino acid). It is blind to the translation ramp, 5' mRNA structure, GC extremes, cryptic regulatory elements, and codon-pair bias. It is a starting point to screen, never a final design. Failure is silent: correct protein, high CAI, poor expression.

## CRITICAL API migration (Biopython 1.82)

`Bio.SeqUtils.CodonUsage` and `Bio.SeqUtils.CodonUsageIndices` were **removed in Biopython 1.82**. Any code calling `generate_index()`, `set_cai_index()`, `cai_for_gene()`, `print_index()`, or importing `SharpEcoliIndex` raises ImportError on any modern install. The replacement is a redesigned class imported directly from `Bio.SeqUtils`:

```python
from Bio.SeqUtils import CodonAdaptationIndex  # NOT Bio.SeqUtils.CodonUsage (removed)
```

| Removed (<=1.81) | Replacement (>=1.82) |
|------------------|----------------------|
| `CodonAdaptationIndex()` then `generate_index(fasta)` | `CodonAdaptationIndex(reference_seqs, table=...)` constructor |
| `cai.cai_for_gene(seq)` | `cai.calculate(seq)` |
| `cai.set_cai_index(d)` | `cai.update(d)` (it is a dict subclass) |
| `SharpEcoliIndex` (bundled) | none -- build from a supplied reference CDS set |
| `cai.print_index()` | iterate the object: `for codon, w in cai.items()` |

## Required Imports

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex, GC123
from Bio.Data.CodonTable import standard_dna_table
from Bio.Data import CodonTable
from collections import Counter
```

## Codon Adaptation Index (CAI)

**Goal:** Measure how closely a gene's codon usage matches the highly expressed genes of a host organism.

**Approach:** Build a `CodonAdaptationIndex` from a reference set of highly expressed CDS (the constructor computes per-codon relative adaptiveness w), then score query sequences with `calculate` (0-1, higher = better adapted).

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex
from Bio.Data.CodonTable import standard_dna_table

# Reference = highly expressed genes of the TARGET host (ribosomal proteins, EFs).
# Parse a FASTA yourself; there is no bundled index. Pass str/Seq/SeqRecord.
reference_seqs = list(SeqIO.parse('highly_expressed_genes.fasta', 'fasta'))
cai = CodonAdaptationIndex(reference_seqs, table=standard_dna_table)

query = Seq('ATGAAACGTGCTGAAGCTAAATAA')
score = cai.calculate(query)   # 0-1; the query MUST be in-frame (see governing principle)
print(f'CAI: {score:.3f}')
```

`CodonAdaptationIndex` **is a dict subclass** -- the codon->w mapping is the object itself. Inspect or override weights directly:

```python
print(cai['GCT'])         # relative adaptiveness of Ala codon GCT
cai.update({'GCT': 0.9})  # override a weight (replaces the old set_cai_index)
```

Verified behavior (Biopython >=1.82):
- **ATG (Met) and TGG (Trp) are excluded** from CAI -- single-codon families, w is always 1.
- **Stop codons are excluded.**
- **Unobserved codons get w = 0.5** (Sharp & Li), softly down-weighted; no division-by-zero.
- **Case-insensitive** -- both the constructor and `calculate` uppercase internally.
- An illegal codon (non-ACGT) in a reference raises `ValueError`; an illegal or trailing-partial codon in a query raises `TypeError`. Out-of-frame input does NOT raise -- it is silently mis-scored.

## RSCU = w is built from these ratios

CAI weights come from RSCU. w_ij = RSCU_ij / RSCU_jmax = (codon count) / (count of the most-used synonymous codon in that family); CAI = exp((1/L) * sum ln w). RSCU itself = observed count / expected-if-uniform within a synonymous family (=1 no bias, >1 over-used, <1 under-used). RSCU normalizes away amino-acid composition, which is why w is built from RSCU ratios rather than raw frequencies.

**Goal:** Quantify synonymous codon bias to detect translational selection or mutational pressure.

**Approach:** Group codons by amino acid via the codon table, then divide each codon's observed count by the family mean.

```python
from Bio.Data import CodonTable
from collections import Counter

def count_codons(seq):
    s = str(seq).upper()
    return Counter(s[i:i+3] for i in range(0, len(s) - 2, 3))

def calculate_rscu(seq, table_id=1):
    '''RSCU per codon: observed / expected-if-uniform within its synonymous family'''
    table = CodonTable.unambiguous_dna_by_id[table_id]
    counts = count_codons(seq)
    back_table = {}
    for codon, aa in table.forward_table.items():
        back_table.setdefault(aa, []).append(codon)
    rscu = {}
    for aa, codons in back_table.items():
        total = sum(counts.get(c, 0) for c in codons)
        expected = total / len(codons) if codons else 0
        for codon in codons:
            rscu[codon] = counts.get(codon, 0) / expected if expected > 0 else 0
    return rscu
```

## Codon optimization for expression

**Goal:** Generate a host-adapted CDS that preserves the protein.

**Approach:** `optimize()` swaps each amino acid for the host's single most-preferred synonymous codon (max-CAI). Always confirm the protein is unchanged, then screen the design against the tradeoffs below.

```python
opt = cai.optimize(query, seq_type='DNA', strict=True)
assert opt.translate() == query.translate()   # protein must be identical
```

`optimize(sequence, seq_type='DNA'|'RNA'|'protein', strict=True)`: `strict=True` **raises ValueError on a tie** (two equally-preferred codons, e.g. `'TTT and TTC are equally preferred.'`); `strict=False` warns and picks one.

### Why naive max-CAI can HURT expression

`optimize()` is blind to everything except single-codon frequency. Screen the output for:
- **The translation ramp** (Tuller et al. 2010): a conserved profile of slow (rare) codons over the first ~30-50 codons spaces ribosomes; flattening it can lower yield and increase misfolding.
- **5' mRNA secondary structure:** strong folding near the start codon impedes initiation. Minimize 5' free energy, sometimes against CAI.
- **GC extremes:** swaps that push GC very high create stable hairpins; very low GC destabilizes.
- **Cryptic elements created by swaps:** splice sites, internal Shine-Dalgarno/RBS, polyadenylation signals, restriction sites, AU-rich destabilizing elements -- silent in protein, corrupting in expression.
- **Codon-pair bias:** decoding efficiency depends on adjacent codon pairs; CAI scores single codons only.

### tAI -- the supply-side alternative

The tRNA Adaptation Index (dos Reis et al. 2004) weights each codon by **tRNA gene copy number** (a proxy for tRNA abundance) scaled by wobble-pairing efficiency at the third position. Where tRNA copy number is a good abundance proxy, tAI tracks expression and elongation speed better than CAI. tAI is **not in Biopython** -- use the R `tAI` package or reimplement.

## Synonymous bias by other metrics

### Effective Number of Codons (Nc)

A reference-free bias measure (lower = more biased; range ~20 fully biased to 61 unbiased). The helper below is a simplified per-amino-acid approximation: Wright's published estimator averages the homozygosity F WITHIN each degeneracy class (2-, 3-, 4-, 6-fold) before combining as `Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6`. The endpoints agree, but intermediate values will not match codonW/standard Nc when families in a class have unequal F. For comparable Nc, average F by class per Wright (1990) or use codonW.

```python
import math
from Bio.Data import CodonTable

def effective_nc(seq, table_id=1):
    table = CodonTable.unambiguous_dna_by_id[table_id]
    counts = count_codons(seq)
    aa_groups = {}
    for codon, aa in table.forward_table.items():
        aa_groups.setdefault(aa, []).append(codon)
    nc_sum = 0
    for aa, codons in aa_groups.items():
        n = sum(counts.get(c, 0) for c in codons)
        if n <= 1:
            continue
        pi_sq = sum((counts.get(c, 0) / n) ** 2 for c in codons)
        F = (n * pi_sq - 1) / (n - 1)
        nc_sum += 1 / F if F > 0 else len(codons)
    return nc_sum if nc_sum > 0 else 61
```

### GC at codon positions (GC123)

```python
from Bio.SeqUtils import GC123

gc_total, gc_pos1, gc_pos2, gc_pos3 = GC123(seq)  # four PERCENTAGES (0-100)
print(f'GC3 (wobble): {gc_pos3:.1f}%')             # correlates with genome GC
```

GC123 returns percentages (0-100), unlike `gc_fraction` which returns 0-1. GC3 at the wobble position usually tracks overall genome GC content.

## Codon tables

```python
from Bio.Data import CodonTable

table = CodonTable.unambiguous_dna_by_id[1]   # standard genetic code
print(table.start_codons, table.stop_codons)
print(table.forward_table['ATG'])             # 'M'
```

| ID | Name | Organism |
|----|------|----------|
| 1 | Standard | Most nuclear genomes |
| 2 | Vertebrate Mitochondrial | Human/mouse mito |
| 4 | Mold/Protozoan Mitochondrial | Fungi, protozoa mito |
| 5 | Invertebrate Mitochondrial | Insects, worms mito |
| 11 | Bacterial/Plastid | E. coli, chloroplasts |

Pass the matching `table=` to `CodonAdaptationIndex` when scoring mitochondrial or bacterial genes.

## Metric reference

| Metric | Range | Reference needed | Interpretation |
|--------|-------|------------------|----------------|
| CAI | 0-1 | Highly expressed genes of the host | Higher = better adapted |
| RSCU | 0-N | None (within-sequence) | 1 = no bias, >1 over-used |
| Nc | ~20-61 | None | Lower = more biased |
| GC3 | 0-100% | None | GC at wobble position |
| tAI | 0-1 | tRNA gene copy numbers | Higher = better tRNA supply |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `ImportError: cannot import name 'CodonUsage'` | `Bio.SeqUtils.CodonUsage` removed in 1.82 | `from Bio.SeqUtils import CodonAdaptationIndex` |
| `AttributeError: 'CodonAdaptationIndex' object has no attribute 'generate_index'` | Old API on new class | Build in the constructor; score with `calculate` |
| Plausible CAI from a frame-shifted CDS | Out-of-frame input silently mis-scored | Confirm frame: length divisible by 3, starts at codon 1 |
| CAI near 1 for every gene | Reference set is whole-genome, not expression-biased | Use only highly expressed genes of the target host |
| `ValueError: ... equally preferred` | `optimize(strict=True)` hit a tie | Pass `strict=False`, or curate weights with `update` |
| High CAI but poor expression in the lab | Max-CAI ignores ramp / 5' structure / cryptic sites | Screen `optimize()` output; treat it as a draft |

## Related Skills

- transcription-translation - Translate CDS and select the correct codon table
- sequence-properties - GC123 and per-position GC content
- sequence-io/read-sequences - Parse reference CDS from FASTA/GenBank for CAI training
- database-access/entrez-fetch - Fetch highly expressed gene sets from NCBI for CAI references

## References

Sharp PM, Li WH (1987) The codon adaptation index -- a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Res 15(3):1281-1295.

dos Reis M, Savva R, Wernisch L (2004) Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res 32(17):5036-5044.

Tuller T, Carmi A, Vestsigian K, Navon S, Dorfan Y, Zaborske J, Pan T, Dahan O, Furman I, Pilpel Y (2010) An evolutionarily conserved mechanism for controlling the efficiency of protein translation. Cell 141(2):344-354.

Wright F (1990) The 'effective number of codons' used in a gene. Gene 87(1):23-29.
