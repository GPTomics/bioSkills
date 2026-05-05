---
name: bio-alignment-msa-statistics
description: Calculate alignment statistics including sequence identity, conservation scores, substitution matrices, and similarity metrics. Use when comparing alignment quality, measuring sequence divergence, and analyzing evolutionary patterns.
tool_type: python
primary_tool: Bio.Align
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# MSA Statistics

Calculate sequence identity, conservation scores, substitution counts, and other alignment metrics.

## Required Import

**Goal:** Load modules for alignment I/O, substitution scoring, and statistical calculations.

**Approach:** Import AlignIO for reading alignments, Counter for column analysis, numpy for matrix operations, and math for entropy calculations.

```python
from Bio import AlignIO
from Bio.Align import substitution_matrices
from collections import Counter
import numpy as np
import math
```

## Pairwise Identity

**"Calculate percent identity"** → Compute the fraction of identical aligned residues between sequence pairs.

**Goal:** Measure sequence similarity as percent identity for individual pairs or across all sequences in an alignment.

**Approach:** Count matching non-gap positions divided by total aligned positions; optionally compute a full N-by-N identity matrix.

### Percent Identity Definitions

There are four common denominators, producing **up to 11.5% difference** on the same alignment. Combined with different alignment algorithms, variation reaches 22%. Always report which method was used.

| Method | Denominator | Code |
|--------|-------------|------|
| PID1 | Aligned positions including internal gaps | `sum(a != '-' or b != '-' for a, b in zip(s1, s2))` |
| PID2 | Aligned residue pairs only (no gaps) | `sum(a != '-' and b != '-' for a, b in zip(s1, s2))` |
| PID3 | Shorter sequence length (ungapped) | `min(len(s1.replace('-', '')), len(s2.replace('-', '')))` |
| PID4 | Mean sequence length (ungapped) | `(len(s1.replace('-', '')) + len(s2.replace('-', ''))) / 2` |

PID2 always gives the highest value; PID4 correlates best with structural similarity (Raghava & Barton 2006 BMC Bioinf, r=0.86 with Q-score) and is recommended for evolutionary analyses.

**Length-asymmetry pathology:** When sequences differ greatly in length, PID4 and PID2 diverge sharply. Example: 80 matches between a 500-residue protein and a 100-residue domain fragment yields PID2 ~84% (matches over aligned residue pairs) but PID4 ~27% (matches over mean ungapped length). Neither is wrong; they answer different questions:
- PID2 -> "how similar is the aligned region?" (motif/domain detection)
- PID4 -> "how similar are the full sequences?" (structural similarity benchmarks)

For ortholog identification at the protein level (full-length, similar size), PID4 is recommended. For domain detection or fragment-vs-genome alignment, PID2 with explicit length annotation is more interpretable. Always report alignment length alongside any percent identity to disambiguate.

### Calculate Identity Between Two Sequences
```python
def pairwise_identity(seq1, seq2, method='pid1'):
    matches = sum(a == b and a != '-' for a, b in zip(seq1, seq2))
    if method == 'pid1':
        denom = sum(a != '-' or b != '-' for a, b in zip(seq1, seq2))
    elif method == 'pid2':
        denom = sum(a != '-' and b != '-' for a, b in zip(seq1, seq2))
    elif method == 'pid3':
        denom = min(len(seq1.replace('-', '')), len(seq2.replace('-', '')))
    elif method == 'pid4':
        denom = (len(seq1.replace('-', '')) + len(seq2.replace('-', ''))) / 2
    return matches / denom if denom > 0 else 0

alignment = AlignIO.read('alignment.fasta', 'fasta')
seq1, seq2 = str(alignment[0].seq), str(alignment[1].seq)
for method in ['pid1', 'pid2', 'pid3', 'pid4']:
    print(f'{method}: {pairwise_identity(seq1, seq2, method) * 100:.1f}%')
```

### Identity Matrix for All Sequences

The double-loop is O(N^2 * L) and fine for hundreds of sequences; for thousands, vectorize via numpy broadcasting:

```python
def identity_matrix_vectorized(alignment):
    # Build N x L character array; for each row, broadcast equality and validity masks against all rows
    ...
```

Full implementation: `examples/identity_matrix.py`. For very large alignments (>10k sequences), switch to k-mer-based distance estimation (e.g. mash) -- exact pairwise identity becomes prohibitive.

## Conservation Scoring Methods

Pick a conservation score by what the downstream task needs:

| Pick this | When |
|-----------|------|
| Majority fraction or Shannon entropy | Quick screening; DNA/RNA logos; coarse column QC |
| Capra-Singh JSD (modern default) | Catalytic-residue / functional-site prediction on protein MSAs |
| ConSurf rate4site | PDB surface mapping when a phylogenetic tree is available |

## Conservation Score

**Goal:** Quantify per-column and overall alignment conservation to identify conserved and variable regions.

**Approach:** Calculate the fraction of the most common residue at each column, optionally ignoring gaps, and smooth with a sliding window.

### Per-Column Conservation
```python
def column_conservation(alignment, col_idx, ignore_gaps=True):
    column = alignment[:, col_idx]
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    most_common_count = counts.most_common(1)[0][1]
    return most_common_count / len(column)

alignment = AlignIO.read('alignment.fasta', 'fasta')
for i in range(min(20, alignment.get_alignment_length())):
    cons = column_conservation(alignment, i)
    print(f'Column {i}: {cons*100:.0f}% conserved')
```

### Average Conservation Across Alignment
```python
def average_conservation(alignment, ignore_gaps=True):
    scores = []
    for col_idx in range(alignment.get_alignment_length()):
        scores.append(column_conservation(alignment, col_idx, ignore_gaps))
    return sum(scores) / len(scores)

avg_cons = average_conservation(alignment)
print(f'Average conservation: {avg_cons*100:.1f}%')
```

### Conservation Profile
```python
def conservation_profile(alignment, window=10):
    profile = []
    for i in range(alignment.get_alignment_length()):
        start = max(0, i - window // 2)
        end = min(alignment.get_alignment_length(), i + window // 2)
        scores = [column_conservation(alignment, j) for j in range(start, end)]
        profile.append(sum(scores) / len(scores))
    return profile

profile = conservation_profile(alignment, window=10)
```

### Capra-Singh Jensen-Shannon Divergence

**Goal:** Score columns by divergence from a residue-frequency background, with a window-smoothed neighbour penalty for catalytic residue prediction.

**Approach:** Compute JSD between the column distribution and a residue-frequency background (Capra & Singh 2007 Bioinf used BLOSUM62-derived; the example uses Robinson & Robinson 1991 PNAS, which gives effectively-equivalent column ranking), then mix with the windowed-neighbour mean. Defaults `window=3`, `lambda_window=0.5` track catalytic-residue annotation in the Catalytic Site Atlas. Full implementation: `examples/capra_singh_jsd.py`.

```python
def capra_singh_score(alignment, background=None, window=3, lambda_window=0.5):
    # raw[i] = JSD(column_i, background) * (1 - gap_fraction_i)   # gap-penalty per Capra-Singh reference impl
    # smoothed[i] = (1 - lambda) * raw[i] + lambda * mean(raw[i-window:i+window+1] excluding i)
    ...
```

**Threshold for catalytic-residue prediction:** Capra & Singh 2007 (Bioinf 23:1875) report AUC ~0.94 and Top-30 score ~0.75 on the Catalytic Site Atlas using JSD with neighbor mixing (window=3, lambda=0.5); the paper does not prescribe a single threshold. Choose by ROC tradeoff for the specific use case. ConSurf-derived rate4site (rate-of-evolution) is competitive but requires a phylogenetic tree; Capra-Singh JSD is the alignment-only equivalent.

## Substitution Counts

**Goal:** Tabulate observed substitution frequencies from the alignment for evolutionary analysis or custom scoring matrices.

**Approach:** Enumerate all pairwise non-gap character comparisons at each column and tally substitution pairs.

### Count Substitutions from Alignment

```python
def substitution_counts(alignment):
    # Tally pairwise non-gap residue mismatches across all column-pair comparisons
    ...
```

Full implementation: `examples/substitution_counts.py`. The script also reports Ti/Tv ratio for DNA alignments.

### Built-in Pairwise Substitutions

For pairwise alignments created with `PairwiseAligner`, use the `.substitutions` property:

```python
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=-1)
print(aligner.align(seq1, seq2)[0].substitutions)
```

### BLOSUM62 Lambda Is Not a Single Number

Different tools calibrate Karlin-Altschul lambda differently (NCBI BLAST tabulates 0.3176; FASTA/SSEARCH recomputes per query; HMMER `phmmer` derives it via Forward calibration; Bio.Align does not implement Karlin-Altschul). Expect ~2% bit-score variation across tools on the same alignment. Always record tool and version when citing bit scores; for borderline (~30 bit) hits this matters.

## Information Content

**Goal:** Measure column variability using Shannon entropy and derive information content for identifying functionally important positions.

**Approach:** Compute Shannon entropy from character frequencies per column; information content is max entropy minus observed entropy.

### Shannon Entropy Per Column
```python
import math

def shannon_entropy(column, ignore_gaps=True):
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    total = len(column)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy

alignment = AlignIO.read('alignment.fasta', 'fasta')
for i in range(min(20, alignment.get_alignment_length())):
    column = alignment[:, i]
    ent = shannon_entropy(column)
    print(f'Column {i}: entropy = {ent:.2f} bits')
```

### Information Content (Kullback-Leibler Divergence)

The classic uniform-background formulation `IC = log2(alphabet_size) - H` (Schneider & Stephens 1990 NAR) is only valid when the genomic background is uniform. This is approximately true for random DNA but emphatically wrong for protein, where amino acid frequencies range from 1.4% (Trp) to 9.4% (Leu). For amino acids, use Kullback-Leibler divergence `IC = sum_i p_i * log2(p_i / b_i)` against the Robinson & Robinson 1991 PNAS empirical background. Full implementation: `examples/entropy_analysis.py`.

```python
ROBINSON_BACKGROUND = {
    'A': 0.0780, 'R': 0.0512, 'N': 0.0427, 'D': 0.0530, 'C': 0.0193,
    'Q': 0.0419, 'E': 0.0629, 'G': 0.0738, 'H': 0.0224, 'I': 0.0526,
    'L': 0.0922, 'K': 0.0596, 'M': 0.0224, 'F': 0.0399, 'P': 0.0508,
    'S': 0.0712, 'T': 0.0584, 'W': 0.0133, 'Y': 0.0327, 'V': 0.0653,
}
DNA_UNIFORM = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

def information_content(column, background, ignore_gaps=True):
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    total = len(column)
    return sum((c / total) * math.log2((c / total) / background.get(r, 1e-9))
               for r, c in counts.items() if c > 0)
```

For sequence-logo letter heights, use the Schneider-Stephens form (`letter_height = p_i * (log2(alphabet) - H_observed)`, uniform background) when the comparison is "informative vs random"; use the KL form (against an empirical background such as Robinson-Robinson 1991) for protein logos where amino-acid frequencies are non-uniform or when the comparison is "informative vs the proteome". When the background is unknown, default to the empirical alignment composition rather than uniform.

## Gap Statistics

**Goal:** Summarize gap distribution across the alignment to assess alignment quality and identify problematic regions.

**Approach:** Calculate gap fractions per column and aggregate statistics including total gaps, gap-free columns, and gappiest sequence/column.

### Gap Fraction Per Column
```python
def gap_profile(alignment):
    profile = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        gap_fraction = column.count('-') / len(alignment)
        profile.append(gap_fraction)
    return profile

gaps = gap_profile(alignment)
avg_gaps = sum(gaps) / len(gaps)
print(f'Average gap fraction: {avg_gaps*100:.1f}%')
```

### Gap Statistics Summary

```python
def gap_statistics(alignment):
    # total gaps, gap fraction, gappiest sequence and column indices, gap-free column count
    ...
```

Full implementation: `examples/gap_statistics.py`.

## Alignment Quality Metrics

**Goal:** Score alignment quality using sum-of-pairs or simple match/mismatch/gap scoring across all columns.

**Approach:** For each column, score all pairwise residue comparisons and sum across the alignment.

### Overall Alignment Score
```python
def alignment_score(alignment, match=1, mismatch=-1, gap=-2):
    total_score = 0
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        for i, c1 in enumerate(column):
            for c2 in column[i+1:]:
                if c1 == '-' or c2 == '-':
                    total_score += gap
                elif c1 == c2:
                    total_score += match
                else:
                    total_score += mismatch
    return total_score

score = alignment_score(alignment)
print(f'Alignment score: {score}')
```

### Sum of Pairs Score

Biopython's `substitution_matrices.load('BLOSUM62')` returns a `Bio.Align.substitution_matrices.Array` object (a numpy-backed 2D array indexed by residue characters), NOT a dict. The correct accessor is `matrix[c1, c2]`; calling `.get((c1, c2), 0)` on an Array silently always returns 0. Standard BLOSUM62 includes `B`, `Z`, `X`, and `*`; pairs containing residues outside the matrix alphabet (e.g. `U` selenocysteine, `J` Leu/Ile) raise `IndexError` and should be skipped or scored zero.

```python
def sum_of_pairs(alignment, substitution_matrix=None):
    if substitution_matrix is None:
        substitution_matrix = substitution_matrices.load('BLOSUM62')

    total = 0.0
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        for i, c1 in enumerate(column):
            for c2 in column[i+1:]:
                if c1 == '-' or c2 == '-':
                    continue
                try:
                    total += substitution_matrix[c1, c2]
                except (KeyError, IndexError):
                    continue
    return total
```

**SP-score is biased on unbalanced datasets.** The above implementation gives equal weight to every sequence pair. On phylogenetically structured datasets (e.g. 95 mammals + 5 outgroups), 99% of pairs are mammal-mammal and the SP score reports only mammal-internal alignment quality. MUSCLE and T-Coffee internally compute weighted SP using sequence weights (Henikoff or position-based) so pair contributions are downweighted by cluster redundancy. For SP-as-quality-score on real data, multiply each pair contribution by `weight[i] * weight[j]` from the Henikoff weights in `alignment/msa-parsing` (`examples/henikoff_weights.py`).

## Position-Specific Score Matrix (PSSM)

**Goal:** Build a position-specific scoring matrix from the alignment for motif analysis or sequence scoring.

**Approach:** Raw counts give frequencies; without pseudocounts, log-odds against background diverge to negative infinity at any column missing a residue. Henikoff JG & Henikoff S 1996 (Bioinf 12:135-143) introduced data-dependent pseudocount weighting; a simple Laplace add-one is the minimal correct approach for production use. Full implementation: `examples/pssm.py`.

```python
def pssm_with_pseudocounts(alignment, background, pseudocount=1.0):
    # log2((counts[r] + pc * background[r]) / (n + pc) / background[r]) per column
    ...
```

`pseudocount=1.0` is the Laplace prior; HMMER uses Dirichlet mixtures for sophisticated smoothing. For motif scanning, score a candidate site by summing per-position log-odds; sites above a calibrated threshold are predicted hits. Use `ROBINSON_BACKGROUND` (defined in the IC section above) for protein.

## Effective Sequence Number (Neff)

**Goal:** Estimate non-redundant sequence count for MSA-depth metrics.

**Approach:** Cluster at an identity threshold (0.62 protein, 0.80 nucleotide) and weight by inverse cluster size; reference implementation lives in `msa-parsing` (`examples/neff.py`). `Neff/L > 0.5` is the rule-of-thumb for direct-coupling-analysis contact prediction; AlphaFold's MSA-depth scoring uses a closely related metric.

## Mutual Information with APC

**Goal:** Detect coevolving column pairs as a coupling/contact signal.

**Approach:** Pairwise MI minus average-product correction (Dunn, Wahl, Gloor 2008 Bioinf). Reference implementation lives in `msa-parsing` (`examples/mi_apc.py`). For production-grade contact prediction beyond a few hundred columns, switch to plmDCA (Ekeberg et al 2013) or EVcouplings (Hopf et al 2017).

## Distance Correction Models

For publication-grade pairwise distances, do NOT pick a model by rule of thumb. Run ModelTest-NG to select the best-fit substitution model by AIC/BIC, then apply that correction via IQ-TREE2 (`.mldist` output) or EMBOSS `distmat`. Hand-coded JC69 / K80 / blosum62 corrections via `Bio.Phylo.TreeConstruction.DistanceCalculator` are appropriate only for exploratory work.

```bash
modeltest-ng -i alignment.fasta -d nt -t ml
modeltest-ng -i alignment.fasta -d aa -t ml -p 4
```

```python
from Bio.Phylo.TreeConstruction import DistanceCalculator

calculator = DistanceCalculator('blosum62')
distance_matrix = calculator.get_distance(alignment)
```

For substitution-model selection in the context of full phylogenetic inference, see `phylogenetics/modern-tree-inference`.

## Alignment Quality Assessment

### When to Worry About Alignment Quality

| Warning Sign | Implication | Action |
|-------------|-------------|--------|
| Average pairwise identity <25% (protein) | Twilight zone; alignment may be unreliable | Use GUIDANCE2 to assess; consider structural alignment |
| >30% of columns have >50% gaps | Possible non-homologous sequences or misalignment | Remove outlier sequences and re-align |
| Identity varies dramatically across regions | Domain architecture mismatch | Align domains separately |
| Conservation pattern absent in expected functional regions | Alignment error or non-homology | Verify with BLAST that sequences are truly homologous |

### Quantifying Alignment Uncertainty

Alignment uncertainty propagates directly into evolutionary inference; for selection analysis (dN/dS), misaligned codons create artificial nonsynonymous differences and false positive signals. For column-confidence quantification (GUIDANCE2, MUSCLE5 ensemble, T-Coffee TCS), see the Confidence Assessment section in `alignment/multiple-alignment`.

## Note on Bio.Align.AlignInfo

The `AlignInfo.SummaryInfo` class is **deprecated** in recent Biopython versions. Use the custom functions in this skill instead:
- For PSSM: use `pssm_with_pseudocounts()` above
- For information content: use `information_content()` function earlier in this skill
- For consensus: see msa-parsing skill

## Quick Reference: Metrics

| Metric | Description | Range |
|--------|-------------|-------|
| Identity | Fraction of identical residues | 0-1 |
| Conservation | Most common residue frequency | 0-1 |
| Shannon Entropy | Variability measure | 0 to log2(alphabet) |
| Information Content | Max entropy - observed entropy | 0 to log2(alphabet) |
| Gap Fraction | Proportion of gaps | 0-1 |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `ZeroDivisionError` | Empty column after gap removal | Check for gap-only columns |
| `KeyError` | Character not in substitution matrix | Handle gaps separately |
| Negative IC | Wrong alphabet size | Use 4 for DNA, 20 for protein |

## Related Skills

- alignment/multiple-alignment - Run MSA tools and quantify alignment confidence with ensembles
- alignment/msa-parsing - Parse, filter, trim, and assess alignment quality (Henikoff weights, Neff, MI-APC live there)
- alignment/alignment-io - Read/write alignment files
- alignment/pairwise-alignment - Create and score pairwise alignments
- alignment/alignment-trimming - Column trimming before downstream statistics
- alignment/structural-alignment - Twilight-zone alternative when sequence MSA is unreliable
- phylogenetics/distance-calculations - Distance models and tree building from corrected distances
- sequence-manipulation/sequence-properties - Sequence-level statistics

## References

- Schneider TD, Stephens RM. 1990. Sequence logos: a new way to display consensus sequences. NAR 18:6097-6100.
- Robinson AB, Robinson LR. 1991. Distribution of glutamine and asparagine residues and their near neighbors in peptides and proteins. PNAS 88:8880-8884.
- Capra JA, Singh M. 2007. Predicting functionally important residues from sequence conservation. Bioinf 23:1875-1882.
- Henikoff JG, Henikoff S. 1996. Using substitution probabilities to improve position-specific scoring matrices. Bioinf 12:135-143.
- Pei J, Grishin NV. 2001. AL2CO: calculation of positional conservation in a protein sequence alignment. Bioinf 17:700-712.
- Mayrose I, Graur D, Ben-Tal N, Pupko T. 2004. Comparison of site-specific rate-inference methods for protein sequences: empirical Bayesian methods are superior. MBE 21:1781-1791.
- Valdar WSJ. 2002. Scoring residue conservation. Proteins 48:227-241.
- Raghava GPS, Barton GJ. 2006. Quantification of the variation in percentage identity for protein sequence alignments. BMC Bioinf 7:415.
- Dunn SD, Wahl LM, Gloor GB. 2008. Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. Bioinf 24:333-340.
- Altschul SF et al. 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. NAR 25:3389-3402.
- Pearson WR. 2013. An introduction to sequence similarity ("homology") searching. Curr Protoc Bioinf 3.1.
- Eddy SR. 2008. A probabilistic model of local sequence alignment that simplifies statistical significance estimation. PLOS CB 4:e1000069.
- Edgar RC. 2004. MUSCLE: multiple sequence alignment with high accuracy and high throughput. NAR 32:1792-1797.
- Darriba D et al. 2020. ModelTest-NG: a new and scalable tool for the selection of DNA and protein evolutionary models. MBE 37:291-294.
