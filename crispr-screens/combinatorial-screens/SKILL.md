---
name: bio-crispr-screens-combinatorial-screens
description: Designs and analyzes combinatorial CRISPR screens covering paired-Cas9 (Big Papi, Najm 2018), enhanced AsCas12a multiplex (enCas12a, DeWeirdt 2021), in4mer 4-guide-array Cas12a (Esmaeili Anvar N et al 2024 Nat Commun 15:3577) and the Inzolia paralog-pair library, paralog-buffering detection (Dede 2020 Genome Biol; Thompson 2021 Cell Reports 36:109597), genetic-interaction (GI) scoring as observed_double_LFC minus expected_additive_double_LFC, synthetic-lethal and synthetic-rescue interaction interpretation, the half-of-essentiality buffered by paralogs phenomenon, multiplex screen statistical analysis with MAGeCK MLE interaction terms, and the relationship to single-cell combinatorial Perturb-seq. Use when designing a paralog or pathway-pair screen, choosing between paired-Cas9 (Big Papi) and Cas12a multiplex (Inzolia), interpreting genetic interaction scores, identifying synthetic-lethal targets for drug development, or scaling beyond single-gene CRISPR screens.
tool_type: mixed
primary_tool: enCas12a
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5.9+ (for MLE with interaction terms), Inzolia library annotation (Bayle 2024), pandas 2.2+, numpy 1.26+, scipy 1.12+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version`; `mageck mle --help`
- For Cas12a libraries: verify against published Inzolia / in4mer / Big Papi annotations

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Combinatorial CRISPR Screen Analysis

**"Run a combinatorial CRISPR screen to find synthetic-lethal interactions"** -> Design a paired or multiplex library, screen for double-knockout fitness, score per-pair genetic interaction (GI = observed_double - expected_additive), and identify synthetic-lethal (negative GI) and synthetic-rescue (positive GI) interactions.

- CLI: `mageck mle` with explicit interaction terms for paired-Cas9 (Big Papi-style)
- Python: custom GI scoring for Cas12a multiplex (in4mer / Inzolia)
- Modality: enCas12a / LbCas12a single-array multiplex (preferred for paralog screens)

## Combinatorial Architecture Decision Tree

| Goal | Architecture | Library | Why |
|------|--------------|---------|-----|
| Paralog buffering, identify synthetic lethal paralog pairs | enCas12a single-array 4-guide multiplex | Inzolia (Bayle 2024) | Cas9 single-KO misses ~half of paralog-buffered essentials |
| Test specific pathway pair (e.g., DNA repair branches) | Big Papi (paired-Cas9 dual sgRNA cassette) | Custom | Mature methodology; SpCas9 well-characterized |
| Combinatorial 3-way / 4-way knockout | in4mer (4-guide single Cas12a array) | Custom (in4mer) | Single transcript processed by Cas12a; multi-gene |
| Single-cell Perturb-seq with multi-pert per cell | Combinatorial Perturb-seq + Cas9 multiplex | Custom | Single-cell readout of multi-perturbation effects |
| Drug-modifier + KO interaction | Cas9 KO + drug treatment | Standard libraries | Drug as second "perturbation" |

**Fails when:**
- Big Papi without proper paired-guide cloning: silent fusion of two guides into one cassette gives single perturbation
- Cas12a screens analyzed as Cas9 screens: MAGeCK normalization fails because Cas12a has different cut profile
- in4mer 4-guide arrays without all-singleton controls: GI scoring requires single-gene baselines

## Cas9 vs Cas12a for Multiplex

| Property | Cas9 paired (Big Papi) | Cas12a multiplex (in4mer / Inzolia) |
|----------|------------------------|--------------------------------------|
| Multiplex capacity per cassette | 2 sgRNAs (paired) | 4 (in4mer); 2 (standard Cas12a) |
| sgRNA processing | Two separate U6 promoters | Single transcript processed by Cas12a itself |
| sgRNA inhibition with multiple targets | None | None (Cas12a's intrinsic processing handles all) |
| Library size for 1,000 pairs | ~2,000 paired cassettes | ~250 4-guide cassettes (in4mer) |
| Validated libraries | Limited (mostly custom) | Inzolia: 18k 4-guide arrays for 4k paralog pairs |
| Per-perturbation editing efficiency | High (each sgRNA independently) | Variable (Cas12a less efficient on some targets) |
| Best for | Pairwise GI of specific interest | Genome-scale paralog buffering; multi-gene perturbation |

**Recommendation:** For modern paralog screens, use Cas12a multiplex with the Inzolia library. The 30% library-size reduction vs paired-Cas9 makes it more cost-effective for genome-scale.

## The Paralog Buffering Phenomenon

**Dede et al 2020 *Genome Biol* 21:262; Thompson et al 2021** demonstrated that approximately half of constitutively-expressed essential genes are never detected in Cas9 single-KO screens. The reason: gene paralogs perform redundant essential functions. Loss of one paralog is buffered by the other; only loss of both creates the essentiality phenotype.

**Quantified impact:** ~24+ synthetic-lethal paralog pairs identified in Dede 2020 across 3 cell lines; 79% reproduce in ≥2 lines, 58% in all 3. These pairs were not findable by single-gene Cas9 screens, requiring combinatorial methodology.

**Examples:**
- **MAPK1 (ERK2) + MAPK3 (ERK1):** ERK family redundancy in proliferation
- **PIK3CA + PIK3CB:** PI3K alpha/beta redundancy
- **AKT1 + AKT2:** AKT family redundancy
- **HSP90AA1 + HSP90AB1:** HSP90 alpha/beta redundancy
- **STAG1 + STAG2:** Cohesion complex paralogs

Each is buffered: loss of one is tolerated; loss of both is lethal.

## Genetic Interaction (GI) Scoring

**Goal:** Identify pairs where the double-knockout fitness differs from the additive expectation.

**Approach:** From per-pair and per-singleton fitness data, compute GI = observed_double_LFC - (single_A_LFC + single_B_LFC). Synthetic lethal: GI < threshold (more depleted than additive). Synthetic rescue: GI > threshold (less depleted than additive).

```python
import pandas as pd
import numpy as np
from scipy.stats import zscore

def gi_score(paired_lfc_df, single_lfc_df):
    '''Score genetic interactions from paired vs single LFCs.

    paired_lfc_df: rows = paired-KO; columns = ['gene_A', 'gene_B', 'paired_lfc']
    single_lfc_df: rows = single-KO; columns = ['gene', 'single_lfc']
    '''
    single = dict(zip(single_lfc_df['gene'], single_lfc_df['single_lfc']))
    df = paired_lfc_df.copy()
    df['single_A_lfc'] = df['gene_A'].map(single)
    df['single_B_lfc'] = df['gene_B'].map(single)
    df['expected_additive'] = df['single_A_lfc'] + df['single_B_lfc']
    df['gi_score'] = df['paired_lfc'] - df['expected_additive']
    df['gi_z'] = zscore(df['gi_score'])
    df['gi_class'] = np.where(df['gi_z'] < -2, 'synthetic_lethal',
                                np.where(df['gi_z'] > 2, 'synthetic_rescue', 'no_interaction'))
    return df.sort_values('gi_z')
```

**Interpretation:**
- GI z-score < -2: Synthetic lethal (double-KO more lethal than expected) -- candidate drug target combinations
- GI z-score > 2: Synthetic rescue (double-KO less lethal than expected) -- compensatory pathway / paradoxical hit
- GI z-score -1 to 1: No interaction; effects are additive

## Run Combinatorial Screen Analysis (MAGeCK MLE with Interaction Indicator)

**Goal:** Use MAGeCK MLE to estimate the effect of each gene independently and the additional effect when both genes are simultaneously perturbed.

**Approach:** Design matrix encodes single-A, single-B, double-AB conditions; the `interaction` column is set to 1 only for double-KO samples. The resulting beta for that column captures the extra effect beyond the sum of single-gene betas. Note: MAGeCK MLE does not natively perform a formal interaction-significance test, but the `interaction|beta` and `|fdr` columns serve as the GI estimate; for formal interaction testing, compute GI = observed_double_lfc - (single_A_lfc + single_B_lfc) explicitly (see GI scoring section below).

```bash
# Design matrix encoding double-KO as a separate "interaction" indicator
# Conditions: NT (control), A_KO, B_KO, A_B_KO
cat > combo_design.txt <<EOF
Samples         baseline    geneA       geneB       interaction
NT_r1           1           0           0           0
NT_r2           1           0           0           0
A_r1            1           1           0           0
A_r2            1           1           0           0
B_r1            1           0           1           0
B_r2            1           0           1           0
AB_r1           1           1           1           1
AB_r2           1           1           1           1
EOF

mageck mle \
    --count-table combo_counts.txt \
    --design-matrix combo_design.txt \
    --output-prefix combo_mle

# Output: per-gene beta scores per design column
# The "interaction" column beta captures additional joint effect beyond additive
```

**Interpretation of MAGeCK MLE output:**

| Column | Meaning |
|--------|---------|
| `geneA|beta` | Single-A effect |
| `geneB|beta` | Single-B effect |
| `interaction|beta` | Additional effect under joint perturbation beyond sum of singles |
| `interaction|p-value`, `|fdr` | Significance vs zero |

A significantly negative `interaction|beta` is synthetic lethal; positive is synthetic rescue. For formal GI hypothesis testing, prefer the explicit GI scoring approach (next section) over MAGeCK MLE interpretation, since MAGeCK MLE does not validate the additive null.

## Inzolia / in4mer 4-Guide Array Analysis

**Bayle 2024 *Nat Commun* 15:3375** introduced in4mer, a Cas12a multiplex library where each cassette contains 4 guides processed by Cas12a's intrinsic crRNA-processing activity. The Inzolia library is the canonical implementation covering ~4,000 paralog pairs.

**Library design:**
- 4 guides per cassette (Cas12a single-transcript array)
- Per pair: 4 cassettes with different guide combinations targeting the same paralog pair
- Includes singleton controls: each gene targeted alone in 4-guide arrays with placeholder guides
- ~18,000 total cassettes covering ~4,000 paralog pairs

```python
# Per-pair analysis from in4mer screen
def in4mer_pair_analysis(paired_counts_df, gene_pairs):
    '''Aggregate cassette-level counts to per-pair LFCs.
    paired_counts_df: rows = cassettes; columns = sample counts.
    gene_pairs: DataFrame with cassette_id -> [gene_A, gene_B, gene_C, gene_D]
    '''
    # Aggregate by (gene_A, gene_B) pair if 2-gene pair, or by all 4 genes if multi-gene
    pair_lfc = paired_counts_df.merge(gene_pairs, on='cassette_id').groupby(['gene_A', 'gene_B'])
    return pair_lfc.agg(['mean', 'std', 'count'])
```

## Failure Modes

### Big Papi cassette fuses sgRNAs

**Trigger:** Improperly designed dual sgRNA cassette where two sgRNAs are read as one fused sequence.
**Mechanism:** Without spacer or terminator between sgRNAs, transcription doesn't terminate correctly.
**Symptom:** Cassette appears as single perturbation (the first sgRNA dominates); GI scoring fails.
**Fix:** Use validated Big Papi or paired-sgRNA cloning protocols; verify by amplicon sequencing of clones.

### Cas12a screen with low editing efficiency

**Trigger:** Cas12a less efficient than Cas9 at some loci; some guides in the 4-guide array don't cut.
**Mechanism:** Cas12a editing rate varies by sequence context; some loci edit at <30%.
**Symptom:** Specific pairs missing expected effects despite cassette presence.
**Fix:** Pilot Cas12a efficiency at the loci before full screen; use enCas12a (enhanced) variant; for known low-efficiency loci, supplement with Cas9.

### GI scoring without singletons

**Trigger:** Library lacks single-gene controls (only paired knockouts).
**Mechanism:** GI = paired - expected_additive requires single-gene LFC; without them, expected cannot be computed.
**Symptom:** Cannot score GI; only paired LFCs available.
**Fix:** Design library to include singletons (place gene A with 3 placeholder guides; gene B with 3 placeholders); re-run with full design.

### Single-gene LFCs from different cell line

**Trigger:** Using public single-gene LFCs (e.g., DepMap) as the baseline for paired-screen GI scoring.
**Mechanism:** Single-gene effects are cell-line specific; using HCT116 single-gene LFCs to score K562 paired-screen GIs is invalid.
**Symptom:** GI scores look noisy; many false positives.
**Fix:** Include singleton controls in the screen; or use cell-line-matched DepMap data.

### Confounding cell-cycle / proliferation in GI scoring

**Trigger:** Paired KO of two cell-cycle-impacting genes; the double-effect saturates cell cycle.
**Mechanism:** If A_KO causes 50% growth arrest and B_KO causes 50%, the combined 75% arrest is already saturating proliferation; additive expectation overestimates double-effect, generating false "synthetic-rescue."
**Symptom:** GI scores positive for pairs of essential cell-cycle genes; biologically unexpected.
**Fix:** Use log-space (LFC) GI scoring rather than linear; saturation is less severe in log-space. Alternative: model with logistic / saturable response curve.

### Library skew amplifying noise

**Trigger:** Inzolia library has uneven cassette representation; some pairs at 10x lower coverage than others.
**Mechanism:** Standard library QC (Gini, skew) applies; low-coverage cassettes yield noisier LFCs.
**Symptom:** GI z-scores vary 2-3x across cassettes targeting the same pair.
**Fix:** Standard library QC; for low-coverage pairs, aggregate fewer cassettes but with more sequencing depth; or drop low-coverage pairs from analysis.

## Cross-Modality Validation

For high-stakes synthetic-lethal hits (drug-target nomination), validate by:

1. **Orthogonal chemistry:** Re-validate with Cas9 if Cas12a, or vice versa
2. **Arrayed validation:** Single-knock-out arrayed setup with same cell line; quantify proliferation
3. **CRISPRi orthogonal:** Use dCas9-KRAB to confirm knockdown phenotype (no DNA damage)
4. **Pharmacological:** Inhibit paralog with drug; confirms target accessibility for drug development

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Synthetic lethal GI z-score | <-2 | Standard convention |
| Synthetic rescue GI z-score | >2 | Standard convention |
| No interaction | -1 to +1 | Within additive expectation |
| Cas9 paired-screen cassette count per pair | 4-6 | Standard library convention |
| Cas12a 4-guide arrays per pair (Inzolia) | 4 | Bayle 2024 |
| Singletons in combinatorial library | At least 4-6 per single gene | For stable expected_additive |
| Cells per cassette for stable GI | 500+ at infection | Standard pooled-screen coverage |
| Cas12a editing efficiency for inclusion | >50% | Below = unreliable signal |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Big Papi cassette acts as single | Fused sgRNAs in cassette | Re-derive validated paired-sgRNA protocol |
| Cas12a low editing | Locus-specific inefficiency | Pilot loci first; use enCas12a |
| Cannot compute GI | No singletons in library | Re-design to include all-singletons |
| GI scores noisy | Library skew | Standard library QC; aggregate cassettes |
| Many false "rescue" GIs | Saturation in linear-space | Use log-space (LFC) GI scoring |
| Drug-target paralog shows no GI in screen | Cell-line-specific buffering | Cross-validate with multiple lines |

## References

- Najm FJ et al. 2018. *Nat Biotechnol* 36:179. Big Papi paired-Cas9 platform.
- DeWeirdt PC et al. 2021. *Nat Biotechnol* 39:94. enAsCas12a multiplex.
- Esmaeili Anvar N et al. 2024. *Nat Commun* 15:3577. in4mer / Inzolia paralog library.
- Dede M et al. 2020. *Genome Biol* 21:262. Paralog buffering in Cas9 screens.
- Thompson NA et al. 2021. *Cell Reports* 36:109597. Paralog synthetic lethals validated. (Earlier "Cell Syst 12:1132" attribution conflated with Pacini 2021.)
- Boettcher M et al. 2018. *Nat Biotechnol* 36:170. Original paired-Cas9 combinatorial screen.
- Horlbeck MA et al. 2018. *Cell* 174:953-967. CRISPRi/a / combinatorial work; paralog buffering as a co-essentiality pattern was characterized more directly in Dede 2020 *Genome Biol* 21:262 and Gonatopoulos-Pournatzis 2020 *Nat Biotechnol* 38:638.

## Related Skills

- crispr-screens/library-design - Inzolia / in4mer / Big Papi library design
- crispr-screens/screen-qc - Library QC including cassette skew
- crispr-screens/mageck-analysis - MAGeCK MLE with interaction terms
- crispr-screens/hit-calling - Cross-method analysis of combinatorial data
- crispr-screens/perturb-seq-analysis - Combinatorial Perturb-seq
- crispr-screens/copy-number-correction - Pre-correction for cancer-line combinatorial screens
- crispr-screens/in-vivo-screens - In-vivo paralog screens
- pathway-analysis/go-enrichment - Functional analysis of GI clusters
