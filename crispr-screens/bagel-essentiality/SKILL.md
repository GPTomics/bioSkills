---
name: bio-crispr-screens-bagel-essentiality
description: Identifies essential genes from CRISPR-Cas9 fitness screens using BAGEL2 (Kim & Hart 2021 Genome Med), a Bayesian classifier scoring per-gene Bayes Factors via log-likelihood ratios over per-sgRNA fold changes, calibrated against CEGv2 core-essentials (Hart 2017 G3, ~684 genes) and NEGv1 non-essentials (Hart 2014, ~927 genes). Covers the fc + bf + pr workflow, the linear-extrapolation improvement over BAGEL1 truncation, multi-target off-target correction, tumor-suppressor sensitivity (BAGEL2 detects enrichment), and BF-to-FDR calibration (BF >6 ≈ FDR 0.05 from Hart 2017). Use when classifying essential vs non-essential genes, calibrating BAGEL2 thresholds against PR curves, identifying tumor suppressors alongside essentials, comparing BAGEL2 hits to MAGeCK / drugZ, or generating publication-quality essentiality calls.
tool_type: cli
primary_tool: BAGEL2
---

## Version Compatibility

Reference examples tested with: BAGEL2 1.0.5+ (hart-lab/bagel), pandas 2.2+, numpy 1.26+, scipy 1.12+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `BAGEL.py fc --help`; `BAGEL.py bf --help`; `BAGEL.py pr --help`
- Python: BAGEL2 is distributed via `git clone` (no canonical PyPI release); confirm `python BAGEL.py --version` after checkout.

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## BAGEL2 Essentiality Analysis

**"Identify essential genes from my CRISPR fitness screen using BAGEL2"** -> Compute per-sgRNA fold changes from counts, derive per-gene log-likelihood ratios against reference essential and non-essential gene sets, sum to Bayes Factor, and apply BF threshold calibrated by precision-recall against the reference.

- CLI: `BAGEL.py fc` to compute fold changes
- CLI: `BAGEL.py bf` to compute Bayes Factors
- CLI: `BAGEL.py pr` for precision-recall curves
- Reference sets: CEGv2 (essentials) and NEGv1 (non-essentials); both at https://github.com/hart-lab/bagel

## The BAGEL2 Bayesian Framework (under the hood)

**Why this matters for postdoc-level use:** BAGEL2 uses a Bayes-factor classifier trained on known essential and non-essential genes. The chain:

1. For each sgRNA, compute log-fold-change (LFC) treatment vs control.
2. For each gene, look up per-sgRNA LFCs.
3. For each sgRNA, compute the log-likelihood ratio: `log( P(LFC | gene is essential) / P(LFC | gene is non-essential) )`. The numerator and denominator are KDEs (kernel density estimates) of LFC distributions from CEGv2 and NEGv1 reference sgRNAs.
4. Sum per-gene log-likelihood ratios across all sgRNAs targeting the gene -> per-gene Bayes Factor.
5. Bootstrap (default 1000 iterations) for confidence interval; BF >6 corresponds to ~FDR 0.05 (Hart 2017 G3 calibration).

**Critical BAGEL2 improvements over BAGEL1:**

- **Linear extrapolation**: BAGEL1 truncated the LLR at the edges of its KDE; BAGEL2 fits a linear regression in the stable region and extrapolates, giving wider dynamic range. This recovers tumor suppressors (highly positive LFC) that BAGEL1 missed.
- **Multi-target correction**: For sgRNAs targeting multiple genomic loci (off-targets), BAGEL2 down-weights their contribution. The original BAGEL counted off-target hits as essentiality signal.
- **Tumor suppressor sensitivity**: BAGEL2 correctly identifies positive selection (enrichment) genes -- not possible in BAGEL1.

## Calibration to CEGv2 / NEGv1

**Why these reference sets matter:** BAGEL2's discriminative power depends on KDEs of LFCs from known essential vs known non-essential genes. CEGv2 (Hart 2017) is 684 core essential genes shared across cell lines; NEGv1 (Hart 2014) is 927 non-essential genes verified across multiple screens. These act as positive and negative controls within every screen.

Reference set integrity:
- CEGv2: pan-cancer essentials -- common dropouts across most cancer cell lines
- NEGv1: confidently non-essential -- genes without expression or genes with verified neutral status

**Critical pitfall:** Using a custom essentiality reference (e.g., a single-cell-line CRISPR screen) instead of CEGv2 biases the BAGEL2 model toward that line's specific biology. Always use the standardized references unless there is a specific reason for custom training.

## Compute Per-Sample Fold Changes

**Goal:** Generate per-sgRNA fold-change matrix as input for Bayes-factor calculation.

**Approach:** Take normalized counts, compute log-fold-change vs a control (Day 0 or plasmid baseline) per sgRNA.

```bash
# BAGEL2 installation: distributed via git clone (no canonical PyPI release).
git clone https://github.com/hart-lab/bagel
cd bagel
# Some forks publish to PyPI (e.g. `bagel-cas9`) but the official distribution is the GitHub repo.

# Inputs:
# counts.txt: tab-separated with columns: sgRNA, GENE, Sample1, Sample2, ...
# Control column(s): typically Day 0 or plasmid sample(s)
# Treatment column(s): screen endpoint

BAGEL.py fc \
    -i counts.txt \
    -o foldchange.txt \
    -c Plasmid \                           # control sample (or Day 0)
    --min-reads 30                         # minimum reads/sgRNA in control
# Output: foldchange.txt - per-sgRNA LFCs
```

## Compute Bayes Factors

**Goal:** Score per-gene essentiality as a Bayes Factor.

**Approach:** Run `BAGEL.py bf` with the fold-change matrix and reference gene sets; specify number of bootstrap iterations.

```bash
BAGEL.py bf \
    -i foldchange.txt \
    -o bayes_factor.txt \
    -e CEGv2.txt \                         # essentials reference (CEGv2)
    -n NEGv1.txt \                          # non-essentials reference
    -c Sample1,Sample2,Sample3 \            # treatment samples to score
    -k 1000                                # bootstrap iterations (1000 default)
# Output: bayes_factor.txt - per-gene Bayes Factor + CI
```

**Output columns:**

| Column | Meaning |
|--------|---------|
| `GENE` | Gene symbol |
| `BF` | Per-gene Bayes Factor (log-likelihood ratio summed across sgRNAs) |
| `STD` | Standard deviation from bootstrap |
| `NumObs` | Number of sgRNAs contributing |

**Interpretation rule:** BF >6 corresponds to FDR ~0.05 against CEGv2; BF >12 corresponds to FDR ~0.005. Higher BF = stronger evidence the gene is essential. BAGEL2 also reports negative BFs which can indicate tumor suppressors (positive selection).

## Precision-Recall Curve

**Goal:** Empirically select BF threshold for a given precision/recall tradeoff.

**Approach:** Run `BAGEL.py pr` to compute precision and recall at every BF level against CEGv2; pick the BF that gives desired precision.

```bash
BAGEL.py pr \
    -i bayes_factor.txt \
    -o precision_recall.txt \
    -e CEGv2.txt \
    -n NEGv1.txt
# Output: precision_recall.txt - precision/recall at each BF threshold
```

**Practical thresholds (Kim & Hart 2021):**

| BF threshold | Precision | Recall | Use case |
|--------------|-----------|--------|----------|
| 0 | 0.85 | 0.95 | Exploratory; high recall |
| 6 | 0.95 | 0.85 | Standard; corresponds to FDR 0.05 |
| 12 | 0.99 | 0.65 | High-confidence; corresponds to FDR 0.005 |
| 30 | 1.00 | 0.20 | Ultra-stringent; near-certain essentials |

**Pick threshold based on application:** For exploratory hit calling, BF >0 with low precision is acceptable; for clinical-grade essentiality calls, BF >12 or higher.

## Interpret BAGEL2 Results

**Goal:** Stratify genes into essential, non-essential, and tumor-suppressor categories.

**Approach:** Apply BF threshold to classify; flag negative BF as candidate tumor suppressors.

```python
import pandas as pd

def interpret_bagel(bf_path, bf_essential=6, bf_tumor_suppressor=-6):
    '''Classify genes from BAGEL2 BF output.'''
    df = pd.read_csv(bf_path, sep='\t')
    df['call'] = 'neutral'
    df.loc[df['BF'] > bf_essential, 'call'] = 'essential'
    df.loc[df['BF'] < bf_tumor_suppressor, 'call'] = 'tumor_suppressor'
    return df.sort_values('BF', ascending=False)
```

**Tumor suppressor identification:** Genes with significantly negative BF (e.g., <-6) are enriched in the screen, indicating fitness advantage from their loss. This is biologically distinct from "non-essential" and may indicate tumor-suppressor function. BAGEL1 could not detect this; BAGEL2's linear extrapolation enables it.

## Bayesian Reasoning Per Sgrna

**Why this matters:** BAGEL2 computes per-sgRNA contributions; a gene with 4 sgRNAs each contributing +5 to BF gets +20 total. A gene with 3 sgRNAs contributing +5 and 1 sgRNA contributing -3 (off-target or low-efficacy) gets +12 net.

```python
# Per-sgRNA contributions for diagnosis
# Output table: each sgRNA's LLR contribution to gene-level BF
# Useful for identifying low-efficacy guides
```

**Critical:** When per-sgRNA contributions are very heterogeneous (one sgRNA dominates BF), the gene is "guide-of-one"; verify with JACKS efficiency analysis or apply the second-best-sgRNA rule from [[hit-calling]].

## Comparing BAGEL2, MAGeCK, drugZ

| Property | BAGEL2 | MAGeCK | drugZ |
|----------|--------|--------|-------|
| Statistical framework | Bayes factor with reference sets | NB GLM | Bidirectional Z-score |
| Calibrated against | CEGv2 / NEGv1 | Internal null | Vehicle distribution |
| Tumor suppressor detection | YES | Limited (RRA positive-selection score) | YES |
| Best for | Essentiality classification | General hit calling | Chemogenomic drug screens |
| Output | Bayes factor + CI | FDR + LFC | Z-score + FDR per direction |
| Hit threshold | BF >6 (≈FDR 0.05) | FDR <0.05 | FDR <0.05 |
| Library calibration | Indirect (reference set) | None | None |

**Reconciliation:** BF >6 ≈ MAGeCK FDR 0.05 (Hart 2017 G3 calibration). BAGEL2 hits absent from MAGeCK suggest weak signal that BAGEL2's reference anchoring detects but MAGeCK's null-based test misses; verify by inspecting per-sgRNA contributions.

## Failure Modes

### BAGEL2 returns no hits despite known essentials

**Trigger:** Wrong reference gene set file; CEGv2 or NEGv1 file may have wrong format or be missing genes.
**Mechanism:** BAGEL2 trains KDEs on the reference; if references are not representative, KDE separation is poor and no gene has BF >6.
**Symptom:** Median BF near zero; no genes >6 even at low FDR.
**Fix:** Re-download CEGv2 / NEGv1 from https://github.com/hart-lab/bagel. Verify gene symbols match the screen's annotation.

### BAGEL2 calls negative-LFC genes "tumor suppressors"

**Trigger:** Heavy dropout screen where many genes drop out; the dropout signal is captured as positive BF but the *enriched* genes (negative BF) are noise.
**Mechanism:** BAGEL2's symmetric distribution treats deeply enriched genes as significant; in a dropout-only screen, the enrichment signal is purely noise.
**Symptom:** Many genes with negative BF; these don't validate as tumor suppressors.
**Fix:** Restrict tumor-suppressor calling to screens specifically expecting enrichment (e.g., drug-resistance, GoF screens); for dropout screens, only interpret positive BF.

### Bootstrap CI is wide; BF estimates unstable

**Trigger:** Per-gene number of sgRNAs too low (e.g., <4 in some libraries).
**Mechanism:** Bootstrap of LLR over very few sgRNAs creates wide CI.
**Symptom:** STD column larger than BF; many genes have CI spanning zero.
**Fix:** Use a library with at least 4-6 sgRNAs/gene; or increase bootstrap iterations to 5000+; or filter out genes with <3 sgRNAs.

### Low BF for known essential despite high LFC

**Trigger:** One sgRNA per gene is contributing very low LLR (off-target or low-efficacy).
**Mechanism:** BAGEL2 sums LLR; one weak guide drags total down.
**Symptom:** Known essential like RPS3 has BF <6 despite 3 of 4 guides showing -5 LFC.
**Fix:** Inspect per-sgRNA LLR; identify the dragging guide; verify whether to exclude or to use JACKS for efficacy-aware analysis.

### Non-cancer cell-line screen with custom essentials

**Trigger:** Iurine embryonic kidney HEK293T or iPSC-derived neurons where standard essentials may not be essential.
**Mechanism:** CEGv2 is calibrated for cancer cell lines; some essentials in tumor cells are not essential in iPSC.
**Symptom:** PR curve against CEGv2 shows poor separation; many CEGv2 essentials don't drop out.
**Fix:** Use cell-type-specific essentialome (e.g., Dempster 2019 *Nat Commun* defined essentials in various cell types); or use MAGeCK / Chronos which doesn't depend on reference sets.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| BF for FDR ~0.05 | >6 | Hart 2017 G3 calibration; Kim & Hart 2021 |
| BF for FDR ~0.005 | >12 | Empirical from PR curve |
| BF for FDR ~0.001 | >30 | Empirical |
| BF for tumor-suppressor candidate | <-6 | Empirical; verify with orthogonal screen |
| Bootstrap iterations | 1000 default; 5000+ for tight CI | Hart-lab convention |
| Min reads per sgRNA in control | 30 | Joung 2017; BAGEL2 default |
| Min sgRNAs per gene for stable BF | 4-6 | Wider with library convention |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| No hits despite essentials present | Wrong reference set | Re-verify CEGv2 / NEGv1 files |
| Wide bootstrap CI | Too few sgRNAs/gene | Increase library coverage; more iterations |
| Negative BF for known essentials | Confounding factor (e.g., CN amplification) | Pre-correct with CRISPRcleanR / Chronos |
| Tumor suppressor calls don't validate | Pure dropout screen; enrichment is noise | Restrict tumor suppressor calls to expected design |
| Per-sgRNA LLR dominated by one guide | Outlier or off-target | Apply second-best-sgRNA rule |

## References

- Kim E & Hart T. 2021. *Genome Medicine* 13:2. BAGEL2 algorithm and improvements.
- Hart T & Moffat J. 2016. *BMC Bioinformatics* 17:164. BAGEL Bayes factor framework.
- Hart T et al. 2017. *G3* 7:2719. CEGv2 / NEGv1 calibration; FDR-BF relationship.
- Hart T et al. 2014. *Mol Syst Biol* 10:733. Original essential gene reference set.
- Dempster JM et al. 2019. *Nat Commun* 10:5817. Cell-type-specific essentialomes.
- Pacini C et al. 2021. *Cell Syst* 12:1132. Reference essentiality benchmarks.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK RRA/MLE alternative
- crispr-screens/jacks-analysis - JACKS for per-guide efficacy
- crispr-screens/drugz-chemogenomic - drugZ for drug screens
- crispr-screens/hit-calling - Cross-method decision tree
- crispr-screens/screen-qc - Pre-BAGEL QC including CEGv2 PR-AUC
- crispr-screens/library-design - 4-6 sgRNAs/gene library standard
- crispr-screens/copy-number-correction - Pre-correction for cancer-line screens
- pathway-analysis/go-enrichment - Downstream functional analysis
