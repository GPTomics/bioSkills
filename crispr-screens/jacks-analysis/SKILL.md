---
name: bio-crispr-screens-jacks-analysis
description: Runs JACKS (Joint Analysis of CRISPR/Cas9 Knockout Screens; Allen et al 2019 Genome Research) which models per-sgRNA log-fold-change as the product of a treatment-dependent gene-essentiality term and a treatment-independent guide-efficacy term. Covers the Bayesian decomposition math, the hierarchical efficacy prior shared across screens performed with the same library, when JACKS outperforms MAGeCK (multi-screen joint analysis, libraries with broad efficacy variance) and when it does not (single screen, novel libraries with no prior efficacy), library-reuse efficacy transfer, downstream essentiality interpretation, and the 2.5x sample-size reduction enabled by efficacy-aware testing. Use when running multiple screens with the same library, when guide-level noise is suspected to dominate per-gene signal, when reusing published essentiality reference screens for efficacy priors, or when comparing screens performed across cell lines that share library but differ biologically.
tool_type: python
primary_tool: JACKS
---

## Version Compatibility

Reference examples tested with: JACKS 0.2.0+ (felicityallen/JACKS), pandas 2.2+, numpy 1.26+, scipy 1.12+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `python run_JACKS.py --help` (run_JACKS.py at the JACKS repo root after clone)
- Python: from jacks.jacks_io import runJACKS; help(runJACKS)
- GitHub: install via `git clone https://github.com/felicityallen/JACKS && cd JACKS && pip install .`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## JACKS CRISPR Screen Analysis

**"Analyze CRISPR screens with guide-level efficacy modeling"** -> Jointly model per-sgRNA log-fold-change across one or more screens as the product of gene essentiality and guide efficacy, sharing efficacy across screens with the same library so that low-quality guides are down-weighted automatically.

- CLI: `python run_JACKS.py countfile replicatefile guidemappingfile [options]` (script at JACKS repo root)
- Python: `from jacks.jacks_io import runJACKS` for programmatic use; lower-level `from jacks.infer import inferJACKS`
- Output: per-gene essentiality (`X1`), per-sgRNA efficacy (`X1`), log-likelihood ratio per gene

## The JACKS Model (under the hood)

**Why this matters for postdoc-level use:** JACKS decomposes the observed per-sgRNA log-fold-change as:

```
LFC[i, c] = gene_effect[g(i), c] * guide_efficacy[i] + noise
```

where `i` is sgRNA index, `c` is screen condition, `g(i)` is the gene targeted by sgRNA i. Gene effect varies by condition (different cell lines, different treatments) but guide efficacy is intrinsic to the sgRNA sequence and is treated as constant across screens. The model fits both parameters via variational Bayes with hierarchical priors:

- `guide_efficacy[i] ~ Beta(alpha, beta)` shared across all sgRNAs (hyperparameters fit empirically)
- `gene_effect[g, c] ~ Normal(0, sigma_c^2)` per condition

The variational posterior gives expected guide efficacy and gene effect; log-likelihood-ratio tests against a null (zero gene effect) provide gene-level significance.

**Critical assumption:** Guide efficacy is treated as cell-line independent within the same chemistry. Allen 2019 reports per-sgRNA Cas9 KO efficacy correlates ~0.7 across cell lines (within-chemistry), supporting library-shared efficacy. **However**, efficacy is NOT shareable across chemistries: Cas9 KO efficacy != CRISPRi knockdown efficiency != CRISPRa activation efficiency. JACKS must be run separately per chemistry; use only within the same chemistry on the same library.

## When JACKS Outperforms MAGeCK and BAGEL2

| Scenario | Advantage | Quantified gain (Allen 2019) |
|----------|-----------|------------------------------|
| Multi-screen joint analysis (>=3 screens with same library) | Efficacy shared; noise averaged | 12-21% lower error vs MAGeCK; 9% vs BAGEL2; 97-99% of cell lines improved |
| Reusing public reference screens (DepMap, Sanger Score) as efficacy prior | Transfer learning | New screens can be smaller; 2.5x sample-size reduction |
| Libraries with broad efficacy variance (e.g. older GeCKOv2) | Down-weights known weak guides | Larger gain than on Brunello (already efficacy-filtered) |
| Heterogeneous quality (mixed plasmid quality across screens) | Per-screen noise estimation | Cleaner per-condition gene effects |

## When JACKS Is Not the Right Tool

- **Single screen, no prior efficacy:** JACKS has nothing to leverage; MAGeCK or BAGEL2 work as well.
- **Single timepoint / two-condition essentiality:** RRA or BAGEL2 simpler and equivalent.
- **Heavy-selection drug screens:** drugZ explicit for chemogenomic; JACKS less sensitive.
- **Cancer-cell-line copy-number screens:** Chronos preferred; jointly models CN bias + screen quality; JACKS does neither.
- **Cross-chemistry sharing (e.g. CRISPRi + Cas9):** Efficacy is chemistry-specific; do not share.

## Run JACKS Joint Analysis

**Goal:** Jointly analyze multiple CRISPR screens performed with the same library and chemistry.

**Approach:** Provide a count matrix with all samples across all screens, a replicate map identifying which samples belong to which screen and condition, and a sgRNA-to-gene map. JACKS learns guide efficacy shared across screens and gene effects per screen.

```python
# Programmatic invocation
from jacks.jacks_io import runJACKS

# Input file paths
counts_path = 'counts.txt'                    # rows=sgRNA; first cols 'sgRNA' (or custom), then sample counts
replicate_map_path = 'replicatemap.txt'       # tab-separated with header: Replicate, Sample, Control
guide_map_path = 'guidemap.txt'               # tab-separated with header: sgRNA, Gene

# Replicate map format (tab-separated WITH header; column names match flags below)
# Replicate                Sample          Control
# Screen1_T1               Screen1_T       Screen1_C
# Screen1_T2               Screen1_T       Screen1_C
# Screen1_C1               Screen1_C       Screen1_C
# Screen2_T1               Screen2_T       Screen2_C
# Screen2_T2               Screen2_T       Screen2_C
# Screen2_C1               Screen2_C       Screen2_C

runJACKS(
    countfile=counts_path,
    replicatefile=replicate_map_path,
    guidemappingfile=guide_map_path,
    rep_hdr='Replicate',
    sample_hdr='Sample',
    ctrl_sample_hdr='Control',                # per-sample control specification
    sgrna_hdr='sgRNA',
    gene_hdr='Gene',
    outprefix='jacks_out',
    apply_w_hp=True,                          # hierarchical prior on efficacy
)
```

```bash
# Equivalent CLI run (run_JACKS.py is at the JACKS repo root after clone)
python run_JACKS.py \
    counts.txt \
    replicatemap.txt \
    guidemap.txt \
    --rep_hdr Replicate \
    --sample_hdr Sample \
    --ctrl_sample_hdr Control \              # per-sample control (or --common_ctrl_sample <name>)
    --sgrna_hdr sgRNA \
    --gene_hdr Gene \
    --outprefix jacks_out \
    --apply_w_hp                              # hierarchical prior on efficacy
# Outputs:
#   jacks_out_gene_JACKS_results.txt    gene-level: X1 (effect), X2 (std), p_neg, p_pos
#   jacks_out_grna_JACKS_results.txt    sgRNA-level: X1 (efficacy 0-1), X2 (std)
#   jacks_out_JACKS_results_full.pickle  full posterior for downstream
```

## Output Interpretation

| Column | Meaning | Direction |
|--------|---------|-----------|
| `X1` (gene file) | Posterior mean of gene_effect | Negative = essential (depleted); positive = enriched |
| `X2` (gene file) | Posterior std of gene_effect | Lower = more confident; combine with X1/X2 ratio for z-like scoring |
| `p_neg` / `p_pos` | One-sided posterior probability | Hits at p_neg <0.05 or p_pos <0.05 |
| `fdr_log10` | log10(FDR) of gene effect | Use <-1 for FDR<0.1; <-2 for FDR<0.01 |
| `X1` (sgRNA file) | Posterior mean of guide efficacy | 0 = ineffective, 1 = maximum efficacy |
| `X2` (sgRNA file) | Posterior std of efficacy | Confidence in efficacy estimate |

**Interpretation rule:** A gene is essential if X1 < 0 AND fdr_log10 < -1. The X1/X2 ratio gives a z-like statistic; |X1/X2| > 2 corresponds to ~95% credible deviation from zero. Sort by X1 (most negative first) for essentiality rank.

## Build Library-Wide Efficacy Prior from Reference Screens

**Goal:** Transfer learned efficacy from a large public screen panel to a new small screen.

**Approach:** Run JACKS on the reference panel (e.g. DepMap CRISPR screens with TKOv3 or Brunello), extract per-sgRNA efficacy posterior, and supply it as the prior for a new screen.

```python
def extract_efficacy_prior(reference_jacks_results):
    '''Build per-sgRNA efficacy prior (mean + std) from a large reference screen.'''
    df = pd.read_csv(reference_jacks_results, sep='\t')
    prior = df[['sgRNA', 'X1', 'X2']].rename(columns={'X1': 'eff_mean', 'X2': 'eff_std'})
    return prior

# Use in new JACKS run via --reffile <path>
# Reference: Allen 2019 Genome Research 29:464; 2.5x sample-size reduction documented
```

## Per-sgRNA Efficacy Diagnostics

**Goal:** Identify low-efficacy guides for library refinement.

**Approach:** Examine the distribution of inferred efficacies; guides below 0.3 are likely non-functional and should be excluded from re-designed libraries.

```python
import matplotlib.pyplot as plt

def efficacy_summary(grna_results_path, low_threshold=0.3):
    df = pd.read_csv(grna_results_path, sep='\t')
    df['low_eff'] = df['X1'] < low_threshold
    summary = {
        'total_guides': len(df),
        'low_efficacy_count': df['low_eff'].sum(),
        'low_efficacy_pct': df['low_eff'].mean() * 100,
        'median_efficacy': df['X1'].median(),
        'q25_q75': (df['X1'].quantile(0.25), df['X1'].quantile(0.75)),
    }
    # Per-gene proportion of low-efficacy guides
    by_gene = df.groupby('Gene')['low_eff'].mean().sort_values(ascending=False)
    summary['genes_with_all_low_eff'] = (by_gene == 1).sum()  # genes where every guide is weak
    return summary, by_gene
```

**Critical:** Genes where every guide is low-efficacy will show no signal regardless of biology. Filter from interpretation; flag for re-design with updated rules (Brunello / TKOv3).

## Comparing JACKS, MAGeCK, BAGEL2

| Property | JACKS | MAGeCK | BAGEL2 |
|----------|-------|--------|--------|
| Statistical framework | Variational Bayes | NB GLM + alpha-RRA / MLE | Bayes factor on per-sgRNA fold change |
| Models guide efficacy | Yes (jointly) | No (optional fixed input) | No |
| Multi-screen joint | Yes (native) | Limited (MLE design matrix) | No (per-screen) |
| Speed | Slow (variational inference) | Fast | Fast |
| Output | gene effect + sgRNA efficacy | beta or RRA score | Bayes Factor |
| Best for | Multi-screen joint analyses, library calibration | General-purpose, single screen | Essentiality classification |
| Quantified accuracy gain (Allen 2019) | 12-21% lower error vs MAGeCK | Reference | 9% lower error vs JACKS for essentiality alone |

**Reconciliation:** Hits identified by JACKS AND MAGeCK are high confidence. JACKS-only hits typically reflect strong gene signals where one or two guides were dragging down MAGeCK; verify the up-weighted high-efficacy guides have the expected sign. MAGeCK-only hits at FDR <0.05 may be single-guide outliers; check sgrna_summary for guide-level dispersion.

## Failure Modes

### Efficacy collapsed near zero for all guides

**Trigger:** Screen used a chemistry the model doesn't support (e.g., CRISPRi screen analyzed with JACKS defaults).
**Mechanism:** CRISPRi efficacy is fundamentally different from Cas9-KO efficacy; the Beta-prior hyperparameters fit on Cas9 data don't transfer.
**Symptom:** Median efficacy <0.2; almost no significant gene effects.
**Fix:** Train per-chemistry priors separately; for CRISPRi/a, current JACKS recommends `--apply_w_hp` with manually set hyperparameters from a CRISPRi reference dataset.

### Cross-cell-line efficacy disagreement

**Trigger:** Pooling screens across cell lines with very different Cas9 expression / chromatin / fitness baselines.
**Mechanism:** Efficacy depends on Cas9 expression and chromatin accessibility; sharing across lines averages real per-line differences.
**Symptom:** Per-line gene effects look noisier than per-line MAGeCK results.
**Fix:** Use Chronos for multi-cell-line screens with screen-quality modeling; reserve JACKS for screens with matched chemistry + cell type / culture conditions.

### MCMC / variational convergence failure

**Trigger:** Too few iterations relative to library size (10k iters for 100k-guide library is sometimes insufficient).
**Mechanism:** Variational lower bound has not plateaued; estimates noisy.
**Symptom:** Repeated runs produce different gene effects.
**Fix:** Increase `n_iterations` (default 1000) to 5000-10000 for publication-grade output; verify ELBO plateaus; use `--seed` for reproducibility.

### sgRNA-to-gene map mismatch

**Trigger:** Guide map and count matrix use different sgRNA naming conventions (e.g. `BRCA1_1` vs `BRCA1.1`).
**Mechanism:** JACKS reads the map as a join; mismatched rows give NaN gene effects.
**Symptom:** Many genes missing from output.
**Fix:** Standardize naming; sanity check `len(jacks_output) == n_genes_expected`.

### Reference efficacy prior from wrong library

**Trigger:** Using DepMap Brunello efficacy as prior for a screen with a custom TKOv3-style library.
**Mechanism:** Per-sgRNA efficacy is sequence-specific; sgRNAs in one library map to different gene contexts than another.
**Symptom:** Worse gene-effect estimation than no prior.
**Fix:** Match library exactly; if no matched reference exists, run without prior.

## Reconciliation: When JACKS and Other Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| JACKS significant, MAGeCK not | One low-efficacy guide dragged MAGeCK; JACKS down-weighted it | Trust JACKS if 3+ high-efficacy guides agree |
| MAGeCK significant, JACKS not | All guides have similar efficacy; JACKS prior shrinks signal | Verify per-guide LFC consistency in MAGeCK sgrna_summary |
| JACKS efficacy ~0.5 for all guides | Hierarchical prior over-shrinkage | Run with `--apply_w_hp` false; refit hyperparameters |
| Gene effect different sign from MAGeCK | Multi-screen pooling created mean effect different from single-screen | Run per-screen separately to confirm |

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Hit fdr_log10 | <-1 (FDR <0.1); <-2 (FDR <0.01) | Allen 2019 Genome Research convention |
| Effective gene signal | X1 < 0 AND abs(X1/X2) > 2 | Bayesian z-equivalent |
| Low-efficacy guide flag | X1 (sgRNA) <0.3 | Allen 2019: below this, guide is non-functional |
| Reference for prior reuse | DepMap or Sanger Score panel | Established efficacy distribution |
| Minimum screens for joint efficacy benefit | 3+ | Below this, single-screen tools (MAGeCK/BAGEL2) equivalent |
| Iterations for variational inference | 5000+ publication; 1000 default | Verify ELBO plateaus |
| Cross-library efficacy transfer | Not supported | Different libraries -> different sequences -> different efficacies |
| Cross-chemistry efficacy transfer | Not supported | Cas9 efficacy != CRISPRi efficacy |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Many NaN gene effects | sgRNA-to-gene map mismatch | Verify naming consistency between count matrix and map |
| Median efficacy <0.2 | Wrong chemistry assumed by prior | Disable `--apply_w_hp` or use matched prior |
| ELBO not plateaued | Too few iterations | Increase iterations to 5000+ |
| Inconsistent gene effects between runs | No `--seed` set | Set seed for reproducibility |
| Library-reuse prior doesn't help | Wrong library reference | Match library exactly |

## References

- Allen F et al. 2019. *Genome Research* 29:464. JACKS; original Bayesian joint analysis paper.
- Yusa K, Allen F (Sanger Wellcome Institute). https://github.com/felicityallen/JACKS. Official repository.
- Behan FM et al. 2019. *Nature* 568:511. Sanger Score CRISPR Panel; library-wide reference screen data.
- Tsherniak A et al. 2017. *Cell* 170:564. DepMap CRISPR; alternative reference for efficacy transfer.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK RRA/MLE comparison
- crispr-screens/bagel-essentiality - Alternative for essentiality without efficacy modeling
- crispr-screens/library-design - sgRNA design rules informed by JACKS efficacy output
- crispr-screens/copy-number-correction - Chronos preferred for cancer-line multi-screen analyses
- crispr-screens/screen-qc - Pre-JACKS QC; replicate Pearson must pass before joint analysis
- crispr-screens/hit-calling - Cross-method decision tree
- crispr-screens/batch-correction - JACKS does not adjust for batch; pre-correct if necessary
