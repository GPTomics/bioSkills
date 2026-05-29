---
name: bio-gene-regulatory-networks-differential-networks
description: Compare gene co-expression and regulatory networks between biological conditions to find rewired relationships using DiffCorr, DiffCoEx, DINGO/iDINGO, and CoDiNA. Covers the differential-connectivity-is-not-differential-expression distinction, the pairwise multiple-testing explosion, marginal vs partial (direct) rewiring, and the underpowered-rewiring failure mode. Use when comparing co-expression networks between disease vs control, treatment, or developmental stages, or finding hub genes that rewire without changing mean expression. For single-condition modules see coexpression-networks; for differential expression of means see differential-expression/de-results.
tool_type: r
primary_tool: DiffCorr
---

## Version Compatibility

Reference examples tested with: DiffCorr 0.4.1+, DINGO/iDINGO 1.0.4+, CoDiNA 1.1.2+; Python path uses scipy 1.12+, statsmodels 0.14+, networkx 3.0+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

In statsmodels, `multipletests()` defaults to method `'hs'` (Holm-Sidak), NOT Benjamini-Hochberg. Always pass `method='fdr_bh'` explicitly for differential-correlation FDR.

# Differential Networks

**"Compare gene co-expression networks between my disease and control groups"** -> Test whether gene-gene relationships differ between two conditions, identifying gained, lost, and reversed edges and the genes that rewire.
- R: `DiffCorr::comp.2.cc.fdr()` (pairwise Fisher z); `DiffCoEx` (module-level); `iDINGO::dingo()` (partial-correlation)
- Python: Fisher z-test with `scipy.stats` + `statsmodels` FDR

## The Single Most Important Modern Insight -- Differential Connectivity Is Not Differential Expression

A gene can have identical mean expression in two conditions yet a completely rewired set of correlation partners -- and that rewiring, not the mean shift, can be the disease signal. The classic demonstration is Hudson, Reverter & Dalrymple 2009 (*PLoS Comput Biol* 5:e1000382): myostatin received the top Regulatory Impact Factor despite **not being differentially expressed**, correctly fingering the gene carrying the causal mutation purely from the change in its correlation wiring to differentially-expressed targets. So differential expression (a shift in means) and differential connectivity (a shift in the correlation structure) are **orthogonal questions**, and the most differentially-connected hub is often not differentially expressed. Three distinct analyses are routinely conflated and must be kept separate: differential **expression** (mean shift), differential **co-expression** (pairwise correlation shift, DiffCorr/DiffCoEx), and differential **connectivity/rewiring** at the conditional-independence level (DINGO).

The dominant practical failure is **statistical power**. The variance of a *difference* of two correlations is large, so rewiring detection needs many samples per group -- far more than differential expression. Worse, pairwise differential-correlation testing has a multiple-testing explosion: p genes produce ~p^2/2 edge tests, so without aggressive FDR (or a module-level method that sidesteps per-edge testing) the results are dominated by false positives. Most "rewired hub" findings in small cohorts are underpowered noise.

## Differential-Network Method Taxonomy

| Method | Citation | Level | Edge type | Note |
|--------|----------|-------|-----------|------|
| DiffCorr | Fukushima 2013 *Gene* | per-edge | marginal | simple Fisher z; p^2/2 tests -> aggressive FDR needed |
| DiffCoEx | Tesson 2010 *BMC Bioinformatics* | module | marginal | WGCNA-based; tests modules, sidesteps per-edge multiplicity |
| DINGO | Ha 2015 *Bioinformatics* | per-edge | **partial** | group-specific GGM; bootstrap differential score (direct rewiring) |
| iDINGO | Class 2018 *Bioinformatics* | per-edge | partial, multi-omics | chain-graph across data types (e.g. miRNA->mRNA->protein) |
| CoDiNA | Gysi 2020 *PLoS ONE* | per-edge | on supplied nets | compares >=2 networks; common/specific/different edge classes |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Two conditions, quick pairwise rewiring | DiffCorr (Fisher z) + strict FDR | simplest; report gained/lost/reversed |
| Want modules that rewire, not edges | DiffCoEx | module-level testing avoids the p^2/2 explosion |
| Need direct (not indirect) rewiring | DINGO/iDINGO | partial correlation removes confounded indirect changes |
| More than two conditions | CoDiNA | n-way comparison with edge classification |
| Multi-omics rewiring | iDINGO | chain-graph respects the biological hierarchy |
| Just want mean-expression changes | -> differential-expression/de-results | that is DE, not rewiring |
| Build the per-condition networks first | -> coexpression-networks | rewiring compares already-built networks |

## DiffCorr: Pairwise Differential Correlation (R)

**Goal:** Find gene pairs whose correlation differs significantly between two conditions.

**Approach:** Fisher z-transform each correlation per condition and test the z-difference with FDR; classify surviving edges as gained, lost, or reversed.

```r
library(DiffCorr)

expr_all <- read.csv('normalized_counts.csv', row.names = 1)        # genes x samples
info <- read.csv('sample_info.csv', row.names = 1)
# Filter to top variable genes first: p^2/2 edge tests make the full matrix intractable.
gene_vars <- apply(expr_all, 1, var)
top <- names(sort(gene_vars, decreasing = TRUE))[1:3000]
d1 <- expr_all[top, info$condition == 'control']
d2 <- expr_all[top, info$condition == 'disease']

# Returns the differential correlations directly (threshold filters exported pairs by lfdr).
# It only writes the file when save = TRUE, so use the returned data.frame. Columns carry
# spaces ('molecule X', 'molecule Y', 'r1', 'r2', 'lfdr (difference)') -- index with [[ ]].
res <- comp.2.cc.fdr(data1 = d1, data2 = d2, threshold = 0.05, save = TRUE,
                     output.file = 'diffcorr.txt')
```

## DINGO: Direct Differential Rewiring (R)

**Goal:** Detect rewiring at the conditional-independence (direct edge) level rather than marginal correlation.

**Approach:** Estimate a group-specific Gaussian graphical model and bootstrap an edge-wise differential score.

```r
library(iDINGO)
# dingo(dat, x, ...): dat = samples x genes; x = the binary group covariate (length n).
fit <- dingo(dat = expr_mat, x = group, B = 100, cores = 8)   # B = bootstrap reps
# fit$diff.score / fit$p.val give edge-wise differential connectivity (direct edges).
```

## Python: Fisher z Differential Network

**Goal:** Compare correlation networks between two conditions in a Python-native workflow.

**Approach:** Compute per-condition correlation matrices, test each pair with Fisher's z, apply BH FDR (explicitly), and classify edges.

```python
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

def fisher_z(r1, n1, r2, n2):
    z1, z2 = np.arctanh(np.clip([r1, r2], -0.9999, 0.9999))
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z = (z1 - z2) / se
    return z, 2 * stats.norm.sf(abs(z))

def differential_network(e1, e2, fdr=0.05):
    genes = e1.columns.tolist()
    n1, n2 = len(e1), len(e2)
    c1, c2 = e1.corr().values, e2.corr().values
    rows = []
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            z, p = fisher_z(c1[i, j], n1, c2[i, j], n2)
            rows.append((genes[i], genes[j], c1[i, j], c2[i, j], z, p))
    df = pd.DataFrame(rows, columns=['g1', 'g2', 'r1', 'r2', 'z', 'p'])
    # statsmodels default is Holm-Sidak ('hs'); BH must be requested explicitly.
    df['padj'] = multipletests(df['p'], method='fdr_bh')[1]
    return df
```

## Per-Method Failure Modes

### Underpowered rewiring claims
**Trigger:** declaring rewired hubs from a small cohort. **Mechanism:** the variance of a correlation difference is large; rewiring needs more samples than DE. **Symptom:** few or no edges survive FDR, or unstable results across resampling. **Fix:** require adequate n per group; treat low-power results as exploratory.

### Pairwise multiple-testing explosion
**Trigger:** testing all gene pairs with weak/no FDR. **Mechanism:** p genes -> ~p^2/2 tests. **Symptom:** thousands of "significant" edges, irreproducible. **Fix:** pre-filter to variable genes, apply strict FDR, or use a module-level method (DiffCoEx).

### Conflating DE with rewiring
**Trigger:** interpreting differentially-connected genes as differentially expressed (or vice versa). **Mechanism:** they are orthogonal. **Symptom:** a rewired hub dismissed because it is not DE. **Fix:** report DE and differential connectivity separately; a non-DE gene can be the key rewired hub.

### Marginal rewiring read as direct
**Trigger:** interpreting a DiffCorr gained edge as a direct regulatory change. **Mechanism:** marginal correlation mixes direct and indirect edges; a changed edge may reflect a shifted common driver. **Symptom:** mechanistic claims from marginal rewiring. **Fix:** use DINGO (partial correlation) when directness matters.

### Holm-Sidak instead of BH
**Trigger:** `multipletests(p)` without `method=`. **Mechanism:** statsmodels defaults to `'hs'`, more conservative than intended. **Symptom:** unexpectedly few hits. **Fix:** pass `method='fdr_bh'`.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| >= 15-20 samples per group | correlation-stability convention | rewiring is lower-powered than DE; small n gives noise |
| Pre-filter to top ~2000-5000 variable genes | practical | bounds the p^2/2 test count |
| BH FDR < 0.05 | standard | controls the false-discovery rate across many edge tests |
| effect-size filter abs(delta r) > 0.3 | convention | avoid reporting trivially different correlations |
| DINGO bootstrap B = 100 | iDINGO default-scale | stabilizes the differential score |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| millions of edge tests / out of memory | full gene matrix | pre-filter to variable genes |
| far fewer hits than expected | statsmodels Holm-Sidak default | use `method='fdr_bh'` |
| rewired hub "should be DE" objection | conflating connectivity with expression | report them as separate, orthogonal results |
| DGCA not installable from CRAN | archived May 2024 | install from GitHub (`andymckenzie/DGCA`) |
| reversed edges look like noise | no effect-size filter | require abs(delta r) above a threshold |

## References

- Hudson NJ, Reverter A, Dalrymple BP. 2009. A differential wiring analysis... correctly identifies the gene containing the causal mutation. *PLoS Comput Biol* 5(5):e1000382.
- de la Fuente A. 2010. From 'differential expression' to 'differential networking'. *Trends Genet* 26(7):326-333.
- Fukushima A. 2013. DiffCorr: analyze and visualize differential correlations in biological networks. *Gene* 518(1):209-214.
- Tesson BM, Breitling R, Jansen RC. 2010. DiffCoEx: differentially coexpressed gene modules. *BMC Bioinformatics* 11:497.
- Ha MJ, Baladandayuthapani V, Do KA. 2015. DINGO: differential network analysis in genomics. *Bioinformatics* 31(21):3413-3420.
- Class CA, Ha MJ, Baladandayuthapani V, Do KA. 2018. iDINGO: integrative differential network analysis in genomics. *Bioinformatics* 34(7):1243-1245.
- Gysi DM, et al. 2020. Co-expression differential network analysis (CoDiNA). *PLoS ONE* 15(10):e0240523.

## Related Skills

- coexpression-networks - build the per-condition co-expression networks being compared
- scenic-regulons - TF regulon activity differences as a complementary rewiring readout
- grn-inference - VIPER differential protein activity between conditions
- differential-expression/de-results - differential expression of means (the orthogonal question)
- temporal-genomics/temporal-grn - time-resolved network change across stages
