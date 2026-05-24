---
name: bio-crispr-screens-batch-correction
description: Batch effect correction for CRISPR screens covering ComBat empirical-Bayes, RUV, SVA, control-sgRNA normalization, and the model-based alternative of including batch as a covariate in MAGeCK MLE or Chronos. Covers screen-specific batch sources (passage cohort, library lot, infection day, sequencing run, Cas9 lot, FBS lot), PCA + variance-decomposition diagnostic to decide if correction is needed, when correction harms biology by over-correcting condition into batch, limma removeBatchEffect for visualization-only correction, and relationship to multi-condition design matrices. Use when combining screens for joint analysis, when passage cohort confounds biology, when DepMap-style panels need Chronos with batch covariates, when picking ComBat vs RUV, or when correction harms biology and should be replaced with explicit covariate modeling.
tool_type: mixed
primary_tool: pyComBat
---

## Version Compatibility

Reference examples tested with: pyComBat 0.3.3+ (epigenelabs/pyComBat), MAGeCK 0.5.9+, R/limma 3.58+, sva 3.50+, RUVSeq 1.36+, pandas 2.2+, numpy 1.26+, scikit-learn 1.4+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show pycombat`; `from combat.pycombat import pycombat`
- R: `packageVersion('sva')`; `?ComBat`; `packageVersion('RUVSeq')`; `?RUVg`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Batch Correction for CRISPR Screens

**"Correct batch effects in my CRISPR screens"** -> Diagnose the batch source, decide whether to remove via empirical-Bayes (ComBat), explicit covariate modeling (MAGeCK MLE / Chronos design matrix), control-guide-anchored normalization, or unwanted-variation decomposition (RUV, SVA), then apply only the correction that preserves biological condition signal.

- Python: `pyComBat.pycombat` for empirical-Bayes correction
- Python: explicit batch covariates in `mageck mle --design-matrix`
- R: `sva::ComBat`, `RUVSeq::RUVg`, `limma::removeBatchEffect`
- CLI: `chronos-cn` natively handles screen-batch covariates

## Batch Sources in CRISPR Screens

| Source | Mechanism | Detectable by |
|--------|-----------|---------------|
| Library lot | Different aliquots or PCR amplifications | Gini shift; plasmid-pool sequencing |
| Cell passage cohort | Cells passaged through different periods | PCA Day-0 samples clustering by passage |
| Infection day | Lentivirus titer drifts; FBS lot changes | PCA Day-0 samples cluster by day |
| Cas9 enzyme lot | Cas9 expression heterogeneity | PR-AUC drift across screens |
| Sequencing run | Lane bias, flowcell variant, machine | Per-sample read-count distribution |
| FBS / culture lot | Fetal bovine serum lot variations confound proliferation | Day-0 vs endpoint differential not present in vehicle |
| Tissue-prep batch | In-vivo: animal cohort, surgical day, organ-prep tech | In-vivo screens (see [[in-vivo-screens]]) |

**Critical:** Batch effects in CRISPR screens often correlate with biology (e.g., the drug arm was processed in batch 2 because that's when the drug arrived). This confounds correction. Always check for confounding before applying ComBat.

## Batch Effect Decision Tree

| Diagnostic finding | Recommended correction |
|--------------------|------------------------|
| PCA shows samples cluster by condition, not batch | No correction needed; biology dominates |
| PCA PC1 separates batches, PC2 separates conditions | Apply ComBat with condition as biological_covariate |
| Batch fully confounded with condition (e.g. all drug in batch 2, all vehicle in batch 1) | Correction will destroy biology; instead redesign next screen with cross-batch balance OR re-analyze with batch in MAGeCK MLE design matrix |
| Day-0 (pre-perturbation) samples cluster by batch | Strong batch effect; ComBat needed |
| Endpoint samples cluster by batch but not Day-0 | Selection-driven artifact (FBS lot etc); correct or include batch as covariate |
| Replicates within a batch are tight; across-batch much wider | Classic batch effect; ComBat |
| Each replicate scatters randomly across PCs | Sample-level noise; no batch correction will help |
| Cancer-line panel with multiple batches | Use Chronos (built-in batch modeling) |

## Diagnose: PCA + Variance Decomposition

**Goal:** Quantify what fraction of variance is batch vs condition before correcting.

**Approach:** Run PCA on log10(counts+1); fit ANOVA decomposing variance into batch and condition components; report variance explained.

```python
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy import stats

def batch_diagnostic(counts_df, metadata_df, batch_col='batch', condition_col='condition'):
    '''Variance decomposition: report fraction of PC1/PC2 variance attributable to batch vs condition.'''
    log_counts = np.log10(counts_df + 1).T  # samples as rows
    pca = PCA(n_components=5)
    pcs = pca.fit_transform(log_counts)
    out = pd.DataFrame({
        'PC': range(1, 6),
        'var_explained': pca.explained_variance_ratio_,
    })
    pc_df = pd.DataFrame(pcs, columns=[f'PC{i+1}' for i in range(5)], index=counts_df.columns).join(metadata_df)
    for i in range(5):
        pc = pc_df[f'PC{i+1}']
        f_b, p_b = stats.f_oneway(*[pc[pc_df[batch_col] == b] for b in pc_df[batch_col].unique()])
        f_c, p_c = stats.f_oneway(*[pc[pc_df[condition_col] == c] for c in pc_df[condition_col].unique()])
        out.loc[i, 'batch_F'] = f_b
        out.loc[i, 'batch_p'] = p_b
        out.loc[i, 'cond_F'] = f_c
        out.loc[i, 'cond_p'] = p_c
    return out
```

**Interpretation:** If PC1 has batch F-stat > condition F-stat by 10x, batch is dominating and correction is warranted. If condition dominates PC1, no correction needed.

## ComBat Empirical-Bayes Correction

**Goal:** Remove batch-specific location and scale shifts while preserving biological condition signal.

**Approach:** Log-transform counts, fit ComBat with explicit `biological_covariate` indicating condition (so the model knows which signal to preserve), back-transform.

```python
import numpy as np
from combat.pycombat import pycombat

def combat_correct(counts_df, batch_vector, condition_vector=None):
    '''ComBat on log-counts with optional biological covariate (condition).
    Preserves condition signal while removing batch shifts.'''
    data = np.log2(counts_df.values + 1)
    if condition_vector is not None:
        mod = pd.get_dummies(condition_vector).values.astype(float)
        corrected = pycombat(data, np.array(batch_vector), mod=mod)
    else:
        corrected = pycombat(data, np.array(batch_vector))
    return pd.DataFrame(np.power(2, corrected) - 1,
                         index=counts_df.index, columns=counts_df.columns).clip(lower=0)
```

**Critical caveat:** ComBat assumes batch effects are linear shifts of mean and variance in log space. Non-linear effects (e.g., gene-specific batch sensitivity) remain. Always re-check PCA after correction to confirm batches now overlap.

## RUV (Remove Unwanted Variation)

**Goal:** Identify hidden batch sources via control sgRNAs whose true signal is known.

**Approach:** Designate non-targeting controls as "negative controls" (assumed unchanged); RUV decomposes their variance into unwanted factors, then subtracts these from all data.

```r
library(RUVSeq)
# counts_df: rows = sgRNAs, columns = samples
ntc_indices <- which(rownames(counts_df) %in% ntc_sgrna_names)
seqset <- newSeqExpressionSet(counts = as.matrix(counts_df))
ruv_corrected <- RUVg(seqset, cIdx = ntc_indices, k = 2)  # k = 2 unwanted factors
# Access corrected data
corrected_counts <- normCounts(ruv_corrected)
```

**When to use:** RUV preferred over ComBat when batches are not annotated (e.g., unknown technical confounders). Worse than ComBat when batch is known and well-annotated; ComBat is more direct.

## SVA (Surrogate Variable Analysis)

**Goal:** Estimate unknown latent factors that may confound the screen.

**Approach:** SVA computes surrogate variables that capture variance not explained by known biological factors; these can then be added to the MAGeCK MLE design matrix as covariates.

```r
library(sva)
# counts_df: rows = sgRNAs, columns = samples
mod <- model.matrix(~ condition, data = metadata)
mod0 <- model.matrix(~ 1, data = metadata)
sv_obj <- sva(as.matrix(counts_df), mod, mod0)
n_sv <- sv_obj$n.sv  # number of surrogate variables
# Add to design matrix for MAGeCK MLE
design_mat <- cbind(mod, sv_obj$sv)
```

**Use case:** When the screen has clear biological signal (e.g. essentiality recovery passes) but small effect sizes are hidden by noise; SVA-discovered latent factors as covariates can recover them.

## Batch as Explicit Covariate (Preferred for MAGeCK MLE / Chronos)

**Goal:** Model batch and biology in the same regression instead of pre-correcting.

**Approach:** Add batch indicator columns to the MLE design matrix. The fitted beta for condition is the effect after accounting for batch; no pre-correction needed.

```bash
# Design matrix for a screen with 2 batches and 2 conditions
cat > design.txt <<EOF
Samples         baseline    batch2    treatment
Veh_b1_r1       1           0         0
Veh_b1_r2       1           0         0
Drug_b1_r1      1           0         1
Drug_b1_r2      1           0         1
Veh_b2_r1       1           1         0
Veh_b2_r2       1           1         0
Drug_b2_r1      1           1         1
Drug_b2_r2      1           1         1
EOF

mageck mle \
    --count-table counts.txt \
    --design-matrix design.txt \
    --output-prefix batch_aware_mle
```

**Why this is preferred:** ComBat shifts counts before testing; the MLE-with-covariates approach correctly propagates uncertainty from the batch term into the condition beta's standard error. ComBat-then-test pretends the corrected counts are noise-free, biasing FDR.

## Control-Sgrna Anchored Normalization

**Goal:** Use non-targeting controls as the per-sample reference so batch shifts cancel.

**Approach:** Scale each sample so its NTC sgRNAs have a constant median. Subsequent fold changes are relative to NTCs in each sample, automatically batch-controlling.

```python
def ntc_anchored_normalize(counts_df, ntc_sgrna_names, target_median=1000):
    '''Scale each sample so its NTC median is target_median. Subsequent LFC is NTC-anchored.'''
    is_ntc = counts_df.index.isin(ntc_sgrna_names)
    ntc_medians = counts_df.loc[is_ntc].median(axis=0)
    scale_factors = target_median / ntc_medians.replace(0, np.nan)
    return counts_df * scale_factors, scale_factors
```

**Critical:** Requires ≥500 NTCs in the library (see [[library-design]]). With fewer, the NTC median is unstable and amplifies noise rather than removing batch.

## When NOT to Correct

| Situation | Why correction hurts |
|-----------|----------------------|
| Batch is fully confounded with condition | Correction destroys biology along with batch; redesign or accept |
| Batch effect is smaller than between-replicate noise | Correction adds noise without removing meaningful variance |
| Replicates already correlate >0.95 within and across batches | No batch effect to correct |
| Single-screen analysis | No "batch" to correct; only replicate noise |
| Per-batch sample size <3 | Cannot estimate batch shift reliably; correction is harmful |

## Failure Modes

### ComBat eliminates biological signal

**Trigger:** Batch is correlated with condition (e.g., all drug-arm samples were processed week 2; all vehicle-arm samples week 1).
**Mechanism:** ComBat without a `mod` covariate treats condition variance as batch variance; corrects it away.
**Symptom:** PR-AUC against CEGv2 drops after ComBat correction.
**Fix:** Always supply `mod` covariate matrix indicating condition; verify by comparing PR-AUC before and after.

### RUV adds noise instead of removing it

**Trigger:** k (number of unwanted factors) set too high.
**Mechanism:** RUV's least-squares decomposition over-fits; "removed" variance includes biology.
**Symptom:** Hits decrease and replicate Pearson drops after correction.
**Fix:** Choose k via cross-validation; default k=1 or 2 for most screens.

### Batch-aware MLE collinear design matrix

**Trigger:** Adding a batch indicator that is fully collinear with another design column (e.g., all of batch 2 is also Day 21).
**Mechanism:** MLE design matrix is singular; betas not estimable.
**Symptom:** MAGeCK MLE errors out or produces NaN betas.
**Fix:** Drop the collinear column; re-design experiment with cross-batch balance.

### ComBat after RUV double-corrects

**Trigger:** Applying multiple corrections sequentially.
**Mechanism:** Both methods remove variance; sequential application removes biology twice.
**Symptom:** All signal gone; counts look uniformly noisy.
**Fix:** Pick one method based on diagnostic; never combine.

### Per-batch sample size too small

**Trigger:** 2 replicates per batch with 3 batches; ComBat estimates batch shift from 2 samples.
**Mechanism:** Insufficient data to estimate batch parameters; high-variance estimates.
**Symptom:** Correction makes some batches worse than uncorrected.
**Fix:** Need ≥3 (preferably 4-6) samples per batch; below this, use covariate modeling instead.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| PC1 batch F vs condition F | F_batch > 10x F_cond -> apply correction | Standard variance-decomposition diagnostic |
| ComBat min samples per batch | ≥3, ideally 4-6 | Empirical Bayes prior estimation |
| RUV `k` (unwanted factors) | k=1 default; k=2 if multiple known batch sources | Risso 2014; cross-validate |
| NTCs needed for NTC-anchored norm | ≥500 in library | Stable median |
| Post-correction PCA check | Batches must overlap in PC1/PC2 plot | Visual sanity check |
| Post-correction PR-AUC | Should be same or higher than pre | If lower, correction destroyed biology |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| PR-AUC drops after ComBat | Batch confounded with condition | Add `mod` covariate; or redesign |
| MAGeCK MLE NaN beta after adding batch column | Collinear design matrix | Drop collinear column |
| Replicates still cluster by batch after RUV | k too low | Increase k; cross-validate |
| Replicates lose internal cohesion after correction | Over-correction | Reduce k or revert |
| NTC-anchored norm worse than median | Too few NTCs | Use median; add NTCs to next library |
| Sequencing-run-level batch survives ComBat | Non-linear sequencing effect | Pre-normalize with `mageck count --norm-method control` first |

## References

- Johnson WE et al. 2007. *Biostatistics* 8:118. Original ComBat algorithm.
- Leek JT et al. 2012. *Bioinformatics* 28:882. SVA package.
- Risso D et al. 2014. *Nat Biotechnol* 32:896. RUVSeq.
- Pacini C et al. 2021. *Cell Syst* 12:1132. CRISPR-screen batch effects analysis.
- Hayer KE et al. 2023. *Genome Biol* 24:1. Benchmark of normalization methods for CRISPR screens.

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK MLE with explicit batch covariates
- crispr-screens/screen-qc - Pre-correction PCA diagnostic
- crispr-screens/copy-number-correction - Chronos handles batch + CN jointly
- crispr-screens/library-design - NTC composition for NTC-anchored normalization
- crispr-screens/jacks-analysis - Joint analysis across batches with shared efficacy
- crispr-screens/hit-calling - Post-correction hit calling
- crispr-screens/in-vivo-screens - In-vivo-specific batch sources (animal cohort, tissue prep)
