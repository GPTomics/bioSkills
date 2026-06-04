---
name: bio-imaging-mass-cytometry-differential-analysis
description: Compare cell-type composition and spatial features across conditions in IMC/MIBI cohorts with the patient as the experimental unit, covering pseudoreplication, per-patient aggregation, mixed models, compositional (Dirichlet/scCODA) differential abundance, diffcyt, per-image-to-patient spatial differential testing (SpaceANOVA), batch covariates, and FDR. Use when testing whether a cell type or spatial niche differs between groups, avoiding cell-level pseudoreplication, choosing a differential-abundance method, or correctly powering an IMC cohort comparison.
tool_type: mixed
primary_tool: diffcyt
---

## Version Compatibility

Reference examples tested with: diffcyt 1.22+ (R), lme4 1.1+ (R), statsmodels 0.14+, scanpy 1.10+, sccoda 0.1.9+, numpy 1.26+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: cell-type proportions are compositional (they sum to 1), so a real increase in one type forces apparent decreases in others -- a Dirichlet/CLR-aware method (scCODA) or a reference cell type is needed, not independent per-type tests. statsmodels `mixedlm` fits patient as a random effect. diffcyt operates on per-sample cluster counts (DA) and per-sample median marker expression (DS).

# IMC Differential Analysis

**"Compare cell types and spatial structure between my conditions"** -> Aggregate to the patient, then test across patients -- never across cells.
- Python: `statsmodels.formula.api.mixedlm`, `sccoda`, `scanpy`
- R: `diffcyt`, `lme4::lmer`, SpaceANOVA for spatial differential testing

## The Single Most Important Modern Insight -- the replicate is the patient, not the cell, and the million-cell count is a red herring

After phenotyping and spatial analysis, every interesting claim ("disease has more Tregs", "responders have more CD8-tumor contact") is a comparison BETWEEN groups, and the single most common fatal error is testing it at the cell level. Hundreds of thousands of cells from one patient are not independent replicates -- they are correlated reads of one biological sample, and cells within an image are massively spatially autocorrelated. Testing at the cell level inflates n by orders of magnitude and manufactures significance: a per-cell test over 50,000 cells reports p~0 for trivial effects because the effective sample size is the number of PATIENTS (often 10-40), not cells (Squair 2021 *Nat Commun* 12:5692). In imaging this is worse than in scRNA because slide and ROI add nesting levels and ROIs are not random samples of the tissue. The correct spine is invariant across every differential question: compute a per-image (or per-ROI) summary, aggregate to ONE value per patient, then test across patients with the patient as the unit -- a mixed model with patient as a random effect (image nested within patient), a pseudobulk-style per-patient summary, or a cell-count-weighted average of per-image statistics (Samorodnitsky and Wu 2024 *Brief Bioinform* 25:bbae522). Two riders complete the picture: cell-type proportions are COMPOSITIONAL (they are constrained to sum to 1, so a real rise in one type mechanically depresses the others, and independent per-type tests double-count this), and acquisition BATCH drifts by day/run and can align with clinical group, so batch must be a covariate and acquisition order randomized against condition. Phenotyping-method choice and statistical-unit choice are orthogonal: getting the cell types right does not excuse testing them wrong.

## Differential Question Taxonomy

| Question | Per-image summary | Patient-level test |
|----------|-------------------|--------------------|
| Cell-type abundance differs between groups | per-image cell-type proportions | mixed model on proportions; scCODA (compositional); diffcyt-DA |
| A functional/state marker differs within a type | per-image median marker per type | pseudobulk per patient + limma/edgeR; diffcyt-DS |
| A spatial interaction/niche differs between groups | per-image enrichment z / Ripley's K / CN abundance | mixed model / cell-count-weighted; SpaceANOVA (FANOVA on cross-K) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Compare cell-type proportions, several types shift | scCODA (or CLR + mixed model) | handles the compositional constraint and the reference-type problem |
| Standard cytometry-style DA on clusters | diffcyt-DA (edgeR on per-sample counts) | designed for per-sample cluster counts; established |
| Multiple ROIs per patient | mixed model, patient random effect (image nested) | respects nesting; or cell-count-weighted aggregate |
| Few patients (n < ~10) | simple per-patient test; report low power honestly | do NOT rescue power with cell count |
| Differential functional-marker expression within a type | pseudobulk per patient + limma/edgeR (diffcyt-DS) | aggregates out cell-level pseudoreplication |
| Spatial interaction/niche across groups | per-image spatial stat -> patient aggregate; SpaceANOVA | the spatial statistic is the summary; the unit is still the patient |
| Acquisition batch aligns with group | include batch covariate; if confounded, the contrast is unrescuable | randomize acquisition order against condition |

## Aggregate to the Patient (the spine)

**Goal:** Collapse millions of cells to one summary per patient before any test.

**Approach:** Compute per-image cell-type proportions, then aggregate images to their patient. Every downstream test consumes this patient-level table, not the cell table.

```python
import pandas as pd

# obs has one row per cell with image_id, patient, condition, cell_type
counts = obs.groupby(['patient', 'condition', 'image_id', 'cell_type']).size().unstack(fill_value=0)
image_prop = counts.div(counts.sum(axis=1), axis=0)            # per-image proportions
patient_prop = image_prop.groupby(['patient', 'condition']).mean()   # one row per patient
```

## Differential Abundance with a Mixed Model

**Goal:** Test a cell type's proportion across groups while respecting patient/ROI nesting.

**Approach:** Fit a mixed model with patient as a random effect when multiple ROIs per patient exist; this absorbs within-patient correlation that a fixed-effect test would treat as independent replication.

```python
import statsmodels.formula.api as smf

# one row per image; proportion of the target type; patient random intercept
df = image_prop.reset_index().rename(columns={'Treg': 'prop'})   # 'Treg' = an actual cell_type column
model = smf.mixedlm('prop ~ condition + batch', df, groups=df['patient'])   # batch as covariate
res = model.fit()
print(res.summary())   # the condition coefficient is tested with patient as the unit
```

## Compositional Differential Abundance (scCODA)

**Goal:** Avoid the false "everything changed" artifact when proportions are constrained to sum to 1.

**Approach:** Model the counts as compositional against a reference cell type; a change is interpreted relative to that reference rather than as an independent per-type shift.

```python
import sccoda.util.cell_composition_data as dat
from sccoda.util import comp_ana as mod

# patient-level cell-type COUNTS (not proportions); pick a biologically stable reference type
data = dat.from_pandas(patient_counts, covariate_columns=['condition'])
analysis = mod.CompositionalAnalysis(data, formula='condition', reference_cell_type='Epithelial')
result = analysis.sample_hmc()
result.summary()
```

## Differential Spatial Feature

**Goal:** Test whether a spatial interaction or niche differs between groups, at the patient unit.

**Approach:** Treat the per-image spatial statistic (a neighborhood-enrichment z, a Ripley's cross-K curve, a CN abundance) as the summary, aggregate to patient, and test across patients with FDR over cell-type pairs. SpaceANOVA does this as a functional ANOVA on per-image cross-K with subject structure.

```python
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# per_image_pair: rows = (image_id, patient, condition, batch, enrichment_z) for ONE type pair
pvals = {}
for pair, sub in per_image_pair.groupby('pair'):
    res = smf.mixedlm('enrichment_z ~ condition + batch', sub, groups=sub['patient']).fit()
    pvals[pair] = res.pvalues['condition[T.responder]']
padj = dict(zip(pvals, multipletests(list(pvals.values()), method='fdr_bh')[1]))   # FDR across pairs
```

## Differential State Within a Type

**Goal:** Test whether a functional/state marker (Ki67, PD-1) differs within a cell type between groups.

**Approach:** Pseudobulk to one value per patient per cell type (median marker expression among that type's cells), then test across patients. diffcyt-DS (R) formalizes this on per-sample medians; the Python pseudobulk path is below.

```python
import numpy as np

t = adata[adata.obs['cell_type'] == 'T cell']
ki67 = t[:, 'Ki67'].X
ki67 = ki67.toarray().ravel() if hasattr(ki67, 'toarray') else np.asarray(ki67).ravel()
pb = t.obs.assign(ki67=ki67).groupby(['patient', 'condition'])['ki67'].median().reset_index()
print(smf.ols('ki67 ~ condition', pb).fit().pvalues['condition[T.responder]'])   # one value per patient
```

## Per-Method Failure Modes

### Cell-level testing
**Trigger:** a t-test/Wilcoxon/regression over individual cells. **Mechanism:** correlated cells from one sample are pseudoreplicates; effective n is the patient count. **Symptom:** p~0 for trivial effects; "significant" findings that do not replicate. **Fix:** aggregate to per-patient summaries; test across patients.

### Independent per-type proportion tests
**Trigger:** a separate test per cell type on proportions. **Mechanism:** proportions sum to 1, so a real rise in one type depresses others mechanically. **Symptom:** many types appear to change in opposite directions. **Fix:** compositional model (scCODA) or CLR transform with a reference type.

### Over-correcting batch
**Trigger:** aggressive integration to make clusters patient-agnostic, then testing on corrected data. **Mechanism:** correction can treat real between-patient biology as batch. **Symptom:** the disease signal disappears. **Fix:** correct minimally, validate that invariant types align while variable types stay separate, and keep integration out of the across-patient inference path.

### Rescuing power with cell count
**Trigger:** claiming significance from n=4 patients because millions of cells were imaged. **Mechanism:** cell count is not replication. **Symptom:** confident claims from few patients. **Fix:** report the patient n and the true power honestly; collect more patients.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Unit of replication = patient count (often 10-40) | Squair 2021 *Nat Commun* 12:5692 | cells/ROIs are pseudoreplicates |
| Cell-count-weighted per-image aggregation | Samorodnitsky and Wu 2024 *Brief Bioinform* 25:bbae522 | controls type-I error with high power (vs unweighted ROI averaging) |
| Reference cell type for compositional DA | Buttner 2021 *Nat Commun* 12:6876 | proportions are not independent |
| BH-FDR across cell-type pairs and radii | multiplicity | ~200 pairs x radii guarantees false positives |
| Batch covariate; randomize acquisition order | spatial dossier | run drift can align with clinical group |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| p~0 with a trivial effect | cell-level test | per-patient aggregation; mixed model |
| All cell types "changed" | compositional constraint ignored | scCODA / CLR with reference type |
| Disease signal vanished after integration | batch over-correction | minimal correction; keep it out of the test |
| Significant with n=4 patients | power rescued by cell count | report patient n; do not over-claim |
| Many significant pairs | no FDR across pairs/radii | BH-FDR over the full grid |

## References

- Squair JW, Gautier M, Kathe C, et al. 2021. Confronting false discoveries in single-cell differential expression. *Nat Commun* 12:5692. — pseudoreplication; aggregate to sample.
- Samorodnitsky S, Wu MC. 2024. Statistical analysis of multiple regions-of-interest in multiplexed spatial proteomics data. *Brief Bioinform* 25(6):bbae522. — ROI aggregation, cell-count-weighted averaging, SPOT omnibus.
- Seal S, Neelon B, Angel PM, et al. 2024. SpaceANOVA: Spatial Co-occurrence Analysis of Cell Types in Multiplex Imaging Data Using Point Process and Functional ANOVA. *J Proteome Res* 23(4):1131-1143. — FANOVA on per-image cross-K with subject structure.
- Weber LM, Nowicka M, Soneson C, Robinson MD. 2019. diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. *Commun Biol* 2:183. — diffcyt-DA/DS.
- Buttner M, Ostner J, Muller CL, Theis FJ, Schubert B. 2021. scCODA is a Bayesian model for compositional single-cell data analysis. *Nat Commun* 12:6876. — compositional differential abundance.
- Schurch CM, Bhate SS, Barlow GL, et al. 2020. Coordinated Cellular Neighborhoods Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive Front. *Cell* 182(5):1341-1359.e19. — cellular neighborhoods compared across patients.

## Related Skills

- phenotyping - supplies the cell-type labels whose proportions are compared
- spatial-analysis - supplies the per-image spatial statistics that become patient-level summaries
- quality-metrics - batch must be diagnosed and entered as a covariate
- experimental-design/randomization-blocking - the experimental-unit and pseudoreplication foundation
- clinical-biostatistics/subgroup-analysis - multiplicity and effect estimation in clinical cohorts
- flow-cytometry/differential-analysis - diffcyt-DA/DS for suspension cytometry
