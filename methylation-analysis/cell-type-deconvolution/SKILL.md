---
name: bio-methylation-cell-type-deconvolution
description: Estimates cell-type composition from bulk DNA methylation and uses it to defuse the single biggest EWAS confounder. Covers reference-based deconvolution (Houseman constrained-projection, minfi estimateCellCounts2 with FlowSorted.Blood.EPIC + IDOL-optimized libraries, EpiDISH RPC/CBS/CP, 12-cell extended, cord-blood nRBC references, EpiSCORE/hepidish for solid tissue), reference-free correction (ReFACTor, RefFreeEWAS, SVA), using fractions as covariates vs the compositionality/collinearity trap, and cell-type-resolved EWAS (CellDMC, TCA, TOAST, omicwas, HIRE). Use when estimating blood/tissue cell fractions, adjusting an EWAS for composition, choosing a deconvolution reference, or attributing a methylation signal to a cell type. For the EWAS confounder-vs-mediator decision see ewas-design; for the IEAA cell-count adjustment of DNAm age see epigenetic-clocks; for clean beta input see array-preprocessing.
tool_type: r
primary_tool: EpiDISH
---

## Version Compatibility

Reference examples tested with: EpiDISH 2.18+, minfi 1.48+, FlowSorted.Blood.EPIC 2.0+, FlowSorted.CordBloodCombined.450k 1.20+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The REFERENCE package is the version that matters most. A reference package is platform-, tissue-, and age-specific: `FlowSorted.Blood.EPIC` ships `IDOLOptimizedCpGs` (EPIC) and `IDOLOptimizedCpGs450klegacy` (450K) as distinct libraries, the 12-cell library lives in a separate `FlowSorted.BloodExtended.EPIC` package, and cord blood needs `FlowSorted.CordBloodCombined.450k` (it carries nucleated red blood cells). `estimateCellCounts2` returns a `Neu` (neutrophil) column where the older `minfi::estimateCellCounts` returns `Gran` - the label changes the downstream column names. Record the reference package and version alongside the array build.

# Cell-Type Deconvolution

**"How much of my methylation signal is just cell composition?"** -> Project the bulk beta matrix onto a purified-cell reference to estimate per-sample fractions, then carry those fractions forward as covariates - because a bulk methylome is a composition-weighted average and a composition difference is a methylation difference.
- R: `epidish(beta.m, ref.m = centDHSbloodDMC.m, method = 'RPC')$estF`

Scope: estimate cell-type fractions from a clean bulk beta matrix and use them downstream. Clean beta read-in and EPICv2 replicate-probe collapse -> array-preprocessing. The confounder-vs-mediator decision and the EWAS regression itself -> ewas-design. Adjusting a DNAm clock for cell counts (IEAA) -> epigenetic-clocks. Single-cell/sorted atlases for reference building and validation ground truth -> single-cell/preprocessing. Predictive-model training/leakage -> machine-learning/biomarker-discovery.

## The Single Most Important Modern Insight -- A Cell-Fraction Estimate Is a Projection, Not a Measurement

A reference-based fraction is not a measurement of a sample's composition; it is a projection of that sample onto cell types someone else purified, on someone else's platform, in someone else's tissue. Three corollaries each common misuse violates:

1. **A composition difference IS a methylation difference.** Bulk DNAm is the fraction-weighted average of its constituent cell-type methylomes, and most CpG variance is between cell types, not between conditions. If cases and controls differ in composition - which they almost always do (age, sex, infection, smoking shift the neutrophil-to-lymphocyte ratio) - the EWAS reports cell-count differences as if they were disease methylation. This is the #1 EWAS confounder (Jaffe & Irizarry 2014 *Genome Biol* 15:R31).
2. **The reference defines the answer.** A cell type present in the sample but absent from the reference is silently redistributed onto the nearest reference types - no error, the fractions still sum to ~1. Cord blood without nRBC, a solid tissue against a blood reference, EPICv2 data against a 450K library: all return confident, wrong proportions.
3. **Fractions are compositional.** They live on a simplex (sum ~1), so they are not independent: one going up forces others down. Naively co-regressing or correlating all K fractions manufactures spurious negative associations.

Organize the analysis around matching the reference and handling compositionality, not around picking an algorithm. Deconvolution turns an uncontrollable confounder into a measurable covariate - but only as accurately as the reference matches the sample.

## Reference-Based: The Houseman Constrained-Projection Foundation

Houseman 2012 (*BMC Bioinformatics* 13:86) is the origin. From a matrix of FACS/MACS-purified cell-type mean methylation at discriminating CpGs (L-DMRs), solve for each sample a constrained quadratic program: the non-negative fraction vector w (w_i >= 0, sum ~1) minimizing the squared distance between observed beta and reference x w over the L-DMR CpGs. This "CP" (constrained projection) is what every later method is measured against. The L-DMR selection is itself a tuning choice that the IDOL work (Koestler 2016 *BMC Bioinformatics* 17:120) optimized into a fixed, benchmarked library.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| EpiDISH (RPC) | Teschendorff 2017 *BMC Bioinformatics* 18:105 | robust partial correlation; downweights noisy CpGs | general; the robust default across tissues/noise |
| EpiDISH (CP) | Houseman 2012 *BMC Bioinformatics* 13:86 | constrained quadratic projection | reproduce the classic Houseman estimate |
| EpiDISH (CBS) | Newman 2015 *Nat Methods* 12:453 | CIBERSORT nu-SVR | borrowed from expression; an alternative |
| minfi estimateCellCounts2 | Salas 2018 *Genome Biol* 19:64 | Houseman projection on the IDOL-optimized EPIC/450K library | from an RGChannelSet; modern 6-cell blood (Neu) |
| FlowSorted.BloodExtended.EPIC | Salas 2022 *Nat Commun* 13:761 | 12-cell IDOL library | naive/memory T, Treg, eosinophil/basophil resolution |
| hepidish | Teschendorff 2017 *BMC Bioinformatics* 18:105 | hierarchical Epi/Fib/Immune then immune subtypes | solid tissue with immune infiltration |
| EpiSCORE | Teschendorff 2020 *Genome Biol* 21:221 | scRNA-seq-imputed DNAm reference | solid tissues with no sorted reference |
| ReFACTor | Rahmani 2016 *Nat Methods* 13:443 | sparse-PCA components as covariates | reference-free; no matched reference exists |
| RefFreeEWAS | Houseman 2014 *Bioinformatics* 30:1431 | NMF/SVD-style latent cell-mixture | reference-free; unlabeled components |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Adult whole blood | estimateCellCounts2 IDOL (6: Neu/CD4T/CD8T/NK/Bcell/Mono) or EpiDISH RPC (centDHSbloodDMC.m gives 7, adds Eos) | benchmarked blood references |
| Need naive/memory T, Treg, Eos, Bas | 12-cell FlowSorted.BloodExtended.EPIC | the 6-cell library cannot resolve these |
| Cord blood / newborn | FlowSorted.CordBloodCombined.450k | adds nRBC; an adult reference is silently wrong |
| Saliva / buccal | epithelial + immune reference (hepidish) | saliva is not blood; epithelial fraction dominates |
| Solid tissue / tumor with infiltration | hepidish or EpiSCORE | flat blood reference on solid tissue is meaningless |
| 450K data | EpiDISH cent*450k.m / IDOLOptimizedCpGs450klegacy | platform-matched CpGs; EPIC library drops CpGs |
| EPICv2 data | collapse replicate probes first -> array-preprocessing | suffixed replicate beads hide the reference CpGs |
| No matched reference (novel tissue) | ReFACTor / RefFreeEWAS + sensitivity | reference-free fallback; components are unlabeled |
| Which cell type drives a signal | CellDMC / TCA / TOAST (below) | model composition, do not just regress it out |
| Confounder-vs-mediator decision | -> ewas-design | upstream: adjust out, or resolve cell-specific? |

## Estimate Blood Fractions with EpiDISH (RPC)

**Goal:** Get per-sample fractions of the major immune cell types from a clean beta matrix to use as EWAS covariates.

**Approach:** Pass the beta matrix and a tissue-matched reference centroid to `epidish` with `method='RPC'` (the robust option), then read the sample-by-cell-type matrix from `$estF`.

```r
library(EpiDISH)
data(centDHSbloodDMC.m)    # 7 immune cell types, adult whole blood

out <- epidish(beta.m = beta_matrix, ref.m = centDHSbloodDMC.m, method = 'RPC')
fractions <- out$estF       # samples x cell types; rows sum to ~1
```

## Estimate Blood Fractions from an RGChannelSet (minfi + IDOL)

**Goal:** Estimate the modern 6-cell IDOL blood composition straight from raw IDAT-derived data.

**Approach:** Run `estimateCellCounts2` on the RGChannelSet with the IDOL probe selection and the platform-matched reference; for 450K data switch the reference library so the same cell types are estimated cross-platform.

```r
library(FlowSorted.Blood.EPIC)

counts <- estimateCellCounts2(
  rgSet,
  compositeCellType = 'Blood',
  processMethod = 'preprocessNoob',
  probeSelect = 'IDOL',
  cellTypes = c('CD8T', 'CD4T', 'NK', 'Bcell', 'Mono', 'Neu'),   # Neu, not Gran
  referencePlatform = 'IlluminaHumanMethylationEPIC'
)$counts
```

## Solid Tissue: Hierarchical Deconvolution

**Goal:** Deconvolve a solid tissue (epithelial + fibroblast + infiltrating immune) rather than forcing a blood reference onto it.

**Approach:** Use `hepidish` to first split Epithelial/Fibroblast/total-Immune, then deconvolve the immune fraction into subtypes and multiply through. For tissues with no sorted reference at all, EpiSCORE builds an imputed DNAm reference from a single-cell RNA atlas.

```r
library(EpiDISH)
data(centEpiFibIC.m)       # Epithelial / Fibroblast / Immune-Cell
data(centBloodSub.m)       # immune subtypes for the second level

frac <- hepidish(beta.m = beta_matrix, ref1.m = centEpiFibIC.m,
                 ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
# h.CT.idx = 3 = the Immune column in ref1 to expand with ref2
```

## Reference-Free Correction

**Goal:** Capture composition structure when no matched reference exists, accepting unlabeled components.

**Approach:** ReFACTor selects the most composition-informative CpGs and runs sparse-PCA; use the top components as EWAS covariates. RefFreeEWAS decomposes the matrix into a latent cell-mixture term. Both correct without naming the cell types, so check that genuine top hits survive (they can absorb real signal).

```r
library(TCA)
ref <- refactor(beta_matrix, k = 6)   # k = expected number of cell types
covariates <- ref$scores              # top sparse-PC components as EWAS covariates
```

## Using the Fractions: Covariate vs Cell-Type-Resolved

There are two distinct moves once fractions exist, and they answer different questions.

**As covariates (the standard EWAS defense).** Include the fractions in the per-CpG design matrix so composition is regressed out. Because fractions are compositional (sum ~1), do NOT enter all K - drop one reference cell type (or use a compositional transform) to avoid perfect collinearity. The confounder-vs-mediator decision (regress out, or treat composition as the mechanism) belongs to ewas-design; execution belongs to differential-cpg-testing.

**Cell-type-resolved EWAS (which cell type drives the signal).** Instead of regressing composition away, model a phenotype x cell-fraction INTERACTION per CpG to ask which cell type carries the differential methylation and in which direction. CellDMC (Zheng 2018 *Nat Methods* 15:1059) is the simplest member; a family generalizes it:

| Method | Citation | Adds beyond the interaction | Output |
|--------|----------|-----------------------------|--------|
| CellDMC | Zheng 2018 *Nat Methods* 15:1059 | per-CpG linear pheno x fraction interaction | which cell type is DM + direction (a test) |
| TCA | Rahmani 2019 *Nat Commun* 10:3417 | tensor model; per-sample per-cell-type levels | cell-type-specific methylation + association test |
| TOAST | Li & Wu 2019 *Genome Biol* 20:190 | iterative csDM; improves reference-free composition | csTest per cell type; runs reference-free |
| omicwas | Takeuchi & Kato 2021 *BMC Bioinformatics* 22:141 | nonlinear ridge for the logit scale + fraction collinearity | cell-type-specific association statistics |
| HIRE | Luo 2019 *Nat Commun* 10:3113 | joint multiplicative-composition hierarchical model | risk-CpG sites per cell type |

```r
library(EpiDISH)
res <- CellDMC(beta.m = beta_matrix, pheno.v = phenotype, frac.m = fractions)
# res$dmct: per-CpG, which cell type is differentially methylated (-1/0/1)
```

A cell-type-resolved call is an ill-posed inverse problem regularized by an assumed reference: rare cell types (2-5% of the mixture) are badly underpowered, fraction collinearity destabilizes the interactions, and deconvolution error propagates straight into the attribution (HIRE's argument for estimating composition jointly). Validation is hard without sorted/single-cell ground truth - method papers lean on simulations and reconstructed mixtures, which are circular. Treat an in-silico cell-type-specific hit as a HYPOTHESIS about what to sort next, not a finding; confirm load-bearing attributions in sorted or single-cell DNAm from independent samples (Walker 2025 *Brief Bioinform* 26:bbaf427).

## The IEAA Link to Clocks

Intrinsic epigenetic age acceleration (IEAA) is DNAm age residualized on chronological age AND estimated blood cell counts - so deconvolution is the prerequisite step: estimate fractions here, then hand them to epigenetic-clocks as the cell-count covariates that distinguish cell-intrinsic aging from a composition shift. Do not teach the clock here; compute the fractions and route the IEAA adjustment to epigenetic-clocks.

## Per-Method Failure Modes

### Missing cell type silently redistributed
**Trigger:** a sample contains a cell type absent from the reference (cord-blood nRBC, a rare infiltrate, a granulocyte subtype collapsed to Gran). **Mechanism:** the constrained projection has no column for it, so its signal lands on the nearest present types. **Symptom:** plausible-looking fractions that sum to ~1 with no warning. **Fix:** match the reference to tissue+age (FlowSorted.CordBloodCombined.450k for newborns; hepidish/EpiSCORE for solid tissue).

### Platform-mismatched reference library
**Trigger:** 450K data with the EPIC IDOL library, or EPICv2 with either. **Mechanism:** reference CpGs are partly absent on the other platform, shrinking the L-DMR set used for the projection. **Symptom:** biased fractions, no error. **Fix:** `IDOLOptimizedCpGs450klegacy` / `cent*450k.m` for 450K; collapse EPICv2 replicate probes first (-> array-preprocessing).

### Collinear cell-fraction covariates
**Trigger:** entering all K fractions (sum ~1) into a design matrix. **Mechanism:** the simplex constraint makes the K-th fraction a linear function of the others. **Symptom:** rank-deficient design, dropped coefficient, or spurious negative fraction-fraction correlations. **Fix:** drop one reference cell type or use a compositional (CLR/ILR) transform.

### Reference-free over-correction
**Trigger:** including too many ReFACTor/RefFreeEWAS components, or using them when a reference exists. **Mechanism:** unlabeled latent components can absorb true biological signal alongside composition. **Symptom:** top EWAS hits vanish; false negatives. **Fix:** prefer reference-based when a reference exists; use reference-free as a fallback/sensitivity check and confirm hits survive.

### Cell-type attribution from rare cells
**Trigger:** reading a CellDMC/TCA call for a 2-5% cell type. **Mechanism:** a rare cell contributes a fraction-attenuated slice of bulk variance, so its interaction estimate is dominated by deconvolution noise. **Symptom:** confident-looking csDM in basophils/eosinophils; nulls misread as "no effect." **Fix:** report each cell type's mean fraction; distrust specific calls for low-abundance types; never infer absence of effect from an underpowered null.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| IDOL EPIC 6-cell library ~450 CpGs | Salas 2018 *Genome Biol* 19:64 | benchmarked L-DMR set; R^2 ~0.992 on reconstructed mixtures |
| method = 'RPC' for EpiDISH | Teschendorff 2017 *BMC Bioinformatics* 18:105 | robust to outlier/noisy CpGs; more stable than CP across tissues |
| drop 1 of K fractions as covariates | compositional constraint | fractions sum to ~1, so all K are perfectly collinear |
| cord blood reference must carry nRBC | Gervin 2019 / CordBloodCombined | nRBC abundant in cord blood, absent from adult references |
| ReFACTor k = expected cell-type count | Rahmani 2016 *Nat Methods* 13:443 | k sets the rank; too high over-corrects, too low under-corrects |
| csDM credible only for abundant types | Walker 2025 *Brief Bioinform* 26:bbaf427 | rare cells are fraction-attenuated and underpowered |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Fractions look fine but EWAS still inflated | unmodeled cell type / wrong reference | match reference to tissue+age+platform |
| Design matrix rank-deficient | all K fractions entered as covariates | drop one cell type or CLR-transform |
| estimateCellCounts2 returns Neu, code expects Gran | minfi vs FlowSorted label difference | use Neu (estimateCellCounts2) consistently |
| Reference CpGs not found on EPICv2 | replicate probes not collapsed | collapse to one value per CpG first |
| Negative or all-zero fraction for a type | platform mismatch / absent in sample | check platform-matched library; inspect mean fraction |
| csDM hit in a rare cell type | underpowered interaction | report the fraction; validate by sorting/single-cell |

## References

- Houseman EA, Accomando WP, Koestler DC, et al. 2012. DNA methylation arrays as surrogate measures of cell mixture distribution. *BMC Bioinformatics* 13:86.
- Jaffe AE, Irizarry RA. 2014. Accounting for cellular heterogeneity is critical in epigenome-wide association studies. *Genome Biol* 15:R31.
- Koestler DC, Jones MJ, Usset J, et al. 2016. Improving cell mixture deconvolution by identifying optimal DNA methylation libraries (IDOL). *BMC Bioinformatics* 17:120.
- Salas LA, Koestler DC, Butler RA, et al. 2018. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. *Genome Biol* 19:64.
- Salas LA, Zhang Z, Koestler DC, et al. 2022. Enhanced cell deconvolution of peripheral blood using DNA methylation for high-resolution immune profiling. *Nat Commun* 13:761.
- Teschendorff AE, Breeze CE, Zheng SC, Beck S. 2017. A comparison of reference-based algorithms for correcting cell-type heterogeneity in epigenome-wide association studies. *BMC Bioinformatics* 18:105.
- Teschendorff AE, Zhu T, Breeze CE, Beck S. 2020. EPISCORE: cell type deconvolution of bulk tissue DNA methylomes from single-cell RNA-Seq data. *Genome Biol* 21:221.
- Houseman EA, Molitor J, Marsit CJ. 2014. Reference-free cell mixture adjustments in analysis of DNA methylation data. *Bioinformatics* 30:1431-1439.
- Rahmani E, Zaitlen N, Baran Y, et al. 2016. Sparse PCA corrects for cell type heterogeneity in epigenome-wide association studies. *Nat Methods* 13:443-445.
- Zheng SC, Breeze CE, Beck S, Teschendorff AE. 2018. Identification of differentially methylated cell types in epigenome-wide association studies. *Nat Methods* 15:1059-1066.
- Rahmani E, Schweiger R, Rhead B, et al. 2019. Cell-type-specific resolution epigenetics without the need for cell sorting or single-cell biology. *Nat Commun* 10:3417.
- Li Z, Wu H. 2019. TOAST: improving reference-free cell composition estimation by cross-cell type differential analysis. *Genome Biol* 20:190.
- Takeuchi F, Kato N. 2021. Nonlinear ridge regression improves cell-type-specific differential expression analysis. *BMC Bioinformatics* 22:141.
- Luo X, Yang C, Wei Y. 2019. Detection of cell-type-specific risk-CpG sites in epigenome-wide association studies. *Nat Commun* 10:3113.
- Walker EM, Dempster EL, Franklin A, et al. 2025. Guidance for the design and analysis of cell-type-specific DNA methylation epidemiology studies. *Brief Bioinform* 26:bbaf427.

## Related Skills

- array-preprocessing - Provides the clean beta matrix deconvolution consumes
- ewas-design - Cell-fraction covariate strategy (confounder vs mediator)
- epigenetic-clocks - IEAA: adjust the clock for estimated cell composition
- differential-cpg-testing - Uses cell fractions as design-matrix covariates
- single-cell/preprocessing - scRNA atlases for reference building (EpiSCORE) and ground truth
- machine-learning/biomarker-discovery - Predictive-model boundary
- workflows/methylation-pipeline - End-to-end pipeline
