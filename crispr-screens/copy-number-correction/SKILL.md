---
name: bio-crispr-screens-copy-number-correction
description: Corrects the gene-independent copy-number artifact in CRISPR-Cas9 screens (Aguirre 2016 / Munoz 2016 Cancer Discov) where amplified loci appear essential from DNA-damage burden of simultaneous cuts. Covers the p53-dependent G2-arrest mechanism, CRISPRcleanR (Iorio 2018) unsupervised pre-hoc correction, CERES (Meyers 2017) joint CN + gene-effect model, Chronos (Dempster 2021) DepMap-standard population-dynamics + CN model with lowest residual bias, the decision tree by data availability, the Spearman LFC-vs-CN diagnostic, focal-amplification examples (ERBB2 in HER2+, MYC in colorectal, FGFR1 in head and neck), and CRISPRi/a alternatives that bypass the artifact. Use when screening cancer cell lines, diagnosing essentiality at amplified loci, choosing CRISPRcleanR / CERES / Chronos, deciding whether CN correction is needed before MAGeCK / BAGEL2 / drugZ, or switching from Cas9 to CRISPRi.
tool_type: mixed
primary_tool: CRISPRcleanR
---

## Version Compatibility

Reference examples tested with: CRISPRcleanR 3.0+ (R/Bioconductor), Chronos 2.0+ (https://github.com/broadinstitute/chronos), CERES (legacy, superseded by Chronos), pandas 2.2+, numpy 1.26+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('CRISPRcleanR')`; `?ccr.CleanCN`
- Python: `pip show chronos-cn`; `chronos --help`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## Copy-Number Bias Correction in CRISPR Screens

**"Correct copy-number artifacts in my cancer-cell-line screen"** -> Identify gene-independent depletion at amplified loci, apply CRISPRcleanR (pre-hoc, unsupervised, position-based) or Chronos (joint model, supervised with CN profile) to remove the artifact, then proceed to hit calling on corrected data.

- R: `CRISPRcleanR::ccr.CleanCN()` for unsupervised pre-hoc correction (no CN profile required)
- Python: Chronos (`chronos-cn`) for joint cell-population dynamics + CN modeling
- Python: CERES (legacy, superseded by Chronos)

## The Copy-Number Artifact (Mechanism)

**Aguirre AJ et al 2016 *Cancer Discov* 6:914** and **Munoz DM et al 2016 *Cancer Discov* 6:900** demonstrated that focal amplification regions in cancer cell lines appear systematically "essential" in CRISPR-Cas9 screens, independent of the gene's actual biology. The mechanism:

1. A focal amplification creates 4-50+ copies of a genomic region.
2. Each sgRNA targeting a gene in that region cuts at all copies simultaneously.
3. Multiple cuts trigger a p53-dependent DNA-damage response.
4. Cells arrest in G2 phase; the sgRNA appears depleted because its bearer cells don't proliferate.
5. The depletion is proportional to the number of simultaneous cuts, not the gene's essentiality.

**Consequence:** ERBB2 appears essential in HER2-amplified SK-BR-3 (24 copies). MYC appears essential in MYC-amplified colorectal lines (10+ copies). FGFR1 appears essential in FGFR1-amplified head-and-neck lines. These are all false positives.

**Affects:** All Cas9-KO screens in cancer cell lines. Universal, not conditional. Cannot be remediated by sequencing depth, library size, or replicate count. Requires explicit correction.

**Bypassed by:**
- CRISPRi (catalytically dead Cas9, no DNA damage) -> no artifact
- CRISPRa (catalytically dead Cas9) -> no artifact
- Base editing (single-strand nick + deaminase) -> reduced artifact (Hess 2016 demonstration)
- Prime editing (nick + RT) -> reduced artifact

## Correction Method Decision Tree

| Available data | Recommended method | Why |
|----------------|---------------------|-----|
| Cell-line panel without matched CN profile | CRISPRcleanR | Unsupervised; uses genomic position only |
| Single cell line with matched WGS/SNP-array CN | CRISPRcleanR or Chronos | Either works; Chronos more rigorous |
| DepMap-scale (1000+ cell lines, longitudinal) | Chronos | Population-dynamics + screen quality + CN; DepMap quarterly standard |
| Single cell line, multi-timepoint | Chronos | Leverages longitudinal counts |
| Need to integrate with downstream MAGeCK | CRISPRcleanR (pre-hoc) | Outputs corrected counts for any downstream tool |
| Multiple cell lines + multiple batches | Chronos | Joint modeling of all dimensions |

## CRISPRcleanR (Iorio 2018) - Unsupervised Pre-Hoc

**Goal:** Correct copy-number bias without requiring matched CN profile by detecting position-based systematic enrichment / depletion patterns.

**Approach:** Order sgRNAs by chromosomal coordinate; detect segments where sgRNAs show systematic depletion (or enrichment) inconsistent with single-gene biology; shift these segments toward the global mean. The intuition: focal amplifications create depletion bands extending tens to hundreds of kb; non-amplified essential genes are punctate.

```r
library(CRISPRcleanR)

# Load library annotation (sgRNA -> chromosomal coordinates)
data(KY_Library_v1.0)   # KY library; replace with your library annotation
# OR use ccr.PrepareAnnotations() to make custom

# Load count data with first 2 cols: sgRNA, gene, then sample counts
counts <- read.table('counts.txt', header=TRUE, sep='\t')

# 1. Normalize and compute logFC
norm_counts <- ccr.NormfoldChanges(counts, min_reads=30,
                                     EXPname='my_screen',
                                     libraryAnnotation=KY_Library_v1.0)

# 2. Compute genome-sorted sgRNA fold changes
gw_log_fc <- ccr.logFCs2chromPos(norm_counts$norm_fold_changes,
                                   KY_Library_v1.0)

# 3. Apply CRISPRcleanR correction
corrected <- ccr.GWclean(gw_log_fc, display=TRUE, label='my_screen')
# Output: corrected$corrected_logFCs and corrected$segments
# corrected_logFCs can replace LFCs downstream

# 4. Re-derive corrected counts for downstream MAGeCK
corrected_counts <- ccr.correctCounts(my_screen=norm_counts,
                                        correction=corrected,
                                        outprefix='cleanr_corrected',
                                        libraryAnnotation=KY_Library_v1.0)
```

**Key parameter:** `min_reads=30` is the lower-count threshold for inclusion. This must match the library-coverage strategy; too high removes legitimate guides, too low keeps noisy guides.

**Output:** Pre-corrected LFCs and counts that can be fed into MAGeCK / BAGEL2 / drugZ as if they were the original screen data. The correction is independent of CN profile (unsupervised) and works on cell lines without matched WGS.

## Chronos (Dempster 2021) - Joint Population-Dynamics + CN Model

**Goal:** Estimate gene fitness while jointly accounting for copy-number-driven depletion, screen quality, and longitudinal cell-population dynamics.

**Approach:** Model the cell population over time as an ODE driven by per-gene fitness effects; add a separate term for copy-number-driven depletion; estimate all parameters via maximum-likelihood with regularization. Outputs a "gene effect score" normalized against the empirical distributions of essential and non-essential reference genes.

```python
# Chronos via the chronos-cn package
import chronos

# Inputs
# 1. Counts: rows = sgRNA, columns = samples (per-timepoint per-cell-line)
# 2. Sequence map: sgRNA -> cell line -> sample timepoint
# 3. Copy-number profile per cell line per genomic region
# 4. Guide-gene map

model = chronos.Chronos(
    sequence_map=sequence_map,           # which samples are from which cell line + condition
    guide_gene_map=guide_gene_map,
    reads=counts_df,
    copy_number=copy_number_df,           # per-cell-line CN profile
    pretrained_offset=None,
)
model.train(
    n_steps=2000,
    learning_rate=0.1,
    verify_normalize=True,
)
gene_effects = model.gene_effect()        # cells x genes; standardized score
gene_probabilities = model.gene_probability()  # probability of being essential
```

**DepMap convention:** A gene-effect score <-1 corresponds to "essential" in that cell line; <-0.5 is "depleting." Each DepMap release (quarterly) provides Chronos gene effects and probabilities.

**Critical:** Chronos benefits most from longitudinal data (multiple timepoints per cell line) but can run with multiple cell lines at single timepoint; it cannot run with a single screen (one line, one timepoint) without matched CN profile. For single-cell-line, single-timepoint screens without matched CN, use CRISPRcleanR instead.

## CERES (Legacy, Superseded by Chronos)

**Meyers RM et al 2017 *Nat Genet* 49:1779** introduced the first formal CN-correction method at DepMap scale. CERES decomposes per-sgRNA LFC as `sgRNA_efficacy * gene_effect - CN_term(copy_number)`, fitting jointly. Superseded by Chronos at DepMap in 2021 due to Chronos' better handling of screen quality and longitudinal data. CERES remains useful for cross-validation.

## Detect Uncorrected CN Bias

**Goal:** Verify that copy-number bias is corrected (or detect it in raw data).

**Approach:** For genes with matched CN profile, compute Spearman ρ between gene-level LFC and copy number. A negative correlation (-ρ) indicates amplified genes are depleted, i.e., CN artifact.

```python
import pandas as pd
from scipy.stats import spearmanr

def detect_cn_bias(gene_lfc_df, cn_df):
    '''Test whether gene-level LFC negatively correlates with copy number.
    A bias-free screen has Spearman rho near zero between CN and LFC.'''
    merged = gene_lfc_df.merge(cn_df, on='gene')
    rho, p = spearmanr(merged['copy_number'], merged['lfc'])
    return {
        'cn_lfc_rho': rho,
        'p_value': p,
        'amplified_mean_lfc': merged[merged['copy_number'] > 4]['lfc'].mean(),
        'diploid_mean_lfc': merged[(merged['copy_number'] >= 1.5) & (merged['copy_number'] <= 2.5)]['lfc'].mean(),
        'bias_present': rho < -0.1 and p < 0.01,
    }
```

**Threshold (Aguirre 2016):** Spearman ρ <-0.10 between LFC and CN indicates significant CN bias. Even modest amplifications (4+ copies) generate detectable artifact. Run this diagnostic before AND after correction.

## Reconciliation: When CN Correction Fails

If post-CRISPRcleanR or post-Chronos the CN-LFC Spearman is still significantly negative, the correction is incomplete. Possible causes:

1. **Insufficient CN resolution:** A specific 4-copy region went undetected. Refine CN profile with deeper WGS.
2. **CRISPRcleanR position-based correction missed it:** The amplification is small relative to the segmentation algorithm's resolution. Use Chronos with matched CN profile.
3. **Genomic rearrangement creates a "ghost" amplification:** A complex rearrangement appears as normal CN but Cas9 cuts at multiple sites due to translocation breakpoints. Combine WGS structural variants with the analysis.
4. **Cell line has p53 wild-type with extreme cut-toxicity sensitivity:** The artifact may persist; use CRISPRi screens for that line.

## Apply CN Correction to Pipeline

**Workflow:**

```
1. mageck count (raw counts)
2. screen-qc verification
3. CN diagnostic: Spearman of LFC vs CN (if CN profile available)
4. If bias detected:
   a. CRISPRcleanR (pre-hoc) -> corrected counts -> MAGeCK / BAGEL2 / drugZ
   OR
   b. Chronos (joint model with CN profile) -> gene effects directly
5. Re-diagnose: Spearman of CORRECTED LFC vs CN should be near zero
6. Hit calling
```

For DepMap-style large panels:
```
Chronos handles batch + CN + screen quality in one step; no pre-correction needed.
```

For Sanger Score-style panel (Behan 2019):
```
CRISPRcleanR was used historically; cross-check with Chronos when CN profile available.
```

## Failure Modes

### CRISPRcleanR removes legitimate essential signal

**Trigger:** A genuine essential gene happens to lie in a region with adjacent uncorrected non-essential signal; the segment-based correction includes the essential.
**Mechanism:** CRISPRcleanR's `ccr.GWclean()` segments sgRNAs by position; segments containing multiple genes with directional consistency are corrected as a unit.
**Symptom:** A known essential drops out of post-correction hit list.
**Fix:** Inspect segments manually; if a known essential was within a corrected segment, investigate. Cross-check with non-CN-corrected MAGeCK + BAGEL2 to see if essential was a hit pre-correction.

### Chronos fails on single-timepoint or single-cell-line data

**Trigger:** Chronos requires multiple timepoints (or multiple cell lines) for population-dynamics estimation.
**Mechanism:** Single observation per condition leaves model under-determined.
**Symptom:** Chronos errors out or produces flat gene-effect distributions.
**Fix:** Use CRISPRcleanR (which handles single-timepoint single-line); collect multi-timepoint data for Chronos.

### Spearman ρ still negative after CRISPRcleanR

**Trigger:** Amplification is too small or complex for the segment-based approach.
**Mechanism:** CRISPRcleanR detects systematic spatial patterns; isolated 4-copy regions can slip through.
**Symptom:** Post-correction Spearman ρ -0.05 to -0.10 between LFC and CN.
**Fix:** Refine CN profile (deeper WGS); apply Chronos with matched CN as alternative; or supplement with focal-amplification-aware methods.

### Cell line lacks matched CN profile

**Trigger:** Newly characterized line or rare patient-derived line; WGS not done.
**Mechanism:** Chronos requires CN as input; CRISPRcleanR doesn't but works better with it.
**Symptom:** Cannot apply Chronos; CRISPRcleanR less precise without supervised CN.
**Fix:** Run SNP-array (cheap, fast) or low-coverage WGS to obtain CN profile; in interim, use CRISPRcleanR unsupervised mode.

### CN amplification at non-coding region drives apparent essentiality

**Trigger:** Amplification at a gene-poor region; sgRNAs at edge genes get artifactually depleted.
**Mechanism:** Even non-essential genes adjacent to amplifications are depleted because the Cas9 cuts are at the amplified loci.
**Symptom:** Non-essential genes near amplification show LFC <0.
**Fix:** Inspect chromosomal position of "essential" hits; flag genes within 100 kb of known amplifications for orthogonal validation. This is the classic Aguirre 2016 observation.

## CRISPRi/a Alternative

**For variant-function or non-cancer-line essentiality screens, switching to CRISPRi (catalytically dead dCas9-KRAB) avoids the artifact entirely.** No DNA double-strand breaks = no p53 G2 arrest = no copy-number-driven depletion.

| Approach | CN artifact | When to use |
|----------|--------------|-------------|
| Cas9 KO | YES; requires correction | Loss-of-function essentiality, traditional screens |
| CRISPRi | NO | Cancer lines with focal amps; knockdown of cuttable-toxic genes |
| CRISPRa | NO | Gain-of-function; activation screens |
| Base editing | Reduced (single-strand nick) | Variant function |
| Prime editing | Reduced | Precise edits |

See [[library-design]] for CRISPRi (Dolcetto) and CRISPRa (Calabrese) library options.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Spearman ρ (CN vs LFC) | <-0.10 -> bias present | Aguirre 2016 *Cancer Discov* 6:914 |
| Copies for detectable artifact | 4+ | Aguirre 2016: even 4-copy regions generate measurable bias |
| CRISPRcleanR `min_reads` | 30 (default) | Iorio 2018; lower thresholds in low-coverage screens |
| Chronos gene-effect threshold for "essential" | <-1 (cancer line) | DepMap convention |
| Chronos gene-probability for "essential" | >0.7 | DepMap convention |
| Post-correction Spearman ρ | abs(ρ) <0.05 | Acceptable correction quality |
| Cell-line CN profile resolution | ≥SNP-array level | Below this, CRISPRcleanR unsupervised |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Chronos errors on single-timepoint screen | Insufficient longitudinal data | Use CRISPRcleanR instead |
| CRISPRcleanR removes a known essential | Segment-based over-correction | Manually inspect segments; cross-check with non-corrected |
| Spearman ρ still -0.15 after correction | Method too coarse for the amp | Refine CN profile; use Chronos |
| ERBB2 listed as essential in SK-BR-3 | Uncorrected HER2 amplification | Always apply correction before hit calling |
| CN profile missing for newly characterized line | Profile not generated | Run SNP-array / low-coverage WGS |
| Hits restricted to non-amplified regions only | Over-correction | Reduce CRISPRcleanR aggressiveness; check known biology |

## References

- Aguirre AJ et al. 2016. *Cancer Discov* 6:914. Copy-number gene-independent toxicity.
- Munoz DM et al. 2016. *Cancer Discov* 6:900. CN amplification CRISPR artifacts.
- Meyers RM et al. 2017. *Nat Genet* 49:1779. CERES; first formal CN correction at DepMap scale.
- Iorio F et al. 2018. *BMC Genomics* 19:604. CRISPRcleanR.
- Dempster JM et al. 2021. *Genome Biol* 22:343. Chronos.
- Hess GT et al. 2016. *Nat Methods* 13:1036. Reduced toxicity of base editor screens (less DNA damage).
- Behan FM et al. 2019. *Nature* 568:511. Sanger Score with CRISPRcleanR-corrected data.
- Pacini C et al. 2024. *Nat Commun* 15:1230. DepMap quality scoring and benchmarking.
- DepMap Q4 2024+ data releases. https://depmap.org/portal/

## Related Skills

- crispr-screens/screen-qc - CN-LFC Spearman diagnostic; pre-correction QC
- crispr-screens/library-design - Switch to Dolcetto (CRISPRi) to bypass artifact
- crispr-screens/mageck-analysis - MAGeCK on CRISPRcleanR-corrected counts
- crispr-screens/bagel-essentiality - BAGEL2 on CRISPRcleanR-corrected counts
- crispr-screens/hit-calling - Cancer-line hit calling with Chronos
- crispr-screens/batch-correction - Chronos handles batch + CN jointly
- crispr-screens/jacks-analysis - JACKS does not handle CN bias
- clinical-databases/clinvar-lookup - Variant annotation downstream
- copy-number/copy-ratio-segmentation - CN profile derivation upstream
