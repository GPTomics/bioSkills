---
name: bio-hi-c-analysis-hic-differential
description: Compares Hi-C contact maps between conditions across the right scale -- differential bin-pair contacts (multiHiCcompare, diffHic), differential A/B compartments (dcHiC), differential TAD boundaries (delta insulation), and differential loops (diffloop, DiffHiChIP) -- with distance-stratified between-sample normalization, replicate-aware NB-GLM FDR, HiCRep SCC reproducibility gating, and CNV correction for cancer/aneuploid samples. Use when comparing Hi-C between treatment and control, finding differential contacts/compartments/boundaries/loops, normalizing two maps of unequal depth, choosing a replicate-aware test, gating replicates with SCC, or correcting copy-number artifacts before a tumor-vs-normal comparison.
tool_type: python
primary_tool: cooltools
---

## Version Compatibility

Reference examples tested with: cooler 0.10+, cooltools 0.7+, bioframe 0.7+, multiHiCcompare 1.20+, diffHic 1.34+, dcHiC (2022 release), hicrep (Bioconductor) 1.26+, edgeR 4.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: cooltools provides the Python feature-extraction parts (expected, eigenvectors, insulation, pileups) but NO turnkey two-condition test -- the differential statistics live in R/Bioconductor (multiHiCcompare, diffHic, dcHiC, hicrep). The `get.scc` signature changed between the TaoYang-dev and current Bioconductor/qunhualilab releases; verify with `?get.scc` before calling. A `.mcool` is multi-resolution: pass a single-resolution URI (`file.mcool::/resolutions/10000`), not the bare file.

# Hi-C Differential Analysis

**"What changed in 3D genome organization between my conditions?"** -> Equalize the two maps with a distance-stratified between-sample normalization, then test at the SCALE of the question (compartment, TAD boundary, loop, or bin-pair) with a replicate-aware method, not pixel-wise log2 subtraction.
- Python (features): `cooltools.expected_cis(clr)`, `cooltools.eigs_cis(clr)`, `cooltools.insulation(clr)`
- R (bin-pair test): `make_hicexp(...) |> cyclic_loess() |> hic_exactTest() |> results()` (multiHiCcompare)
- R (compartments): `Rscript dchicf.r --pcatype cis|select|analyze` (dcHiC)

## The Single Most Important Modern Insight -- balancing makes a matrix self-consistent, NOT cross-comparable

ICE/KR/SCALE balancing equalizes the marginals WITHIN one map (it removes per-bin visibility bias). It says nothing about whether map A and map B are on the same footing. Two balanced maps still differ in (a) total sequencing depth and (b) cis/trans ratio, and a naive `log2(A/B)` is dominated by those two nuisances plus the shared distance-decay P(s) -- with biology buried underneath. "I balanced both, so I can subtract them" is the single most common error in differential Hi-C. The fix is a BETWEEN-sample, DISTANCE-STRATIFIED normalization (multiHiCcompare's cyclic loess on the M-D plot, or diffHic's trended loess offsets) before any difference is interpretable.

The M-D plot is the RNA-seq MA-plot's distance-aware cousin: M = log2(IF1/IF2), but plotted against genomic DISTANCE D instead of mean abundance. The loess fit is done PER distance stratum because both bias and variance depend on distance -- a sparser library loses long-range pairs faster, so the depth bias is itself distance-dependent and a single global size-factor cannot fix it. After normalization M should center on 0 at every D; a residual M-trend at large D means normalization failed at long range. Inspect `MD_hicexp()` -- do not trust the call set blind.

## Differential-Method Taxonomy (scale-matched -- one tool cannot do all four)

| Method | Scale / object | Mechanism | Replicates | When |
|--------|----------------|-----------|------------|------|
| multiHiCcompare | bin-pair contacts | cyclic-loess M-D normalization + edgeR exactTest/GLM | >=2 for FDR | multi-group, joint normalization, covariates |
| diffHic | bin-pair contacts | squareCounts -> trended/CNV loess offsets -> edgeR NB-GLM | >=2 (required) | replicate-rich, CNV correction, full edgeR machinery |
| dcHiC | A/B compartment (Mb) | sign/PC-consistent eigenvectors -> quantile-norm -> Mahalanobis | works at 1, better with reps | quantitative compartment SHIFTS across samples |
| delta insulation | TAD boundary (sub-Mb) | difference of Crane-style insulation scores per condition | reps for significance | boundary strengthening/loss, not bin-pair counts |
| diffloop / DiffHiChIP | loops (anchored) | edgeR NB on anchored counts; IHW by distance (DiffHiChIP) | YES for FDR | Hi-C/HiChIP/ChIA-PET loop sets, long-range power |
| Selfish | regions (n=1) | Gaussian-pyramid self-similarity, distance-controlled | none (descriptive) | replicate-poor 2-map ranking, no honest FDR |
| FIND | bin-pair (n=1) | spatial Poisson process over the 2D neighborhood | none (descriptive) | neighborhood-aware 2-map ranking |
| HiCRep SCC | whole-matrix similarity | stratum-adjusted correlation (NOT a per-feature test) | n/a (pairwise) | replicate QC gate + coarse condition distance |

Using a bin-pair tool (multiHiCcompare/diffHic) to "find compartment changes" is a category error: compartments are an eigenvector property of the whole chromosome, not a sum of independent bin-pairs, so the result is a noisy bin-pair list, not coherent compartment calls. Match the tool to the scale of the question.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| First question for any comparison | HiCRep SCC (within- vs between-condition) | gate reproducibility before any call |
| Differential contacts, >=2 reps/condition | multiHiCcompare (cyclic_loess + edgeR) or diffHic | distance-stratified norm + replicate-aware FDR |
| Differential contacts, n=1 vs n=1 | Selfish or FIND, DESCRIPTIVE only | no replicates -> no honest FDR |
| Differential A/B compartments | dcHiC (cis -> select -> analyze) | sign-consistent eigenvectors + Mahalanobis shifts |
| Differential TAD boundaries | delta insulation (cooltools insulation per condition) -> stats | boundaries are local insulation depth, not counts |
| Differential loops (Hi-C) | diffloop or DiffHiChIP on per-condition loop calls | anchored count test with distance modeling |
| Differential loops (HiChIP / PLAC) | -> hichip-plac-loops then DiffHiChIP | peak-biased data needs anchored null |
| Tumor vs normal / aneuploid | CNV-correct FIRST (diffHic normalizeCNV / OneD) | CNV masquerades as differential contacts |
| Two maps, unequal depth, quick look | downsample to equal valid pairs, then M-D normalize | depth fix only; still distance-stratify |
| Call the per-condition features in Python | -> compartment-analysis, tad-detection, loop-calling | cooltools extracts; bring to R for the test |
| Annotate the differential anchors/boundaries | -> chip-seq/peak-annotation, genome-intervals/overlap-significance | overlap with TF peaks / enrichment test |

## Replicate QC Gate with HiCRep SCC (do this FIRST)

**Goal:** Decide whether replicates are reproducible enough that a differential call is meaningful at all.

**Approach:** Plain Pearson on Hi-C always looks reproducible -- the shared P(s) decay alone drives r > 0.9 even between unrelated maps. HiCRep's SCC stratifies by distance (removing the decay) and smooths for sparsity, then variance-weights the strata. Compute SCC for within-condition replicate pairs and between-condition pairs; within must clearly exceed between or there is nothing to call.

```r
library(hicrep)

# Current Bioconductor/qunhualilab interface: dat is a 4-column table
# (mid1, mid2, IF_A, IF_B); resol = bin size; max = max distance considered.
# Verify with ?get.scc -- the older TaoYang-dev interface is get.scc(mat1, mat2, resol, h, lbr, ubr).
scc_out <- get.scc(dat_repA_vs_repB, resol = 50000, max = 5000000)
scc_out$scc   # stratum-adjusted correlation coefficient in [-1, 1]
```

A differential claim is only meaningful when within-condition SCC clearly exceeds between-condition SCC. If they overlap, the "differential" signal is replicate noise.

## Differential Bin-Pair Contacts with multiHiCcompare

**Goal:** Find individual bin-pairs whose contact frequency changes between conditions, with a calibrated FDR from replicate variance.

**Approach:** Build a Hi-C experiment from per-replicate sparse upper-triangular tables, normalize jointly across all samples with cyclic loess on the M-D plot, then run edgeR's exact test (2 groups) or GLM (covariates). Filtering on mean abundance happens at `make_hicexp` time -- it is independent of the contrast, so it shrinks the multiple-testing burden without inflating FDR.

```r
library(multiHiCcompare)

# Each replicate is a 4-column sparse table: chr, region1(bp), region2(bp), IF
# chr coded 1-22, 23=X, 24=Y.
hicexp <- make_hicexp(c1_r1, c1_r2, c2_r1, c2_r2,
                      groups = c(0, 0, 1, 1),
                      zero.p = 0.8,                 # drop bin-pairs >80% zero across samples
                      A.min  = 5,                   # drop bin-pairs with low mean IF (independent filter)
                      filter = TRUE,
                      remove.regions = hg19_cyto)   # blacklist centromeres/telomeres
hicexp <- cyclic_loess(hicexp, span = NA)           # span=NA -> GCV chooses the loess span
hicexp <- hic_exactTest(hicexp)                     # 2-group; use hic_glm(hicexp, design) for covariates
res <- results(hicexp)                              # chr, region1, region2, D, logFC, logCPM, p.value, p.adj
sig <- topDirs(hicexp, logfc_cutoff = 1, logcpm_cutoff = 1, p.adj_cutoff = 0.1, return_df = 'pairedbed')
MD_hicexp(hicexp)                                   # diagnostic: M should center on 0 at every D
```

diffHic is the alternative when full edgeR control is wanted: `squareCounts` -> `filterDirect` -> `normOffsets(type='loess')` -> `asDGEList` -> `estimateDisp` -> `glmQLFit` -> `glmQLFTest`. Same NB-GLM engine as edgeR/csaw; requires biological replicates for dispersion.

## Differential A/B Compartments with dcHiC

**Goal:** Detect compartment changes between conditions, including graded shifts that never cross the A/B boundary, with cross-sample-comparable eigenvectors.

**Approach:** Everyone computes PC1, but PC1's sign is arbitrary per chromosome per sample and sometimes the compartment signal is in PC2 -- naive multi-sample comparison silently compares flipped or mismatched axes. dcHiC anchors the sign (GC/gene-density) and selects the correct PC, quantile-normalizes the scores, then uses a multivariate Mahalanobis distance per bin to flag outliers across all samples. Run as a staged CLI; the input file lists `<matrix> <bed> <replicate_prefix> <experiment_prefix>` per replicate (no dashes/dots in prefixes).

```bash
# Staged CLI (dchicf.r). cis = per-sample compartments; select = pick PC per chr;
# analyze = differential PCA (Mahalanobis); viz = IGV-style browser.
Rscript dchicf.r --file input.txt --pcatype cis    --dirovwt T --cthread 2 --pthread 4
Rscript dchicf.r --file input.txt --pcatype select --dirovwt T --genome hg38
Rscript dchicf.r --file input.txt --pcatype analyze --dirovwt T --diffdir cond1_vs_cond2
Rscript dchicf.r --file input.txt --pcatype viz    --diffdir cond1_vs_cond2 --genome hg38
# Differential calls land in DifferentialResult/<diffdir>/fdr_result/ ; optional: --pcatype subcomp (HMM) and dloop.
```

For a quick Python eyeball of compartment switching (NOT a replicate-aware test -- use dcHiC for that), difference the phased E1 from cooltools per condition and flag sign flips. State it as exploratory.

## Differential TAD Boundaries via Delta Insulation

**Goal:** Find boundaries that strengthen, weaken, or appear/disappear between conditions.

**Approach:** Boundary changes are about local insulation DEPTH, not bin-pair counts, so compute the Crane-style insulation score per condition with cooltools at a fixed window, then difference the scores; attach significance with replicate insulation profiles. The window must match the bin size (typically 5-25x the bin).

```python
import cooltools

ins1 = cooltools.insulation(clr1, window_bp=[200000], ignore_diags=2)   # window ~5-25x bin size
ins2 = cooltools.insulation(clr2, window_bp=[200000], ignore_diags=2)
merged = ins1.merge(ins2, on=['chrom', 'start', 'end'], suffixes=('_1', '_2'))
merged['delta_insulation'] = merged['log2_insulation_score_200000_2'] - merged['log2_insulation_score_200000_1']
```

GENOVA and FAN-C also compute the insulation score for differencing; significance is added separately (paired test across replicate insulation tracks).

## CNV Correction Before Cancer Comparisons

**Goal:** Avoid calling copy-number differences as 3D-structure differences in tumor-vs-normal data.

**Approach:** Contact count scales with copy number, and balancing assumes equal visibility -- so ICE on aneuploid data SHIFTS contacts between amplified and deleted regions instead of removing the artifact, lighting up enormous spurious "differential interaction" blocks. Either regress out the marginal (1D) coverage log-ratio as a covariate (diffHic `marginCounts` + `normalizeCNV`, or OneD's 1D GAM), OR call CNV first (HiNT/HiCnv) and interpret rearrangement blocks separately. Never run vanilla balanced-matrix differential on cancer data.

```r
library(diffHic)

margins <- marginCounts(data)              # 1D marginal coverage per bin (a RangedSummarizedExperiment)
nb.off  <- normalizeCNV(data, margins)     # 2D loess on (abundance, marginal log-ratio) -> GLM offsets; matches margins internally
y <- asDGEList(data)
y$offset <- nb.off                         # attach the CNV offsets; asDGEList does not carry them automatically
# then estimateDisp -> glmQLFit -> glmQLFTest as usual
```

## Per-Method Failure Modes

### Naive log2 of two balanced maps
**Trigger:** `log2((mat2+1)/(mat1+1))` on two balanced coolers of different depth. **Mechanism:** balancing is within-map only; depth + cis/trans + P(s) differences dominate. **Symptom:** a smooth distance-dependent gradient in the "differential" map, strongest off-diagonal. **Fix:** distance-stratified between-sample normalization (multiHiCcompare cyclic_loess / diffHic loess offsets) before differencing.

### Pooling distances into one FDR
**Trigger:** one BH correction over all bin-pairs regardless of distance. **Mechanism:** counts and variance span 2-3 orders of magnitude across distance; short-range is high-count/low-variance, long-range is sparse. **Symptom:** almost all hits are short-range; >500 kb changes vanish. **Fix:** distance-stratified testing; IHW weighting by distance (DiffHiChIP) recovers long-range loops.

### No-replicate FDR
**Trigger:** n=1 vs n=1 reported with a p-value/FDR. **Mechanism:** dispersion needs within-condition variability; with n=1 there is none, so any FDR is fabricated. **Symptom:** thousands of "significant" hits that do not replicate. **Fix:** n>=2 (ideally 3) + diffHic/multiHiCcompare; for n=1 use Selfish/FIND DESCRIPTIVELY only.

### Plain Pearson "reproducibility"
**Trigger:** reporting raw matrix Pearson/Spearman as a QC number. **Mechanism:** shared P(s) decay inflates r > 0.9 even between unrelated maps. **Symptom:** everything looks reproducible, including failed libraries. **Fix:** HiCRep SCC (stratum-adjusted); compare within- vs between-condition.

### CNV read as structure
**Trigger:** tumor-vs-normal on vanilla balanced matrices. **Mechanism:** count scales with copy number; balancing worsens it by shifting contacts. **Symptom:** huge block-shaped "differential" regions aligned to known CNVs. **Fix:** CNV-correct (diffHic normalizeCNV / OneD) or call CNV first and interpret separately.

### Wrong scale for the question
**Trigger:** a bin-pair tool used to find "compartment changes". **Mechanism:** compartments are a whole-chromosome eigenvector property, not a sum of bin-pairs. **Symptom:** a scattered bin-pair list with no coherent A/B structure. **Fix:** dcHiC for compartments; delta insulation for boundaries; diffloop for loops.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Replicates n>=2 (ideally 3) per condition | NB dispersion estimation | within-condition variance is required for an honest FDR |
| Within-condition SCC > between-condition SCC | replicate QC gate | if they overlap, "differential" signal is replicate noise |
| `zero.p = 0.8` (drop >80% zero) | multiHiCcompare default | sparse long-range pairs break NB and waste FDR budget |
| `A.min = 5` mean-IF filter | multiHiCcompare independent filter | filter independent of contrast preserves FDR validity |
| Bin-pair / loop FDR <= 0.1 | genome-wide multiple testing | millions of bin-pairs need FDR control, not raw p |
| Long-range threshold ~500 kb | DiffHiChIP 2025 benchmark | distance-aware IHW recovers >500 kb loops flat tests miss |
| Insulation window 5-25x bin size | insulation-score scale | window much smaller than this is noisy; much larger blurs boundaries |
| Compartment resolution 100kb-1Mb | compartment scale | A/B is chromosome-scale; finer bins mix in TAD/loop structure |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Differential map is a smooth distance gradient | naive log2 of balanced maps, no between-sample norm | cyclic_loess / diffHic loess offsets first |
| All hits short-range, no long-range | distances pooled into one FDR | distance-stratified test; IHW by distance |
| Thousands of "significant" n=1 hits | no replicates -> no real dispersion | n>=2; Selfish/FIND descriptive for n=1 |
| Block-shaped diff regions over known CNVs | CNV not corrected on cancer data | diffHic normalizeCNV / OneD before testing |
| `get.scc` argument error | TaoYang-dev vs Bioconductor signature skew | `?get.scc`; 4-col `dat` + `resol` + `max` (current) |
| dcHiC fails on prefixes | dashes/dots in replicate/experiment names | use underscores only in prefix columns |
| Compartment call is sign-scrambled | comparing raw PC1 across samples | dcHiC anchors sign + selects the right PC |
| Empty cooltools result | chrom naming mismatch (`chr1` vs `1`) | harmonize names across cooler, fasta, tracks |

## References

- Stansfield JC, Cresswell KG, Dozmorov MG. 2019. multiHiCcompare: joint normalization and comparative analysis of complex Hi-C experiments. *Bioinformatics* 35(17):2916-2923. doi:10.1093/bioinformatics/bty950
- Stansfield JC, Cresswell KG, Vladimirov VI, Dozmorov MG. 2018. HiCcompare: an R-package for joint normalization and comparison of HI-C datasets. *BMC Bioinformatics* 19:279. doi:10.1186/s12859-018-2288-x
- Lun ATL, Smyth GK. 2015. diffHic: a Bioconductor package to detect differential genomic interactions in Hi-C data. *BMC Bioinformatics* 16:258. doi:10.1186/s12859-015-0683-0
- Chakraborty A, Wang JG, Ay F. 2022. dcHiC detects differential compartments across multiple Hi-C datasets. *Nat Commun* 13:6827. doi:10.1038/s41467-022-34626-6
- Roayaei Ardakany A, Ay F, Lonardi S. 2019. Selfish: discovery of differential chromatin interactions via a self-similarity measure. *Bioinformatics* 35(14):i145-i153. doi:10.1093/bioinformatics/btz362
- Djekidel MN, Chen Y, Zhang MQ. 2018. FIND: difFerential chromatin INteractions Detection using a spatial Poisson process. *Genome Res* 28(3):412-422. doi:10.1101/gr.212266.116
- Lareau CA, Aryee MJ. 2018. diffloop: a computational framework for identifying and analyzing differential DNA loops from sequencing data. *Bioinformatics* 34(4):672-674. doi:10.1093/bioinformatics/btx623
- Bhattacharyya S, Salgado Figueroa D, Georgopoulos K, Ay F. 2025. DiffHiChIP: identifying differential chromatin contacts from HiChIP data. *Cell Rep Methods* 5(11):101214. doi:10.1016/j.crmeth.2025.101214
- Yang T, Zhang F, Yardimci GG, Song F, Hardison RC, Noble WS, Yue F, Li Q. 2017. HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient. *Genome Res* 27(11):1939-1949. doi:10.1101/gr.220640.117
- Vidal E, le Dily F, Quilez J, Stadhouders R, Cuartero Y, Graf T, Marti-Renom MA, Beato M, Filion GJ. 2018. OneD: increasing reproducibility of Hi-C samples with abnormal karyotypes. *Nucleic Acids Res* 46(8):e49. doi:10.1093/nar/gky064
- Imakaev M, Fudenberg G, McCord RP, Naumova N, Goloborodko A, Lajoie BR, Dekker J, Mirny LA. 2012. Iterative correction of Hi-C data reveals hallmarks of chromosome organization. *Nat Methods* 9:999-1003. doi:10.1038/nmeth.2148
- Open2C, Abdennur N, et al. 2024. Cooltools: enabling high-resolution Hi-C analysis in Python. *PLoS Comput Biol* 20(5):e1012067. doi:10.1371/journal.pcbi.1012067

## Related Skills

- compartment-analysis - Per-condition A/B eigenvectors that dcHiC differences
- tad-detection - Per-condition insulation scores for delta-insulation boundary tests
- loop-calling - Per-condition loop calls fed to diffloop/DiffHiChIP
- matrix-operations - Balancing and expected/O-E that precede any comparison
- hichip-plac-loops - Peak-anchored loop calls and DiffHiChIP for HiChIP/PLAC-seq
- hic-data-io - Load and convert the cooler files this skill compares
- hic-visualization - Render differential maps and split-view comparisons
- chip-seq/peak-annotation - Annotate differential anchors/boundaries with TF peaks
- genome-intervals/overlap-significance - Permutation test for differential-feature enrichment
- differential-expression/de-results - The edgeR/FDR mental model reused here (MA-plot, independent filtering, IHW)
