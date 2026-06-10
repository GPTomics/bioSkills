---
name: bio-methylation-array-qc-filtering
description: Performs probe filtering and sample-level QC on Illumina Infinium methylation arrays (450K / EPIC / EPICv2) to decide which probes and samples to trust. Drops detection-p-failed and low-bead-count probes, removes cross-reactive/non-specific probes (Chen 2013 / Pidsley 2016 lists via maxprobes), excludes SNP-overlapping probes with dropLociWithSnps, and handles sex-chromosome probes. Collapses EPICv2 replicate probes with betasCollapseToPfx and harmonizes across array versions (EPICv2 hg38 vs 450K/EPIC hg19, intersect plus mLiftOver). Runs sample-identity QC: getSex sex prediction vs sample sheet for swap detection, rs-SNP fingerprint clustering for duplicates/swaps, and Sentrix chip/array-position batch diagnosis. Use when filtering methylation array probes, detecting sample swaps or mislabels, collapsing EPICv2 replicates, or merging 450K/EPIC/EPICv2 cohorts. For IDAT-to-corrected-beta normalization see array-preprocessing; for batch correction and study design see ewas-design.
tool_type: r
primary_tool: minfi
---

## Version Compatibility

Reference examples tested with: minfi 1.48+, sesame 1.20+, maxprobes 0.0.2+, ChAMP 2.32+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The ARRAY VERSION is the version that matters most here. Cross-reactive probe lists, SNP-overlap annotation, the manifest, and the genome build are all array-version-specific. Record whether the data is 450K (hg19), EPIC v1 (hg19), or EPIC v2 (hg38), and which annotation package supplies the probe metadata (`IlluminaHumanMethylation450kanno.ilmn12.hg19`, `...EPICanno.ilm10b4.hg19`, `...EPICv2anno.20a1.hg38`). EPICv2 carries ~5,100 replicate probes (2-10 designs per locus) and is hg38-native; both facts break naive cross-version merges if ignored.

# Array QC and Filtering

**"Which probes and samples can I trust on my methylation array?"** -> Mask the failed/cross-reactive/SNP-overlapping probes, collapse EPICv2 replicates, and check every sample's predicted sex and rs-SNP fingerprint against the sample sheet - because a raw array beta is uninterpretable until it is detection-masked and probe-filtered, and a sample swap is the failure no downstream model can rescue.
- R: `dropLociWithSnps(gset, snps=c('CpG','SBE'), maf=0)` then `getSex()` and `getSnpBeta()` for identity QC

Scope: probe filtering + sample-level QC + EPICv2 replicate collapse + cross-version harmonization for Infinium arrays. IDAT reading, background/dye correction, detection-p masking, and normalization -> array-preprocessing (it produces the matrix QC'd here). Explicit chip/position batch CORRECTION and study design -> ewas-design. Per-CpG testing on the filtered matrix -> differential-cpg-testing. Cohort cell-composition QC -> cell-type-deconvolution. Short-read bisulfite or long-read MM/ML calling are different modalities (see the bisulfite skills and long-read-sequencing/nanopore-methylation).

## The Single Most Important Modern Insight -- A Raw Array Beta Is Uninterpretable Until Detection-Masked and Probe-Filtered, and a Sample Swap Is the Failure Nothing Downstream Can Rescue

A beta value of 0.5 from a failed probe, a cross-reactive probe, or a SNP-overlapping probe looks exactly like a real intermediate methylation call - confident, reproducible, and wrong. Deciding which probes and samples to TRUST is the measurement's integrity layer, not optional cleanup. Three corollaries every misuse violates:

1. **A confident beta can be pure noise or pure genotype.** A detection-p-failed probe returns a number with no signal behind it. A cross-reactive probe sums fluorescence from multiple genomic locations. A CpG-SNP or SBE-SNP probe reports the donor's GENOTYPE, not methylation - producing reproducible-but-genetic "associations." None of these are visible in the beta value itself.
2. **The most common and most embarrassing failure is a sample swap or chip-confounded batch.** A `getSex()` prediction that disagrees with the sample sheet, or rs-SNP fingerprints that cluster two "different" samples together, reveals a mislabel that no model corrects after the fact. If Sentrix chip or array-position is confounded with the biological group, the technical and biological signals are mathematically inseparable - randomize at design (-> ewas-design), do not try to rescue it.
3. **Merging across array versions silently misaligns loci.** EPICv2 measures ~5,100 loci with 2-10 replicate probes and is annotated on hg38; 450K/EPIC are hg19 and have unique probe IDs. Collapse replicates and intersect/liftover BEFORE merging, or the same locus is counted multiple times (inflating its weight and breaking per-CpG FDR) and coordinates clash across builds.

Organize the work around delivering a trustworthy, merge-safe matrix: filter probes, collapse replicates, verify sample identity. Over-filtering is its own error - dropping every flagged probe discards real signal, so the maf cutoff and the cross-reactive list are calibrated decisions, not a fixed recipe.

## Filtering Taxonomy

| Filter | Tool / function | Citation | What it removes |
|--------|-----------------|----------|-----------------|
| Failed detection-p | `detectionP()` (minfi) / `pOOBAH` (sesame) | Aryee 2014 *Bioinformatics* 30:1363; Zhou 2018 *NAR* 46:e123 | probes with signal indistinguishable from background (per applied in array-preprocessing) |
| Low bead count | `getNBeads()` (minfi) / sesame bead data | Aryee 2014 *Bioinformatics* 30:1363 | probes built from too few beads (<3), unreliable |
| Cross-reactive / non-specific | `dropXreactiveLoci()` / `xreactive_probes()` (maxprobes) | Chen 2013 *Epigenetics* 8:203; Pidsley 2016 *Genome Biol* 17:208 | probes co-hybridizing to multiple loci; list is ARRAY-VERSION-specific |
| SNP-overlapping | `dropLociWithSnps()`, `getSnpInfo()` (minfi) | Aryee 2014 *Bioinformatics* 30:1363 | CpG-SNP / SBE-SNP probes that report genotype not methylation |
| Sex-chromosome | annotation `chr` (minfi) | Aryee 2014 *Bioinformatics* 30:1363 | chrX/chrY probes (sex-confounded; drop or sex-stratify) |
| EPICv2 replicates | `betasCollapseToPfx()` (sesame) | Kaur 2023 *Epigenetics Commun* 3:6 | extra designs per locus; collapse to one value per cg core ID |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| 450K cohort, EWAS-grade filter | maxprobes 450K list + `dropLociWithSnps` + drop chrX/Y | Chen 2013 list is 450K-specific; standard EWAS attrition |
| EPIC v1 cohort | maxprobes EPIC list + `dropLociWithSnps` + sex-chr decision | Pidsley 2016 list is EPIC-specific; do not reuse the 450K list |
| EPIC v2 cohort | sesame `betasCollapseToPfx` FIRST, then filter on hg38 anno | collapse replicates before any per-locus filter or FDR |
| Sex is the phenotype | keep chrX/chrY, analyze sex-stratified | dropping sex probes throws away the signal of interest |
| Suspected mislabels / replicates | `getSex()` vs sheet + `getSnpBeta()` rs-fingerprint clustering | swaps and duplicates are invisible in methylation alone |
| Merge 450K + EPIC + EPICv2 | intersect probe IDs after collapse; `mLiftOver` for coordinates | EPICv2 is hg38, others hg19; counts/coords clash otherwise |
| Chip/position confounded with group | -> ewas-design | unrecoverable by filtering; a design problem |
| Need the corrected matrix to filter | -> array-preprocessing | this skill QC's a matrix it does not produce |

## Probe Filtering on a GenomicRatioSet

**Goal:** Reduce a corrected `GenomicRatioSet` to the probes whose beta values reflect methylation rather than noise, genotype, or cross-hybridization.

**Approach:** Drop SNP-overlapping probes with minfi, remove the array-version-matched cross-reactive list with maxprobes, optionally drop sex-chromosome probes, and record the attrition at each step for the methods section.

```r
library(minfi)
library(maxprobes)

# gset is a corrected GenomicRatioSet from array-preprocessing (IDAT -> noob/funnorm -> ratios)
start_n <- nrow(gset)

# SNP at the CpG interrogation or single-base-extension site reports genotype, not methylation.
# maf=0 drops ANY annotated SNP (conservative EWAS default); raise maf to keep rare variants.
gset <- dropLociWithSnps(gset, snps = c('CpG', 'SBE'), maf = 0)

# Cross-reactive list is ARRAY-VERSION-specific: 'EPIC' (Pidsley 2016) vs '450K' (Chen 2013).
gset <- maxprobes::dropXreactiveLoci(gset)

# Sex-chromosome probes are sex-confounded; drop for autosomal EWAS or analyze sex-stratified.
anno <- getAnnotation(gset)
autosomal <- !(anno$chr %in% c('chrX', 'chrY'))
gset <- gset[autosomal, ]

attrition <- c(start = start_n, after_snp_xreact_sex = nrow(gset))
attrition
```

## Detection-p and Low-Bead Masking (boundary with array-preprocessing)

Per-probe detection-p and low-bead masking depend on the raw two-channel signal and control probes, which only exist at the `RGChannelSet`/`SigDF` stage handled in array-preprocessing. That skill applies `detectionP()` (minfi) or `pOOBAH` (sesame, which also catches deletion-driven false-intermediate calls) and `getNBeads()` before producing the corrected matrix. This skill assumes that masking is already done; if a supplied beta matrix has NOT been detection-masked, route back to array-preprocessing rather than trusting the betas. The thresholds (detection-p, fraction-of-samples-failed) live with the masking step, not here.

## EPICv2 Replicate Collapse

**Goal:** Reduce EPICv2's multiple probe designs per locus to one value per legacy CpG before any per-CpG analysis or cross-version merge.

**Approach:** Collapse replicate betas by probe-ID prefix with sesame, choosing mean (default) or the minimum-detection-p replicate, which also strips the design suffix so IDs revert to the classic cg form.

```r
library(sesame)

# EPICv2 IDs carry a design/replicate suffix (e.g. cg00000029_TC21); ~5,100 loci have 2-10 designs.
# Leaving replicates uncollapsed counts a locus multiple times: inflates its weight, makes
# correlated duplicate "tests" break per-CpG FDR, and corrupts any cross-version merge.
betas_collapsed <- betasCollapseToPfx(betas_epicv2)   # averages the replicate designs to one value per cg core ID

# betasCollapseToPfx only AVERAGES (it takes betas and nothing else). To keep the best-detection
# replicate instead, request collapse at the SigDF stage from the IDATs (a beta matrix has already
# discarded the per-probe detection p that minPval needs):
# betas <- openSesame(idat_prefixes, func = getBetas, collapseToPfx = TRUE, collapseMethod = 'minPval')
```

## Cross-Version Harmonization

**Goal:** Merge 450K, EPIC, and EPICv2 cohorts (or apply a 450K-trained clock/EWAS signature to EPICv2) without double-counting loci or clashing genome builds.

**Approach:** Collapse EPICv2 replicates first, intersect on the shared cg core IDs, then liftover coordinates because EPICv2 is hg38 while 450K/EPIC are hg19.

```r
library(sesame)

# 1. Collapse EPICv2 to cg core IDs (above), then intersect probe sets across versions.
shared <- Reduce(intersect, list(rownames(betas_450k), rownames(betas_epic), rownames(betas_collapsed)))

# 2. Coordinates differ by build: EPICv2 is hg38, 450K/EPIC are hg19. mLiftOver harmonizes
#    probe-level data across platforms/builds; intersect IDs first, lift coordinates before merging.
# betas_v2_hg19 <- mLiftOver(betas_collapsed, target_platform = 'HM450')

merged <- cbind(betas_450k[shared, ], betas_epic[shared, ], betas_collapsed[shared, ])
dim(merged)   # a 450K-trained clock/EWAS does not transfer to EPICv2 without this intersection
```

## Sample-Level Identity QC

**Goal:** Catch sample swaps, mislabels, and unintended duplicates before any analysis - the single most common data-integrity failure.

**Approach:** Predict sex from chrX/chrY intensity and compare to the sample sheet, then cluster samples on the rs-SNP genotyping probes (65 on 450K, ~59 on EPIC) to find duplicates and swaps independent of methylation.

```r
library(minfi)

# Sex from log2(median chrY intensity) - log2(median chrX intensity); two clusters = M/F.
# A predicted sex that disagrees with the sample sheet is the canonical sample-swap flag.
predicted <- getSex(gmset)              # gmset = mapped MethylSet/GenomicMethylSet
mismatch <- predicted$predictedSex != sample_sheet$Sex
sample_sheet$Basename[mismatch]

# rs-SNP fingerprint: ~59 explicit rs genotyping probes. Clustering on these betas (each ~0/0.5/1)
# reveals duplicate individuals and swaps regardless of methylation - genotype is identity.
snp_betas <- getSnpBeta(rgset)          # rgset = the raw RGChannelSet from array-preprocessing
identity_clusters <- hclust(dist(t(snp_betas)))
plot(identity_clusters)                 # technical replicates of one person cluster tightly
```

## Chip / Array-Position Batch Diagnosis

Sentrix chip (BeadChip barcode) and array position (`Sentrix_Position`, the row/column on the chip) are the dominant technical axes in Infinium data. This skill DIAGNOSES whether they associate with top variance components; it does NOT correct them. ChAMP's `champ.SVD()` regresses the leading singular vectors of the beta matrix against chip, position, plate, and the biological factors, flagging which technical axis loads on real variance. If chip or position is confounded with the biological group, it is mathematically unrecoverable - hand the explicit correction (ComBat/SVA, or chip/position as covariates/random effects) and the design fix to ewas-design.

## Per-Method Failure Modes

### SNP-overlapping probes left in
**Trigger:** running a per-CpG test without `dropLociWithSnps`. **Mechanism:** a SNP at the CpG or SBE site makes the probe report genotype, not methylation. **Symptom:** reproducible "associations" that are actually genetic (often mQTL-driven, trimodal beta). **Fix:** `dropLociWithSnps(snps=c('CpG','SBE'), maf=0)`; raise maf only to deliberately keep rare variants.

### Wrong cross-reactive list for the array
**Trigger:** applying the Chen 2013 450K list to EPIC/EPICv2 data (or vice versa). **Mechanism:** the cross-reactive probe set is array-version-specific. **Symptom:** wrong probes dropped, real cross-reactive probes retained. **Fix:** use the array-matched list (`xreactive_probes(array_type='EPIC')` vs `'450K'`); maxprobes maps EPICv2 via the collapsed EPIC core IDs.

### EPICv2 replicates not collapsed
**Trigger:** treating EPICv2 betas as if probe IDs were unique. **Mechanism:** ~5,100 loci have 2-10 designs; the same locus appears multiple times. **Symptom:** duplicated rownames, inflated locus weight, broken per-CpG FDR, corrupted cross-version merge. **Fix:** `betasCollapseToPfx()` first; strip the suffix back to the cg core ID before anything downstream.

### Build mismatch on merge
**Trigger:** merging EPICv2 (hg38) coordinates with 450K/EPIC (hg19). **Mechanism:** EPICv2 annotation is hg38-native. **Symptom:** loci silently misaligned by the hg19/hg38 offset. **Fix:** intersect on cg IDs and `mLiftOver` (or restrict to shared IDs and track the build per version).

### Sample swap not checked
**Trigger:** analyzing without the sex/identity QC. **Mechanism:** a mislabeled IDAT carries the wrong phenotype. **Symptom:** weakened or spurious associations; `getSex()` disagrees with the sheet; rs-fingerprints cluster two "different" samples. **Fix:** run `getSex()` vs sample sheet and `getSnpBeta()` fingerprint clustering as mandatory pre-analysis QC.

### Over-filtering
**Trigger:** dropping every flagged probe reflexively. **Mechanism:** some "cross-reactive" probes are fine for the specific locus of interest; maf=0 removes any-SNP probes including innocuous ones. **Symptom:** real signal discarded; clock/signature CpGs lost. **Fix:** treat the maf cutoff and cross-reactive list as calibrated to the question; report attrition and check that target CpGs survive.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| detection-p > 0.01 = failed | Aryee 2014 *Bioinformatics* 30:1363 | signal indistinguishable from background; applied in array-preprocessing |
| bead count < 3 = unreliable | minfi docs | too few beads per probe to trust the intensity |
| `dropLociWithSnps(maf=0)` | minfi docs | maf=0 drops any annotated SNP; raise to keep rare variants (calibrated) |
| ~6% of 450K probes cross-reactive | Chen 2013 *Epigenetics* 8:203 | ~29-39K loci co-hybridize; array-version-specific list |
| EPICv2 ~5,100 replicate loci (2-10 designs) | Kaur 2023 *Epigenetics Commun* 3:6 | collapse to one cg core ID before per-locus FDR |
| 65 rs-SNP probes on 450K (~59 on EPIC) | minfi annotation | enough genotype to fingerprint identity and catch swaps |
| getSex on log2 medY - log2 medX | Aryee 2014 *Bioinformatics* 30:1363 | X/Y intensity clusters by sex; mismatch = swap flag |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Duplicated rownames in EPICv2 beta matrix | replicate probes not collapsed | `betasCollapseToPfx()` before merge/test |
| Reproducible genetic-looking hits | SNP-overlap probes retained | `dropLociWithSnps(snps=c('CpG','SBE'), maf=0)` |
| Coordinates off when merging cohorts | EPICv2 hg38 vs 450K/EPIC hg19 | intersect cg IDs; `mLiftOver` before merge |
| `getSex` disagrees with sample sheet | sample swap/mislabel | trace the IDAT; rs-SNP fingerprint to confirm |
| Cross-reactive filter drops too few/many | wrong array_type list | match the list to the array version |
| `dropXreactiveLoci` errors on EPICv2 object | maxprobes keys on EPIC core IDs | collapse EPICv2 to cg core IDs first |

## References

- Aryee MJ, Jaffe AE, Corrada-Bravo H, et al. 2014. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. *Bioinformatics* 30:1363-1369.
- Zhou W, Triche TJ Jr, Laird PW, Shen H. 2018. SeSAMe: reducing artifactual detection of DNA methylation by Infinium BeadChips in genomic deletions. *Nucleic Acids Res* 46:e123.
- Chen YA, Lemire M, Choufani S, et al. 2013. Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray. *Epigenetics* 8:203-209.
- Pidsley R, Zotenko E, Peters TJ, et al. 2016. Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling. *Genome Biol* 17:208.
- Kaur D, Lee SM, Goldberg D, et al. 2023. Comprehensive evaluation of the Infinium human MethylationEPIC v2 BeadChip. *Epigenetics Commun* 3:6.

## Related Skills

- array-preprocessing - Produces the corrected beta/M matrix being QC'd and filtered
- ewas-design - Chip/position batch correction and study design
- cell-type-deconvolution - Cohort composition QC
- differential-cpg-testing - Downstream per-CpG testing on the filtered matrix
- workflows/methylation-pipeline - End-to-end pipeline
