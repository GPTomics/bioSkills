---
name: bio-copy-number-hrd-scoring
description: Quantify homologous recombination deficiency (HRD) from tumor copy number using the three genomic-scar metrics — loss of heterozygosity (LOH), large-scale state transitions (LST), and telomeric allelic imbalance (TAI) — with scarHRD, and via the whole-genome HRDetect and CHORD models. Covers the genomic instability score, the PARP-inhibitor clinical context, whole-genome-doubling correction, and the scar-versus-state distinction. Use when computing an HRD score for PARP-inhibitor eligibility, deriving LOH/LST/TAI scars from allele-specific copy number, deciding between scar-based and mutational-signature HRD methods, or interpreting an HRD result in a BRCA-reverted or low-purity tumor.
tool_type: mixed
primary_tool: scarHRD
---

## Version Compatibility

Reference examples tested with: R 4.3+ with scarHRD 0.1.1+, sequenza 3.0+ (allele-specific input); HRDetect / CHORD as their respective R packages where whole-genome data is available.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('scarHRD')` then `?scar_score` to confirm arguments
- scarHRD is GitHub-only (`sztup/scarHRD`); install with `remotes::install_github`

scarHRD consumes allele-specific copy number — a Sequenza `.seqz` file or an ASCAT/allele-specific segment table. It cannot run on relative log2 copy ratio.

# HRD Scoring

**"Is this tumor homologous-recombination deficient"** -> HRD leaves characteristic copy-number scars. Three are quantified and summed into an HRD score: loss of heterozygosity (LOH), large-scale state transitions (LST), and telomeric allelic imbalance (TAI). A high score predicts response to platinum chemotherapy and PARP inhibitors. The scar score is a *consequence* of past HR deficiency — which is both its strength (it integrates over tumor history) and its key limitation.

- R: `scarHRD` — the three genomic scars and their sum
- Whole-genome: `HRDetect` (weighted multi-signature model), `CHORD` (random forest)
- Input: allele-specific copy number from Sequenza or ASCAT (see allele-specific-copy-number)

## The Three Genomic Scars

| Scar | Definition | Captures |
|------|------------|----------|
| HRD-LOH | Number of LOH segments > 15 Mb but shorter than a whole chromosome | Large interstitial allelic loss |
| LST | Chromosomal breaks between adjacent segments each >= 10 Mb, separated by < 3 Mb | Large-scale rearrangement burden |
| TAI | Number of subtelomeric regions with allelic imbalance not crossing the centromere | Telomere-bounded allelic imbalance |

The HRD score is the sum of the three (the "genomic instability score", GIS). Each component has a precise size rule — these thresholds (15 Mb, 10 Mb, 3 Mb) are not arbitrary; they were selected to correlate with BRCA1/BRCA2/RAD51C deficiency (Abkevich 2012, Popova 2012, Birkbak 2012).

## Method Selection

| Method | Input | Strength | Fails when |
|--------|-------|----------|------------|
| scarHRD (LOH+LST+TAI) | Allele-specific CN (panel/WES/WGS) | Works on panels; the clinical-assay basis | Low purity; LST not WGD-corrected; relative CN input |
| HRDetect | Whole-genome (SNV sig 3, SV signatures, HRD index, indel microhomology) | Most accurate; integrates substitution + rearrangement signatures | Needs WGS; not applicable to panels/WES |
| CHORD | Whole-genome somatic mutation contexts | Distinguishes BRCA1- vs BRCA2-type deficiency | Needs WGS; somatic calls required |

Decision: for a targeted panel or WES the genomic-scar score (scarHRD-style) is the only option and is the basis of approved companion diagnostics; for whole-genome data, HRDetect or CHORD are more accurate because they add mutational-signature evidence.

## Computing Genomic Scars with scarHRD

**Goal:** Compute LOH, LST, TAI, and the HRD sum from allele-specific copy number.

**Approach:** Run scarHRD on a Sequenza `.seqz` file (or an allele-specific segment table); supply the genome build and ploidy so LST is correctly normalized.

```r
library(scarHRD)

# From a Sequenza .seqz file (allele-specific copy number, with BAF).
hrd <- scar_score('sample.small.seqz.gz',
                  reference = 'grch38',
                  seqz = TRUE)
# hrd is a one-row data frame with columns 'HRD' (LOH), 'Telomeric AI', 'LST', 'HRD-sum'.

# From a pre-computed allele-specific segment table (ASCAT-style: SampleID, Chromosome,
# Start_position, End_position, total_cn, A_cn, B_cn, ploidy):
hrd_seg <- scar_score('sample_allele_specific.txt',
                      reference = 'grch38', seqz = FALSE)
print(hrd_seg)
```

## The Postdoc-Level Caveats

Three points separate a correct HRD interpretation from a naive one:

1. **HRD is a scar, not a current state.** The score reflects HR deficiency that *occurred* during tumor evolution. A tumor that has acquired a BRCA reversion mutation — a real platinum/PARP-inhibitor resistance mechanism — still carries the scars and still scores HRD-high. A high score is not a guarantee of current HR deficiency or of drug response.
2. **LST is ploidy-dependent.** Whole-genome doubling adds breakpoints and inflates the LST count independently of HR status. The score must be computed with the correct ploidy so LST is normalized; an uncorrected WGD tumor can score falsely high.
3. **The score needs allele-specific input.** LOH and TAI are allelic-imbalance metrics — they cannot be derived from total copy number or relative log2. Garbage allele-specific input (low purity, sparse hets) gives a garbage score.

## Failure Modes

### Relative copy number used as input

**Trigger:** Feeding log2 copy ratio or total-CN segments to a scar calculator.

**Mechanism:** LOH and TAI require the minor allele copy number; relative or total CN has no allelic information.

**Symptom:** LOH and TAI near zero regardless of true HRD; nonsensical score.

**Fix:** Use allele-specific copy number from Sequenza or ASCAT (allele-specific-copy-number). The `.seqz` file or an A/B-allele segment table is the correct input.

### LST inflated by uncorrected whole-genome doubling

**Trigger:** Running the scar score without supplying the tumor's ploidy, on a WGD tumor.

**Mechanism:** WGD multiplies segments and breakpoints; LST counts breaks and rises with ploidy independent of HR deficiency.

**Symptom:** A WGD tumor with no BRCA/HR pathway lesion scores HRD-high, driven by LST.

**Fix:** Compute the score with the correct ploidy so LST is normalized. Cross-check a high LST-driven score against HR-pathway gene status and against mutational signature 3.

### Treating a high score as proof of drug response

**Trigger:** Equating HRD-high with current HR deficiency and predicted PARP-inhibitor benefit.

**Mechanism:** The scar persists after HR function is restored (BRCA reversion, other resistance mechanisms); the score integrates over history.

**Symptom:** An HRD-high tumor fails to respond; the score was correct but the tumor is no longer HR-deficient.

**Fix:** Interpret the score as evidence of past HRD. Where possible, integrate current HR-pathway status (BRCA1/2 reversion screening, RAD51 foci assays) before predicting response.

### Low tumor purity

**Trigger:** Computing HRD on a low-purity sample (< ~30-40%).

**Mechanism:** Allele-specific calling fails at low purity (see allele-specific-copy-number); scar counts then derive from an unreliable profile.

**Symptom:** Score unstable across reruns; LOH/TAI near zero on a genome with visible imbalance.

**Fix:** Confirm purity is adequate before scoring; report indeterminate below ~30%.

### Panel HRD score read as a whole-genome score

**Trigger:** Comparing a targeted-panel HRD score directly to a WGS-derived score or to a companion-diagnostic cutoff.

**Mechanism:** Genomic coverage and segment resolution differ; scar counts are not numerically interchangeable across assays.

**Symptom:** A panel score compared to the GIS >= 42 cutoff gives the wrong call.

**Fix:** Use the cutoff validated for the specific assay. Companion-diagnostic thresholds (e.g. Myriad myChoice GIS >= 42) are validated for that assay's design, not portable.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| scarHRD high, HRDetect low | LST-driven score from WGD, not true HRD | Check ploidy correction and signature 3 |
| HRD-high tumor, BRCA wild-type | Other HR lesion, or false-high from WGD/quality | Check RAD51C/PALB2, methylation; verify input |
| HRD-high tumor fails PARP-inhibitor | Scar persists after BRCA reversion | Screen for reversion mutations |
| Panel and WGS scores disagree | Different assay resolution | Use the assay-validated cutoff for each |

**Operational rule:** An HRD score is interpretable only when (1) the input is allele-specific copy number from an adequately pure sample, (2) LST is computed with the correct ploidy, (3) the assay-validated cutoff is used, and (4) the score is read as evidence of *past* HR deficiency, integrated with current HR-pathway status before predicting therapy response.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| HRD-LOH segment size | > 15 Mb, < whole chromosome | Abkevich 2012; correlates with BRCA1/2/RAD51C deficiency |
| LST adjacent-segment size | each >= 10 Mb, gap < 3 Mb | Popova 2012 |
| TAI | subtelomeric allelic imbalance not crossing the centromere | Birkbak 2012 |
| Genomic instability score (GIS) cutoff | >= 42 (Myriad myChoice) | Telli 2016; assay-specific, not portable |
| Purity floor for scoring | ~30-40% | Below this, allele-specific input is unreliable |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| LOH/TAI ~0 on an imbalanced genome | Relative/total CN used as input | Use allele-specific CN (Sequenza/ASCAT) |
| BRCA-wild-type tumor scores HRD-high | LST inflated by uncorrected WGD | Supply correct ploidy; check signature 3 |
| HRD-high tumor does not respond | Scar persists after BRCA reversion | Screen for reversion; assay current HR status |
| Score unstable across reruns | Low purity | Confirm purity; report indeterminate if low |
| Panel score fails the GIS >= 42 call | Cross-assay cutoff misuse | Use the assay-validated threshold |
| scarHRD install fails | GitHub-only package | `remotes::install_github('sztup/scarHRD')` |

## References

- Abkevich V et al 2012. Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. Br J Cancer 107:1776
- Popova T et al 2012. Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation. Cancer Res 72:5454
- Birkbak NJ et al 2012. Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents. Cancer Discov 2:366
- Telli ML et al 2016. Homologous recombination deficiency (HRD) score predicts response to platinum-containing neoadjuvant chemotherapy. Clin Cancer Res 22:3764
- Davies H et al 2017. HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures. Nat Med 23:517
- Sztupinszki Z et al 2018. Migrating the SNP array-based homologous recombination deficiency measures to next generation sequencing data (scarHRD). NPJ Breast Cancer 4:16

## Related Skills

- copy-number/allele-specific-copy-number - Allele-specific copy number input for the scars
- copy-number/subclonal-copy-number - Whole-genome-doubling detection for LST correction
- copy-number/recurrent-cnv - Copy-number signatures, including the HRD-associated signature
- copy-number/cnv-annotation - Annotating HR-pathway gene copy-number status
- clinical-databases/somatic-signatures - SNV mutational signature 3 (HRD substitution signature)
- clinical-databases/variant-prioritization - BRCA1/2 and HR-pathway variant interpretation
