# Copy-Number Correction - Usage Guide

## Overview

Decision-grade copy-number-bias correction for cancer-cell-line CRISPR screens. Addresses the Aguirre 2016 / Munoz 2016 p53-dependent DNA-damage artifact where amplified loci appear "essential" purely from cut-count toxicity. Covers CRISPRcleanR (unsupervised pre-hoc), CERES (legacy joint model), Chronos (DepMap-standard joint cell-population-dynamics + CN), the decision tree among them by data availability, the diagnostic correlation of gene LFC with copy number, and the alternative of switching from Cas9 to CRISPRi to bypass the artifact entirely.

## Prerequisites

```bash
# CRISPRcleanR (R)
R -e "BiocManager::install('CRISPRcleanR')"
# Chronos (Python)
pip install chronos-cn
# or from source: git clone https://github.com/broadinstitute/chronos
# Helpers
pip install pandas numpy scipy
```

Required inputs:
- Raw or normalized counts from mageck count
- Copy-number profile per cell line (from WGS / SNP-array / ASCAT / matched cell-line database). Required for Chronos; optional for CRISPRcleanR.
- Library annotation: sgRNA -> gene -> chromosomal coordinates (for CRISPRcleanR)

## Quick Start

Tell the AI agent what to do:
- "Diagnose copy-number bias in my SK-BR-3 (HER2+) screen: Spearman ρ between gene LFC and CN profile. If ρ <-0.1, apply correction"
- "Apply CRISPRcleanR to my screen (no CN profile available) and output corrected counts ready for MAGeCK"
- "Run Chronos on my multi-cell-line cancer panel with longitudinal counts and matched CN profiles; output gene-effect scores"
- "Compare CRISPRcleanR-corrected vs Chronos-corrected hit lists for the same screen"
- "Decide: switch from Cas9 to CRISPRi (Dolcetto) to bypass the artifact entirely in our HER2-amplified line"

## Example Prompts

### Diagnostic

> "Detect copy-number bias in my cancer-cell-line screen. Compute Spearman ρ between gene-level LFC and copy number profile. Flag bias if ρ <-0.10 and p <0.01. Report mean LFC for amplified (CN >4) and diploid (CN 1.5-2.5) genes."

> "After correction, re-run the diagnostic. Confirm post-correction Spearman ρ is near zero (abs(ρ) <0.05)."

### CRISPRcleanR Application

> "Apply CRISPRcleanR to my single-cell-line Brunello screen without CN profile. Output corrected logFCs and corrected counts compatible with MAGeCK input."

> "After CRISPRcleanR correction, run MAGeCK test on the corrected counts. Compare hit list to pre-correction hit list; identify which hits were artifacts."

### Chronos Application

> "Run Chronos on my 5-cell-line panel with 14-day longitudinal sampling and matched WGS CN profiles. Output gene-effect scores per cell line and gene probability of essentiality."

> "Match the DepMap quarterly Chronos workflow on my panel. Output standardized gene effects (essential <-1, non-essential >0)."

### Comparison

> "Apply CRISPRcleanR vs Chronos to the same screen; compare hit lists; identify cases where they disagree."

> "Reconcile a hit at MYC in my MYC-amplified colorectal line: is it the gene effect or the CN artifact?"

### Bypass via CRISPRi

> "My SK-BR-3 (HER2-amplified) screen has severe CN artifact at ERBB2. Switch from Cas9 to CRISPRi (Dolcetto library) to bypass the artifact. Quantify the gain in interpretable hits."

## What the Agent Will Do

1. Identify if screen is in cancer cell line; if yes, CN correction needed
2. Check if CN profile is available; if not, default to CRISPRcleanR (unsupervised)
3. If multi-cell-line + longitudinal + CN available, default to Chronos
4. Run pre-correction diagnostic: Spearman ρ between gene LFC and CN profile
5. Apply correction method:
   - CRISPRcleanR: `ccr.GWclean()` -> `ccr.correctCounts()` for MAGeCK input
   - Chronos: train model with `sequence_map`, `guide_gene_map`, `reads`, `copy_number`; extract gene_effect()
6. Run post-correction diagnostic; confirm Spearman ρ near zero
7. Pass corrected counts to downstream hit calling (see [[hit-calling]])
8. Cross-check hits against known amplifications in the cell-line (e.g., COSMIC, DepMap)
9. Report correction effectiveness and per-amplification false-positive rates

## Tips

- The CN artifact is universal in cancer cell lines, not conditional. Always check; always correct. ERBB2 in HER2+, MYC in MYC-amplified, FGFR1 in head-and-neck are textbook cases.
- CRISPRcleanR is unsupervised (no CN profile needed) and works on any screen. Chronos requires matched CN AND longitudinal data but is more rigorous.
- For DepMap-style large panels, Chronos is the standard; no need to use CRISPRcleanR. For Sanger Score-style or single-cell-line screens without DepMap-grade data, CRISPRcleanR is the workhorse.
- The artifact applies to Cas9-KO screens only; CRISPRi/a (catalytically dead Cas9), base editor, and prime editor screens are largely free of it.
- Switching from Cas9 to CRISPRi (Dolcetto library) is the cleanest way to bypass the artifact in cancer-line essentiality screens. The tradeoff: CRISPRi knockdown is less complete than Cas9 KO, so some essentials are missed.
- After correction, re-running the diagnostic (Spearman) is mandatory. If correction was insufficient, refine CN profile, try Chronos with the new CN, or accept the bias for specific amplified loci.
- For non-cancer cell lines (HEK293T, U2OS, iPSC), the artifact is much less pronounced because there are few focal amplifications. Correction is still good practice but rarely critical.
- For variant-function screens in cancer lines, use base editor or prime editor (which create fewer DSBs) rather than Cas9 KO.

## Decision Cheat Sheet

| Have CN profile? | Multi-cell-line + multi-timepoint? | Use |
|-------------------|-------------------------------------|-----|
| No | N/A | CRISPRcleanR |
| Yes | No | CRISPRcleanR or Chronos |
| Yes | Yes | Chronos (DepMap standard) |

## Cancer-Line CN Artifacts to Watch For

| Cell line | Amplification | Gene that looks essential |
|-----------|---------------|----------------------------|
| SK-BR-3 | HER2 (ERBB2) ~24 copies | ERBB2 |
| MCF7 | MYC moderate amp | MYC |
| K562 | BCR-ABL | ABL1 |
| HepG2 | MYC amp | MYC |
| A549 | KRAS WT, but some lines have KRAS amp | KRAS |
| Head-and-neck panel | FGFR1 11-20 copies | FGFR1 |
| Colon panel | MYC ~10 copies | MYC |

## Validation Checklist

- [ ] CN diagnostic Spearman ρ computed pre and post correction
- [ ] Post-correction Spearman ρ abs <0.05
- [ ] Known amplification "hits" are gone post-correction
- [ ] Known true essentials (CEGv2) preserved post-correction
- [ ] If switching to CRISPRi, library validates against Dolcetto reference
- [ ] Hit list cross-checked against COSMIC / DepMap amplification database

## Related Skills

- crispr-screens/screen-qc - CN-LFC Spearman diagnostic
- crispr-screens/library-design - Switch to Dolcetto (CRISPRi) to bypass artifact
- crispr-screens/mageck-analysis - MAGeCK on CRISPRcleanR-corrected counts
- crispr-screens/bagel-essentiality - BAGEL2 on corrected counts
- crispr-screens/hit-calling - Cancer-line hit calling
- crispr-screens/batch-correction - Chronos handles batch + CN jointly
- copy-number/copy-ratio-segmentation - CN profile derivation upstream
