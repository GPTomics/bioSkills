# Base Editing Analysis - Usage Guide

## Overview

Decision-grade analysis of base-editor variant-function screens. Covers GRACE-style library design (Sanson 2020 *Nat Commun*), the Hanna 2021 BRCA1/2 SNV-scanning methodology (*Cell* 184:1066), Cuella-Martin 2021 corroboration (*Mol Cell* 81:1107), CBE vs ABE chemistry selection (BE3/BE4 vs ABE7.10/ABE8.20/ABE8e), editing window math (positions 4-8 from PAM-distal end), bystander attribution strategies, editing-efficiency filtering before hit calling, indel-byproduct interpretation, and the Broad be-validation-pipeline end-to-end workflow.

## Prerequisites

```bash
pip install CRISPResso2 pandas biopython numpy scipy
# BE-Hive for editing efficiency prediction
pip install be-hive
# Broad end-to-end pipeline
git clone https://github.com/broadinstitute/be-validation-pipeline
cd be-validation-pipeline && conda env create -f environment.yaml
```

Required inputs:
- Amplicon sequencing FASTQ (per sgRNA or per pool)
- Library file: per-sgRNA spacer, target base, target amino acid, predicted bystander pattern
- BE chemistry vector (CBE or ABE)
- (For variant attribution) ClinVar / COSMIC variant annotation

## Quick Start

Tell the AI agent what to do:
- "Design a CBE library tiling BRCA1 exons 1-23 with 10-15 sgRNAs per amino acid; flag intended target Cs and bystander Cs in the editing window"
- "Run editing-efficiency filtering on pilot library counts; drop sgRNAs <30% target editing"
- "Score variant function per amino acid; deconvolute bystander confounding via allele-frequency tables"
- "Reconcile a screen hit between Hanna's CBE methodology and orthogonal prime-editor follow-up"
- "Choose BE vs PE for installing a C->T variant: BE3 if target at pos 5 with no bystanders; PE2 if zero-bystander requirement"
- "Diagnose: my BE sample has 60% editing but 35% are indels -- is this Cas9 contamination?"

## Example Prompts

### Library Design

> "Design a CBE saturation library tiling BRCA1 RING domain (amino acids 1-100). 10-15 sgRNAs per amino acid where at least one C in the editing window (positions 4-8) hits the target codon. Annotate each sgRNA with predicted target + bystander variants. Output library.tsv with sgRNA, target_aa, target_variant, bystander_variants columns."

> "I need to install MLH1 c.677A>G as a single intended variant. CBE won't work (need ABE). Find ABE7.10 or ABE8e sgRNAs that place A at position 5 with no bystanders. If no zero-bystander spacer exists, list candidates sorted by bystander_count and recommend prime editor as alternative."

### Editing Efficiency Filtering

> "Run CRISPResso2 on my pilot timepoint samples. Compute target editing % per sgRNA. Apply Hanna 2021 threshold: keep sgRNAs >30% target editing for primary screen; flag >50% as validation-grade. Output filtered library and editing efficiency report."

> "My library shows median 22% editing across guides. Diagnose: cell-line BE activity issue, vector mismatch, or BE chemistry choice (CBE vs ABE)?"

### Hit Calling for Variant Function

> "Apply MAGeCK MLE to the BE screen counts (vehicle vs PARPi drug). Use only editing-efficient sgRNAs. Aggregate per-sgRNA LFC to per-variant scores; deconvolute bystanders via allele tables; output per-variant fitness with confidence based on the number of sgRNAs hitting each variant."

> "Compare drugZ vs MAGeCK MLE on the same BE drug-modifier screen; identify novel resistance/sensitivity variants in MCL1, BCL2L1, PARP1 (Hanna 2021 method)."

### Bystander Deconvolution

> "From CRISPResso2 allele tables, separate reads by edit pattern: target only, target+bystander_1, target+bystander_2. Compute per-pattern fitness contribution; flag variants where bystander dominates."

> "For BRCA1 R71 -> R71X variant, I have 5 sgRNAs with different bystander patterns. Identify variants attributable to R71X alone by finding consistent signal across diverse bystander backgrounds."

### Method Comparison

> "Compare CBE BE3 vs BE4max vs eA3A-BE3 vs evoCDA-BE for tiling BRCA1. Score: editing window width, bystander rate, indel byproduct, target-base preference (TC vs other)."

> "Cross-validate Hanna 2021 BRCA1 BE screen hits with parallel PE screen at the same variants. Report concordant variants and BE-only (likely bystander-confounded) hits."

### Diagnostics

> "My BE sample has 70% editing but 25% indels. Compute substitution-vs-indel ratio; if <3, diagnose as Cas9 contamination."

> "Allele table shows 35% target+bystander double-edit. Deconvolute: how much is target-attributable vs bystander-driven?"

## What the Agent Will Do

1. Design BE library: for each target amino acid, find NGG-PAM spacers with target base at editing window positions 4-8; minimize bystanders
2. Annotate each sgRNA: target_variant, bystander_pattern (other edits in window)
3. Order library; receive plasmid pool; sequence to verify Gini <0.1
4. Lentiviral package and infect at MOI 0.3 in BE-validated cell line (HEK293T, U2OS, K562)
5. Pilot timepoint amplicon sequencing of ~20 representative loci
6. Run CRISPResso2 with `--base_editor_output`, `--conversion_nuc_from C --conversion_nuc_to T` (or A->G for ABE)
7. Compute editing efficiency per sgRNA; filter to >30% (primary) or >50% (validation)
8. Run main screen with vehicle vs drug treatment
9. Per-sgRNA hit calling via MAGeCK MLE or drugZ
10. Aggregate to per-variant scores; deconvolute bystander via allele tables
11. Cross-check substitution-vs-indel ratio per sgRNA (>10 for clean BE)
12. Validate top hits via orthogonal prime-editor or arrayed BE
13. Annotate against ClinVar / COSMIC for clinical interpretation

## Tips

- Editing efficiency filtering is non-negotiable. Unedited reads carry no biological perturbation; including <30% efficiency sgRNAs adds noise and dilutes signal. Hanna 2021 set the field-standard 30% primary / 50% validation threshold.
- Bystander confounding is the central interpretation challenge in BE screens. Plan the library with multiple sgRNAs per target (10-15 in Hanna 2021) so bystander patterns vary; consistent signal across diverse bystanders attributes to target.
- For pristine variant-function calls (clinical-grade), validate with prime editor (zero bystanders). Convergent BE + PE signal is the gold standard for variant pathogenicity in pooled screens.
- Cell-line BE activity varies; pilot in your target line before designing the full library. Median library editing <30% is a cell-line issue, not a library issue.
- ABE7.10 has lower indel rate (<2%) than CBE (5-10%); when target is A->G, ABE is the cleaner choice.
- For C->T variants with target at position 4-5 (PAM-distal end) and bystander at position 6-7, both will edit; this is the most common interpretation trap.
- Editor chemistries have different windows: BE3/BE4 = 4-8, eA3A = 5-7 (narrower; good for bystander minimization), ABE8e = 4-10 (wider; more bystanders), evoCDA = 1-9 (broadest; high bystander rate).
- Run substitution-vs-indel ratio per sgRNA as a per-sample QC. <3 means Cas9-like activity (vector mismatch or contamination); >10 means clean BE.
- The Broad be-validation-pipeline is the canonical reference workflow; reuse it before writing custom scripts.
- For drug-modifier BE screens (like Hanna 2021 PARPi), drugZ is more sensitive than MAGeCK for chemogenomic interactions.

## Chemistry Cheat Sheet

| Need | Editor | Why |
|------|--------|-----|
| C->T at TC context | BE3, BE4max, eA3A-BE3 | eA3A-BE3 narrowest window |
| C->T at any context | BE4max | Standard CBE |
| A->G | ABE7.10, ABE8.20, ABE8e | ABE8e fastest editing |
| Multi-base / no bystander | Prime editor | See [[prime-editing-screens]] |
| Transversions (C->G, C->A) | CGBE1, GBE | Rare use; less mature |
| Saturation mutagenesis | BE3 + PE for orthogonal validation | Combined coverage |

## Validation Strategy

| Tier | Validation requirement |
|------|-------------------------|
| Tier 1 (high confidence) | BE + PE concordant at same variant + arrayed confirmation |
| Tier 2 (medium) | BE alone, multiple sgRNAs converge despite bystander differences |
| Tier 3 (exploratory) | Single sgRNA hit; bystander confounded; not interpretable |

## Related Skills

- crispr-screens/crispresso-editing - CRISPResso2 modes and allele tables
- crispr-screens/library-design - GRACE-style BE library design
- crispr-screens/prime-editing-screens - Orthogonal PE for variant attribution
- crispr-screens/hit-calling - Variant-level hit aggregation
- crispr-screens/screen-qc - Editing-efficiency QC
- crispr-screens/drugz-chemogenomic - drugZ for BE drug-modifier screens
- crispr-screens/mageck-analysis - MAGeCK MLE for BE screen sgRNA-level analysis
- clinical-databases/clinvar-lookup - Pathogenicity annotation
- variant-calling/variant-annotation - VEP for predicted protein effects
