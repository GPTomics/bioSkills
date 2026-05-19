# Pharmacogenomics - Usage Guide

## Overview

Call CYP2D6, CYP2C9, CYP2C19, DPYD, TPMT, NUDT15, UGT1A1, SLCO1B1, CYP3A5, VKORC1, and other CPIC Level A pharmacogenes from sequencing data; translate diplotype to activity score (Caudle 2020) and phenotype; apply CPIC and DPWG dosing guidance. Covers HLA-B*57:01 / B*15:02 / B*58:01 / A*31:01 / B*35:02 screening for ICI / antiepileptics / abacavir / minocycline, with explicit recognition of the CYP2D6 SV footgun (*4xN clinical silence), the 2024 DPYD activity-score framework, the 2025 TPMT+NUDT15 compound-IM update, and the COAG paradigmatic ancestry-algorithm failure.

## Prerequisites

```bash
# PharmCAT (CPIC-recommended multi-gene)
java -jar pharmcat.jar --help
# Download from https://pharmcat.org

# Cyrius (CYP2D6 SV-aware; mandatory for CYP2D6)
pip install cyrius
# or: conda install -c bioconda cyrius

# Aldy (alternative multi-gene)
pip install aldy

# HLA imputation from SNP arrays
# R: install.packages('HIBAG')
# Download ancestry-matched reference panel from https://hibag.s3.amazonaws.com/

# PharmGKB API
pip install requests pandas
```

## Quick Start

Tell the agent what to do:
- "Run PharmCAT on this VCF; supplement CYP2D6 with Cyrius; cross-reference HLA from typing"
- "Calculate CYP2D6 activity score for *4xN/*10 (this is the clinical-silence footgun)"
- "For this East-Asian patient, screen NUDT15 *3 before mercaptopurine; do NOT rely on TPMT alone"
- "Apply CPIC 2024 DPYD activity-score framework to determine fluoropyrimidine dose"
- "Pre-emptive PGx panel from VCF: report all CPIC Level A genes with phenotype + dosing"
- "Screen HLA-B*57:01 before abacavir; require 4-field specificity to distinguish from *57:03"

## Example Prompts

### Multi-Gene Pipeline

> "Run PharmCAT on my cohort VCF. For samples with CYP2D6 suspected SVs (CN!=2 in BAM coverage), run Cyrius and merge outside calls. Generate per-sample CPIC report."

> "Apply Swen 2023 PREPARE-style 12-gene panel to this cohort; flag patients with actionable interactions for their prescribed drugs."

### CYP2D6-Specific

> "For sample with CYP2D6 *4xN/*10, compute activity score using Caudle 2020 values. Confirm *4xN is clinically silent (AS contribution = 0)."

> "Run Cyrius on this WGS BAM; report CYP2D6 diplotype + activity + phenotype + copy number. Cross-check against PharmCAT prediction."

> "Patient is *2xN/*4. Is this an Ultra-rapid Metabolizer? Show activity score calculation."

### DPYD 2024 Framework

> "Apply CPIC 2024 DPYD activity-score framework to my fluoropyrimidine-treated cancer cohort. Identify AS 0 (avoid), AS 1.0-1.5 (50% start + TDM), AS 2.0 (full dose). Note: c.85T>C is NOT in the CPIC 2024 actionable set."

> "For African-ancestry patients, supplement the 4-variant CPIC core panel with extended DPYD coverage."

### TPMT + NUDT15

> "Apply Maillard 2026 update: compound TPMT IM + NUDT15 IM requires more aggressive dose reduction than single-gene IM. Report dose recommendation."

> "Screen East-Asian thiopurine cohort for NUDT15 *3 (c.415C>T, rs116855232). Do NOT rely on TPMT-only testing."

### HLA Screening

> "Pre-prescription HLA-B*57:01 screen for abacavir; require 4-field. Distinguish *57:01 (risk) from *57:02 / *57:03 (no risk)."

> "Screen this Han Chinese patient for B*15:02 (carbamazepine SJS), B*58:01 (allopurinol), B*13:01 (dapsone) before any prescribing."

> "Impute HLA from SNP-array genotypes using HIBAG with the ancestry-matched reference panel."

### Clopidogrel and CYP2C19

> "For this post-PCI ACS patient with CYP2C19 *2/*17, recommend P2Y12 selection per CPIC. Cite TAILOR-PCI sensitivity analyses + Pereira 2021 meta-analysis."

### Warfarin Cross-Ancestry

> "Apply IWPC warfarin algorithm explicitly INCLUDING CYP2C9 *5/*6/*8/*11 for this African-ancestry patient. Document COAG paradigmatic ancestry-algorithm failure rationale."

## What the Agent Will Do

1. Choose appropriate tool stack: PharmCAT for the multi-gene panel + Cyrius for CYP2D6 (SV-aware) + dedicated HLA typer for HLA loci.
2. Apply Caudle 2020 activity values for CYP2D6 (note *10 = 0.25, not 0.5).
3. Handle CYP2D6 SV correctly: *4xN is clinically silent; *1xN/*2xN multiply functional activity.
4. Apply CPIC 2024 DPYD activity-score framework (not single-variant logic).
5. Pair TPMT + NUDT15 for thiopurines; apply Maillard 2026 compound-IM dosing.
6. Require 4-field HLA resolution for PGx screening; distinguish risk-allele from non-risk family members.
7. For African-ancestry warfarin, explicitly include CYP2C9 *5/*6/*8/*11; document COAG failure rationale.
8. Cite both CPIC and DPWG when they disagree.

## Tips

- PharmCAT is the CPIC-recommended multi-gene tool; pair with Cyrius for CYP2D6 SVs.
- Cyrius is 99.3% accurate on GeT-RM reference samples vs Aldy 82-87% vs Stargazer 84%.
- CYP2D6 *4xN is clinically silent (no-function * N = 0). This is the most common reportable error in clinical PGx.
- Caudle 2020 reset CYP2D6 *10 activity from 0.5 to 0.25; reclassified large fractions of East-Asian populations to IM.
- CPIC 2024 DPYD update moves to activity-score framework; c.85T>C is NOT in the CPIC 2024 actionable set despite frequent commercial reporting.
- TPMT-only testing misses NUDT15 *3 (9.8% frequency in Han Chinese vs <1% EUR). Always pair TPMT + NUDT15.
- HLA-B*57:01 vs *57:02 / *57:03 distinguishes abacavir risk vs no-risk; require 4-field.
- HLA-B*35:02 (minocycline DILI) vs *35:01 (TMP-SMX DILI); different drugs, different alleles.
- COAG (Kimmel 2013 *NEJM*) failed in African Americans because the algorithm omitted CYP2C9 *5/*6/*8/*11; paradigmatic ancestry-algorithm failure.
- TAILOR-PCI primary endpoint was negative (HR 0.66, p=0.06), but Pereira 2021 meta-analysis of 7 RCTs shows ~30% MACE reduction.
- PREPARE (Swen 2023) showed 30% reduction in actionable ADRs with pre-emptive 12-gene panel.
- CPIC vs DPWG: CPIC tells you what to do once you have a result; DPWG also advises whether to test. They sometimes disagree.
- HLA imputation from SNP arrays drops 10-20 percentage points without ancestry-matched reference panels.
- CYP3A5 *3/*3 is the COMMON state in non-AFR (non-expressers); expressers (any *1) need 1.5-2x higher tacrolimus dose.

## Related Skills

- clinical-databases/hla-typing - HLA-B*57:01, B*15:02, B*58:01, A*31:01 typing
- clinical-databases/clinvar-lookup - Non-PGx variant pathogenicity
- clinical-databases/variant-prioritization - Rare-disease pipeline
- clinical-databases/myvariant-queries - Aggregated PGx variant annotation
- chemoinformatics/admet-prediction - Drug metabolism prediction
