# Somatic Mutational Signatures - Usage Guide

## Overview

Extract and assign COSMIC v3.4 mutational signatures (84 SBS / 11 DBS / 18 ID / 24 CN / 16 SV) from somatic VCFs using SigProfilerSuite, MutationalPatterns, MuSiCal (mvNMF), SigNet (deep learning for low mutation counts), or HRDetect (BRCA1/2 deficiency classifier). Covers de novo extraction vs refit-to-COSMIC choice by cohort size, FFPE-artifact-as-SBS30-not-SBS33 correction, Petljak 2022 APOBEC3A-dominance subtyping, and clinical actionability for PARP inhibitor (HRD) / ICI (MMR-D, POLE) / therapy-induced (5-FU SBS17b, platinum SBS31/35).

## Prerequisites

```bash
# SigProfilerSuite (recommended)
pip install SigProfilerMatrixGenerator SigProfilerExtractor SigProfilerAssignment
# One-time reference download: ~3 GB
python -c "from SigProfilerMatrixGenerator import install; install.install('GRCh38')"

# MuSiCal (mvNMF for non-uniqueness)
pip install musical

# SigNet (deep learning, low mutation count)
# git clone https://github.com/weghornlab/SigNet

# R: MutationalPatterns + BSgenome
# BiocManager::install(c('MutationalPatterns', 'BSgenome.Hsapiens.UCSC.hg38'))

# HRDetect
# devtools::install_github('Nik-Zainal-Group/signature.tools.lib')
```

## Quick Start

Tell the agent what to do:
- "Extract de novo SBS signatures from this 80-sample WGS cohort with SigProfilerExtractor; report stability gates"
- "Refit my tumor cohort to COSMIC v3.4 with SigProfilerAssignment; identify dominant signatures"
- "Run HRDetect on these BRCA1/2-status-unknown breast cancer samples; flag HRD-positive for PARP eligibility"
- "Differentiate APOBEC3A vs APOBEC3B activity in this cohort using YTCA/RTCA tetranucleotide ratios"
- "For my 5-FU-treated CRC samples, identify SBS17b therapy-induced contribution"

## Example Prompts

### De Novo vs Refit

> "For my 100-sample WGS cohort with putative novel mutagen exposure, run SigProfilerExtractor de novo with nmf_replicates=100, validate stability gates (min >= 0.2; avg >= 0.8), and report any signatures not in COSMIC v3.4."

> "For my 30-sample cohort (below de novo threshold), use SigProfilerAssignment to refit to COSMIC v3.4 with strict forward-backward selection."

### HRD Classification

> "Run HRDetect on this breast cancer cohort. Compute SBS3 + ID6 + RS3 + RS5 + HRD-LOH; threshold BRCA_prob >= 0.7. Cross-check against germline BRCA1/2 + BRCA1 methylation status."

> "For HR-deficient tumors without BRCA1/2 germline mutations, identify PALB2 / FBXW7 / CDK12 alterations and BRCA1 promoter hypermethylation."

### APOBEC Mechanism

> "Distinguish APOBEC3A vs APOBEC3B activity in this bladder cancer cohort using YTCA vs RTCA 5' tetranucleotide ratios (per Petljak 2022). Co-locate kataegis loci."

### Therapy-Induced

> "For this CRC cohort treated with 5-FU, separate SBS17b therapy-induced from SBS17a unknown-etiology contributions."

> "For platinum-treated tumors, identify SBS31 + SBS35 + SBS86 + SBS87 contributions; report total platinum-attributable mutation burden."

### FFPE Artifact Correction

> "These are FFPE samples; identify FFPE-artifact signal (resembling SBS30, NOT SBS33 as legacy literature claims). Recommend running matched fresh-frozen controls or enzymatic uracil pretreatment."

### Tumor Evolution

> "Run MutationTimer on this WGS cohort with Battenberg CN to time mutations relative to chromosomal-instability events. Identify clonal vs subclonal signature activity."

## What the Agent Will Do

1. Generate 96-context (SBS) / 78-DBS / 83-ID / 48-CN matrix from VCFs via SigProfilerMatrixGenerator (set `exome=True` for WES).
2. Choose de novo vs refit by cohort size: >= 50 WGS with possible novel etiology -> de novo; < 50 or single sample -> refit.
3. Run with stability gates: nmf_replicates=100, min_stab >= 0.2, avg_stab >= 0.8.
4. Decompose dominant signatures into etiology categories: HRD (SBS3 + ID6 + CN17), MMR-D (SBS6/14/15/20/21/26/44 + ID1/2), POLE (SBS10a/10b/28), APOBEC (SBS2/13), UV (SBS7a-d), tobacco (SBS4), aflatoxin (SBS24), 5-FU (SBS17b), platinum (SBS31/35).
5. For HRD classification, run HRDetect 6-feature lasso (SBS3 + ID6 + RS3 + RS5 + HRD-LOH + indel-microhomology).
6. Flag FFPE-artifact-as-SBS30 (NOT SBS33); recommend matched fresh-frozen controls or enzymatic uracil pretreatment.
7. Use SigProfilerTopography for spatial preferences (replication timing, transcribed-strand, nucleosome occupancy); detect aristolochic-acid (SBS22) transcribed-strand bias.
8. Pin COSMIC version; cross-version comparison requires re-extraction.

## Tips

- COSMIC v3.4 (Sept 2024) is current; signatures split since v3.3 include SBS40 -> 40a/b/c (Senkin 2024), SBS17 -> 17a/b (5-FU; Christensen 2019), SBS10 -> 10a/b/c/d (POLE/POLD1).
- Single-sample de novo is unstable; require >= 200 mutations per sample AND cohort >= 50.
- SigProfilerExtractor stability gates are non-negotiable: nmf_replicates=100, min >= 0.2, avg >= 0.8.
- deconstructSigs is **operationally deprecated**; NNLS overfits onto reference set. Use SigProfilerAssignment or MutationalPatterns strict refit.
- FFPE artifact resembles SBS30 (NOT SBS33 as commonly cited); enzymatic uracil pretreatment shifts to SBS1-like.
- WES requires trinucleotide-context correction (`exome=True`); WES-derived signatures NOT directly comparable to WGS without correction.
- Aristolochic-acid SBS22 shows strong transcribed-strand bias; use SigProfiler `tsb_stat=True` or MutationalPatterns.
- APOBEC: A3A vs A3B distinguished by YTCA vs RTCA tetranucleotide ratio (Petljak 2022 *Nature*); not all tools surface this.
- HRDetect requires 6 features (SBS3 + ID6 + RS3 + RS5 + HRD-LOH + indel-microhomology) for 98.7% sensitivity; single-feature (SBS3 alone) is insufficient.
- SBS5 etiology is contested; report as "unknown, clock-like"; NOT polymerase fidelity errors.
- POLE-exo + MMR-D produces ultra-hypermutator (>500 mut/Mb) phenotype with SBS14 + SBS20.
- For low-mutation-count single samples, SigNet (Serrano 2023) outperforms NMF-based approaches.
- Cross-cancer HRD generalization: HRDetect was trained on breast cancer; cross-cancer requires revalidation with HRD-LOH.

## Related Skills

- clinical-databases/tumor-mutational-burden - TMB and ICI biomarker
- clinical-databases/msi-detection - MSI as MMR-D biomarker (paired with SBS6/15/26/44)
- variant-calling/variant-calling - Somatic VCF input
- variant-calling/variant-calling - Mutect2 / Strelka2 somatic upstream
- data-visualization/heatmaps-clustering - Signature contribution visualization
