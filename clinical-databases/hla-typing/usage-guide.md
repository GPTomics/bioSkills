# HLA Typing - Usage Guide

## Overview

Call HLA class I (A, B, C) and class II (DRB1, DRB3/4/5, DQA1, DQB1, DPA1, DPB1) alleles at the resolution required for HSCT, solid-organ transplant, neoantigen prediction, pharmacogenomic screening, and disease-association studies. Covers short-read tools (OptiType, HLA-LA, T1K, Polysolver, HLA-HD, arcasHLA), long-read (StarPhase, FuFiHLA), and SNP-array imputation (HIBAG, HLA-TAPAS) with explicit reconciliation when tools disagree.

## Prerequisites

```bash
# Short-read all-rounder
conda install -c bioconda t1k

# Class I gold standard (TCGA-compatible)
conda install -c bioconda optitype

# Class II reference grade (needs ~30-100 GB scratch)
conda install -c bioconda hla-la
# Build PRG (one-time, ~12 hours): HLA-LA.pl --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT

# RNA-seq
pip install arcas-hla
arcasHLA reference --update  # Pull current IPD-IMGT/HLA release

# Long-read (PacBio HiFi)
# StarPhase: install from PacBio Bioconda channel or follow PacBio docs

# SNP-array imputation (R)
# install.packages('HIBAG')
# Download ancestry-matched reference panel from https://hibag.s3.amazonaws.com/
```

## Quick Start

Tell the agent what to do:
- "Run T1K on my WGS BAM to call HLA-A, B, C, DRB1, DRB3/4/5, DQB1, DPB1 at 4-field"
- "Type HLA from my RNA-seq with arcasHLA for neoantigen prediction"
- "Screen this patient for HLA-B*57:01 before abacavir; require 4-field specificity"
- "Impute HLA from SNP-array genotypes using HIBAG with the African-ancestry panel"
- "For my HSCT donor search, type all class I + II at 6-field on PacBio HiFi using StarPhase"
- "Compare T1K vs HLA-LA outputs on these 20 samples and reconcile discrepancies"

## Example Prompts

### Standard Typing

> "Run T1K on these 50 WES samples. Output 4-field calls for all class I + II loci. Verify the DRB1+DRB3/4/5 linkage rule on each sample and flag any violations."

> "Type HLA-A, B, C with OptiType on this TCGA-compatible workflow. Compare to T1K and report concordance."

> "For my tumor + normal RNA-seq pairs, run arcasHLA on both. Flag any tumors with HLA LOH (allele loss vs normal)."

### Pharmacogenomic Screening

> "Screen these patients for HLA-B*57:01 before abacavir; require 4-field resolution. Distinguish *57:01 from *57:02 and *57:03 explicitly."

> "Pre-prescription HLA screen for carbamazepine: check B*15:02 (Han Chinese SJS risk) and A*31:01 (European DRESS). Report ancestry-specific risk."

> "For this Han Chinese cohort, screen B*58:01 before allopurinol and B*13:01 before dapsone."

### Transplant Matching

> "Type all class I + II at 6-field for this HSCT donor with PacBio HiFi using StarPhase. Output null alleles (N-suffix) explicitly; do not strip suffixes for clinical reporting."

> "For this DPB1 mismatch, classify as TCE3 core (GvHD-reduction permissive) vs non-core (relapse-protection) per Meurer 2024."

> "Calculate eplet mismatch score (HLAMatchmaker) for this donor-recipient pair at HLA-A, B, DRB1, DQB1."

### SNP-Array Imputation

> "From this UK Biobank-style array, impute HLA-A, B, C, DRB1 at 4-field using HIBAG. Use ancestry-matched panels (EUR for NFE samples, AFR for AFR samples)."

> "For multi-ancestry GWAS, use HLA-TAPAS multi-ancestry reference panel for HLA imputation."

### Cross-Tool Reconciliation

> "Run both T1K and HLA-LA on these 20 samples; reconcile class II disagreements and flag samples needing long-read confirmation."

> "Tumor HLA differs from germline HLA in this patient; run LOHHLA to confirm somatic LOH and report germline class I as the actionable genotype."

## What the Agent Will Do

1. Choose the right tool: T1K for general WGS/WES (class I + II + KIR), OptiType for class I only or TCGA convention, HLA-LA for class II reference grade, arcasHLA for RNA-seq, StarPhase for long-read transplant grade, HIBAG for SNP arrays.
2. Verify alt-aware alignment was used; extract chr6 + HLA alt contigs explicitly.
3. Update the tool's IPD-IMGT/HLA reference bundle if stale (>6 months old).
4. Run typing with appropriate threads / RAM (HLA-LA needs 30-100 GB scratch).
5. Verify DRB1+DRB3/4/5 linkage rule as routine QC.
6. Preserve N/L/S/Q/A expression suffixes for clinical reports.
7. Reconcile cross-tool disagreements; for PGx, require 4-field specificity (B*57:01 not B*57).
8. For tumor samples, check HLA LOH via LOHHLA.

## Tool Selection Reference

| Scenario | Tool | Resolution |
|----------|------|-----------|
| WGS/WES general purpose | T1K | 4-field, class I + II + KIR |
| WGS/WES class I only | OptiType | 4-field |
| WGS/WES class II reference | HLA-LA | 4-field |
| RNA-seq | arcasHLA | 4-field |
| TCGA cancer cohort | Polysolver | 4-field, class I |
| Long-read transplant | StarPhase (PacBio HiFi) | 8-field |
| Long-read platform-agnostic | FuFiHLA (PacBio + ONT) | 8-field |
| SNP array | HIBAG (ancestry-matched) | 4-field |
| Multi-ancestry GWAS | HLA-TAPAS | 4-field |
| Cost-effective targeted typing | pbaa + StarPhase | 8-field |

## Tips

- Tool reference-bundle vintage matters more than algorithm choice for non-European samples. Update IPD-IMGT/HLA bundles quarterly.
- 4-field resolution is the clinical standard for PGx screening; 2-field (e.g., "B*57") is insufficient for distinguishing *57:01 (abacavir risk) from *57:03 (no risk).
- DRB1+DRB3/4/5 haplotype linkage is fixed; use as routine QC check.
- N-suffix alleles (null, no protein expression) must be preserved in transplant typing reports.
- Pre-2021 WES capture kits (SureSelect v5, Nextera Rapid Capture) under-cover DPB1 exon 2; switch to WGS or amplicon for class II certainty.
- DRA is monomorphic; do not waste effort typing it.
- KIR locus is functionally paired with HLA (KIR3DL1 binds HLA-Bw4) but is on chr19; T1K co-types both.
- Tumor HLA can differ from germline due to LOH (~40% of NSCLC); run LOHHLA on tumor-normal pairs.
- HLA-G, HLA-E, HLA-F are non-classical (NK-related); not typed by standard tools.
- MICA/MICB are HLA-region but NK ligands; useful for transplant but typed separately.
- HLA-Bw4 vs Bw6 motif (residues 77-83) determines KIR3DL1 ligand status; callable directly from any 2-field B genotype.
- For African-ancestry transplant cases, long-read (PacBio HiFi or ONT R10.4 duplex) is recommended; IPD-IMGT/HLA still under-represents AFR alleles ~30-40%.
- alt-aware alignment is essential: bwa-mem-alt with the GRCh38 HLA alt contigs; non-alt-aware drops accuracy 5-10 percentage points.

## Related Skills

- clinical-databases/pharmacogenomics - HLA-drug interactions
- immunoinformatics/mhc-binding-prediction - Downstream neoantigen prediction
- workflows/neoantigen-pipeline - End-to-end neoantigen workflow
- clinical-databases/clinvar-lookup - HLA disease associations
- population-genetics/population-structure - Ancestry-aware imputation context
