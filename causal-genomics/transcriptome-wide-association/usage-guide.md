# Transcriptome-Wide Association (TWAS) Skill Usage Guide

## Overview

Transcriptome-wide association studies (TWAS) test gene-level association with a GWAS trait via genetically predicted tissue expression. The agent runs TWAS from GWAS summary statistics using pre-trained eQTL prediction models (PrediXcan / FUSION / MetaXcan family), aggregates multi-tissue signals (S-MultiXcan / UTMOST), and probabilistically fine-maps co-significant gene clusters (FOCUS / MA-FOCUS).

TWAS is the natural cousin of Mendelian randomization and colocalization in the causal-genomics workflow. TWAS asks "is gene expression in tissue T associated with the trait?"; cis-eQTL MR asks "does that association reflect a causal effect under IV assumptions?"; coloc asks "do the GWAS and eQTL share a causal variant?". A single TWAS hit is associational, not causal -- the dominant failure modes are LD-induced co-significance at gene-dense loci, tissue mis-specification, and ancestry mismatch between GWAS and prediction weights. Strong causal claims require triangulation across at least 3 of {TWAS, FOCUS PIP, coloc PP.H4, cis-eQTL MR}.

## Prerequisites

- GWAS summary statistics with standard columns (SNP, A1/A2, BETA or Z, SE, P, N or use_n)
- Pre-trained prediction weights: FUSION panels from gusevlab.org/projects/fusion (.pos + per-gene RData) OR PredictDB models from predictdb.org (GTEx v8 elastic-net or MASHR, .db + covariance)
- Ancestry-matched LD reference panel (e.g. 1000 Genomes EUR; PLINK .bim/.bed/.fam by chromosome)
- Python 3.9-3.11 + R 4.3+ + PLINK 1.9 + PLINK 2.0
- Disk space: GTEx v8 PredictDB MASHR-EUR is approximately 8 GB; FUSION GTEx weights per tissue are 1-3 GB

Install Python tools:

```bash
pip install metaxcan pyfocus
```

Install R / FUSION:

```bash
git clone https://github.com/gusevlab/fusion_twas
R -e "install.packages(c('plink2R', 'optparse', 'glmnet', 'methods', 'RColorBrewer', 'GBJ'))"
```

For MA-FOCUS and multi-ancestry analyses, additionally pull `mancusolab/ma-focus`. For UTMOST, `Joker-Jerome/UTMOST`. For TIGAR-V2, `yanglab-emory/TIGAR`.

## Quick Start

Tell the AI agent what to do:
- "Run TWAS on my GWAS sumstats using FUSION with GTEx whole blood weights and fine-map hits with FOCUS"
- "Use S-PrediXcan across all 49 GTEx tissues and combine with S-MultiXcan for a joint multi-tissue test"
- "I have a list of TWAS hits at a gene-dense locus -- run FOCUS to identify the probable causal gene"
- "Triangulate TWAS, coloc, and cis-eQTL MR to nominate a candidate causal gene at my GWAS lead locus"
- "Run a multi-ancestry TWAS combining EUR and AFR sumstats with MA-FOCUS fine-mapping"
- "Restrict TWAS to brain tissues for a psychiatric trait using PsychENCODE weights"

## Example Prompts

### Single-tissue TWAS from GWAS summary statistics
> "Run S-PrediXcan on my LDL-cholesterol GWAS using the GTEx v8 MASHR liver model. Filter to genes with significant association at p < 2.3e-6 (Bonferroni for 22k genes) and report the top 20 by Z magnitude. Note the assumption that EUR-trained weights apply to an EUR-only GWAS."

### Multi-tissue prioritisation
> "Run S-PrediXcan across all 49 GTEx v8 MASHR tissues for my schizophrenia GWAS. Combine the per-tissue outputs with S-MultiXcan for a joint multi-tissue test. Compare the joint-significant gene set against the per-tissue Bonferroni-significant lists, and flag any genes where the joint test gains power over the best per-tissue test."

### TWAS fine-mapping at a gene-dense locus
> "At the chr19:11p15.5 locus my TWAS reports 7 genes passing significance. Run FOCUS using the GTEx v8 whole blood DB and 1000G EUR LD reference to compute per-gene PIPs and report the credible gene set. Genes with PIP >= 0.8 are causal candidates; co-significant genes with PIP < 0.5 are LD-tagged."

### Multi-ancestry TWAS
> "Run MA-FOCUS on my T2D GWAS combining EUR (BBJ + UKB), EAS (BBJ), and AFR (AAGILE) sumstats using ancestry-matched 1000 Genomes LD references. Use the MA-FOCUS DBs for each ancestry's PredictDB v8 MASHR adipose tissue. Report joint PIPs and per-ancestry contribution."

### Triangulation with MR and coloc
> "I have a TWAS hit for SORT1 in liver from S-PrediXcan. Triangulate this with (a) a cis-eQTL MR using SORT1 cis-eQTLs from GTEx liver as the exposure and my LDL GWAS as the outcome with TwoSampleMR, and (b) coloc.abf between the GTEx liver cis-eQTL for SORT1 and the GWAS at the SORT1 locus. Report a 4-way concordance summary across TWAS, FOCUS PIP, cis-MR, and coloc PP.H4."

### Splice-TWAS for a neuropsychiatric trait
> "Run TWAS using GTEx splicing models (sQTL-based PredictDB) instead of expression for my major depressive disorder GWAS in brain frontal cortex. Identify splice-mediated gene hits and contrast with the expression-TWAS hit list at the same loci."

### FUSION conditional analysis
> "Run FUSION on my CAD GWAS with GTEx artery coronary weights for all 22 autosomes. Then run FUSION.post_process.R --joint at every significant locus to identify conditionally independent genes; report the joint-Z table per locus."

## What the Agent Will Do

1. **Format inputs:** Harmonise GWAS sumstats to the expected schema for the chosen tool (FUSION wants SNP A1 A2 Z; S-PrediXcan wants explicit column mapping flags). Flip alleles if needed.
2. **Run per-tissue TWAS:** FUSION across chromosomes 1-22 OR S-PrediXcan per tissue.
3. **Aggregate / joint-test:** Combine multi-tissue outputs via S-MultiXcan or UTMOST if tissue prior is uncertain.
4. **Apply FOCUS / MA-FOCUS:** Probabilistic gene-level fine-mapping at each significant TWAS locus.
5. **Triangulate:** Cross-reference each candidate causal gene with cis-eQTL MR (causal-genomics/mendelian-randomization) and coloc (causal-genomics/colocalization-analysis). Report concordance.
6. **Flag caveats:** Tissue mismatch, ancestry mismatch, HLA region, low-N tissue panels, gene-dense LD-tied loci, FOCUS panel coverage gaps.

## Ancestry-Matched Prediction Panels

Use the panel matching the GWAS ancestry; document any mismatch and apply MA-FOCUS for cross-ancestry credible sets. Patel 2022 Genome Med 14:31 quantifies the ancestry-transfer power loss.

| Panel | N donors | Tissues / cells | Ancestry | Use case |
|-------|----------|-----------------|----------|----------|
| GTEx v8 (full) | 838 | 49 tissues | EUR (~85%) | Default standard for general TWAS |
| GTEx v8 MASHR-EUR | -- | 49 tissues | EUR | Primary; sparser SNP set per gene |
| eQTLGen | 31,684 | Whole blood | EUR (>95%) | Highest blood power; cis + trans available |
| MESA Whole Blood (Mogil 2018) | 1,163 | Blood | Multi-ancestry | AFR / HIS-relevant analyses |
| MESA Monocytes-AFA | 183 | Monocytes | AFR | AFR-specific immune traits |
| AFGR | ~2,000 | Whole blood | AFR | AFR (emerging; release-dependent) |
| OneK1K (Yazar 2022) | 982 | PBMC, 14 cell types | EUR | sc-TWAS in immune cell types |
| PsychENCODE (Wang 2018) | ~1,300 | Prefrontal cortex | EUR | Neuropsychiatric traits |
| BrainSeq Phase 2 | ~350 | DLPFC, hippocampus | EUR + AFR | Neuropsychiatric replication |

## GTEx v8 Tissues with N < 100 (Skip Unless Essential)

Below ~100 donors the per-gene elastic-net cross-validation R^2 is unstable and the FUSION heritability filter drops many genes. Skip these unless biologically required:

| Tissue | GTEx v8 N |
|--------|-----------|
| Kidney Medulla | 4 |
| Cervix - Endocervix | 10 |
| Cervix - Ectocervix | 9 |
| Fallopian Tube | 9 |
| Bladder | 21 |
| Brain - Substantia nigra | 139 (marginal; use with caution) |
| Brain - Spinal cord (cervical c-1) | 159 (marginal) |

For brain analyses requiring sub-region resolution at low N, substitute PsychENCODE (~1,300) or BrainSeq.

## Tips

- **Tissue selection:** Pick the tissue based on biology (a priori) and verify with LDSC-SEG / CELL-TYPE-SPECIFIC LDSC on the GWAS sumstats. Do not let the TWAS pick its tissue post-hoc -- this inflates false positives.
- **FOCUS is mandatory at gene-dense loci:** Always run FOCUS after FUSION/PrediXcan. Reporting raw TWAS hits without FOCUS PIPs misleads readers about which co-significant gene is causal.
- **Ancestry-match weights:** Use GTEx EUR weights for EUR GWAS; MESA / eQTLGen-Asian / AFGR for non-EUR. If unavoidable, document the mismatch and use MA-FOCUS for cross-ancestry credible sets.
- **HLA caveat:** chr6:25-35 Mb (hg38) violates the gene-by-gene LD assumption. Exclude from genome-wide TWAS summaries by default; analyse HLA classical alleles separately.
- **Low-N tissue panels:** Skip GTEx tissues with N < 100 donors unless biologically essential; per-gene weights overfit and produce inflated TWAS Z under the null.
- **Triangulation rule:** A "strong candidate causal gene" requires 3-of-4 concordance across {TWAS, FOCUS PIP >= 0.8, coloc PP.H4 >= 0.7, cis-eQTL MR p < threshold}. 2-of-4 is suggestive; 1-of-4 is associational only.
- **Drug-target framing:** For drug-target prioritisation, cis-MR with a colocalised eQTL is the canonical evidence; TWAS adds the gene-by-tissue power but should not replace MR + coloc.
- **Splice-TWAS is underutilised:** sQTL-based PredictDB models exist for GTEx v8 and recover splicing-mediated GWAS hits missed by expression-only TWAS, especially in neuropsychiatric and immune traits.
- **Conditional analysis:** FUSION's `--joint` is a lighter-weight alternative to FOCUS when a full DB build is unavailable; it does not return PIPs but identifies conditionally independent gene signals.
- **Document version pinning:** PredictDB model versions and FUSION weight panel dates matter; cite the exact version (e.g. "GTEx v8 MASHR-EUR, PredictDB release 2022-01") in methods.

## Related Skills

- causal-genomics/fine-mapping - Variant-level credible sets are the precursor to FOCUS gene-level fine-mapping
- causal-genomics/colocalization-analysis - Coloc PP.H4 triangulation with TWAS hits
- causal-genomics/mendelian-randomization - cis-eQTL MR triangulation; drug-target prioritisation
- causal-genomics/mediation-analysis - Downstream gene-mediated trait effects given TWAS hits
- population-genetics/association-testing - Upstream GWAS summary statistic generation
- population-genetics/linkage-disequilibrium - LD reference panel construction
- differential-expression/deseq2-basics - eQTL count data for custom prediction-weight training
- single-cell/preprocessing - Cell-type-resolved eQTL panels for sc-TWAS
- workflows/gwas-pipeline - Upstream GWAS pipeline producing TWAS input
- variant-calling/variant-annotation - Functional annotation of TWAS / FOCUS top variants
