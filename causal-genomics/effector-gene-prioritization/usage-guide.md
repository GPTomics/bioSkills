# Effector Gene Prioritization Skill Usage Guide

## Overview

Effector gene prioritization is the central interpretation step that bridges variant-level statistical fine-mapping and gene-level biological hypothesis at GWAS loci. The agent takes a GWAS lead variant (or credible set), interrogates a portfolio of variant-to-gene (V2G) methods (Open Targets L2G, MAGMA, FUMA, cS2G, PoPS, FUSION/S-PrediXcan TWAS, ABC, ENCODE-rE2G), and integrates orthogonal evidence streams (fine-mapping PIP, colocalization PP.H4, distance, chromatin enhancer-gene linking, polygenic similarity priors) into a concordance-based confidence tier.

The dominant failure mode is the nearest-gene assumption: roughly 30-50% of fine-mapped GWAS variants regulate a gene that is NOT the closest TSS. The L2G + PoPS combination plus enhancer-gene evidence (ABC / ENCODE-rE2G) plus eQTL/pQTL coloc is the modern publication-grade triangulation standard. No single tool is sufficient; concordance across >= 3 of 6 orthogonal evidence streams is the operational rule for a "high-confidence causal effector" call.

## Prerequisites

- GWAS summary statistics (SNP, CHR, BP, A1, A2, BETA or Z, SE, P, N)
- Fine-mapping output (SuSiE / FINEMAP credible sets) from causal-genomics/fine-mapping
- Colocalization output (coloc.abf or coloc.susie PP.H4) from causal-genomics/colocalization-analysis
- Ancestry-matched LD reference panel (e.g. 1000 Genomes EUR PLINK bfile)
- Pre-computed integrative scores when available: Open Targets L2G, cS2G, PoPS feature matrix
- For matched-tissue enhancer-gene linking: ATAC-seq + H3K27ac ChIP + Hi-C/Micro-C in the candidate causal tissue (cross-reference atac-seq/enhancer-gene-linking)
- Python 3.9-3.11 + R 4.3+ + PLINK 1.9 / 2.0 + MAGMA 1.10+

Install:

```bash
# MAGMA ships as a pre-compiled zip, not a git repo. Download it manually from the official CNCR page
# (https://cncr.nl/research/magma/) -- the old ctg.cncr.nl and vu.data.surfsara.nl direct links now
# redirect to an HTML landing page, so a hard-coded wget would save that page AS the zip and unzip
# would fail. Once the real zip is in the working directory, unzip it (filename tracks the version):
unzip magma_v*.zip

# PoPS
git clone https://github.com/FinucaneLab/pops
pip install -r pops/requirements.txt

# Open Targets Genetics (community wrapper)
pip install opentargets-genetics
# OR query GraphQL directly via requests
```

Pre-computed scores:

```bash
# cS2G
wget https://alkesgroup.broadinstitute.org/cS2G/cS2G_UKBB/cS2G.tsv.gz

# PoPS pre-built features
# See FinucaneLab/pops releases for PoPS_features_full.txt

# Open Targets L2G: queried via API; no local download
```

### eQTL / pQTL Panel Selection

The panel choice dominates coloc and TWAS sensitivity. Match to the causal tissue identified via LDSC-SEG or stratified LDSC before locking in.

| Panel | N donors | Tissues / cells | Use case |
|-------|----------|-----------------|----------|
| GTEx v8 | 838 | 49 tissues | PredictDB standard for 2026 TWAS pipelines |
| GTEx v10 (2024) | ~1,000 | 49 tissues | Limited PredictDB coverage; transitional |
| eQTLGen | 31,684 | Whole blood | Maximum statistical power for blood-relevant traits |
| eQTLGen Phase 2 (sc) | ~30,000 | sc immune | Cell-type-resolved blood eQTLs |
| OneK1K | 982 | PBMC, 14 cell types | sc-eQTL discovery for immune traits |
| INTERVAL pQTL | 3,301 | Plasma proteome | Drug-target cross-reference via pQTL-MR |

## Quick Start

Tell the AI agent what to do:
- "Prioritise effector genes at my GWAS lead locus using Open Targets L2G + PoPS + coloc concordance"
- "Run MAGMA gene-based and gene-set enrichment on my GWAS sumstats, then layer PoPS for polygenic priority"
- "I have a fine-mapped credible set; which gene is the most likely causal effector?"
- "Triangulate L2G, PoPS, coloc PP.H4, ABC, and distance to nominate a candidate causal gene at the SORT1 locus"
- "Run FUMA SNP2GENE on my GWAS and compare with Open Targets L2G calls for the top 20 hits"
- "Use ABC enhancer-gene predictions to assign a distal regulatory variant to its likely target gene"

## Example Prompts

### Open Targets covers the trait
> "My GWAS is on type 2 diabetes (UKB-2020). Query Open Targets Genetics for the L2G top-ranked genes at all genome-wide-significant loci. Cross-check each with V2G sub-scores and flag genes where yProbaDistance dominates yProbaMolecularQTL (distance-only candidates). Report the top 50 with L2G score >= 0.5."

### Custom trait from scratch
> "I have GWAS summary statistics for a custom rare-disease phenotype not in Open Targets. Run MAGMA gene-based with a 35kb upstream + 10kb downstream window using the 1000G EUR LD reference, then run PoPS on the MAGMA Z output with the pre-built feature matrix. Cross-reference with SuSiE fine-mapping credible sets and coloc PP.H4 against eQTLGen whole blood. Report concordance per gene."

### Tissue-known prioritisation
> "For my LDL-cholesterol GWAS, the causal tissue is liver. Run S-PrediXcan with GTEx liver MASHR weights and coloc.susie with GTEx liver eQTL panel. Layer ABC enhancer-gene predictions from HepG2 (ENCODE) for distal regulation. Integrate with Open Targets L2G and report effector-gene candidates with >= 3 concordant evidence streams."

### Tissue-unknown prioritisation
> "I have a schizophrenia GWAS but the causal cell type within brain is unclear. Run LDSC-SEG to prioritise brain tissues, then S-MultiXcan across all GTEx brain sub-regions, then PoPS for polygenic priority. Report per-gene concordance across L2G, PoPS, and tissue-prioritised coloc."

### Distal regulation suspected
> "At my chr8:9p21 GWAS lead, the nearest gene is CDKN2A but I suspect long-range regulation. Run ABC and ENCODE-rE2G in matched cell type (CMK or vascular smooth muscle), cross-check with Open Targets L2G, and report whether the fine-mapped credible variant maps to CDKN2A or a distal gene (e.g. MTAP, ANRIL)."

### Publication-grade triangulation
> "I am writing up an effector-gene nomination for ANGPTL4 at a triglycerides GWAS lead. Triangulate (a) fine-mapping (SuSiE PIP), (b) coloc PP.H4 with GTEx subcutaneous adipose eQTL, (c) Open Targets L2G score, (d) PoPS score, (e) ABC enhancer-gene linking in adipocytes (Engreitz portal), (f) distance to TSS. Report concordance per locus and flag whether ANGPTL4 meets the >= 3 of 6 high-confidence threshold."

### Multi-effector locus
> "At my chr11p15.5 lipid GWAS locus, three genes (CLU, NCAM1, MTHFD1L) all show modest L2G scores. Run conditional analysis (FUSION --joint or GCTA-COJO) to test independence, then per-gene coloc.susie at each independent signal, and report whether the locus is multi-effector or LD-tied."

## What the Agent Will Do

1. **Gather upstream evidence**: Fetch fine-mapping credible sets (causal-genomics/fine-mapping), coloc PP.H4 (causal-genomics/colocalization-analysis), and TWAS hits (causal-genomics/transcriptome-wide-association) for each GWAS lead.
2. **Run MAGMA**: Gene-based annotation with FUMA-default window (35kb upstream + 10kb downstream); gene-set enrichment against MSigDB or curated gene sets.
3. **Query L2G**: Open Targets GraphQL for studyLocusL2G; record per-gene `yProbaModel` and per-component sub-scores.
4. **Run PoPS**: LASSO on MAGMA gene Z with curated gene-feature matrix; record per-gene polygenic priority score and top-decile flags.
5. **Cross-check cS2G**: Lookup pre-computed cS2G per-SNP gene allocations for credible-set variants.
6. **Layer enhancer-gene** (if matched epigenome): ABC and ENCODE-rE2G predictions for fine-mapped credible variants (cross-reference atac-seq/enhancer-gene-linking).
7. **Score concordance**: Per (locus, gene) candidate, count passing evidence streams from {fine-mapping, coloc, distance, PoPS, L2G, ABC/rE2G}. Tier as high-confidence (>= 3), strong (>= 4), or near-certain (>= 5).
8. **Reconcile disagreement**: Flag loci where L2G and PoPS top genes diverge; report both with caveats; recommend follow-up (CRISPRi-FlowFISH at top candidates).
9. **Report caveats**: Tissue choice, ancestry of weights, HLA exclusion, panel coverage limits, single-locus multi-effector signatures.

## Tips

- **Nearest-gene is wrong 30-50% of the time**: Never report the nearest gene as causal without checking distal regulation. L2G and PoPS are designed to overcome this; use them.
- **Multi-evidence concordance is the operational standard**: Require >= 3 of 6 streams. Two streams is suggestive; one stream is associational only. The 6 streams are fine-mapping, coloc, distance, PoPS, L2G, and ABC/rE2G (chromatin enhancer-gene).
- **L2G and PoPS are orthogonal**: L2G uses per-locus features; PoPS uses similarity-based features. Concordance across both is a strong signal; disagreement is informative, not failure.
- **Tissue prioritisation precedes effector-gene assignment**: Identify the causal tissue with LDSC-SEG / CELLEX / EWCE before locking on a single eQTL panel. Wrong-tissue eQTL coloc produces spurious "causal genes" that are LD-tagged in unrelated tissues.
- **MAGMA window choice is a methodological lever**: 35kb upstream + 10kb downstream is FUMA default; 0kb+0kb is conservative; 50kb+50kb is liberal. Run sensitivity over windows for high-stakes reports.
- **HLA is a no-fly zone**: chr6:25-35 Mb (hg38) breaks every gene-by-gene method. Exclude from genome-wide effector-gene summaries; analyse HLA classical alleles separately.
- **Multi-effector loci are real**: 5-10% of GWAS loci have multiple causal genes. Allow multi-gene reporting; test independence with conditional analysis (FUSION --joint, GCTA-COJO).
- **cS2G is a per-SNP heritability-calibrated aggregator**: useful as a cross-check on L2G but not a replacement; cS2G calibrates against heritability enrichment, L2G calibrates against curated gold-standard genes.
- **Functional validation is the gold standard**: CRISPRi-FlowFISH (Fulco 2019) is the experimental ground truth. For drug-target nominations, computational concordance is necessary but not sufficient.
- **Document panel versions**: Open Targets release, PredictDB version, ABC release, ENCODE-rE2G model cell type. Effector-gene calls change across releases; cite the exact versions.

## Related Skills

- causal-genomics/fine-mapping - Variant-level credible sets feeding L2G and concordance scoring
- causal-genomics/colocalization-analysis - coloc PP.H4 as one of the six evidence streams
- causal-genomics/transcriptome-wide-association - Gene-level association and FOCUS gene fine-mapping
- causal-genomics/mendelian-randomization - cis-eQTL MR as orthogonal causal evidence
- causal-genomics/mediation-analysis - Downstream gene-mediated trait effects
- causal-genomics/proteome-mr-drug-target - pQTL-based drug-target nomination
- atac-seq/enhancer-gene-linking - ABC and ENCODE-rE2G enhancer-gene predictions
- atac-seq/atac-peak-calling - Enhancer candidates in matched tissue
- gene-regulatory-networks/coexpression-networks - Gene co-expression features feeding PoPS
- gene-regulatory-networks/scenic-regulons - TF-target regulons as supporting evidence
- pathway-analysis/go-enrichment - Pathway context for candidate effector genes
- population-genetics/association-testing - Upstream GWAS summary-statistic generation
- variant-calling/variant-annotation - Coding-consequence annotation
- workflows/gwas-pipeline - End-to-end GWAS pipeline producing input
