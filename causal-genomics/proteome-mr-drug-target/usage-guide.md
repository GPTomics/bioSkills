# Proteome-Wide Drug-Target Mendelian Randomization - Usage Guide

## Overview

Cis-pQTL Mendelian randomization is the genetic-instrument analogue of a placebo-controlled drug trial: a genetic variant near the gene encoding a target protein perturbs that protein's plasma level from conception, and the variant's effect on a downstream phenotype estimates the causal effect of pharmacological inhibition of that target. The PCSK9-LDL-cholesterol-CAD story is the canonical positive example: cis-MR at the PCSK9 locus correctly predicted the direction and magnitude of clinical effect of evolocumab and alirocumab years before the FOURIER and ODYSSEY trials. The PCSK9 -> type 2 diabetes signal (Schmidt 2017 Lancet Diabetes Endocrinol) was the canonical on-target adverse-effect discovery.

This skill operationalises the Schmidt 2020 cis-MR framework (Nat Commun 11:3255) for plasma-proteome drug-target validation. It handles cis-pQTL instrument selection from UKB-PPP (Olink antibody), deCODE / Fenland / INTERVAL / ARIC (SomaScan aptamer), and FinnGen-PPP; cis-IVW with correlated instruments; colocalization triangulation; phenome-wide cis-MR for on-target adverse-effect scanning; cross-platform replication; and PAV (protein-altering-variant) sensitivity. The agent enforces the triangulation rule that a single significant cis-MR is necessary but never sufficient for a drug-target nomination.

## Prerequisites

```r
install.packages(c('remotes', 'MendelianRandomization', 'coloc', 'susieR', 'dplyr'))
remotes::install_github('MRCIEU/TwoSampleMR')
remotes::install_github('MRCIEU/ieugwasr')
remotes::install_github('MRCIEU/genetics.binaRies')   # bundles plink2 binary
remotes::install_github('rondolab/MR-PRESSO')
```

```bash
# Ensembl VEP for PAV annotation
conda install -c bioconda ensembl-vep
vep_install -a cf -s homo_sapiens -y GRCh38 -c $HOME/.vep

# 1000 Genomes EUR plink reference (for local clumping and LD matrix)
# https://mrcieu.github.io/ieugwasr/articles/local_ld.html
```

pQTL summary statistics:
- UKB-PPP Olink (Sun 2023): ukb-ppp.gwas.eu portal (download by protein UniProt ID)
- deCODE SomaScan (Ferkingstad 2021): https://www.decode.com/summarydata/
- Fenland SomaScan (Pietzner 2021): EBI GWAS catalog
- FinnGen-PPP Olink: FinnGen DF12 release portal

## Quick Start

Tell the AI agent what to do:
- "Run a cis-MR of PCSK9 on coronary artery disease using UKB-PPP cis-pQTLs"
- "Test IL6R inhibition on rheumatoid arthritis with cis-pQTL instruments from deCODE plus coloc triangulation"
- "Run a phenome-wide cis-MR of PCSK9 across OpenGWAS to find on-target adverse effects"
- "Replicate this cis-MR target across UKB-PPP Olink and deCODE SomaScan and flag PAV-confounded instruments"
- "Test a target with two independent cis-signals using coloc.susie per credible set"
- "Apply the Burgess 2016 sample-overlap correction to a UKB-PPP exposure with a UKB phenotype outcome"

## Example Prompts

### Single Drug Target, Single Outcome

> "I have UKB-PPP cis-pQTLs for PCSK9 within +/-500 kb of the gene and CARDIoGRAMplusC4D CAD summary stats. Run cis-IVW with weak-IV filtering, harmonise with action=2, triangulate with coloc.abf at p12=5e-6, and report PP.H4 plus PAV-excluded sensitivity. Annotate every cis-pQTL with VEP first."

> "Test IL6R protein -> rheumatoid arthritis using cis-pQTLs from deCODE SomaScan plus colocalization. Flag any cis-pQTL coloc'd with neighbouring genes' eQTLs in GTEx whole blood."

### Phenome-Wide On-Target Adverse-Effect Scan

> "Hold the PCSK9 cis-pQTL instrument set fixed and run cis-MR against all OpenGWAS outcomes with sample size >= 50,000 in European-ancestry cohorts. Bonferroni-correct over outcomes and report a phewas-style forest plot of significant hits."

> "Screen IL23R cis-pQTL effects across all FinnGen DF12 disease endpoints for on-target adverse-effect discovery prior to a clinical-trial protocol."

### Target with Multiple Independent Cis-Signals (Allelic Heterogeneity)

> "ANGPTL3 has two independent cis-pQTLs in low LD. Run coloc.susie on the cis-window with in-sample LD, report per-credible-set PP.H4 against triglycerides GWAS, and run a Wald ratio per credible set."

### Cross-Platform Replication

> "Run the same cis-MR of target X on disease Y independently in UKB-PPP (Olink) and deCODE (SomaScan). Report direction agreement, magnitude ratio, and flag the protein as platform-discordant if direction disagrees."

### Sample-Overlap Correction

> "Both my exposure GWAS (UKB-PPP cis-pQTL for protein X) and outcome GWAS (UKB HES-derived phenotype) are from UK Biobank. Apply MR-RAPS with the one-sample-equivalent treatment, or switch the outcome to FinnGen for an independent-cohort test."

### Drug Repurposing / Target Nomination

> "Cross-reference cis-MR estimates with Open Targets L2G scores and the Open Targets Drug platform to nominate druggable proteins for an autoimmune indication. Require cis-MR P < 1.7e-5, coloc PP.H4 >= 0.7, and cross-platform replication."

## What the Agent Will Do

1. Identify the target gene's coordinates (hg38), define the +/-500 kb cis-window
2. Download or load cis-pQTL summary stats for the target protein from UKB-PPP (Olink) and ideally also deCODE / Fenland (SomaScan)
3. Filter to genome-wide-significant cis-SNPs (P < 5e-8) inside the cis-window
4. Compute per-SNP F-statistic from the EXPOSURE GWAS (not outcome); drop F < 10
5. Annotate every retained cis-SNP with Ensembl VEP; tag PAVs (missense/nonsense/splice/protein_altering)
6. LD-clump within the window at r2 < 0.1 using local plink + 1KG (ancestry-matched)
7. Harmonise with the outcome GWAS via `harmonise_data(action = 2)`; document palindromic drops
8. Run cis-MR panel: Wald ratio (sentinel) + cis-IVW (random) + Egger (if >=10 SNPs) + weighted median + MR-PRESSO outlier test
9. Run cis-IVW with `correl = TRUE` if instruments are in moderate LD; supply the in-window LD matrix from plink
10. Run coloc.abf (or coloc.susie if multi-causal evidence) on the full cis-window with p12 = 5e-6; report PP.H4 + sensitivity
11. For each cis-pQTL, run coloc against non-target genes' eQTL/pQTL in the same window (drop if PP.H4 >= 0.5 with neighbour)
12. Run PAV-excluded sensitivity: re-run cis-MR after dropping all PAV cis-SNPs; report both
13. If cross-platform replication: repeat steps 2-11 on the alternate platform; report agreement
14. If sample overlap (UKB-on-UKB): apply Burgess 2016 correction or switch outcome to FinnGen
15. For pheWAS: loop step 7-8 over the OpenGWAS outcome catalogue; Bonferroni-correct
16. Produce STROBE-MR-compliant report listing: cis-MR estimate + 95% CI, PP.H4, PAV-excluded estimate, platform agreement, neighbour-gene coloc results, sample-overlap statement, claim-strength ladder

## Tips

- The exclusion-restriction in cis-MR is relaxed (the protein product directly mediates) but NOT eliminated; neighbour-gene mediation is the dominant residual concern
- Always VEP-annotate cis-pQTLs and run a PAV-excluded sensitivity panel; PAVs can produce apparent platform-specific pQTLs without changing actual protein abundance
- Olink (antibody PEA) and SomaScan (aptamer SOMAmer) disagree on ~15-30% of proteins; replicate every clinically-actionable claim on both
- Cross-platform direction disagreement is a strong flag for measurement artifact; do not advance such targets to nomination
- Sample-overlap: treating UKB-PPP exposure + UKB outcome as two-sample is a common methodological error; bias is one-sample-equivalent, toward observational
- Cis-MR clumping uses r2 < 0.1 within window (not r2 < 0.001 polygenic-MR convention); the window itself is the LD-pruning mechanism
- Egger needs >=10 cis-pQTLs to be powered; most cis-windows have fewer, so cis-MR-Egger is usually exploratory only
- Phenome-wide cis-MR is the canonical use case for discovering on-target adverse effects (PCSK9 -> T2D is the textbook example)
- For correlated cis-pQTLs (r2 0.1-0.7), use `MendelianRandomization::mr_ivw(correl = TRUE, correl.x = ld_matrix)`; supply the matrix from local plink against an ancestry-matched reference
- coloc PP.H4 >= 0.7 is the Open Targets clinical bar; PP.H4 >= 0.5 alone is exploratory; require PP.H3 to be substantially less than PP.H4 for confidence
- The Steiger filter has a known caveat under unmeasured confounding (Hemani Tilling 2022); cross-validate direction with bidirectional cis-MR
- L2G score (Mountjoy 2021) plus cis-MR plus coloc plus existing-drug evidence is the modern Open Targets target-prioritization stack
- Trans-pQTLs (outside +/-500 kb) cannot serve as primary instruments; their inclusion violates exclusion-restriction by definition
- Some proteins (e.g. complement factors, immunoglobulins) have many PAVs in cis; the PAV-excluded sensitivity may discard most instruments and the cis-MR claim cannot be made
- The cis-window default is +/-500 kb (Schmidt 2020); some pipelines use 1 Mb. Pre-specify and document the window in methods

## Related Skills

causal-genomics/mendelian-randomization - Parent polygenic-MR framework; cis-MR is the drug-target specialization
causal-genomics/colocalization-analysis - Required PP.H4 triangulation for any cis-MR drug-target claim
causal-genomics/fine-mapping - Credible-set construction prior to coloc.susie at the cis-locus
causal-genomics/pleiotropy-detection - MR-PRESSO / Egger diagnostics adapted to cis-window
causal-genomics/transcriptome-wide-association - eQTL-based parallel evidence for the same target
causal-genomics/mediation-analysis - Step from cis-MR to downstream mediator pathway
population-genetics/association-testing - Source GWAS pipelines for pQTL discovery
population-genetics/linkage-disequilibrium - LD-matrix construction for cis-IVW correl=TRUE
variant-calling/variant-annotation - VEP PAV annotation for sensitivity analysis
clinical-databases/clinvar-lookup - Pathogenic-variant context for nominated targets
