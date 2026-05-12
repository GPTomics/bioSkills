# Colocalization Analysis - Usage Guide

## Overview

Colocalization tests whether two or more association signals at a genomic locus are driven by the same causal variant. It is the standard tool for integrating GWAS hits with molecular QTLs (eQTL, sQTL, pQTL, mQTL), prioritising candidate genes, and distinguishing shared causality from coincidental overlap due to linkage disequilibrium. Methods range from the fast single-causal coloc.abf (Giambartolomei 2014) to the multi-causal coloc.susie (Wallace 2021), the multi-trait HyPrColoc (Foley 2021), and the pleiotropy-vs-linkage SMR/HEIDI (Zhu 2016).

The agent will pick the appropriate method based on the experimental scenario: number of traits, expected number of causal variants per locus, availability of LD reference, ancestry composition, and whether the question is causal-mediation (SMR) or shared-variant (coloc). It will harmonise allele coding across summary statistics, format inputs correctly (`type`, `sdY`, `s`, `N`), run colocalization, and ALWAYS report sensitivity over the p12 prior alongside the headline PP.H4.

## Prerequisites

```r
install.packages(c('coloc', 'susieR', 'ggplot2', 'patchwork', 'data.table'))
install.packages('remotes')
remotes::install_github('jrs95/hyprcoloc')          # GitHub-only, never CRAN

# Optional, for multi-omic
remotes::install_github('clagiamba/moloc')           # 3-5 traits

# CLI tools (outside R)
# SMR: download binary from https://yanglab.westlake.edu.cn/software/smr/
# eCAVIAR: compile from https://github.com/fhormoz/caviar
# PWCoCo: compile from https://github.com/jwr-git/pwcoco
# SharePro_coloc: pip install sharepro-coloc

# LD reference panel (e.g. 1000 Genomes phase 3) via plink2
# plink2 --pfile 1KG_EUR_chr6 --chr 6 --from-bp 30000000 --to-bp 31000000 --r-phased square --out ld
```

## Quick Start

Tell the AI agent what is needed:
- "Test if my GWAS lead SNP at chr6:30500000 colocalizes with the IL2 eQTL in whole blood"
- "Run coloc.susie at this locus; I have two independent GWAS signals after conditional analysis"
- "Check colocalization between my GWAS and all 49 GTEx tissues using HyPrColoc"
- "I have a GWAS + eQTL + sQTL + mQTL at the same locus; run moloc"
- "Harmonise these summary stats and run coloc with sensitivity over p12"
- "Run SMR + HEIDI between my GWAS and eQTLGen blood eQTL"
- "The locus is in the MHC -- what method should be used"
- "My GWAS is in East Asian ancestry but the eQTL is GTEx EUR -- can coloc be trusted"

## Example Prompts

### Single-Tissue GWAS-eQTL

> "Take this 1 Mb window centred on the GWAS lead SNP and run coloc.abf against the eQTL for the nearest gene. Report PP.H4 with p12 sensitivity."

> "I have summary stats with p-values and MAF but no betas -- run coloc using p-value / MAF input format."

### Multi-Causal / Allelic Heterogeneity

> "GCTA-COJO conditional analysis identified two independent GWAS signals at this locus. Run coloc.susie with the LD matrix and report all credible-set pair PPs."

> "PP.H3 is dominating coloc.abf despite obvious overlap in LocusZoom. Switch to coloc.susie or eCAVIAR and re-test."

### Multi-Tissue

> "Run coloc.abf between this GWAS locus and the same gene across all 49 GTEx v8 tissues; rank tissues by PP.H4 to identify causal cell type."

> "Use HyPrColoc across GTEx tissues to cluster tissues that share the causal variant."

### Multi-Trait / Multi-Omic

> "Integrate GWAS + eQTL + sQTL + pQTL at this locus with moloc and report the PPA for the all-share hypothesis."

> "Cluster 12 cardiovascular GWAS traits at the LDLR locus with HyPrColoc."

### SMR / HEIDI

> "Run SMR + HEIDI between my disease GWAS .ma file and eQTLGen .besd. Report SMR p, HEIDI p, and number of HEIDI SNPs."

### MHC / HLA

> "This locus is in the MHC. Flag the long-range-LD problem and recommend HLA-coloc or exclude-and-report-HLA-association strategy."

### Ancestry-Mismatched

> "GWAS is FinnGen (FIN) and eQTL is GTEx (EUR). Compute z-score vs LD diagnostics with estimate_s_rss; if lambda > 0.05, switch to coloc.abf or SharePro_coloc."

### Visualization

> "Make a LocusCompare plot and a stacked regional Manhattan for this colocalization."

## What the Agent Will Do

1. Identify the question type (single-causal, multi-causal, multi-trait, pleiotropy vs linkage) and select method from the decision tree.
2. Extract a +/- 500 kb to 1 Mb window centred on the GWAS lead SNP (or joint top-variant when both traits available).
3. Harmonise allele coding between summary statistics; flip betas where A1/A2 swap; drop palindromic SNPs at high MAF.
4. Format the per-trait coloc input lists (beta, varbeta, snp, position, N, type, sdY for quant, s for cc).
5. Run the chosen colocalization method (coloc.abf, coloc.susie, HyPrColoc, moloc, SMR, eCAVIAR, PWCoCo, or SharePro_coloc).
6. For coloc.susie: run `susieR::estimate_s_rss` to verify z-score vs LD consistency; abort if lambda > 0.05.
7. ALWAYS run `coloc::sensitivity(res, 'H4 > 0.75')` and report the p12 range over which PP.H4 stays above threshold.
8. Interpret PP.H0-H4 against the appropriate threshold (>= 0.75 screening, >= 0.80 published, >= 0.90 stringent).
9. Generate regional association and LocusCompare plots colored by LD to lead.
10. For multi-causal results, report per-(CS1, CS2) PP and lead SNP per credible set.
11. For multi-omic results, report the all-share PPA and per-pair PPs.
12. Flag known failure modes (MHC, ancestry mismatch, low-N eQTL, lead-SNP-swap window bias) explicitly in the report.

## Plain-Language H0-H4 (for Methods Section)

For a methods or supplementary description in plain prose:

- **H0** -- Neither trait has an association signal at this locus. Posterior reflects "no signal anywhere".
- **H1** -- Only trait 1 has a causal variant in the window; trait 2 has no signal.
- **H2** -- Only trait 2 has a causal variant in the window; trait 1 has no signal.
- **H3** -- Both traits have causal variants in the window, but they are different SNPs (linkage / coincidence under LD).
- **H4** -- Both traits share a single causal variant in the window (colocalization).

A high PP.H4 (>= 0.75) supports a shared-causal-variant interpretation. A high PP.H3 supports distinct causal variants in linkage. PP.H0 / PP.H1 / PP.H2 indicate the locus is underpowered for at least one trait. Always report all five posteriors, not PP.H4 alone.

## GWAS Summary-Stats Harmonisation (Worked Example)

```r
harmonise <- function(df1, df2) {
    m <- merge(df1, df2, by='SNP', suffixes=c('.1','.2'))
    same <- m$A1.1 == m$A1.2 & m$A2.1 == m$A2.2
    flip <- m$A1.1 == m$A2.2 & m$A2.1 == m$A1.2
    palindromic <- (m$A1.1 %in% c('A','T') & m$A2.1 %in% c('A','T')) |
                   (m$A1.1 %in% c('C','G') & m$A2.1 %in% c('C','G'))
    m$BETA.2[flip] <- -m$BETA.2[flip]
    m$MAF.2[flip] <- 1 - m$MAF.2[flip]
    keep <- (same | flip) & !(palindromic & m$MAF.1 > 0.42)
    m[keep, ]
}
```

Harmonisation pitfalls to watch for:

- **Allele coding mismatch.** GWAS may report effect allele as A1 while eQTL reports it as A2. Always check both and flip betas where needed.
- **Build mismatch.** hg19 GWAS coords + hg38 eQTL coords silently merge on rsID but break on chr:pos. Lift over with `rtracklayer::liftOver` or CrossMap before merging.
- **Palindromic SNPs at high MAF.** A/T and C/G SNPs at MAF > 0.42 cannot be unambiguously strand-resolved; drop them or resolve with reference-panel MAF.
- **Multi-allelic SNPs.** Many summary stats collapse multi-allelic loci by keeping only the most-frequent alt; if datasets pick different alts, harmonisation drops the SNP. Split on chr:pos:ref:alt as a unique key.
- **rsID dependence.** rsID can be remapped across dbSNP builds (e.g. merge of two rsIDs into one). Prefer chr:pos:ref:alt keys for cross-study merges.

## Worked PWCoCo Conditional Recipe

```bash
# Step 1: GCTA-COJO identifies independent signals at the locus
gcta64 --bfile 1KG_EUR --chr 6 --extract locus.snplist \
       --cojo-file gwas.ma --cojo-slct --out gwas_cojo

# gwas_cojo.jma.cojo lists independent signals (per --cojo-p 5e-8 default)

# Step 2: PWCoCo runs pairwise conditional coloc.abf per signal pair
pwcoco --bfile 1KG_EUR --sum_stats1 gwas.txt --sum_stats2 eqtl.txt \
       --p_cutoff1 5e-8 --p_cutoff2 5e-5 \
       --chr 6 --top_snp rs12345 \
       --out pwcoco_result
```

Output: one coloc.abf result per (conditional signal in trait 1, conditional signal in trait 2) pair. Interpret each row as a separate single-signal coloc test. If COJO finds 2 GWAS + 1 eQTL signals, expect 2 result rows. PWCoCo requires individual-level reference (plink bfile); cannot run on summary stats alone.

## Lead-SNP-Swap Operational Steps

To diagnose window-centring bias, re-run coloc.abf three times with different window centres:

1. **GWAS-centred window:** +/- 500 kb around the GWAS lead SNP.
2. **eQTL-centred window:** +/- 500 kb around the eQTL top SNP for the gene.
3. **Joint top-variant window:** +/- 500 kb around the SNP with the lowest min-p across both traits.

Report all three PP.H4 values. If they agree within 0.1, the result is stable. If they swing > 0.2, the locus is borderline and the report must list all three centrings. For multi-causal loci (allelic heterogeneity), the joint top-variant window typically gives the most-defensible result for coloc.abf; coloc.susie removes the centring sensitivity by construction.

## Tips

- **Default to coloc.abf + sensitivity** for first-pass screening; escalate to coloc.susie only when multi-causal is biologically plausible (conditional analysis identified independent signals).
- **PP.H3 + PP.H4 framing**: report both. PP.H4 = 0.6 with PP.H3 = 0.05 is qualitatively different from PP.H4 = 0.6 with PP.H3 = 0.3.
- **p12 sensitivity is non-negotiable**: every reported PP.H4 should be accompanied by the p12 range over which it stays above threshold. Reviewer-grade reporting expects this.
- **Trans-eQTL needs lower p12** (5e-6 or 1e-6) to avoid over-claiming biologically improbable shared causality.
- **Match LD ancestry to GWAS ancestry**. EUR LD on AFR z-scores produces spurious credible sets in coloc.susie.
- **MHC is special**: never report a standard coloc PP.H4 for MHC without the long-range-LD caveat. Use HLA-coloc on classical alleles instead.
- **N < 200 eQTL is underpowered for coloc**: report H0/H1/H2 dominance honestly; do not claim absence of colocalization.
- **sdY semantics**: provide `sdY=1` ONLY if the trait is standardised (e.g. inverse-normal-transformed eQTL expression). If unsure, omit and let coloc estimate from MAF + varbeta.
- **Palindromic SNPs**: drop A/T and C/G SNPs with MAF > 0.42; their strand cannot be inferred from coding alone.
- **Lead-SNP-swap diagnostic**: re-run with the window centred on the OTHER trait's lead. If PP.H4 swings > 0.2, the locus is borderline; report multiple centrings.
- **SMR vs coloc give different answers** to different questions. Run both for therapeutic-target prioritization; report agreement metrics.
- **Open Targets / FinnGen PP.H4 default is >= 0.75**, not 0.8. Use >= 0.7 as triangulation tier when PP.H4 is one of several lines of evidence (TWAS, cis-MR, effector-gene), >= 0.75 for screening, >= 0.80 for published stand-alone coloc claims, >= 0.90 for stringent clinical follow-up, >= 0.95 for industry / regulatory drug-target submission.

## Related Skills

- causal-genomics/mendelian-randomization - Downstream MR using colocalized SNPs as IVs
- causal-genomics/fine-mapping - SuSiE credible sets that feed coloc.susie
- causal-genomics/mediation-analysis - Causal mediation building on shared causal variants
- causal-genomics/pleiotropy-detection - Distinguishing horizontal pleiotropy from shared causality
- population-genetics/association-testing - GWAS summary stat generation
- population-genetics/linkage-disequilibrium - LD panel construction for coloc.susie / PWCoCo
- variant-calling/variant-annotation - Functional annotation for variant-specific priors
- single-cell/scatac-analysis - Per-cell-type chromatin context for coloc results
- differential-expression/deseq2-basics - eQTL count generation
- workflows/gwas-pipeline - Upstream GWAS producing coloc input
