# Biodiversity Metrics Usage Guide

## Overview

Quantifies biodiversity from species abundance and incidence tables using the Hill-number framework (richness, Shannon, Simpson reported as effective species counts rather than raw entropies), coverage-based rather than sample-size-based rarefaction-extrapolation (the postdoc-grade default per Chao & Jost 2012), asymptotic richness via Chao1/ACE/jackknife with explicit acknowledgement that Chao1 is a LOWER BOUND not a point estimate, beta diversity decomposition into turnover vs nestedness components in the Baselga framework with the Podani alternative as a sensitivity check, mandatory Hellinger transformation prior to ordination (solves the double-zero problem), phylogenetic diversity (Faith PD, MPD, MNTD with SES_MPD/MNTD requiring explicit null-model justification), and functional diversity with Maire 2015 trait-axis dimensionality optimization to avoid FRic inflation from collinear traits.

## Prerequisites

- R with iNEXT (`install.packages('iNEXT')`)
- iNEXT.3D for phylogenetic and functional diversity Hill numbers (`install.packages('iNEXT.3D')`)
- vegan for ordination, Hellinger transform, and classic indices (`install.packages('vegan')`)
- betapart for Baselga turnover/nestedness partitioning (`install.packages('betapart')`)
- picante for SES_MPD, SES_MNTD, Faith's PD (`install.packages('picante')`)
- mFD for trait-space dimensionality optimization (`install.packages('mFD')`)
- ggplot2 for visualization (`install.packages('ggplot2')`)

## Quick Start

Tell your AI agent what you want to do:

- "Compare species diversity across sites using Hill numbers and coverage-based rarefaction"
- "Build rarefaction/extrapolation curves bounded by the doubling rule"
- "Decompose beta diversity into turnover and nestedness using both Baselga and Podani frameworks"
- "Estimate asymptotic richness with Chao1, but check whether singletons are real biology first"
- "Compute SES_MPD with explicit null-model justification for phylogenetic clustering"
- "Optimize trait-axis dimensionality with Maire 2015 mAD before computing FRic"

## Example Prompts

### Alpha Diversity with Hill Numbers

> "I have species abundance data from 12 sites. Compute Hill numbers q=0,1,2 with iNEXT coverage-based rarefaction-extrapolation, bound extrapolation at 2x the reference sample size, and report results as effective species counts at 95% coverage."

> "Estimate the lower bound on richness for each site using Chao1, but first check that singletons are biologically real (not PCR-error artifacts). If amplicon data, report Chao2 from incidence instead."

### Coverage-Standardized Comparison

> "My eDNA samples vary 10-fold in read depth. Use coverage-based rarefaction in iNEXT to standardize diversity to 95% coverage across all sites. Report the effective species counts as Hill numbers, not raw Shannon."

### Beta Diversity Partition

> "Decompose pairwise Sorensen beta diversity into turnover and nestedness using the Baselga partition in betapart. Then run the Podani alternative partition and report whether the two frameworks agree on the dominant component."

> "Calculate multi-site beta diversity for my sampling region. Report whether turnover or nestedness dominates AND acknowledge the Baselga partition is not mathematically unique."

### Phylogenetic Diversity

> "Calculate Faith's PD for each community against an ultrametric tree. Then compute SES_MPD with the 'taxa.labels' null model and justify the null-model choice."

> "Test phylogenetic clustering versus overdispersion using SES_MPD AND SES_MNTD. Do not interpret 'clustering implies environmental filtering' from sign alone; test trait conservatism (Blomberg's K) first per Mayfield & Levine 2010."

### Functional Diversity

> "I have species traits and a community matrix. Run Maire 2015 mAD optimization to choose the number of PCoA axes, then compute FRic, FEve, FDiv with the optimal k."

### Compositional Data Considerations

> "My data are amplicon ASV read counts. Apply Hill-number diversity OR switch to presence-absence Chao2 metrics that are robust to compositional distortion. Document which choice was made."

## What the Agent Will Do

1. Inspect the input data type (raw abundance, ASV/OTU counts, incidence frequencies, sample-by-species matrix) and decide whether singletons are biologically real or PCR/sequencing-error artifacts
2. Apply Hellinger transformation when ordination follows; use raw counts only when feeding directly into iNEXT
3. Compute Hill numbers at q = 0, 1, 2 with iNEXT coverage-based rarefaction-extrapolation, bound extrapolation at 2x the reference sample size, and report as effective species counts (numbers-equivalent of Shannon/Simpson)
4. Choose the appropriate asymptotic richness estimator: Chao1 if f2 > 0 and singletons are real; Chao1bc or jackknife1 if f2 = 0; ACE if many rare classes; Chao2 from incidence for amplicon data with unreliable singletons
5. Report Good's coverage alongside richness estimators to indicate sampling adequacy
6. Decompose beta diversity using both Baselga (turnover + nestedness) AND a richness-difference alternative (Podani/Carvalho), acknowledging the partition is not unique
7. Compute Faith's PD and SES_MPD/SES_MNTD with explicit, documented null-model choice and avoid the "clustering implies filtering" interpretation trap
8. Optimize trait-axis dimensionality via Maire 2015 mAD before computing FRic/FEve/FDiv
9. Flag the compositional-data caveat if data come from amplicon sequencing and recommend either restricting to presence-absence metrics or applying a CLR-based diversity framework
10. Produce publication-ready rarefaction curves with bootstrap CIs (nboot >= 200)

## Tips

- Coverage-based rarefaction (iNEXT type 3) is the postdoc-grade default; sample-size rarefaction is now considered an anti-pattern for cross-site comparison when assemblages differ in true diversity (Chao & Jost 2012)
- Hill numbers at q = 0 (richness), q = 1 (Shannon-equivalent exp(H)), q = 2 (Simpson-equivalent 1/D) form a unified parametric family; report all three to reveal sensitivity to rare-vs-common species
- Chao1 is a LOWER BOUND under a Gamma-Poisson model, NOT a point estimate; report it with "at least" phrasing and Good's coverage to indicate informativeness
- For amplicon/eDNA data, the standard practice is to either skip Chao1 (Callahan 2017 ASV philosophy) or denoise aggressively and report observed ASV count plus Good's coverage; raw Chao1 on un-denoised data measures PCR error, not biology
- The doubling rule: iNEXT extrapolation is reliable only up to 2x the reference sample size; the default endpoint enforces this, do not override above
- Bootstrap CIs need nboot >= 200 for publication; default nboot = 50 is too few
- Hellinger transformation (`vegan::decostand(matrix, method = 'hellinger')`) before PCA/RDA on community data is mandatory; without it, ordination is dominated by sample-size effects and the double-zero problem
- Baselga and Podani beta-diversity partitions give different ecological interpretations of the same data; report both and acknowledge the partition is not mathematically unique
- For SES_MPD/SES_MNTD, the null-model choice dominates the result: `taxa.labels` randomizes which species are present holding sample richness constant; `independentswap` preserves both row and column marginals; `richness` randomizes within sample; document the choice explicitly
- Mayfield & Levine 2010 showed phylogenetic clustering can result from competition (not just environmental filtering) when traits track phylogeny; do not infer process from SES sign alone, test trait conservatism with Blomberg's K and Pagel's lambda
- Functional richness (FRic) inflates with correlated trait axes; use `mFD::quality.fspaces()` to choose the Maire 2015 mAD-optimal axis count
- Bray-Curtis is not a true metric (violates triangle inequality); use Sorensen instead when downstream methods require metricity (e.g., Ward's clustering, complete linkage on transformed Bray-Curtis)
- Amplicon counts are compositional (Gloor 2017); applying Hill numbers and Bray-Curtis directly assumes absolute abundances and may distort results; consider presence-absence metrics (Chao2, Sorensen) or CLR-based diversity for amplicon data

## Related Skills

- ecological-genomics/edna-metabarcoding - Generate ASV/species tables from raw eDNA reads before diversity analysis
- ecological-genomics/community-ecology - Constrained ordination, indicator species, PERMANOVA on transformed community matrices
- microbiome/diversity-analysis - 16S clinical microbiome diversity metrics with explicit compositional considerations
- data-visualization/ggplot2-fundamentals - Customize diversity plots, rarefaction curves, and beta-diversity triangles
- phylogenetics/tree-io - Prepare ultrametric phylogenetic trees for PD, MPD, and MNTD
