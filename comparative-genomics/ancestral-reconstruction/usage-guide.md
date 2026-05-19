# Ancestral Reconstruction - Usage Guide

## Overview

Ancestral state reconstruction (ASR) infers states at internal phylogenetic nodes for sequences (codon, protein, nucleotide), discrete traits, continuous traits, and gene content. The choice of framework depends on the data class and inference question. For sequences, **PAML codeml `RateAncestor=1`** (Yang 1995 Genetics 141:1641) and **IQ-TREE2 `--ancestral`** provide marginal ML reconstructions; **GRASP** (Foley 2022) handles indels probabilistically -- critical for protein resurrection. For discrete traits, **corHMM** (Boyko & Beaulieu 2021) with hidden Markov rate models supersedes simpler Mk approaches.

The most consequential mistake is reconstructing under a site-independent or trait-stationary model when the biology demands a hidden-rate / non-stationary / epistatic model -- the resulting "ancestor" is mathematically optimal under the wrong model and silently wrong.

## Prerequisites

```bash
# Sequence ASR
conda install -c bioconda paml iqtree

# RevBayes for Bayesian ASR
# https://revbayes.github.io/download

# GRASP for indel-aware protein ASR
# https://github.com/bodenlab/GRASP

# FastML
conda install -c bioconda fastml

# R packages for trait ASR
Rscript -e "install.packages(c('ape', 'phytools', 'geiger', 'corHMM', 'phangorn', 'OUwie', 'bayou', 'RERconverge'))"
Rscript -e "remotes::install_github('thej022214/hisse')"

# Python
pip install biopython ete4
```

## Quick Start

Tell the AI agent what to do:
- "Reconstruct the ancestral protein sequence at the mammalian root for resurrection construct design"
- "Infer ancestral state of trait X using corHMM with hidden rate categories"
- "Perform stochastic character mapping with 1000 simulations to get full posterior over character histories"
- "Reconstruct ancestral continuous trait values along the tree under OU model with regime shifts"

## Example Prompts

### Protein Resurrection

> "Reconstruct the ancestral cytochrome c sequence at the eukaryote root using GRASP for indel-aware ASR. Design 8 alternative constructs varying ambiguous sites (P < 0.8). Validate ML construct against extant homologs as functional baseline."

### Discrete Trait Evolution

> "Reconstruct ancestral states of a binary trait (flight ability) on 50 bird species. Test ER vs SYM vs ARD rate models via AIC. If rates are asymmetric, use ARD. Cross-validate against stochastic mapping with `phytools::make.simmap(nsim=5000, pi='fitzjohn')`."

### State-Dependent Diversification

> "Test whether a discrete trait (sociality) influences speciation/extinction across 200 species using HiSSE as the required null. BiSSE has 40% Type-I false-positive rate (Rabosky-Goldberg 2015); only conclude trait-dependent diversification if HiSSE-null is rejected."

## What the Agent Will Do

1. **Validate inputs**: tree rooted; MSA quality; trait data completeness
2. **Run appropriate framework**:
   - Sequence ASR: PAML codeml + IQ-TREE2 + GRASP for indels
   - Discrete trait: corHMM with hidden rates + ape::ace for AIC model comparison
   - Continuous trait: phytools::fastAnc + geiger::fitContinuous for model selection
   - Bayesian: RevBayes or BayesTraits for full posterior
3. **Apply model selection**: ER/SYM/ARD for discrete; BM/OU/EB/lambda for continuous; site-homogeneous vs CAT-PMSF for deep sequence
4. **Stochastic mapping** for discrete trait posterior; >=1000 simulations
5. **Quantify confidence**: per-site posteriors; phylogenetic signal (Pagel's lambda, Blomberg's K)
6. **Report**: ancestral states with uncertainty, model fit, ambiguous sites for downstream testing
7. **Caveats**: epistasis, LBA at deep nodes, model adequacy, root sensitivity

## Tips

- For protein resurrection, use GRASP (handles indels) + design alternative constructs at ambiguous sites (P < 0.8)
- corHMM with rate.cat=2 is the modern standard for discrete trait ASR; legacy ape::ace lacks hidden rates
- HiSSE is the required null for any trait-diversification correlation test (Rabosky-Goldberg 2015)
- For continuous traits, always run BM vs OU vs EB AIC comparison; report Pagel's lambda
- Stochastic mapping with `pi='fitzjohn'` root prior; >=1000 simulations
- For deep sequence ASR, use site-heterogeneous models (CAT-PMSF, PhyloBayes) to mitigate LBA
- Marginal reconstruction for per-site uncertainty; joint reconstruction for internally consistent ancestor
- Multiple-rooting sensitivity tests; STRIDE / MAD rooting alternatives
- Document model adequacy: AIC + posterior predictive checks
- Codon-level ASR uses PAML M0 (single omega) -- site models are for selection, not ASR
- Bayesian methods (RevBayes, MrBayes) provide honest uncertainty; computationally heavy
- For convergent rate shifts, RERconverge (categorical phenotype) or CSUBST (amino-acid convergence)
- PhyloAcc for noncoding accelerated evolution; not codon-based
- Document phylogenetic signal: lambda > 0.7 strong; 0.3-0.7 moderate; < 0.3 weak

## Related Skills

comparative-genomics/positive-selection - Branch- and site-level selection on ancestral branches
comparative-genomics/ortholog-inference - Orthogroups feed sequence ASR
comparative-genomics/gene-tree-species-tree-reconciliation - DTL-aware ancestral gene content
comparative-genomics/whole-genome-duplication - Time scale from Ks for trait dating
comparative-genomics/comparative-annotation-projection - Project ancestral CDS to descendants
phylogenetics/modern-tree-inference - Rooted ML/Bayesian trees
phylogenetics/bayesian-inference - RevBayes / MrBayes priors
phylogenetics/divergence-dating - Time-calibrated trees for absolute-time ASR
alignment/multiple-alignment - PRANK / MACSE indel-aware alignment
alignment/alignment-trimming - PREQUAL / HmmCleaner filtering before ASR
