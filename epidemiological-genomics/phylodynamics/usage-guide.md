# Phylodynamics Usage Guide

## Overview

This skill drives postdoc-grade phylodynamic inference -- time-scaled phylogenies, R_e estimation, structured-population analysis -- from dated pathogen genomes. It covers ML phylodynamics (TreeTime), full Bayesian (BEAST2 with BICEPS, BDSKY, MASCOT, ORC clock), Bayesian bacterial dating (BactDating after Gubbins), pandemic-scale phylogenetics (UShER + matUtils + RIPPLES), and recombination masking for bacteria (Gubbins, ClonalFrameML). The decision-grade content is everything that surrounds the tool calls: TempEst temporal-signal QC first, recombination masking before any bacterial clock inference, BDSKY origin-vs-rootHeight discipline, multi-chain convergence for any BEAST report, MASCOT over BEAST DTA for source attribution, BICEPS over BSP as the modern skyline default, and reconciliation between phylodynamic R_e and case-based R_t with explicit sampling-bias accounting.

Distinct from `phylogenetics/divergence-dating` (general-purpose dating beyond outbreaks) and from `phylogenetics/bayesian-inference` (BEAST mechanics beyond BDSKY). This skill owns the surveillance application -- R_e estimation, Ne(t), source attribution, recombination-aware bacterial phylodynamics.

## Prerequisites

```bash
conda install -c bioconda beast2 treetime iqtree gubbins clonalframeml snippy snp-dists usher

packagemanager -add BDSKY BEASTLabs feast ORC MASCOT BICEPS

pip install dendropy biopython baltic

R -e "install.packages(c('BactDating', 'bdskytools', 'coda', 'ape', 'ggplot2'), repos='https://cran.r-project.org')"
```

## Quick Start

Tell the AI agent what phylodynamic question to answer:

- "Build a time-scaled tree for these dated SARS-CoV-2 consensus sequences with TreeTime and report root-to-tip R^2"
- "Estimate R_e through time for this Mtb outbreak using BEAST2 BDSKY with 3 independent chains"
- "Mask recombination in this S. pneumoniae core-genome alignment with Gubbins, then run BactDating"
- "Run MASCOT structured-coalescent for migration rates between 5 demes (>=20 sequences per deme)"
- "Place these new SARS-CoV-2 sequences on the daily UShER MAT and extract the subtree for downstream BDSKY"
- "Run a date-randomisation test and report whether the temporal signal is informative"
- "Reconcile phylodynamic R_e from BDSKY with case-based R_t from EpiNow2"

## Example Prompts

### Outbreak time-scaling with TempEst QC

> "Sixty *Mycobacterium tuberculosis* clinical isolates from a hospital cluster investigation spanning 18 months. Snippy against H37Rv, snippy-core for the core alignment, Gubbins on `core.full.aln` for recombination masking (even though Mtb is largely clonal, document recombination check for cross-lineage analysis), IQ-TREE on the masked alignment with `+ASC`, TreeTime for time-scaling with `--coalescent skyline --clock-filter 4`. Inspect `root_to_tip_regression.pdf` and report R^2; if R^2 < 0.3, flag insufficient temporal signal and do not proceed to BDSKY. Run date-randomisation as secondary diagnostic. Document the bundled Coll/Napier barcode version used by TB-Profiler for lineage assignment."

### Full Bayesian R_e estimation for a bacterial outbreak

> "Two hundred *Klebsiella pneumoniae* isolates from a 3-year regional carbapenemase surveillance project. Snippy + snippy-core, Gubbins to mask recombination on `core.full.aln` (Klebsiella is recombining; mandatory), IQ-TREE on the masked alignment with `+ASC`, then BEAST2 with BDSKY in BEAUti template. Set `origin > rootHeight` per Stadler 2013 (init to tMRCA + 0.1*tMRCA). Set sampling proportion per quarter (the actual fraction of clinical isolates sequenced). Become-uninfectious rate from epi knowledge (1/14 days). Run 4 independent BEAST chains from seeds 42, 17, 99, 7; combine via logcombiner only after marginal posteriors overlap; report Gelman-Rubin via coda. Compare phylodynamic R_e to case-based R_t from local surveillance data."

### Recombination masking for S. pneumoniae phylodynamics

> "One thousand *Streptococcus pneumoniae* genomes from a national post-PCV13 carriage surveillance. Build core-genome alignment via snippy-core against ATCC 700669 reference; run Gubbins on `core.full.aln` (S. pneumoniae recombines intensely; unmasked clock rate would inflate 2-5x); rebuild tree with IQ-TREE on the recombination-masked alignment; compare clock rate to Croucher 2013 literature (~1.5e-6 subs/site/year) as a sanity check. If rate is >5e-6, recombination masking is insufficient -- consider fastGEAR (Mostowy 2017) for the highly-recombinogenic lineages."

### Structured-coalescent migration with MASCOT

> "Four hundred SARS-CoV-2 sequences from 5 European countries during Alpha emergence. >=20 sequences per country. BEAUti -> MASCOT template; one trait per country; if sampling intensity covaries with travel volume, switch to MASCOT-GLM. Document why MASCOT not BEAST DTA (Lemey 2009 DTA is sampling-biased per De Maio 2015). Report migration rate matrix with HPD intervals; flag any HPD spanning >2 orders of magnitude as identifiability failure."

### Pandemic-scale UShER + TreeTime + BDSKY workflow

> "Place 5,000 new SARS-CoV-2 sequences on the UCSC daily-updated UShER MAT; use matUtils to extract the JN.1-descendant subtree; re-estimate branch lengths with TreeTime (UShER parsimony branch lengths are biased low); fit BDSKY on the subtree for lineage-specific R_e through time. Document that R_e is lineage-specific phylodynamic, NOT population-wide; reconcile with case-based R_t."

### Date-randomisation and temporal-signal diagnostics

> "These 80 dated influenza H3N2 HA sequences span 6 months of a single season. Run TempEst root-to-tip regression; if R^2 < 0.3, the data do not support time-scaled inference. Run date-randomisation (TreeTime `clock --reassign-dates`) as secondary check; note that date-randomisation can pass with narrow sampling windows (false negative). If temporal signal insufficient, report a tree without time-scaling and explain why R_e estimation is not supported."

## What the Agent Will Do

1. Identify the phylodynamic question (Ne(t), R_e through time, divergence date, migration rates, source attribution) and the data shape (number of tips, sampling window, pathogen).
2. Build a raw topology with IQ-TREE (outbreak scale) or RAxML-NG (larger).
3. For bacteria, mask recombination FIRST: Gubbins on `core.full.aln` (NOT `core.aln`); rebuild tree on the recombination-masked alignment.
4. Run TempEst root-to-tip regression as the first temporal-signal diagnostic; R^2 >= 0.3 minimum.
5. Run date-randomisation as a secondary diagnostic; flag that it can pass with narrow sampling windows.
6. For ML phylodynamics, run TreeTime with `--coalescent skyline --clock-filter 4 --confidence`.
7. For Bayesian phylodynamics, use BEAUti to set up the appropriate template (BICEPS for Ne(t); BDSKY for R_e; MASCOT for structured); set `origin > rootHeight` for BDSKY; set sampling proportion per epoch.
8. Run 3-4 independent BEAST2 chains from different seeds; combine via logcombiner only after marginal posteriors overlap; report ESS and Gelman-Rubin.
9. For source attribution, never use BEAST DTA; use MASCOT or MASCOT-GLM; frame DTA results as "consistent with" only.
10. For pandemic-scale data, use UShER for placement; re-estimate branch lengths via TreeTime before downstream R_e.
11. Reconcile phylodynamic R_e with case-based R_t; report both with sampling-bias assumptions.
12. Pin BEAST version AND every BEAST2 package version installed via `packagemanager`; record random seed.

## Tips

- TempEst R^2 < 0.3 means time-scaling is not supported. Date-randomisation can pass anyway with narrow sampling windows -- false sense of security.
- `run_gubbins.py core.aln` is wrong. The input must be `core.full.aln` (full positions including invariant), or Gubbins cannot estimate background SNP density.
- BDSKY `origin` is the time from the start of the epidemic to the most recent sample. `rootHeight` is the time back to the MRCA of the sample. `origin > rootHeight` always; setting them equal biases R_e upward systematically.
- BSP / Skygrid skyline is superseded by BICEPS (Bouckaert 2022) in BEAST 2.7+. BICEPS uses analytic Ne integration per epoch and new tree-flexing operators; weeks-of-BSP becomes hours.
- BEAST DTA (Lemey 2009) is sampling-biased for source attribution per De Maio 2015. Use MASCOT or MASCOT-GLM for any migration / phylogeography claim.
- MASCOT requires >=20 sequences per deme for migration to be identifiable.
- BEAST ESS > 200 is necessary but not sufficient. Run >=3 chains; check marginal posteriors per chain; combine only after overlap.
- For S. pneumoniae, recombination-masked core clock should be ~1.5e-6 subs/site/year. Rates >5e-6 indicate unmasked recombination.
- TreeTime is 100-1000x faster than BEAST for outbreak-scale data. Reserve BEAST for posteriors on R_e / origin date / migration rate.
- UShER branch lengths are parsimony-derived and systematically shorter than ML; re-estimate before BDSKY.
- R_e (phylodynamic, current immunity context) is not R_0 (basic, fully susceptible). Use R_e (or R_t) consistently.
- BEAST2 XML is fragile across minor releases. Pin BEAST version + every package version + random seed for reproducibility.
- For MASCOT-GLM, provide predictor variables (travel volume, geographic distance) that plausibly drive migration; without covariates, MASCOT-Skyline is the time-varying alternative.
- Featherstone & Duchêne 2023 *Mol Biol Evol* 40:msad132 quantified that for shallow trees with many samples, sampling times dominate over sequence information for R_e inference. Bias in sampling drives bias in R_e regardless of sequence quantity.

## Related Skills

- pathogen-typing - Defines isolates and clusters feeding phylodynamic inference
- transmission-inference - Phylodynamic R_e + transmission tree are complementary
- variant-surveillance - Lineage assignment runs upstream of skygrid / BDSKY by lineage
- phylogenetics/divergence-dating - General-purpose clock-model review beyond outbreaks
- phylogenetics/bayesian-inference - BEAST / RevBayes mechanics beyond BDSKY
- phylogenetics/modern-tree-inference - IQ-TREE / RAxML upstream of time-scaling
- phylogenetics/tree-io - Tree parsing and format conversion
- comparative-genomics/whole-genome-alignment - Core-genome alignment for bacterial phylodynamics
- variant-calling/vcf-basics - Per-isolate variant calls feeding core-SNP alignment
- read-alignment/bwa-alignment - Read mapping upstream of variant calling
- data-visualization/multipanel-figures - Skyline / R_e trajectory plotting
- workflows/somatic-variant-pipeline - End-to-end orchestration patterns
