# Differential Abundance - Usage Guide

## Overview

Differential abundance (DA) testing identifies which individual taxa differ between groups on an amplicon ASV/feature table while respecting the compositional, relative nature of the counts. The headline workflow is a CONSENSUS of two or more compositionally-aware tools: which taxa come out "significant" depends more on the DA tool than on the biology (Nearing 2022, across 38 datasets), so a single tool's hit list is not a defensible deliverable. Run at least two of ALDEx2, ANCOM-BC2, MaAsLin2/3, LinDA, or ZicoSeq, report the intersection as high-confidence and the union as exploratory, name every tool, and disclose disagreement. A relative-abundance increase is also not an absolute increase without an external load anchor (spike-in / flow cytometry / qPCR).

This skill owns per-taxon DA on an amplicon table. Whole-community alpha/beta diversity and PERMANOVA live in diversity-analysis; the same DA math on shotgun profiler tables lives in metagenomics/metagenome-visualization; the shared compositional/closure/CLR/zero theory lives in metagenomics/abundance-estimation.

## Prerequisites

```r
BiocManager::install(c('ALDEx2', 'ANCOMBC', 'Maaslin2'))
install.packages(c('MicrobiomeStat', 'GUniFrac'))   # LinDA, ZicoSeq
```

Conceptual prerequisites:
- A feature table of integer COUNTS (ASVs or taxa collapsed to genus/species), plus sample metadata - typically a phyloseq object. ALDEx2, ANCOM-BC2, LinDA, and ZicoSeq expect counts; MaAsLin2 TSS-normalizes internally and expects features in COLUMNS.
- Answer the whole-community question first (diversity-analysis): knowing whether the communities differ at all frames the per-taxon hits.
- Decide the DA-tool panel (>=2 tools) a priori, before seeing any result, to avoid cherry-picking.
- Decide and declare a prevalence filter; it is a modeling knob that reshapes the FDR landscape.
- Know whether the design has repeated/paired samples; if so, a random effect is required.
- Remove host organelle (Mitochondria/Chloroplast) features and, for low-biomass samples, reagent contaminants upstream (taxonomy-assignment / amplicon-processing) before testing - they otherwise surface as spurious hits or skew the closure.

## Quick Start

Tell your AI agent what you want to do:
- "Find differentially abundant taxa between treatment and control with two methods and give me the consensus"
- "Run ALDEx2 and report effect sizes with BH-adjusted q-values"
- "Run ANCOM-BC2 with age and sex as covariates and only keep hits that pass the sensitivity analysis"
- "Analyze a longitudinal study with a subject random effect using LinDA or MaAsLin2"
- "Filter taxa present in fewer than 10% of samples before testing and check the result is not sensitive to that cutoff"

## Example Prompts

### Consensus across tools
> "I have an ASV table and metadata with two treatment groups. Run at least two compositionally-aware DA methods, report the intersection as high-confidence and the union as exploratory, and tell me which tools found each taxon."

### Conservative two-group test
> "Run ALDEx2 to compare taxon abundance between healthy and diseased samples, gate on both BH-adjusted q below 0.05 and an effect-size floor, and explain the effect-size metric."

### Covariates and sensitivity
> "Run ANCOM-BC2 with age and sex as covariates, set the p-adjust method to BH, and only report hits that pass the pseudo-count sensitivity analysis (passed_ss)."

### Repeated measures
> "My samples are repeated within subjects over time. Use a DA method with a subject random effect so I do not pseudo-replicate, and cross-check the mixed-model hits across tools."

### Relative vs absolute
> "Is this taxon's increase relative or absolute? I do not have load data - explain what I can and cannot claim, and what a spike-in or qPCR anchor would add."

## What the Agent Will Do

1. Confirm the input is integer counts with matching metadata and check feature orientation per tool.
2. Apply and declare a prevalence/abundance filter (e.g. taxa present in >=10% of samples).
3. Run the first compositionally-aware tool (default ALDEx2) and extract BH-adjusted q plus effect size.
4. Run at least one more tool (ANCOM-BC2, LinDA, MaAsLin2/3, or ZicoSeq), adding covariates or a random effect as the design requires.
5. For ANCOM-BC2, set p_adj_method to BH and require passed_ss for confident hits.
6. Intersect the per-tool significant sets: report the intersection as high-confidence, the union as exploratory, and tabulate which tools agree per taxon.
7. State whether claims are relative or absolute and whether a load anchor exists.
8. Produce a results table and an effect-size plot; never pool p-values across tools.

## Tips

- The deliverable is a consensus, not one tool's list. Decide the panel a priori and report all tools, including the ones that disagree.
- Gate on effect size AND q-value, not p alone - large n makes trivially small differences "significant."
- ANCOM-BC2 defaults to Holm, not BH; set p_adj_method to BH deliberately if FDR is wanted, and require passed_ss.
- Repeated/paired samples need a random effect (ANCOM-BC2 rand_formula, MaAsLin2 random_effects, LinDA mixed formula); ignoring it inflates significance.
- The prevalence filter is a modeling choice. Declare the threshold and confirm the headline result survives moving it from 10% to 25%.
- A relative increase is not an absolute increase without a spike-in / flow / qPCR anchor or MaAsLin3's absolute-abundance mode.
- DESeq2/edgeR are RNA-seq-native and misfire on sparse zero-heavy tables (the geometric-mean size factor collapses) - treat as a caveat, not a recipe.
- An uncorrected Wilcoxon/t-test on relative abundances is wrong twice (closure and multiple testing). A BH-corrected simple test honestly labelled as relative is defensible, ideally alongside a compositional tool.
- There is no settled best tool: Nearing favors conservative ALDEx2/ANCOM-II, Yang and Chen favor ZicoSeq/LinDA for power, Pelto finds elementary methods most replicable. Verify current best practice against the latest docs.
- Filter host organelle features and decontaminate low-biomass samples before DA (taxonomy-assignment, amplicon-processing); contaminant ASVs otherwise appear among the hits.

## Related Skills

- diversity-analysis - Whole-community alpha/beta/PERMANOVA; answer "do the communities differ" first
- taxonomy-assignment - Collapse ASVs to genus/species before per-taxon testing
- amplicon-processing - Produces the ASV feature table tested here
- qiime2-workflow - The qiime composition ancombc CLI route
- metagenomics/abundance-estimation - Shared compositional/closure/CLR/zero/load-anchor theory
- metagenomics/metagenome-visualization - The same DA mechanics on shotgun profiler tables
- experimental-design/multiple-testing - FDR control and multiplicity across taxa
