# Diversity Analysis - Usage Guide

## Overview

Diversity analysis summarizes the whole microbial community: alpha diversity (within-sample richness and evenness) and beta diversity (between-sample dissimilarity) of an amplicon ASV/OTU table, with QIIME2 core-metrics-phylogenetic, phyloseq/vegan, and scikit-bio. The central discipline is that a diversity number is the output of three choices made before the number appears - the rarefaction depth, the tree, and the metric - so the agent declares all three and shows the conclusion survives a second reasonable choice rather than reading a single p-value.

The two decisions that dominate everything: the sampling depth in `core-metrics` silently deletes every sample below it (biased toward low-biomass samples), and the tree behind UniFrac/Faith PD is a modeling choice (de novo from short reads loses to SEPP fragment-insertion and Greengenes2 placement). Per-taxon testing belongs to differential-abundance; shotgun profiler tables belong to metagenomics/metagenome-visualization; the shared compositional and rarefaction-debate theory lives in metagenomics/abundance-estimation, and the Hill-number and PERMANOVA-dispersion theory in metagenomics/metagenome-visualization.

## Prerequisites

```r
BiocManager::install(c('phyloseq', 'picante'))
install.packages(c('vegan', 'GUniFrac'))
```

```bash
# QIIME2 installs as its own conda env (the release tag defines the plugin API and .qza format)
conda env create -n qiime2-amplicon-2024.2 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py39-linux-conda.yml
# scikit-bio (Python engine under q2-diversity)
pip install scikit-bio
```

Conceptual prerequisites:
- An ASV/OTU feature table of integer counts (DADA2/Deblur or a QIIME2 `FeatureTable[Frequency]`), plus representative sequences if a tree is needed.
- A phylogenetic tree is REQUIRED for UniFrac and Faith PD; prefer a SEPP-into-reference or Greengenes2 tree over a de novo build from short reads.
- Sample metadata with the grouping/covariate columns.
- Reference packages for SEPP (e.g. Greengenes 13_8 or SILVA) are large downloads; the QIIME2 release tag is recorded with the artifacts.
- The table should already have host mitochondria/chloroplast features removed (taxonomy-assignment) and, for low-biomass samples, reagent contaminants removed with decontam (amplicon-processing) - both distort every diversity metric.

## Quick Start

Tell your AI agent what you want to do:
- "Pick a rarefaction sampling depth from my feature table and tell me which samples it drops"
- "Run core-metrics-phylogenetic and report alpha and beta diversity"
- "Calculate Shannon and Faith PD per sample and compare my groups"
- "Compute weighted and unweighted UniFrac and test the group difference with PERMANOVA and betadisper"
- "Build a SEPP fragment-insertion tree instead of a de novo tree for UniFrac"

## Example Prompts

### Choosing the sampling depth
> "I have a QIIME2 feature table. Summarize the per-sample frequencies, show me an alpha-rarefaction curve, recommend a sampling depth on the plateau, and tell me how many and which samples that depth would drop."

### Alpha diversity
> "Rarefy my ASV table to the depth I chose, calculate observed features, Shannon, and inverse Simpson, report effective species, and test whether diversity differs between control and treatment with a non-parametric test."

### The tree decision
> "My reads are V4 16S. Build a SEPP fragment-insertion tree against a full-length reference instead of a de novo MAFFT+FastTree tree, filter out the ASVs that failed to insert, and explain why this matters for UniFrac."

### Beta diversity
> "Compute weighted and unweighted UniFrac plus generalized UniFrac at alpha 0.5, make PCoA plots, run PERMANOVA on each, and pair every PERMANOVA with betadisper so I know whether a significant result is a location shift or a dispersion difference."

### Compositional ordination
> "Run a robust Aitchison PCA (RPCA) on my ASV table and tell me which ASVs load on the first axis."

## What the Agent Will Do

1. Confirm the input is an amplicon ASV/OTU count table (routes shotgun tables to metagenomics/metagenome-visualization).
2. Read the per-sample frequency distribution and an alpha-rarefaction curve to choose a sampling depth that saturates richness while retaining samples.
3. Report the chosen depth and the list of samples that depth drops, and confirm the conclusion at a nearby depth.
4. Select or build a tree for phylogenetic metrics, preferring SEPP fragment-insertion or Greengenes2 over a de novo build, and filter ASVs that failed to insert.
5. Calculate alpha diversity spanning richness and evenness (Hill q0/q1/q2), report effective species, and test group differences.
6. Calculate beta diversity with both weighted and unweighted UniFrac (and generalized UniFrac alpha=0.5), then ordinate by PCoA.
7. Run PERMANOVA (adonis2) for the group effect and pair it with betadisper to separate location from dispersion.
8. Keep the raw counts and route per-taxon differential-abundance testing to the differential-abundance skill on the unrarefied table.

## Tips

- The sampling depth is the single biggest lever: it silently deletes every sample below it, and the dropped samples are rarely random. Always report the depth AND the dropped-sample list.
- Do not rarefy to `min(sample_sums)` - one tiny library drags every sample down to under-saturated noise. Pick a depth from the rarefaction plateau instead.
- A de novo tree from short reads is the weakest option; prefer SEPP-into-reference or Greengenes2 placement. Treat unweighted UniFrac on a de novo tree as suspect.
- Report weighted AND unweighted UniFrac. Unweighted listens to rare lineages and topology; weighted listens to abundant lineages; they can tell opposite stories.
- Observed features counts ASVs, not species, and moves with DADA2 denoising settings. Prefer Hill q1 (exp Shannon) / q2 (inverse Simpson) as primary richness metrics.
- QIIME2 reports Shannon in log2 (bits) and R reports it in natural log (nats). Report `exp(H')` (effective species) to compare across the two ecosystems.
- Always run betadisper alongside PERMANOVA. A significant PERMANOVA can be a dispersion difference, not a composition shift.
- Rarefy for diversity, never for differential abundance. Keep the raw counts and send DA to the differential-abundance skill.
- Filter host mitochondria/chloroplast and decontaminate low-biomass samples upstream (taxonomy-assignment, amplicon-processing) before diversity; organelle and contaminant reads shift richness and ordination.

## Related Skills

- amplicon-processing - Generate the ASV table and representative sequences upstream
- taxonomy-assignment - Label the ASVs summarized here
- differential-abundance - Per-taxon between-group testing on the unrarefied counts
- qiime2-workflow - The QIIME2 CLI home for core-metrics and tree building
- phylogenetics/tree-io - Read, write, and root the UniFrac/Faith PD tree
- metagenomics/abundance-estimation - Shared CoDA and rarefaction-debate theory
- metagenomics/metagenome-visualization - Shared Hill-number and PERMANOVA-dispersion theory; diversity/ordination on shotgun profiler tables
- data-visualization/ggplot2-fundamentals - Custom ordination and diversity plots
