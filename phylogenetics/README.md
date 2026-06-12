# phylogenetics

## Overview

Phylogenetic tree analysis covering I/O, manipulation, visualization, distance-based methods, maximum likelihood inference, Bayesian analysis, divergence time estimation, and coalescent species tree methods. Provides expert-level decision guidance built on one organizing idea: a tree is an estimate under an assumed model, and support is not accuracy.

**Tool type:** python, cli, mixed | **Primary tools:** Bio.Phylo, IQ-TREE2, RAxML-NG, MrBayes, BEAST2, ASTRAL-III

## Skills

| Skill | Description |
|-------|-------------|
| tree-io | Read, write, convert tree files (Newick, Nexus, NHX, PhyloXML, NeXML) without dropping BEAST/MrBayes annotations |
| tree-visualization | Draw trees with matplotlib; decision guide for ggtree, ETE4, iTOL; a figure is an argument |
| tree-manipulation | Root, prune, ladderize, collapse trees; rooting as a separate inference, soft vs hard polytomies |
| distance-calculations | Model-corrected distance matrices, NJ/BIONJ/FastME/UPGMA trees, saturation testing, bootstrap consensus |
| modern-tree-inference | ML trees with IQ-TREE2/RAxML-NG, ModelFinder, UFBoot2/SH-aLRT, concordance factors, topology tests |
| bayesian-inference | MrBayes, BEAST2, RevBayes, PhyloBayes, MCMC convergence diagnostics, stepping-stone model comparison, CAT models |
| divergence-dating | Molecular clocks, fossil calibration, BEAST2/MCMCTree/TreePL, tip-dating, effective-prior checks |
| species-trees | ASTRAL/coalescent methods, concordance factors, gene tree discordance, ILS and the anomaly zone |

## Example Prompts

- "Which inference method should this dataset use, and can I trust the support values?"
- "Every node has 100% bootstrap but two papers disagree. How do I adjudicate?"
- "Read this Newick tree file and show the taxa"
- "Convert my Nexus tree to Newick format"
- "Draw this tree and save as PDF"
- "Root this tree using Mouse as outgroup. Should I use midpoint or outgroup rooting?"
- "Build a neighbor joining tree from this alignment as a quick sanity check"
- "Run IQ-TREE2 with ModelFinder, ultrafast bootstrap, and SH-aLRT"
- "Find the best substitution model. Should I use MFP or TEST?"
- "Compute concordance factors for my phylogenomic dataset"
- "My tree has low support throughout. What does that mean?"
- "I suspect long branch attraction. How do I detect and fix it?"
- "Run a Bayesian analysis with MrBayes and check MCMC convergence"
- "Compare two substitution models using stepping-stone sampling"
- "Estimate divergence times with BEAST2 using fossil calibrations"
- "Should I use a strict or relaxed molecular clock?"
- "Infer a species tree from my multi-locus dataset using ASTRAL"
- "Should I concatenate or use coalescent methods for my phylogenomic data?"
- "Analyze my multi-gene dataset with partitioned models"

## Requirements

```bash
pip install biopython matplotlib numpy

# ML inference
conda install -c bioconda iqtree raxml-ng

# Bayesian inference
conda install -c bioconda mrbayes beast2

# Species trees
conda install -c bioconda aster  # includes ASTRAL-III, wASTRAL, ASTRAL-Pro
```

## Related Skills

- **alignment** - Prepare MSAs for tree building
- **sequence-io** - Read sequences for alignment and tree building
- **database-access** - Fetch sequences from NCBI for phylogenetic analysis
- **epidemiological-genomics** - Phylodynamics and outbreak investigation
