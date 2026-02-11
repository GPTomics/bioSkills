# Species Delimitation Usage Guide

## Overview

Delimits putative species boundaries from molecular data using complementary approaches: distance-based partitioning (ASAP), tree-based branching rate models (bPTP, GMYC), and full coalescent analysis (BPP). Follows integrative taxonomy best practice by comparing results across multiple methods to identify consensus species. Designed for DNA barcoding datasets (COI, 12S, ITS, matK), cryptic species complexes, and taxonomic validation studies.

## Prerequisites

- ASAP web interface (https://bioinfo.mnhn.fr/abi/public/asap/) or CLI
- PTP-pyqt5 for bPTP (`pip install PTP-pyqt5` from PyPI, or install from GitHub: `pip install git+https://github.com/iTaxoTools/PTP-pyqt5`)
- R with splits for GMYC (`install.packages('splits', repos = 'http://R-Forge.R-project.org')`)
- R with ape for tree manipulation (`install.packages('ape')`)
- BPP v4 (https://github.com/bpp/bpp) for coalescent delimitation
- fossil for partition comparison (`install.packages('fossil')`)

## Quick Start

Tell your AI agent what you want to do:

- "Delimit species in my COI barcoding dataset using ASAP"
- "Run GMYC species delimitation on my ultrametric phylogeny"
- "Compare species delimitation results from ASAP, bPTP, and GMYC"
- "Set up BPP for multi-locus species delimitation"
- "Find the barcode gap in my sequence alignment"

## Example Prompts

### Distance-Based

> "I have an aligned FASTA of 200 COI barcode sequences. Run ASAP with K2P distances and report the top 5 partition rankings."

> "Compute pairwise distances for my barcoding alignment, identify the barcode gap, and cluster sequences at the optimal threshold."

### Tree-Based

> "I have a rooted ML tree from IQ-TREE. Run bPTP with 100,000 MCMC generations and report species partitions with support values."

> "Convert my ML tree to ultrametric with chronos, run single-threshold GMYC, and show the species partition on the tree."

### Coalescent

> "I have alignments for 4 loci and a map of individuals to putative species. Set up and run BPP A10 analysis for joint species delimitation and tree estimation."

### Multi-Method Comparison

> "I ran ASAP, bPTP, and GMYC on my dataset. Compare the three partitions, compute adjusted Rand indices, and identify consensus species boundaries."

### Integrative Taxonomy

> "Delimit species using both distance-based and tree-based methods. For each putative species, report the number of individuals, geographic range, and method agreement."

## What the Agent Will Do

1. Prepare the input data (aligned sequences for ASAP, rooted tree for bPTP, ultrametric tree for GMYC, multi-locus alignments for BPP)
2. Convert ML trees to ultrametric using chronos if needed for GMYC
3. Run the selected delimitation method(s) with appropriate parameters
4. Extract putative species assignments with support values or partition scores
5. Compare results across methods using adjusted Rand index
6. Generate visualizations: color-coded trees, barcode gap histograms, partition comparison tables
7. Report consensus species boundaries supported by multiple methods
8. Provide recommendations for further validation (morphology, ecology, additional loci)

## Tips

- Always run at least two methods and compare results; single-method delimitation is unreliable
- ASAP is the fastest starting point and works directly from aligned sequences without tree building
- bPTP tends to over-split when populations have strong geographic structure; compare with ASAP to detect this
- GMYC requires a strictly ultrametric tree; use chronos() for quick conversion or BEAST2 for rigorous time-calibration
- BPP is the most statistically rigorous method but requires multi-locus data and careful prior specification
- The COI barcode gap is typically at 2-3% for animals, but varies across taxa; examine the distance histogram
- K2P distance is standard for COI barcoding; use p-distance for very closely related taxa
- For GMYC, the multiple-threshold variant handles datasets with variable population sizes better than single-threshold
- Always verify ultrametric property with `is.ultrametric()` before running GMYC
- BPP priors on theta and tau should reflect reasonable biological values; consult published studies for the target group
- Report the number of species found by each method and highlight discordant assignments for further investigation

## Related Skills

- edna-metabarcoding - Generate barcode sequences from eDNA
- conservation-genetics - Population-level genetic assessment
- phylogenetics/tree-io - Tree input/output for delimitation methods
- phylogenetics/modern-tree-inference - Phylogenetic tree construction
- database-access/entrez-fetch - Retrieve barcode sequences from GenBank
