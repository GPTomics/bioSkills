---
name: bio-ecological-genomics-species-delimitation
description: Delimits species boundaries from molecular data using distance-based (ASAP), tree-based (bPTP, GMYC), and coalescent (BPP) methods. Compares multiple delimitation results with delimtools. Use when delineating putative species from DNA barcoding data, resolving cryptic species complexes, or validating taxonomic assignments. Emphasizes multi-method consensus following integrative taxonomy best practice.
tool_type: mixed
primary_tool: ASAP
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, numpy 1.26+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Species Delimitation

**"Delineate species boundaries from my DNA barcoding data"** → Apply multiple delimitation methods (distance-based ASAP, tree-based bPTP/GMYC, coalescent BPP) to molecular data and compare results for multi-method consensus following integrative taxonomy best practice.
- CLI: ASAP web tool or standalone for distance-based partitioning
- Python: bPTP via `PTP-pyqt5` for Bayesian branching-rate analysis
- R: `splits::gmyc()` for coalescent/speciation transition model

Delimits putative species from molecular data using complementary distance-based, tree-based, and coalescent methods.

## Overview of Methods

| Method | Input | Approach | Strengths | Limitations |
|--------|-------|----------|-----------|-------------|
| ASAP | Aligned sequences | Distance-based partitioning | Fast, automatic, ranks partitions | Single locus only |
| bPTP | Rooted phylogeny | Bayesian branching rates | Bayesian support values | Over-splits with structure |
| GMYC | Ultrametric tree | Coalescent/speciation transition | Well-established theory | Requires ultrametric tree |
| BPP | Multi-locus alignment | Full coalescent model | Multi-locus, statistically rigorous | Computationally intensive |

## ASAP (Assemble Species by Automatic Partitioning)

**Goal:** Partition aligned barcode sequences into putative species using pairwise genetic distances.

**Approach:** Run ASAP with a substitution model (K2P, p-distance) and rank candidate partitions by asap-score.

Distance-based method that finds the optimal partition of sequences into putative species:

```bash
# ASAP web interface: https://bioinfo.mnhn.fr/abi/public/asap/asapweb.html
# Upload aligned FASTA, select distance model

# Command-line ASAP (if available)
ASAP -i aligned_sequences.fasta -d K2P -o asap_results/
# -d K2P: Kimura 2-parameter distance (standard for COI barcoding)
# -d p: uncorrected p-distance (for closely related taxa)
# -d JC: Jukes-Cantor (simplest substitution model)
```

### Interpreting ASAP Results

```
# ASAP ranks partitions by asap-score (lower = better)
# Reports: partition rank, number of groups, asap-score, p-value, threshold distance
# Always examine top 5-10 partitions, not just the best one
# Large gaps between consecutive asap-scores indicate robust partitioning
```

### Python ASAP-Style Analysis

```python
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
import numpy as np

alignment = AlignIO.read('aligned_sequences.fasta', 'fasta')

calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

names = [r.id for r in alignment]
n = len(names)
dist_array = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        dist_array[i][j] = dm[names[i], names[j]]

condensed = squareform(dist_array)
Z = linkage(condensed, method='average')

# Scan thresholds to identify barcode gap
# COI barcode gap typically at 2-3% for animals
thresholds = np.arange(0.01, 0.10, 0.005)
for t in thresholds:
    clusters = fcluster(Z, t=t, criterion='distance')
    n_clusters = len(set(clusters))
    print(f'Threshold {t:.3f}: {n_clusters} groups')
```

## bPTP (Bayesian Poisson Tree Processes)

**Goal:** Delimit species by detecting the transition from speciation to coalescent branching rates on a phylogeny.

**Approach:** Run Bayesian PTP MCMC on a rooted tree to obtain posterior support values for each species partition.

Tree-based method that models speciation as a branching rate shift:

```bash
# iTaxoTools PTP (Python 3 compatible)
pip install PTP-pyqt5

# Run bPTP from command line
python -m PTP.PTP -t rooted_tree.nwk -o bptp_results \
    --type bayesian --ngen 100000 --burnin 0.1 --seed 42

# --type bayesian: Bayesian MCMC (preferred over ML for support values)
# --ngen 100000: MCMC generations; increase for large trees
# --burnin 0.1: discard first 10% of samples
```

### Interpreting bPTP Output

```
# Support values for each species partition:
#   > 0.95: strong support
#   0.80-0.95: moderate support
#   < 0.80: weak support, boundary uncertain
# bPTP can over-split when populations have strong geographic structure
# Compare against distance-based methods to identify over-splitting
```

## GMYC (Generalized Mixed Yule Coalescent)

**Goal:** Identify the threshold on an ultrametric tree where branching shifts from interspecific (Yule) to intraspecific (coalescent) processes.

**Approach:** Prepare an ultrametric tree with chronos, run single-threshold GMYC, evaluate with a likelihood ratio test, and extract species assignments.

Models the transition between interspecific (Yule) and intraspecific (coalescent) branching rates on an ultrametric tree:

```r
library(splits)
library(ape)

# Step 1: Prepare ultrametric tree
# GMYC requires an ultrametric (time-calibrated) tree
tree <- read.tree('rooted_tree.nwk')

# chronos() for quick ultrametric conversion (penalized likelihood)
# lambda=1: smoothing parameter; higher = more clock-like
ultrametric_tree <- chronos(tree, lambda = 1)
class(ultrametric_tree) <- 'phylo'

# Verify ultrametric property
is.ultrametric(ultrametric_tree)

# Step 2: Run single-threshold GMYC
gmyc_result <- gmyc(ultrametric_tree, method = 'single')

# Summary: number of entities, confidence interval, likelihood ratio test
summary(gmyc_result)

# Number of putative species (ML entities)
cat('GMYC species:', gmyc_result$entity[1], '\n')

# Likelihood ratio test p-value
# p < 0.05: significant speciation-coalescent transition detected
cat('LR test p-value:', gmyc_result$p.value[1], '\n')

# Step 3: Extract species assignments
species_list <- spec.list(gmyc_result)
cat('Species assignments:\n')
print(species_list)
```

### Multiple-Threshold GMYC

```r
# Allows different coalescent-to-speciation transitions across the tree
# Better for datasets with variable population sizes or divergence times
gmyc_multi <- gmyc(ultrametric_tree, method = 'multiple')
summary(gmyc_multi)
```

### GMYC Visualization

```r
# Color-coded tree by species partition
pdf('gmyc_species_tree.pdf', width = 10, height = 12)
plot(gmyc_result, cex = 0.6)
dev.off()

# Support surface: likelihood across threshold positions
pdf('gmyc_support.pdf', width = 8, height = 5)
plot(gmyc_result, type = 'support')
dev.off()
```

## BPP (Bayesian Phylogenetics and Phylogeography)

**Goal:** Perform rigorous multi-locus species delimitation under the multispecies coalescent model.

**Approach:** Configure BPP A10 analysis (joint delimitation + species tree estimation) with priors on theta and tau, then run MCMC to obtain posterior probabilities for each delimitation model.

Full multispecies coalescent model for rigorous species delimitation:

```
# BPP control file for species delimitation (A10 analysis)
# A10: joint species delimitation and species tree estimation

seed = 12345

seqfile = alignment.phy
Imapfile = imap.txt      # maps individuals to putative species
outfile = bpp_results.txt
mcmcfile = bpp_mcmc.txt

speciesdelimitation = 1   # 1 = delimitation mode
speciestree = 1           # 1 = estimate species tree
speciesmodelprior = 1     # 1 = uniform prior on all rooted trees

# species&tree: species names, sample sizes, starting tree
species&tree = 4  speciesA speciesB speciesC speciesD
                  10 8 12 6
                  ((speciesA, speciesB), (speciesC, speciesD));

# Prior on theta (ancestral population size)
# thetaprior = alpha beta: inverse-gamma(alpha, beta)
# alpha=3, beta=0.004: mean=0.002 (typical for moderate diversity)
thetaprior = 3 0.004

# Prior on tau (divergence time)
# tauprior = alpha beta: inverse-gamma(alpha, beta)
# alpha=3, beta=0.002: mean=0.001 (relatively recent divergence)
tauprior = 3 0.002

# MCMC settings
nsample = 100000
sampfreq = 2
burnin = 50000
```

### Interpreting BPP Results

```
# Posterior probability for each delimitation model
# PP > 0.95: strong support for that delimitation
# PP 0.50-0.95: moderate support
# PP < 0.50: uncertain, consider lumping species
# Always run multiple independent chains to check convergence
```

## Multi-Method Comparison

**Goal:** Evaluate concordance across multiple delimitation methods to identify robust consensus species boundaries.

**Approach:** Compare partition vectors from ASAP, bPTP, and GMYC using the adjusted Rand index and count species per method.

Compare delimitation results across methods to identify consensus species:

```r
library(ape)

# Partition vectors from each method (species assignment per individual)
asap_partition <- c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4)
bptp_partition <- c(1, 1, 1, 2, 2, 3, 3, 4, 5, 5)
gmyc_partition <- c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4)

# Compare with adjusted Rand index
# ARI = 1: perfect agreement; ARI = 0: random; ARI < 0: less than random
library(fossil)
cat('ASAP vs bPTP ARI:', adj.rand.index(asap_partition, bptp_partition), '\n')
cat('ASAP vs GMYC ARI:', adj.rand.index(asap_partition, gmyc_partition), '\n')
cat('bPTP vs GMYC ARI:', adj.rand.index(bptp_partition, gmyc_partition), '\n')

# Consensus: individuals assigned to the same species by >= 2 methods
comparison <- data.frame(
    individual = paste0('ind_', 1:10),
    ASAP = asap_partition,
    bPTP = bptp_partition,
    GMYC = gmyc_partition
)

# Count number of species per method
cat('\nSpecies counts:\n')
cat('  ASAP:', length(unique(asap_partition)), '\n')
cat('  bPTP:', length(unique(bptp_partition)), '\n')
cat('  GMYC:', length(unique(gmyc_partition)), '\n')
```

## Ultrametric Tree Preparation

**Goal:** Convert a maximum-likelihood phylogeny into an ultrametric (time-calibrated) tree suitable for GMYC and other coalescent methods.

**Approach:** Apply penalized likelihood (chronos) with a relaxed clock model, then verify and fix near-ultrametric rounding errors.

Most tree-based methods require ultrametric input:

```r
library(ape)
library(phytools)

tree <- read.tree('ml_tree.nwk')

# Method 1: chronos (penalized likelihood, fastest)
# model='correlated': rates correlated across branches (relaxed clock)
ultra_tree <- chronos(tree, model = 'correlated')

# Method 2: BEAST2 (Bayesian, most rigorous but slow)
# Produces ultrametric tree with credibility intervals
# Run externally, import the MCC tree

# Verify
is.ultrametric(ultra_tree)

# Fix near-ultrametric trees with small rounding errors
ultra_tree <- force.ultrametric(ultra_tree, method = 'extend')
```

## Related Skills

- edna-metabarcoding - Generate barcode sequences from eDNA
- conservation-genetics - Population-level genetic assessment
- phylogenetics/tree-io - Tree input/output for delimitation methods
- phylogenetics/modern-tree-inference - Phylogenetic tree construction
- database-access/entrez-fetch - Retrieve barcode sequences from GenBank
