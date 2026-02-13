# Reference: BioPython 1.83+, numpy 1.26+, scipy 1.12+ | Verify API if version differs
library(splits)
library(ape)
library(phytools)

# --- Step 1: Load rooted phylogenetic tree ---
# Input: rooted ML or Bayesian tree (Newick format)
# Must be rooted; if not, use root() with an outgroup
tree <- read.tree('rooted_tree.nwk')

cat('Tips:', Ntip(tree), '\n')
cat('Rooted:', is.rooted(tree), '\n')

# Root tree if needed (uncomment and specify outgroup)
# tree <- root(tree, outgroup = 'outgroup_taxon', resolve.root = TRUE)

# --- Step 2: Convert to ultrametric ---
# GMYC requires a strictly ultrametric (clock-like) tree
# chronos(): penalized likelihood method (Sanderson 2002)
# model='correlated': rates correlated across branches (relaxed clock)
# lambda=1: smoothing parameter; higher = more clock-like
ultra_tree <- chronos(tree, model = 'correlated', lambda = 1)
class(ultra_tree) <- 'phylo'

cat('Ultrametric:', is.ultrametric(ultra_tree), '\n')

# Fix near-ultrametric trees with tiny rounding errors
if (!is.ultrametric(ultra_tree)) {
    ultra_tree <- force.ultrametric(ultra_tree, method = 'extend')
    cat('Forced ultrametric (extend method)\n')
}

# --- Step 3: Run single-threshold GMYC ---
# Single threshold: one speciation-to-coalescent transition point
# Appropriate when all species have similar effective population sizes
gmyc_single <- gmyc(ultra_tree, method = 'single')

cat('\n--- Single-Threshold GMYC Results ---\n')
summary(gmyc_single)

n_entities <- gmyc_single$entity[1]
cat('ML number of entities (species):', n_entities, '\n')

# Confidence interval for number of entities
cat('CI for entities:', gmyc_single$entity[2], '-', gmyc_single$entity[3], '\n')

# Likelihood ratio test
# p < 0.05: significant transition from speciation to coalescent branching
cat('LR test p-value:', gmyc_single$p.value[1], '\n')
if (gmyc_single$p.value[1] < 0.05) {
    cat('Significant speciation-coalescent transition detected\n')
} else {
    cat('No significant transition; all samples may belong to one species\n')
}

# --- Step 4: Extract species assignments ---
species_assignments <- spec.list(gmyc_single)
cat('\nSpecies assignments:\n')
print(species_assignments)

n_per_species <- table(species_assignments$GMYC_spec)
cat('\nIndividuals per species:\n')
print(n_per_species)

# --- Step 5: Run multiple-threshold GMYC ---
# Multiple thresholds: different transition points across the tree
# Better for datasets with variable Ne or divergence rates
gmyc_multi <- gmyc(ultra_tree, method = 'multiple')

cat('\n--- Multiple-Threshold GMYC Results ---\n')
summary(gmyc_multi)

cat('ML entities (multiple):', gmyc_multi$entity[1], '\n')
cat('LR test p-value:', gmyc_multi$p.value[1], '\n')

# --- Step 6: Compare single vs multiple threshold ---
cat('\n--- Method Comparison ---\n')
cat('Single threshold:', gmyc_single$entity[1], 'species\n')
cat('Multiple threshold:', gmyc_multi$entity[1], 'species\n')

# --- Step 7: Visualization ---
# Color-coded tree by species partition
pdf('gmyc_single_tree.pdf', width = 10, height = max(12, Ntip(tree) * 0.15))
plot(gmyc_single, cex = 0.5)
title('Single-Threshold GMYC Species Delimitation')
dev.off()

# Support surface (likelihood across threshold positions)
pdf('gmyc_support_surface.pdf', width = 8, height = 5)
plot(gmyc_single, type = 'support')
title('GMYC Likelihood Support Surface')
dev.off()

# Multiple-threshold tree
pdf('gmyc_multi_tree.pdf', width = 10, height = max(12, Ntip(tree) * 0.15))
plot(gmyc_multi, cex = 0.5)
title('Multiple-Threshold GMYC Species Delimitation')
dev.off()

# --- Step 8: Export results ---
write.csv(species_assignments, 'gmyc_species_assignments.csv', row.names = FALSE)
cat('\nResults written to gmyc_species_assignments.csv\n')
