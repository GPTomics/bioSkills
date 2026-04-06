---
name: bio-phylo-species-trees
description: Estimate species trees using coalescent methods including ASTRAL-III, wASTRAL, ASTRAL-Pro, SVDQuartets, and BPP. Use when multi-locus data shows gene tree discordance from incomplete lineage sorting, when in the anomaly zone where concatenation is misleading, or when computing concordance factors to assess topological support.
tool_type: mixed
primary_tool: ASTRAL-III
---

## Version Compatibility

Reference examples tested with: ASTER 1.15+ (ASTRAL-III/wASTRAL/ASTRAL-Pro), IQ-TREE 2.2+ (concordance factors), PAUP* 4.0a168+ (SVDQuartets)

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `astral --version` then `astral --help` to confirm flags
- CLI: `iqtree2 --version` then `iqtree2 --help` to confirm flags
- CLI: PAUP* version displayed at startup

If commands fail, introspect the installed version and adapt flags rather than retrying.

# Coalescent-Based Species Tree Estimation

**"Estimate a species tree from multi-locus data"** -> Infer the species-level phylogeny accounting for gene tree discordance caused by incomplete lineage sorting (ILS), using summary coalescent or full-likelihood methods.
- CLI: `astral -i gene_trees.tre -o species_tree.tre` (ASTER/ASTRAL-III)
- CLI: `iqtree2 -t species.tre --gcf gene_trees.tre -s concat.fasta --scfl 100` (concordance factors)
- CLI: `svdquartets nquartets=all bootstrap=standard nreps=100` (PAUP*)

## When Concatenation Fails

The multispecies coalescent (MSC) models gene tree discordance arising from incomplete lineage sorting. Concatenation assumes all genes share the same tree, an assumption violated in three key scenarios:

1. **Short internal branches relative to population size**: When internodes are short (in coalescent units), lineages fail to coalesce within species boundaries, producing ILS. The shorter the internode and the larger the effective population size, the more discordance.

2. **The anomaly zone**: A region of tree space where the most common gene tree topology does NOT match the species tree. Concatenation becomes statistically inconsistent in the anomaly zone; more data makes the wrong answer more confident. This occurs when successive internal branches are both very short.

3. **Rapid radiations**: Short successive speciation events maximize the proportion of tree space in the anomaly zone. Phylogenomic datasets spanning rapid radiations routinely show >50% of gene trees conflicting with the species tree.

**Practical signal**: If gene concordance factor (gCF) values frequently fall below 50% across backbone nodes, coalescent methods should be considered. Concordance factors below 33% at a node indicate that a different resolution has more support than the species tree at that branch.

## Concatenation vs Coalescent Decision

| Condition | Recommended Approach |
|-----------|---------------------|
| All gCF > 70% | Concatenation is reliable; coalescent and concatenation will agree |
| Some gCF 30-70% | Run both; compare topologies and report concordance factors |
| Many gCF < 30% | Coalescent methods required; concatenation is unreliable |
| Known hybridization/introgression | Neither pure method works; use network approaches |

When gCF and species tree disagree at a node, the coalescent analysis should be given priority, as it is statistically consistent under the MSC while concatenation is not.

## Method Comparison

| Tool | Input | Speed | Key Features |
|------|-------|-------|--------------|
| ASTRAL-III | Single-copy gene tree topologies | Fast | Most widely used; handles hundreds of genes |
| wASTRAL | Gene trees with branch lengths/support | Fast | Uses gene tree uncertainty for improved accuracy |
| ASTRAL-Pro 2/3 | Multi-copy gene family trees | Fast | Handles gene duplication and loss |
| CASTER | Multiple sequence alignments | Fast | Bypasses gene tree estimation entirely |
| SVDQuartets | Concatenated alignment (site patterns) | Fast | In PAUP*; works from site patterns without gene trees |
| BPP | Multi-locus alignments | Very slow | Most accurate; jointly estimates species tree + divergences; species delimitation |
| StarBEAST3 | Multi-locus alignments | Very slow | Integrates with BEAST2 for dating |

**ASTER package (2025)**: Consolidates ASTRAL-III, wASTRAL, ASTRAL-Pro 2/3, CASTER, and WASTER into a single C++ package. Much faster and more memory-efficient than the original Java ASTRAL implementations. Install from https://github.com/chaoszhang/ASTER.

## ASTRAL Species Tree Pipeline

### Infer Gene Trees and Estimate Species Tree

**Goal:** Estimate a species tree from multi-locus sequence data by first inferring per-locus ML gene trees, then summarizing them under the coalescent model.

**Approach:** Run IQ-TREE2 independently on each locus alignment with model selection and bootstrap, collect all gene trees, then run ASTRAL to find the species tree that agrees with the largest number of quartet topologies from the gene trees.

```bash
# 1. Infer gene trees (one per locus, single-threaded for parallelism across loci)
for f in loci/*.fasta; do
    iqtree2 -s "$f" -m MFP -B 1000 -bnni -T 1 --prefix "${f%.fasta}"
done

# 2. Collect gene trees into a single file
cat loci/*.treefile > gene_trees.tre

# 3. Run ASTRAL species tree estimation (via ASTER)
# Local posterior probability (pp) replaces bootstrap for branch support
astral -i gene_trees.tre -o species_tree.tre

# 4. wASTRAL for datasets with noisy gene trees (uses branch length information)
wastral -i gene_trees.tre -o species_tree_wastral.tre
```

### Compute Concordance Factors

**Goal:** Quantify how much gene-level and site-level data supports each branch in the species tree.

**Approach:** Map each gene tree and each informative site onto the species tree, counting concordant vs discordant quartets at every internal branch.

```bash
# Gene concordance factors (gCF) + likelihood-based site concordance factors (sCF)
iqtree2 -t species_tree.tre --gcf gene_trees.tre -s concat.fasta --scfl 100

# Output: species_tree.tre.cf.tree  (tree with gCF/sCF annotations)
#         species_tree.tre.cf.stat  (per-branch statistics)
```

## Interpreting Concordance Factors

| Metric | Range | Meaning |
|--------|-------|---------|
| gCF | 0-100% | Percentage of decisive gene trees containing this branch |
| sCF | ~33-100% | Percentage of decisive sites supporting this branch (three possible resolutions) |

### Interpretation Guide

- **gCF = 100%, sCF = 100%**: Complete concordance across all loci and sites (rare in practice).
- **gCF > 50%, sCF > 50%**: Majority support; branch is well-supported under the coalescent.
- **gCF ~ 33%, sCF ~ 33%**: Completely equivocal; three possible resolutions are equally likely. No topological information at this node.
- **gCF ~ 0%**: No gene tree contains this branch. Could indicate severe ILS, introgression, or systematic gene tree estimation error.
- **gCF much lower than sCF**: Gene tree estimation error is inflating discordance. Individual locus alignments lack sufficient signal to resolve gene trees correctly. Consider wASTRAL, which accounts for gene tree uncertainty.
- **sCF < 33%**: A different resolution is better supported than the current species tree at this node. The species tree topology may be incorrect here.

### Distinguishing Quartet Asymmetry

For any internal branch, gene trees partition into three quartet topologies (q1, q2, q3). Under ILS alone, the two minor quartets should be approximately equal (q2 ~ q3). A significant excess of one minor quartet over the other suggests introgression or hybridization rather than ILS.

## Gene Tree Estimation Error vs Biological Discordance

Short alignments produce noisy gene trees, and low gCF may reflect estimation error rather than genuine ILS. Strategies to distinguish:

- **Compare gCF and sCF**: If gCF is much lower than sCF, gene tree estimation error is the primary driver. Sites contain signal that gene tree inference failed to capture.
- **Filter gene trees**: Remove gene trees inferred from very short alignments (<200 informative sites) or with very low bootstrap support.
- **Use wASTRAL**: Weights gene tree quartets by branch support, down-weighting poorly resolved splits.
- **Use CASTER**: Bypasses gene tree estimation entirely, working directly from alignments.

## SVDQuartets (Site-Pattern Method)

When individual locus alignments are too short for reliable gene tree estimation, SVDQuartets works directly from site patterns in the concatenated alignment:

```
# In PAUP* (interactive or via script)
execute concat.nex;
svdquartets evalQuartets=all bootstrap=standard nreps=100;
savetrees file=svdq_species.tre;
```

SVDQuartets is statistically consistent under the MSC and avoids gene tree estimation entirely. Best suited for RADseq/UCE datasets where per-locus alignments are short.

## BPP for Species Delimitation and Species Trees

BPP (Bayesian Phylogenetics and Phylogeography) is the most statistically rigorous coalescent method but computationally expensive:

- **A01 analysis**: Species delimitation on a fixed guide tree. Tests which populations represent distinct species.
- **A11 analysis**: Joint species tree estimation and species delimitation. Simultaneously infers the tree and tests species boundaries.
- Appropriate for small numbers of closely related species/populations (typically <20 species, <50 loci due to computational cost).
- Requires an Imap file mapping individuals to populations/species.

## ASTRAL Branch Support and Branch Lengths

ASTRAL reports local posterior probabilities (localPP) rather than bootstrap values:

| localPP | Interpretation |
|---------|----------------|
| >= 0.95 | Strong support |
| 0.7-0.95 | Moderate support |
| < 0.7 | Weak support |

ASTRAL branch lengths are in coalescent units (not substitutions per site). A branch length of 1 coalescent unit means divergence equals one effective population size generation. Short branches (<0.5 CU) indicate high expected ILS at that node.

## Network Phylogenetics (When Trees Are Insufficient)

When gene tree discordance cannot be explained by ILS alone, reticulate evolution (hybridization, introgression, horizontal gene transfer) may be involved.

**Detection**: If one minor quartet frequency is significantly greater than the other (q2 >> q3 or q3 >> q2) across multiple branches, this violates the ILS-only expectation of q2 ~ q3.

- **D-statistics (ABBA-BABA)**: Quick initial test for introgression between specific taxa. Available in genomics-population or via custom scripts.
- **PhyloNetworks.jl (SNaQ)**: Scalable pseudolikelihood network inference; estimates the number and placement of hybridization events.
- **SplitsTree**: Exploratory visualization of conflicting phylogenetic signal via split networks.

When introgression is detected, report the species tree alongside the network; the tree still represents the dominant vertical signal.

## Practical Workflow Summary

1. Infer per-locus gene trees (IQ-TREE2 with `-m MFP -B 1000 -bnni`)
2. Compute concordance factors to assess discordance levels
3. If gCF > 70% everywhere, concatenation and coalescent methods will agree
4. If gCF < 50% at backbone nodes, run ASTRAL and trust the coalescent result
5. Check quartet asymmetry for introgression signal
6. For short loci (RADseq/UCE), consider SVDQuartets or CASTER instead
7. For species delimitation with few taxa, BPP provides the gold standard

## Related Skills

- modern-tree-inference - ML gene tree inference needed as input for ASTRAL
- bayesian-inference - Bayesian tree inference; StarBEAST3 for co-estimation of species tree and gene trees
- phylogenetics/divergence-dating - Dating divergences after species tree is resolved
- tree-manipulation - Rooting and manipulating species trees
- tree-visualization - Visualizing species trees with concordance factor annotations
- alignment/multiple-alignment - Alignment quality directly affects gene tree accuracy
