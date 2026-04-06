---
name: bio-phylo-modern-tree-inference
description: Build maximum likelihood phylogenetic trees using IQ-TREE2 and RAxML-NG with expert model selection, branch support assessment, and topology testing. Use when inferring publication-quality ML trees, selecting substitution models, interpreting bootstrap and concordance factor support, or running partitioned phylogenomic analyses.
tool_type: cli
primary_tool: IQ-TREE2
---

## Version Compatibility

Reference examples tested with: IQ-TREE 2.2+, RAxML-NG 1.2+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `iqtree2 --version` then `iqtree2 --help` to confirm flags
- CLI: `raxml-ng --version` then `raxml-ng --help` to confirm flags

If commands fail, introspect the installed version and adapt flags rather than retrying.

# Modern ML Tree Inference

**"Infer a maximum likelihood tree from my alignment"** -> Build publication-quality ML trees with automatic substitution model selection, ultrafast bootstrap support, concordance factors, and topology testing.
- CLI: `iqtree2 -s alignment.fa -m MFP -B 1000 -alrt 1000 -bnni` (IQ-TREE2)
- CLI: `raxml-ng --all --msa alignment.fa --model GTR+G --bs-trees 100` (RAxML-NG)

## IQ-TREE2 vs RAxML-NG Decision

| Factor | IQ-TREE2 | RAxML-NG |
|--------|----------|----------|
| Model selection | Built-in ModelFinder | External ModelTest-NG |
| Ultrafast bootstrap | UFBoot2 | No |
| Branch lengths | Good | More accurate |
| Concordance factors | Built-in (gCF/sCF) | No |
| Very large trees (>1000 taxa) | Good | Better |
| Transfer bootstrap | No | Yes |
| Partition models | Extensive | Good |

**Default recommendation:** IQ-TREE2 for most workflows (integrated model selection, UFBoot2, concordance factors). Use RAxML-NG when precise branch lengths matter, for very large trees, or when using transfer bootstrap for rogue-taxon-prone datasets.

**Best practice for important results:** Run both tools and compare topologies.

## Model Selection

### Use ModelFinder (`-m MFP`), Not `-m TEST`

**Goal:** Select the substitution model that best describes the evolutionary process in the alignment.

**Approach:** ModelFinder (`-m MFP`) tests standard models plus FreeRate (+R) models and performs concurrent model-tree search. The older `-m TEST` does not test FreeRate models and tests on a fixed tree.

```bash
# Recommended: ModelFinder Plus (includes FreeRate models)
iqtree2 -s alignment.fasta -m MFP -B 1000 -alrt 1000 -bnni -T AUTO

# Model selection only (no tree inference)
iqtree2 -s alignment.fasta -m MF -T AUTO

# Partition model with automatic merging
iqtree2 -s concat.fasta -p partitions.nex -m MFP+MERGE -B 1000 -bnni -T AUTO
```

### Rate Heterogeneity Models

| Model | Description | When Selected |
|-------|-------------|---------------|
| +G4 | Discrete gamma (4 categories) | Standard default; sufficient for most datasets |
| +I+G4 | Invariant sites + gamma | Often selected despite theoretical identifiability concerns; safe to use |
| +R4/+R5 | FreeRate model | Better fit for large datasets; relaxes gamma assumption |
| +R (auto k) | FreeRate with automatic categories | **Only tested by `-m MFP`**, not `-m TEST` |

FreeRate models can absorb rate variation from long tails that discretized gamma cannot. For large datasets or datasets with complex rate variation, FreeRate often fits better.

### BIC vs AIC for Model Selection

BIC is the IQ-TREE default and recommended for most analyses. BIC penalizes complexity more heavily than AIC, reducing overfitting risk. AIC tends to select overly complex models.

### DNA Model Hierarchy

| Model | Free Parameters | When Appropriate |
|-------|----------------|------------------|
| JC69 | 0 | Almost never in practice; null model |
| K2P/K80 | 1 (kappa) | Very closely related sequences with balanced composition |
| HKY85 | 4 | Moderate divergence, single-gene analyses |
| GTR | 8 | Default; almost always selected by model testing |

The real decision is usually not GTR vs HKY but which rate heterogeneity model (+G vs +I+G vs +R).

### Protein Models

Let ModelFinder choose. For deep phylogenies, profile mixture models (C10-C60 in IQ-TREE, CAT in PhyloBayes) can outperform fixed-matrix models by capturing site-specific amino acid preferences.

```bash
# Protein with ModelFinder
iqtree2 -s protein.fasta -m MFP -B 1000 -bnni -st AA -T AUTO

# Protein with mixture model (for deep phylogenies)
iqtree2 -s protein.fasta -m LG+C60+F+G -B 1000 -bnni -st AA -T AUTO
```

## Branch Support Assessment

### Standard Analysis: UFBoot2 + SH-aLRT

```bash
# Recommended for most analyses
iqtree2 -s alignment.fasta -m MFP -B 1000 -alrt 1000 -bnni -T AUTO --seed 12345
```

The `-bnni` flag is critical: it further optimizes each bootstrap tree with NNI, reducing overestimation from model violations. This flag is default since IQ-TREE 2.2.0 but should be specified explicitly for clarity.

### Interpreting Support Values

| UFBoot2 | SH-aLRT | Interpretation |
|---------|---------|----------------|
| >= 95 | >= 80 | Strong support |
| 80-94 | 70-79 | Moderate support |
| < 80 | < 70 | Weak support |

**Critical nuance:** UFBoot values are NOT comparable to standard bootstrap. UFBoot >= 95 corresponds roughly to standard bootstrap >= 70. Do not apply the traditional >= 70 threshold to UFBoot values.

### When Low Support Matters

- Low support on backbone branches: genuine topological uncertainty. Investigate with concordance factors or coalescent methods
- Low support on recent divergences within a well-sampled clade: may reflect insufficient data rather than genuine uncertainty
- Low support throughout the tree: suspect rapid radiation, incomplete lineage sorting (ILS), hybridization, or inadequate data

### Transfer Bootstrap (RAxML-NG)

For large trees (>1000 taxa), transfer bootstrap expectation (TBE) is less sensitive to rogue taxa than standard bootstrap:

```bash
raxml-ng --all --msa alignment.fasta --model GTR+G --bs-trees autoMRE --bs-metric tbe
```

### Concordance Factors

Concordance factors quantify agreement among loci (gCF) and sites (sCF), complementing bootstrap:

```bash
# After obtaining gene trees and species tree
iqtree2 -t species.treefile --gcf gene_trees.treefile -s concat.fasta --scf 100

# Likelihood-based sCF (more accurate; requires recent IQ-TREE)
iqtree2 -t species.treefile --gcf gene_trees.treefile -s concat.fasta --scfl 100
```

| Metric | Interpretation |
|--------|---------------|
| gCF/sCF > 50% | Majority of loci/sites support this branch |
| gCF/sCF ~ 33% | Completely equivocal (three resolutions equally likely) |
| gCF << sCF | Gene tree estimation error, not genuine discordance |
| sCF < 33% | A different topology is better supported at this node |

For publication, report UFBoot + SH-aLRT at minimum; add concordance factors for phylogenomic datasets.

## Partitioned Analysis

For multi-gene concatenated datasets where genes evolve at different rates:

```bash
# Edge-linked proportional (recommended default)
iqtree2 -s concat.fasta -p partitions.nex -m MFP -B 1000 -bnni -T AUTO

# Edge-unlinked (independent branch lengths per partition)
# Most general but parameter-rich; risk of overfitting with missing data
iqtree2 -s concat.fasta -Q partitions.nex -m MFP -B 1000 -bnni -T AUTO

# With automatic partition merging
iqtree2 -s concat.fasta -p partitions.nex -m MFP+MERGE -B 1000 -bnni -T AUTO
```

| Flag | Model | Recommendation |
|------|-------|----------------|
| `-q` | Edge-linked equal | Unrealistic; not recommended |
| `-p` | Edge-linked proportional | **Recommended default** |
| `-Q` | Edge-unlinked | Justified when different genes have different relative rates across lineages (heterotachy) |

## Topology Testing

### AU Test (Approximately Unbiased)

**Goal:** Test whether alternative tree topologies are significantly worse than the best tree.

**Approach:** Compare candidate topologies using the AU test, which provides proper multiple-testing correction. Preferred over the overly conservative SH test.

```bash
# Compare candidate topologies
iqtree2 -s alignment.fasta -m GTR+F+R3 --trees candidates.treefile --test-au --test 10000 -n 0
```

- p-AU >= 0.05: Tree cannot be rejected
- p-AU < 0.05: Tree is significantly worse

Use topology tests for evaluating specific competing hypotheses, not for fishing through thousands of random topologies.

## Long Branch Attraction (LBA) Awareness

LBA causes distantly related long-branched taxa to group together artifactually due to model misspecification. It affects ML and Bayesian methods, not just parsimony.

**Detection signs:**
- Two long-branched taxa group together when they should not
- Removing one long-branch taxon causes the other to move
- Switching to a more complex model changes the placement

**Mitigation (in order of effectiveness):**
1. Add taxa that break long branches (most effective)
2. Use site-heterogeneous models (C60 in IQ-TREE, CAT-GTR in PhyloBayes)
3. Remove saturated sites (3rd codon positions, hypervariable regions)
4. Use amino acids instead of nucleotides for coding regions
5. RY-coding or Dayhoff-6 recoding for proteins

When LBA is suspected, consider running PhyloBayes with CAT-GTR model (see bayesian-inference skill).

## IQ-TREE2 Output Files

| File | Description |
|------|-------------|
| `.treefile` | Best ML tree (Newick) |
| `.iqtree` | Full report with model parameters |
| `.contree` | Consensus tree with support values |
| `.splits.nex` | Bootstrap splits (Nexus) |
| `.model.gz` | Model parameters |
| `.log` | Run log |
| `.ckp.gz` | Checkpoint for resuming |

## RAxML-NG Usage

```bash
# ML search + bootstrap
raxml-ng --all --msa alignment.fasta --model GTR+G --bs-trees 100

# Thorough search with multiple starting trees
raxml-ng --msa alignment.fasta --model GTR+G --tree pars{10} --prefix ml_search

# Protein models
raxml-ng --msa protein.fasta --model LG+G8+F

# Constrained tree search
raxml-ng --msa alignment.fasta --model GTR+G --tree-constraint constraint.tre

# Check alignment before full run
raxml-ng --check --msa alignment.fasta --model GTR+G
```

## RAxML-NG Output Files

| File | Description |
|------|-------------|
| `.raxml.bestTree` | Best ML tree |
| `.raxml.support` | Tree with bootstrap support |
| `.raxml.bootstraps` | All bootstrap trees |
| `.raxml.mlTrees` | All ML trees from search |
| `.raxml.log` | Analysis log |

## Large Dataset Strategies

```bash
# IQ-TREE2 fast mode for >500 taxa
iqtree2 -s large.fasta -m GTR+G -B 1000 -bnni -T 4 -mem 8G -fast

# RAxML-NG with limited starting trees
raxml-ng --msa large.fasta --model GTR+G --tree pars{5} --threads 8
```

For very large trees (>1000 taxa), consider FastTree 2 for an initial exploratory tree, then refine with RAxML-NG.

## Constrained Analysis

```bash
# Enforce monophyly constraint
iqtree2 -s alignment.fasta -m MFP -g constraint.tre -B 1000 -bnni

# Constraint file: Newick with taxa to constrain
# ((Human,Chimp),Gorilla);
```

## Reproducibility

```bash
# Always set random seed for reproducible results
iqtree2 -s alignment.fasta -m MFP -B 1000 -bnni --seed 12345 -T AUTO
raxml-ng --msa alignment.fasta --model GTR+G --seed 12345 --bs-trees 100
```

## Resuming Interrupted Runs

```bash
iqtree2 -s alignment.fasta -m MFP -B 1000 --redo-tree
raxml-ng --msa alignment.fasta --model GTR+G --redo
```

## Related Skills

- bayesian-inference - Bayesian tree inference with MrBayes, BEAST2, convergence diagnostics
- species-trees - Coalescent methods (ASTRAL) when gene tree discordance is high
- divergence-dating - Molecular clock analysis and divergence time estimation
- tree-io - Read and convert output tree files
- tree-visualization - Visualize trees with support values
- distance-calculations - Compare with distance-based methods
- alignment/alignment-io - Prepare alignments for tree inference
- alignment/multiple-alignment - Alignment quality affects tree inference
