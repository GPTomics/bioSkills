# Gene Tree Species Tree Reconciliation - Usage Guide

## Overview

Reconciliation methods compare a gene tree against a species tree under explicit probabilistic models of duplication (D), horizontal transfer (T), and loss (L). The output is a mapping of D/T/L events to species-tree branches, revealing the evolutionary history of each gene family. ALE / GeneRax / AleRax (Szöllősi 2013; Morel 2020/2024) define the modern probabilistic standard; parsimony-based alternatives (RANGER-DTL, NOTUNG) remain useful for screening.

The most consequential decision in reconciliation is **how to handle gene-tree uncertainty**. Single-ML-tree reconciliation (NOTUNG, RANGER, basic GeneRax) treats input trees as truth and amplifies gene-tree noise as spurious D/T/L events. ALE integrates over a sample of bootstrapped gene trees; AleRax co-estimates everything jointly. For publication-grade DTL inference, ALE or AleRax are required.

## Prerequisites

```bash
# Probabilistic reconciliation
conda install -c bioconda ale generax
git clone --recursive https://github.com/BenoitMorel/AleRax && cd AleRax && ./install.sh

# Parsimony / alternative methods
# RANGER-DTL: https://compbio.engr.uconn.edu/software/RANGER-DTL/
# NOTUNG: https://www.cs.cmu.edu/~durand/Notung/download

# Bayesian + WGD
julia -e 'using Pkg; Pkg.add("Whale")'

# Tree manipulation
conda install -c bioconda newick_utils iqtree
pip install ete4
```

## Quick Start

Tell the AI agent what to do:
- "Reconcile these orthogroup gene trees against the species tree using ALE to identify per-branch D/T/L events"
- "Use AleRax to co-estimate gene trees, species tree, and DTL rates from a set of UFBoot bootstrap samples"
- "Apply Whale.jl to perform Bayesian DTL inference with explicit WGD nodes for a vertebrate species tree"
- "Run RANGER-DTL parsimony sweep across cost weights to identify robust transfer events"

## Example Prompts

### Bacterial Phylogenomics

> "I have UFBoot trees for 5000 orthogroups across 100 bacterial genomes and want to identify horizontal gene transfer events per genome and rates of DTL. Run ALEml_undated on all gene-tree distributions against the rooted species tree, then aggregate the per-branch transfer counts into a summary matrix. Include the optional ALE-rooting (Williams 2017) to validate the species tree root."

### Eukaryote Gene Family with WGD

> "Analyze ohnolog gene-family history across teleost fish, accounting for the 3R WGD. Use Whale.jl with the WGD node placed at the teleost root; report per-WGD retention probabilities and post-WGD duplication / loss rates per family. Cross-validate the WGD detection by inspecting the per-family Bayes factor for the WGD-vs-non-WGD model."

### Refining Gene Trees

> "I have orthogroup gene trees from OrthoFinder but they have low bootstrap support. Use GeneRax with `--strategy SPR` to refine each tree against the species tree under a DL model, producing reconciled gene trees and per-family duplication-loss event histories."

## What the Agent Will Do

1. **Validate input format**: Check species labels match between species tree and gene tree leaves; verify gene IDs encode species via prefix
2. **Build per-orthogroup gene-tree bootstrap distribution** with IQ-TREE2 UFBoot if not provided
3. **Encode gene-tree samples** for ALE (`ALEobserve` produces `.ale` files)
4. **Run reconciliation** with appropriate model:
   - ALEml_undated for posterior over D/T/L
   - GeneRax for ML reconciliation with SPR refinement
   - AleRax for joint co-estimation
   - Whale.jl for Bayesian DTL + WGD
5. **Aggregate per-branch event counts** across families
6. **Cross-validate** with parsimony (RANGER-DTL, NOTUNG) for screening
7. **Report**: per-branch event posteriors, family-level reconciled trees, species-tree refinements
8. **Caveats**: explicit ILS vs HGT note; ghost-lineage caveat; WGD if applicable

## Tips

- ALE undated is the default for bacterial DTL; dated for time-calibrated trees
- For gene-tree-error-resistant inference, use ALE with UFBoot bootstrap samples (>= 100 per family)
- AleRax (Morel 2024) co-estimates gene tree + species tree + DTL rates; preferred for new publication-grade work
- Species labels must match exactly between species tree and gene trees; ALE / GeneRax silently fail on mismatch
- Use `nw_labels -I` (newick utilities) to verify species sets before reconciliation
- For HGT-affected clades (bacteria, archaea), transfers >> duplications; opposite for opisthokonts
- ILS confounded with HGT at short internal branches; use ecceTERA or DLCpar for ILS-aware inference
- Multifurcations in species tree break ALE / GeneRax; resolve via `ape::multi2di()` or randomization
- For WGD lineages, use Whale.jl with explicit WGD branch placement
- RANGER-DTL parsimony is fine for screening; results sensitive to cost weights; sensitivity sweep required
- NOTUNG is duplication-loss only by default; use NOTUNG-HGT extension for HGT-affected; or ALE
- Gene-ID separator convention: `species|gene_id` for ALE / GeneRax; `species_gene_id` for Whale.jl
- AleRax can detect species tree root mis-specification; trust its rooting
- GeneRax MPI parallelization requires careful MPI launcher configuration on HPC

## Related Skills

comparative-genomics/hgt-detection - HGT inference uses DTL reconciliation
comparative-genomics/ortholog-inference - OrthoFinder HOGs feed reconciliation
comparative-genomics/gene-family-evolution - CAFE5 birth-death complementary to per-family reconciliation
comparative-genomics/whole-genome-duplication - Whale.jl native WGD modeling
comparative-genomics/ancestral-reconstruction - DTL informs ancestral gene-content inference
phylogenetics/modern-tree-inference - UFBoot bootstrap gene trees
phylogenetics/bayesian-inference - MrBayes / RevBayes alternative
alignment/multiple-alignment - Quality MSA precedes gene-tree inference
