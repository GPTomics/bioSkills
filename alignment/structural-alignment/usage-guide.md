# Structural Alignment - Usage Guide

## Overview

This skill covers backbone-aware structural alignment: pairwise (TM-align, US-align, Bio.PDB.Superimposer), database search at scale (Foldseek 3Di), structural multiple alignment (Foldmason, PROMALS3D, T-Coffee Expresso), and pLM embedding aligners (TM-Vec, vcMSA). Use it when sequence identity drops below ~25% (twilight zone), when remote-homology detection is required, or when fold-similarity quantification is the goal.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda foldseek tmalign usalign foldmason
# OR direct downloads
# Foldseek: https://github.com/steineggerlab/foldseek
# TM-align: https://zhanggroup.org/TM-align/
# US-align: https://zhanggroup.org/US-align/
# Foldmason: https://github.com/steineggerlab/foldmason

# Python
pip install biopython

# Optional: pLM aligners (heavier installs; GPU strongly recommended)
pip install fair-esm  # ESM-2 embeddings
# vcMSA: https://github.com/clairemcwhite/vcmsa
# DEDAL: https://github.com/google-research/google-research/tree/master/dedal
# TM-Vec: https://github.com/valentynbez/tmvec (active fork; original at tymor22/tm-vec is in limited maintenance)
```

For Foldseek searches against AlphaFoldDB, download the database:

```bash
foldseek databases Alphafold/UniProt /path/to/afdb tmp/
foldseek databases PDB100 /path/to/pdb tmp/
```

## Quick Start

Tell your AI agent what you want to do:
- "Compute the TM-score between these two PDB structures"
- "Search this AlphaFold model against AlphaFoldDB to find structural homologs"
- "Build a structural MSA from these 50 PDB structures with Foldmason"
- "These proteins share <15% identity; align them via predicted structures instead of sequence"

## Example Prompts

### Pairwise Structural Alignment
> "Align these two crystal structures with TM-align and report the TM-score and superposition"

> "I have two homologs at 18% sequence identity. Compute structural alignment with US-align and tell me whether they share a fold."

> "Compare apo and holo conformations of the same protein and quantify the conformational change."

### Database Search at Scale
> "Search my predicted structure against AlphaFoldDB to find evolutionary distant homologs"

> "I have a metagenomic protein with no PDB hit. Predict its structure with ESMFold and search the ESM Atlas."

> "Cluster all PDB structures sharing TM-score > 0.5 with my query."

> "Search this antibody-antigen complex against AFDB-Multimer with Foldseek-Multimer and rank by complex TM-score"

### Structural Multiple Alignment
> "Build an MSA from these 30 PDB structures using Foldmason"

> "Generate a structure-aware MSA combining sequence and PDB templates with PROMALS3D"

> "I need a structural MSA for tree inference from a remote homologue family"

### Hybrid Sequence-Structure
> "Run T-Coffee Expresso to use PDB templates wherever they exist for my sequences"

> "Improve my sequence MSA by adding structural information from AlphaFold predictions"

### pLM Aligners
> "These proteins are 12% identical; align them with a protein language model aligner"

> "Use TM-Vec to predict TM-scores between sequences without computing structure"

## What the Agent Will Do

1. Assess whether structural alignment is the right tool (sequence identity, structure availability, downstream goal)
2. Choose pairwise vs database-search vs MSA based on dataset
3. Select tool: TM-align/US-align (pairwise), Foldseek (search), Foldmason (MSA), PROMALS3D (hybrid)
4. Run with appropriate flags and output formats
5. Interpret TM-score, RMSD, and LDDT in context (TM > 0.5 = same fold; TM > 0.8 = equivalent topology)
6. Optionally derive a sequence MSA from the structural MSA for downstream phylogenetics

## Tips

- TM-score normalised by the shorter chain is the standard fold-similarity metric; report it alongside RMSD and alignment length
- RMSD alone is misleading (length-dependent, outlier-sensitive); never report RMSD without TM-score and aligned residue count
- US-align extends TM-align to multi-chain complexes, RNA, DNA, and protein-NA mixed assemblies; Foldseek-Multimer (Kim et al 2025) is the database-scale successor (10-100x faster pairwise, 1000-10000x faster in DB search)
- DALI Z-score is the right metric for ranking twilight-fold hits (TM 0.3-0.5 regime); TM-score saturates there. Z>20 = definitely homologous, Z 8-19 = probable, Z 2-8 = candidate (verify with biology), Z<2 = noise
- Foldseek's 3Di alphabet provides ~1000x speedup over TM-align with comparable sensitivity in the same-fold regime
- For dark-proteome problems where sequence identity is below 15%, predict structures with ColabFold or ESMFold first, then run Foldseek
- AlphaFoldDB clusters (Barrio-Hernandez et al 2023) are the canonical entry point for AFDB-scale remote-homology search
- pLM aligners (TM-Vec, vcMSA, DEDAL) complement structural alignment when no structures are available; they are NOT a substitute when structures exist
- PyMOL `super` is for distantly related structures (cycle-fitting); `align` is for closer matches; `cealign` runs CE
- TM-score < 0.2 is statistically random; do not interpret low TM as "weak homology"

## Related Skills

- alignment/multiple-alignment - Sequence MSA tools (MAFFT, MUSCLE5, ClustalOmega) for >25% identity
- alignment/pairwise-alignment - Sequence pairwise alignment when twilight zone is not crossed
- alignment/msa-statistics - Per-column statistics on structurally-derived MSAs
- alignment/alignment-trimming - Trim structural MSAs before phylogenetics
- phylogenetics/modern-tree-inference - Build trees from structural MSAs
