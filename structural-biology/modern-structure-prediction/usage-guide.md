# Modern Structure Prediction

## Overview

Predict protein and complex structures with deep-learning models and interpret the confidence they report. The choice of predictor is set by the input and the question, not by novelty: ESMFold is single-sequence, no-MSA, and fast enough for metagenomic scale but lower accuracy and weakest on orphan proteins, while AlphaFold3, Chai-1, and Boltz co-fold complexes, ligands, nucleic acids, ions, and PTMs. MSA depth is the dominant accuracy determinant for the coevolution-based models, so quality tracks how well-represented a sequence's family is. A default prediction is a single dominant conformer with spatially varying reliability - not an ensemble, not a variant-effect or affinity engine, and not an experiment. It is a hypothesis to reconcile with data, read through pLDDT (local), PAE (relative positioning), and pTM/ipTM (fold and interface).

## Prerequisites

```bash
# ESMFold (local, reliable path; needs a CUDA GPU with ~16 GB)
pip install fair-esm

# Confidence parsing and comparison
pip install biopython numpy requests

# Co-folders for complexes/ligands (Linux + CUDA GPU)
pip install chai_lab      # Chai-1
pip install boltz         # Boltz-1/2

# AlphaFold2 with fast MMseqs2 MSA
pip install colabfold
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict this single-chain protein with ESMFold and summarize pLDDT"
- "Choose the right predictor for my protein-ligand complex"
- "Run a co-folder on my heterodimer and tell me if the interface is trustworthy"
- "Compare an ESMFold and an AlphaFold3 model and show where they agree"
- "Explain whether this low-confidence region is disordered or wrong"

## Example Prompts

### Single Protein Prediction
> "Predict the structure of this sequence with the fastest method and flag low-confidence regions"

> "Run ESMFold locally on my sequence and band the pLDDT scores"

> "This protein has almost no homologs - which predictor should I use and why"

### Complex and Ligand Prediction
> "Predict this protein-protein complex and gate the interface on ipTM and inter-chain PAE"

> "Co-fold my protein with this small-molecule ligand and validate the pose"

> "Run Boltz on my heterodimer and tell me whether the two chains are predicted to interact"

### Confidence and Reconciliation
> "Compare ESMFold, AlphaFold3, and Chai-1 predictions and report where they disagree"

> "This region has pLDDT below 50 - is it a modeling error or an intrinsically disordered region"

> "Each domain looks confident but the arrangement seems off - what does the PAE matrix say"

### What NOT to ask these tools
> "I mutated one residue and the structure did not change - is the mutation tolerated" (use a variant-effect tool)

> "Give me the binding affinity from this co-folded pose" (co-fold geometry is not a Kd)

## What the Agent Will Do

1. Pick a predictor from the input and question: ESMFold for single-chain speed and scale, ColabFold/AlphaFold for maximum monomer accuracy, AlphaFold3/Chai-1/Boltz for complexes, ligands, and nucleic acids.
2. Prepare the input (FASTA per chain, YAML for a Boltz affinity request, or an AlphaFold Server JSON with multiple seeds).
3. Run the prediction locally or via server, keeping MSA depth, recycles, and seeds recorded so runs stay comparable.
4. Read the confidence metrics for their distinct questions: pLDDT (local), PAE (relative positioning), pTM (global fold), ipTM (interface).
5. Gate complexes on ipTM plus the inter-chain PAE block rather than per-chain pLDDT.
6. Reconcile multiple predictions, report the aligned selection behind any RMSD, and frame the result as a hypothesis to validate.

## Tips

- Run ESMFold locally; the hosted esmatlas API is intermittently down with SSL/internal-server errors.
- A long pLDDT-below-50 stretch usually marks an intrinsically disordered region, not a modeling failure.
- High per-residue pLDDT with a high inter-domain PAE block means each domain is confident but their relative arrangement is not - do not trust the linker.
- Judge a complex on ipTM (interface) and the inter-chain PAE, never on per-chain pLDDT; ipTM below ~0.6 is unreliable or the chains likely do not interact, 0.6-0.8 is uncertain (let the inter-chain PAE block decide), and above ~0.8 is confident.
- These models are insensitive to single point mutations - do not use them for variant effect, ddG, or pathogenicity; use AlphaMissense, FoldX, or ESM.
- No co-folder gives a trustworthy Kd; Boltz-2's affinity output is a screening prior on a log10(IC50) scale, not a measured constant.
- AF3-class diffusion models can hallucinate order in disordered regions and violate chirality - validate every ligand pose with PoseBusters.
- A single prediction is one dominant conformer; MSA subsampling and AF-Cluster sample alternate states but are unreliable hypotheses.
- Two predictions are not comparable if MSA depth, recycles, seeds, or templates differ - record and report the settings.
- Verify chai-lab and boltz invocation with `--help`; these packages and their output filenames change quickly.

## Related Skills

- alphafold-predictions - Retrieve precomputed AlphaFold DB models and read their pLDDT/PAE
- structure-io - Parse and write predicted PDB/mmCIF files
- geometric-analysis - RMSD, superposition, and TM-score caveats for comparing models
- structure-navigation - Walk chains/residues/atoms in a predicted structure
- alignment/structural-alignment - Structure-based alignment before comparing sequence-different models
- chemoinformatics/virtual-screening - Dock into a predicted pocket (inherits predicted rotamer/backbone error)
- chemoinformatics/ml-docking-rescoring - Rescore co-folded poses; co-fold geometry is not affinity
