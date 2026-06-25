# Structure Probing - Usage Guide

## Overview
Process experimental RNA structure probing data (SHAPE-MaP, DMS-MaPseq) into per-nucleotide reactivity profiles with ShapeMapper2, then use the reactivities as soft restraints on thermodynamic folding. Reactivity reports nucleotide flexibility/accessibility, not base pairing per se, and it is a constraint on folding, not a structure. A single profile is a population average; multi-conformation RNAs need per-read clustering.

## Prerequisites
```bash
# ShapeMapper2 (Linux only; use Docker on macOS)
conda install -c bioconda shapemapper2
docker pull shapemapper2/shapemapper2   # macOS alternative

# DMS-MaPseq and multi-conformation clustering
pip install seismic-rna

# RNA Framework (alternative MaP/RT-stop pipeline)
conda install -c bioconda rnaframework

# ViennaRNA for reactivity-restrained folding
conda install -c bioconda viennarna

# Python dependencies
pip install matplotlib pandas numpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Process my SHAPE-MaP reads into reactivity profiles"
- "Use my SHAPE reactivities to constrain folding"
- "My data is DMS-MaPseq, handle the A/C-only signal correctly"
- "My reactivity profile fits no single structure; find the conformations"
- "Is this protected region base-paired or protein-bound?"

## Example Prompts

### SHAPE-MaP Analysis
> "I have SHAPE-MaP paired-end reads for modified and untreated samples targeting my RNA. Run ShapeMapper2 to get normalized reactivities."

> "Process my amplicon SHAPE-MaP data and plot the reactivity profile."

### Reactivity-Restrained Folding
> "Use my SHAPE reactivities to restrain RNAfold and compare with the unrestrained structure."

> "I have DMS-MaPseq data; fold with the right A/C-only parameters and mask G/U."

### Multiple Conformations
> "My reactivity profile looks inconsistent with one structure. Cluster the reads to find coexisting conformations."

### In-Cell Interpretation
> "I probed in cells; help me tell apart base-pairing from protein protection."

## What the Agent Will Do
1. Identify the readout (MaP vs RT-stop) and reagent (SHAPE vs DMS) and pick the matching tool/scoring
2. Run ShapeMapper2 (or SEISMIC-RNA / RNA Framework) with modified, untreated, and optional denatured samples
3. Assess QC (effective depth, untreated vs modified mutation rates) and carry low-depth positions as no-data
4. Use normalized reactivities to restrain ViennaRNA folding with reagent-correct parameters
5. Flag multi-conformation or in-cell occupancy cases that a single profile cannot resolve

## Tips
- **Reactivity is flexibility, not pairing** - High reactivity = flexible/accessible; low = constrained, which can mean base-paired OR tertiary-contacted, protein-bound, or stacked. Do not equate low reactivity with "base-paired".
- **MaP vs RT-stop** - Mutational profiling (ShapeMapper2, SEISMIC-RNA) encodes modifications as mutations; RT-stop methods (icSHAPE, Structure-seq, DMS-seq) encode them as truncations. They need different scoring; do not mix pipelines.
- **Normalize per transcript** - Raw reactivities are on an arbitrary, experiment-specific scale; the 2-8%/box-plot normalization makes them comparable WITHIN a transcript only. Never pool or compare raw reactivities across transcripts/experiments; compare conditions with delta-SHAPE and standard errors.
- **Parameter attribution** - The m=1.8, b=-0.6 SHAPE pair is Hajdin et al. 2013 (the ViennaRNA default), not Deigan et al. 2009's own m=2.6/b=-0.8. There is no separate DMS-specific standard pair; reuse 1.8/-0.6 (Cordero et al. 2012 showed SHAPE parameters transfer to DMS) or tune, applying the restraint only to A/C and masking G/U to -999.
- **Controls** - The untreated control is mandatory (background subtraction); the denatured control improves sequence-bias normalization.
- **Depth and no-data** - Aim for >=5000 effective depth; carry low-depth nucleotides as -999, never as 0.
- **Multiple conformations** - If a profile fits no single structure, cluster MaP reads with SEISMIC-RNA / DREEM rather than forcing one fold.
- **In-cell** - In cells, protection can mean protein/ligand occupancy, and mRNA is actively unfolded; the in-cell-minus-in-vitro difference maps binding.
- **macOS** - ShapeMapper2 is Linux-only; use Docker or Singularity.

## Related Skills
- secondary-structure-prediction - The folding engine the reactivities restrain
- ncrna-search - Identify the RNA family and a CM consensus structure to probe against
- covariation-analysis - Independent evolutionary evidence for the suggested pairs
- epitranscriptomics/m6a-peak-calling - RNA modifications that confound DMS/SHAPE reactivity
- clip-seq/binding-site-annotation - In-cell protection as an RBP footprint
- read-qc/quality-reports - QC of the underlying sequencing reads
