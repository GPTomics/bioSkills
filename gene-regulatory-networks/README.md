# gene-regulatory-networks

## Overview

Infer and analyze gene regulatory networks from expression and chromatin data. Covers co-expression module discovery, bulk GRN inference and TF protein-activity, single-cell transcription-factor regulon discovery, enhancer-driven multiomic GRN inference, in silico perturbation simulation, and differential network comparison. A unifying theme runs through every skill: an expression-derived network is undirected statistical association, not causal regulation, unless direction is imported from motifs, accessibility, time, or perturbation.

**Tool type:** mixed | **Primary tools:** pySCENIC, SCENIC+, WGCNA, CellOracle, VIPER, DiffCorr

## Skills
| Skill | Description |
|-------|-------------|
| coexpression-networks | Build weighted co-expression networks and identify gene modules with WGCNA/hdWGCNA (marginal vs partial correlation) |
| grn-inference | Bulk GRN inference (ARACNe/GENIE3/GRNBoost2) and TF protein-activity with VIPER (activity, not edges) |
| scenic-regulons | Infer single-cell TF regulons with pySCENIC (motif-pruning as directionality) |
| multiomics-grn | Enhancer-driven eGRNs from paired/unpaired scRNA+scATAC with SCENIC+ |
| perturbation-simulation | Simulate TF perturbation effects on cell state with CellOracle/Dynamo (direction, not magnitude) |
| differential-networks | Compare networks between conditions with DiffCorr/DINGO (connectivity is not expression) |

## Example Prompts
- "Build a signed co-expression network from my RNA-seq data and find hub genes by kME"
- "Infer a regulatory network from bulk RNA-seq and find master regulators with VIPER"
- "Identify transcription factor regulons in my single-cell data"
- "Infer enhancer-driven gene regulatory networks from my 10x multiome data"
- "Simulate what happens if I knock out this transcription factor"
- "Compare co-expression networks between disease and control"

## Requirements
```bash
# Python
pip install pyscenic celloracle loompy arboreto
```

```r
# R
install.packages(c('WGCNA', 'DiffCorr', 'GeneNet'))
BiocManager::install(c('CEMiTool', 'hdWGCNA', 'viper', 'GENIE3', 'iDINGO'))
```

ARACNe-AP (Java) is built from source; SCENIC+ uses its current Snakemake workflow plus cisTarget databases.

## Related Skills

- **single-cell** - Preprocessing, clustering, and doublet removal for scRNA-seq inputs
- **chip-seq** - TF binding data for GRN validation and motif analysis
- **atac-seq** - Chromatin accessibility for regulatory region identification
- **differential-expression** - DE analysis and signatures for network gene prioritization and VIPER
- **pathway-analysis** - Functional enrichment of network modules and regulons
- **causal-genomics** - Population-genetics route to regulatory causality (contrast with perturbation)
