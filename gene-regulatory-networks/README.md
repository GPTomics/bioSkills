# gene-regulatory-networks

## Overview

Infer and analyze gene regulatory networks from expression and chromatin data. Covers co-expression network analysis, transcription factor regulon discovery, multiomics GRN inference, perturbation simulation, and differential network comparison.

**Tool type:** python | **Primary tools:** pySCENIC, SCENIC+, WGCNA, CellOracle, DiffCorr

## Skills
| Skill | Description |
|-------|-------------|
| coexpression-networks | Build weighted co-expression networks and identify gene modules with WGCNA |
| scenic-regulons | Infer TF regulons from scRNA-seq with pySCENIC |
| multiomics-grn | Enhancer-driven GRNs from paired scRNA+scATAC with SCENIC+ |
| perturbation-simulation | Simulate TF perturbation effects on cell state with CellOracle |
| differential-networks | Compare co-expression networks between conditions with DiffCorr |

## Example Prompts
- "Build a co-expression network from my RNA-seq data and find hub genes"
- "Identify transcription factor regulons in my single-cell data"
- "Infer gene regulatory networks from my 10x multiome data"
- "Simulate what happens if I knock out this transcription factor"
- "Compare co-expression networks between disease and control"

## Requirements
```bash
# Python
pip install pyscenic celloracle loompy

# R
install.packages('WGCNA')
install.packages('CEMiTool')
install.packages('DiffCorr')
BiocManager::install('hdWGCNA')
```

## Related Skills

- **single-cell** - Preprocessing and clustering for scRNA-seq inputs
- **chip-seq** - TF binding data for GRN validation
- **atac-seq** - Chromatin accessibility for regulatory region identification
- **differential-expression** - DE analysis for network gene prioritization
- **pathway-analysis** - Functional enrichment of network modules
