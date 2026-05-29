# flow-cytometry

## Overview

Flow, spectral, and mass cytometry (CyTOF) analysis from raw FCS files to differentially abundant cell populations. The category encodes the decisions that separate a defensible analysis from a misleading one: read/compensate/transform correctly, QC on the Time axis in the right order, discriminate doublets by pulse geometry, choose manual gating versus unsupervised clustering, and test between groups with the sample (not the cell) as the experimental unit.

**Tool type:** mixed (R/Bioconductor-centric, with Python readers) | **Primary tools:** flowCore, flowWorkspace, CATALYST, diffcyt, FlowSOM, CytoNorm

## Skills

| Skill | Description |
|-------|-------------|
| fcs-handling | Read/write FCS (conventional, spectral, CyTOF), map channels, choose R vs Python readers |
| compensation-transformation | Spillover compensation, spectral unmixing, logicle/arcsinh; spreading error and cofactor choice |
| cytometry-qc | Time-based cleaning (flowAI/flowCut/PeacoQC), margins, drift, dead cells, calibration |
| doublet-detection | FSC-A vs FSC-H singlet discrimination; CyTOF Gaussian/DNA; residual heterotypic doublets |
| gating-analysis | Manual + automated (openCyto/flowDensity) hierarchies; FMO boundaries; rare-event/MRD |
| clustering-phenotyping | FlowSOM/PhenoGraph, type-vs-state markers, metacluster annotation, UMAP for display |
| differential-analysis | diffcyt DA/DS; sample-as-unit; compositionality; mixed models |
| bead-normalization | EQ-bead drift correction (normCytof) and cross-batch harmonization (CytoNorm) |

## Example Prompts

- "Load my FCS files raw, apply the recorded compensation, and logicle-transform the markers"
- "Run QC in the right order - margins before the density-based cleaning"
- "Gate singlets on the FSC-A vs FSC-H diagonal, then live CD3+ T cells, using my FMO controls for the boundaries"
- "Cluster my CyTOF data with FlowSOM on lineage markers and annotate the metaclusters"
- "Test which populations differ between responders and non-responders, treating donors as the unit"
- "Normalize my multi-batch CyTOF study with EQ beads then CytoNorm using my reference sample"

## Requirements

```r
# R/Bioconductor
BiocManager::install(c("flowCore", "flowWorkspace", "CytoML", "openCyto", "flowDensity",
                       "CATALYST", "FlowSOM", "diffcyt", "flowAI", "PeacoQC", "ggcyto"))
# CytoNorm (cross-batch), Rphenograph (clustering) - GitHub
# remotes::install_github(c("saeyslab/CytoNorm", "JinmiaoChenLab/Rphenograph"))
```

```bash
# Python (optional - FCS reading and the scanpy/AnnData bridge)
pip install flowkit readfcs
```

## Related Skills

- **single-cell** - Shared clustering, doublet, batch-integration, and annotation concepts (FlowSOM <-> Leiden)
- **imaging-mass-cytometry** - Same CATALYST/metal-channel/arcsinh conventions for tissue imaging
- **spatial-transcriptomics** - Spatial proteomics (CODEX/IMC) adjacency to suspension cytometry
- **differential-expression** - diffcyt reuses the edgeR/limma engines and output semantics
- **experimental-design** - Batch layout, anchor-sample design, and FDR for high-dimensional discovery
