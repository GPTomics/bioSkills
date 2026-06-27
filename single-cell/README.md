# single-cell

## Overview

Single-cell genomics with Seurat (R) and Scanpy (Python). Covers the core scRNA-seq workflow (loading, QC, ambient/doublet removal, normalization, batch integration, clustering, marker detection, and cell-type annotation) plus differential abundance, trajectory inference and RNA velocity, cell-cell and metabolite communication, multimodal CITE-seq/Multiome integration, scATAC-seq, Perturb-seq, and lineage tracing.

**Tool type:** mixed | **Primary tools:** Seurat, Scanpy, Pertpy, Cassiopeia, MeboCost

## Skills

| Skill | Description |
|-------|-------------|
| data-io | Load 10X (raw vs filtered), read/write h5ad/RDS/h5mu/zarr, and convert across AnnData/Seurat/SCE without losing layers or transposing the matrix |
| preprocessing | MAD-adaptive QC, ambient-RNA removal, normalization choice (shifted-log/scran/sctransform/Pearson), and HVG selection with tissue-aware guidance |
| doublet-detection | Per-sample doublet calling before integration with scDblFinder/Scrublet/DoubletFinder, rate set from recovered cells, and homotypic/over-removal cautions |
| hashing-demultiplexing | Assign pooled hashed cells to their sample of origin and call cross-sample doublets from HTO counts (HTODemux, hashsolo, demuxEM, GMM-Demux) |
| batch-integration | Integrate batches with Harmony, scVI/scANVI, Seurat, fastMNN, Scanorama, or BBKNN; method choice by design, over-correction diagnosis, and scIB scoring without metric gaming |
| clustering | Graph-based clustering and dimensionality reduction; Leiden vs Louvain, sweeping and validating resolution, and why post-clustering marker p-values are not inference |
| markers-annotation | Cluster marker detection and manual cell-type labeling; marker ranking vs pseudobulk condition DE, with the double-dipping caveat |
| cell-annotation | Automated reference-based label transfer with CellTypist, SingleR, Azimuth, scANVI, scmap; rejection thresholds and novel-vs-artifact triage |
| cnv-inference | Infer large-scale CNVs from tumor scRNA-seq to separate malignant from normal cells and call subclones (inferCNV, copyKAT, Numbat, SCEVAN) |
| differential-abundance | Test cell-type proportion/composition changes between conditions with Milo, scCODA, sccomp, propeller; compositional-simplex and replicate-unit caveats |
| trajectory-inference | Pseudotime, fate probabilities, and RNA velocity with PAGA, Slingshot, Monocle3, Palantir, scVelo, CellRank 2; topology-first method choice and velocity validation |
| cell-communication | Consensus-first ligand-receptor inference (LIANA) with CellPhoneDB, CellChat, NicheNet; every call is a co-expression proxy to validate orthogonally |
| multimodal-integration | Classify multimodal tasks by anchor structure and pick a joint method (WNN, totalVI, MultiVI, MOFA+, GLUE, Seurat v5 bridge), denoising CITE-seq ADT background first |
| scatac-analysis | Process scATAC fragments with Signac/ArchR/SnapATAC2: QC, TF-IDF/LSI with depth-component diagnosis, consensus peaks, chromVAR motifs, doublets, without binarizing |
| perturb-seq | Single-cell CRISPR screens: mixture guide assignment, Mixscape escaper removal, calibrated SCEPTRE testing, E-distance effect size, composition-vs-expression separation |
| lineage-tracing | Lineage trees and clonal dynamics from CRISPR scars, barcodes, or mtDNA with Cassiopeia, Startle, CoSpar; solver choice under homoplasy/dropout, state-vs-fate |
| metabolite-communication | Metabolite-mediated communication by enzyme-to-sensor scoring (MEBOCOST), with scFEA/Compass/NeuronChat; a doubly-inferred, validation-required layer |

## Example Prompts

- "Load my 10X data into Scanpy"
- "Create a Seurat object from this count matrix"
- "Convert h5ad to Seurat format"
- "Run QC and filter cells with >20% mitochondrial"
- "Detect doublets with Scrublet"
- "Remove doublets from my Seurat object"
- "Assign my hashed cells back to their sample of origin and flag cross-sample doublets"
- "Run HTODemux on my Seurat HTO assay"
- "Normalize and find highly variable genes"
- "Preprocess this data using SCTransform"
- "Run PCA and cluster at resolution 0.5"
- "Generate a UMAP colored by cluster"
- "Try different clustering resolutions"
- "Integrate my samples with Harmony to remove batch effects"
- "Find marker genes for each cluster"
- "What cell types are these clusters?"
- "Annotate my cells with CellTypist or a reference atlas"
- "Show a dot plot of canonical PBMC markers"
- "Score cell-cycle phase and check it isn't driving my clusters"
- "Did cell-type proportions change between disease and control?"
- "Test for differential abundance with Milo"
- "Order my cells along a pseudotime trajectory"
- "Run RNA velocity to find the direction of differentiation"
- "Load my CITE-seq data with both RNA and ADT"
- "Run WNN clustering combining RNA and protein"
- "Analyze my 10X Multiome data"
- "Process my scATAC-seq data with Signac"
- "Run chromVAR motif analysis on scATAC"
- "Find differentially accessible peaks"
- "Which cells in my tumor scRNA-seq are malignant?"
- "Infer copy-number and call tumor subclones with Numbat"
- "Analyze my Perturb-seq CRISPR screen"
- "Find genes affected by each perturbation"
- "Build lineage tree from CRISPR barcodes"
- "Track clonal dynamics with CoSpar"
- "Find ligand-receptor communication between my cell types"
- "Which cell types signal to my macrophages, and through what ligands?"
- "Analyze metabolite-receptor signaling between cell types"
- "Find metabolite-mediated communication in my tumor microenvironment"

## Requirements

```bash
# Python (Scanpy + Scrublet)
pip install scanpy anndata leidenalg matplotlib scrublet

# R (Seurat + DoubletFinder)
install.packages('Seurat')
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
BiocManager::install('scDblFinder')

# R (Signac for scATAC)
install.packages('Signac')
BiocManager::install(c('chromVAR', 'motifmatchr', 'JASPAR2020'))

# R (differential abundance)
BiocManager::install('miloR')
install.packages('sccomp')

# R (hashtag demultiplexing + CNV inference from tumor scRNA-seq)
BiocManager::install(c('demuxmix', 'infercnv'))
install.packages('numbat')
remotes::install_github(c('navinlabcode/copykat', 'AntonioDeFalco/SCEVAN'))

# Python (hashtag demultiplexing)
pip install solo-sc pegasuspy demuxEM

# Python (differential abundance, perturbation, lineage, metabolite communication)
pip install sccoda pertpy cassiopeia-lineage cospar mebocost
```

## Related Skills

- **differential-expression** - Bulk RNA-seq DE analysis
- **pathway-analysis** - GO/KEGG enrichment of marker genes
- **rna-quantification** - Count matrix generation
