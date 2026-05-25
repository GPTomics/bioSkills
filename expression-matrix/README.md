# expression-matrix

## Overview

Postdoc-grade ingestion, ID mapping, metadata joining, normalization, and sparse handling for RNA-seq count matrices (bulk and single-cell). Encodes the silent-miscounting traps (featureCounts -p v2.0.2 API break, STAR strandedness column choice, salmon NumReads sum, lengthScaledTPM-is-not-TPM), the HGNC rename and PAR_Y suffix gotchas, the alphabetical reference level / interaction resultsNames pitfalls, the TMM/RLE catastrophic failure modes (MYC, apoptosis, viral host shutoff, prokaryotic stress), and the R/Python sparse interop transpose.

**Tool type:** mixed | **Primary tools:** tximport, pandas, DESeq2, edgeR, anndata, scanpy, biomaRt, AnnotationDbi, scipy.sparse

## Skills

| Skill | Description |
|-------|-------------|
| counts-ingest | featureCounts/HTSeq/STAR/salmon/kallisto/RSEM/10X/H5AD import; tximport countsFromAbundance decision tree; featureCounts -p API break; STAR column choice; tximeta provenance |
| gene-id-mapping | Ensembl/Entrez/HGNC/UniProt/MANE mapping; Ensembl version stripping with _PAR_Y preservation; SEPT/MARCH/MARC HGNC renames; biomaRt version pinning; cross-species orthology |
| metadata-joins | Sample alignment; alphabetical reference level trap; LRT reduced-model rules; interaction resultsNames trap; repeated measures (dream); sample swap detection (XIST/Y, somalier); SABV |
| normalization | TMM/TMMwsp/RLE/upper-quartile; VST/rlog visualization-only; log-CPM; spike-in SBN; scran for single-cell; GC correction (cqn, EDASeq); failure modes under MYC/viral/prokaryotic stress |
| sparse-handling | dgCMatrix/dgRMatrix/dgTMatrix; R/Python interop transpose; HDF5 vs Zarr; HDF5SE + DelayedArray; scanpy backed mode; 10X formats; dense conversion blow-up; ~10-15% density crossover; Dask + Zarr |

## Example Prompts

- "I have featureCounts output for 12 paired-end RNA-seq samples -- get me a clean count matrix"
- "I have STAR-aligned TruSeq Stranded BAMs -- get gene counts"
- "Import my Salmon quantifications for DESeq2"
- "I'm running limma-voom on Salmon output -- import correctly"
- "Convert Ensembl IDs to symbols, including the recently renamed septin and march genes"
- "I need a stable Ensembl version pin for a publication -- make this reproducible"
- "Map my human DE gene list to mouse orthologs for a follow-up experiment"
- "Build the transcript-to-gene table I need for tximport"
- "Align sample metadata with my count matrix and check for swapped samples"
- "Verify the reported sex of each sample matches its expression"
- "Set up a paired-design DE for tumor vs normal from the same patients"
- "Choose the right variance-stabilizing transformation for my 20-sample study with heterogeneous library sizes"
- "My MYC-amplified samples have suspiciously muted fold changes -- diagnose the normalization"
- "Normalize my single-cell matrix properly without breaking the zero-heavy assumption"
- "Convert a 100k-cell AnnData to a SingleCellExperiment without orientation mistakes"
- "I have a 35 GB single-cell h5ad and 16 GB RAM -- work without loading it all"
- "Process the combined recount3 + TCGA matrix without running out of memory"

## Requirements

```bash
pip install pandas numpy scipy anndata scanpy mygene pyensembl pydeseq2 zarr dask

pyensembl install --release 110 --species human
```

```r
BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'tximport', 'tximeta',
                        'biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db',
                        'GenomicFeatures', 'rtracklayer',
                        'HDF5Array', 'DelayedArray', 'DelayedMatrixStats',
                        'SummarizedExperiment', 'SingleCellExperiment',
                        'scran', 'scater', 'EDASeq', 'cqn',
                        'variancePartition', 'zellkonverter'))
install.packages(c('Matrix', 'pheatmap', 'dplyr'))
```

```bash
brew install somalier
pip install ngscheckmate
```

## Related Skills

- **differential-expression** - DESeq2/edgeR/limma-voom DE; tximport feeds DESeq2 directly
- **rna-quantification** - featureCounts, salmon, kallisto, RSEM upstream of ingest
- **single-cell** - Single-cell preprocessing, normalization, integration
- **pathway-analysis** - Entrez IDs for KEGG/Reactome; symbols for MSigDB
- **database-access** - biomaRt, Ensembl REST, UniProt detailed patterns
- **read-qc** - RIN/DV200 covariates; degradation effects
- **read-alignment** - STAR/BWA upstream; strandedness verification
- **alignment-files** - BAM/CRAM upstream of featureCounts
- **workflows** - End-to-end pipelines wrapping these skills
