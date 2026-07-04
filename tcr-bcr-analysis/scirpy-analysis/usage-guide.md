# scirpy Analysis - Usage Guide

## Overview

scirpy is the scanpy-native toolkit for single-cell paired TCR/BCR: it ingests AIRR receptor records, QCs chain pairing, defines clonotypes, and computes repertoire statistics with everything living in the same AnnData/MuData object as the gene expression, so clonality can be overlaid on the transcriptomic UMAP. Two choices dominate downstream validity and must be stated with every result. First, the clonotype definition is receptor-specific: TCR clones share an exact CDR3 nucleotide sequence (identity clonotyping is correct), but B cells somatically hypermutate, so identity clonotyping shatters one BCR lineage into fake singletons and BCR requires nucleotide distance clustering (normalized Hamming, same V and J gene). Second, chain QC filtering is a sampling decision: multichain and TCR+BCR-ambiguous cells are doublets and are dropped, but blanket removal of orphan and extra-chain cells preferentially deletes small clones and biases clonal-expansion and diversity estimates upward. Since scirpy 0.13 the receptor data is an awkward array in obsm['airr'] read via get.airr after pp.index_chains, not legacy per-chain obs columns.

## Prerequisites

```bash
pip install 'scirpy>=0.24' scanpy mudata
```

Upstream contigs come from CellRanger vdj, dandelion, or airrflow. For rigorous BCR, reannotate contigs with IgBLAST (via dandelion or airrflow) before scirpy, since CellRanger BCR contigs are not IMGT-numbered.

## Quick Start

Tell your AI agent what you want to do:
- "Load my 10x VDJ contigs and pair them with my scRNA-seq object"
- "QC chain pairing and remove doublets"
- "Define TCR clonotypes by exact CDR3"
- "Define BCR clonotype clusters accounting for somatic hypermutation"
- "Show clonal expansion overlaid on my UMAP"
- "Compare repertoire diversity and overlap between samples"

## Example Prompts

### Loading and QC

> "Read filtered_contig_annotations.csv, build a MuData with my GEX AnnData, and index the chains"

> "Run chain QC and tell me how many cells are single pair, orphan, extra, or multichain"

> "Filter out doublets but keep orphan cells for the cell-state analysis"

### Clonotype definition

> "Define TCR clonotypes by exact CDR3-nucleotide identity with both arms matching"

> "Define BCR clonotype clusters with normalized Hamming distance on nucleotides, requiring the same V and J gene"

> "Cluster my TCRs by tcrdist similarity to catch antigen-convergent receptors"

### Expansion, diversity, integration

> "Bin cells by clonal expansion and plot the fraction expanded per cell subtype"

> "Compute normalized Shannon diversity per sample and remind me it depends on the clonotype definition"

> "Color my UMAP by clonal expansion and test whether the big clones occupy one transcriptional state"

> "Annotate antigen specificity by matching my TCRs against VDJdb"

## What the Agent Will Do

1. Ingest receptor contigs (read_10x_vdj / read_airr / from_dandelion) and wrap them with the GEX AnnData in a MuData, then run pp.index_chains.
2. Run tl.chain_qc, drop multichain and ambiguous doublets, and decide orphan/extra retention by the downstream question.
3. Cache distances with pp.ir_dist and define clonotypes: define_clonotypes (identity) for TCR, define_clonotype_clusters (normalized_hamming, nt, same V/J) for BCR.
4. Compute clonal expansion, alpha diversity, and repertoire overlap, reporting them alongside the clonotype definition and filter used.
5. Push clonality into the GEX modality and overlay it on the transcriptomic UMAP; optionally test clonotype modularity.
6. Annotate specificity against a reference DB (ir_query) or export AIRR for interchange, handing BCR lineage/SHM work to Immcantation/dandelion.

## Tips

- Access receptor fields with get.airr and get.airr_context, never obs indexing; obsm['airr'] is an awkward array. Migrate pre-0.13 objects with io.upgrade_schema.
- Keep pp.ir_dist metric and sequence identical to the subsequent define_* call, or the cached distances silently mismatch the grouping.
- Never use identity clonotyping for BCR; somatic hypermutation makes lineage members non-identical.
- Report diversity and expansion with the clonotype definition and QC filter that produced them; every knob (receptor_arms, dual_ir, same_v_gene, orphan policy) moves the numbers.
- Extra-VJ cells are often real dual-TCR (allelic inclusion, ~30% of T cells), not junk; set dual_ir deliberately rather than discarding them.
- Filter contigs by duplicate_count/consensus_count in samples dominated by a hyperexpanded clone, where ambient VDJ mRNA manufactures spurious shared chains.
- Cluster on gene expression, then overlay clonality; never cluster cells on receptor sequence and treat it as a cell state.

## Related Skills

- mixcr-analysis - Process raw single-cell VDJ FASTQ
- immcantation-analysis - Proper BCR clonal lineages and SHM downstream
- specificity-annotation - Antigen-specificity clustering on single-cell clonotypes
- single-cell/data-io - Load and manage the GEX AnnData/MuData
- single-cell/clustering - Cell-state clustering to overlay clonality
- single-cell/doublet-detection - Corroborate multichain doublet calls
