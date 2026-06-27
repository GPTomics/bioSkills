# small-rna-seq

## Overview

Analyze small RNA sequencing data across all classes - miRNA, isomiRs, tRNA-derived fragments (tRFs), piRNA, and rRNA/snoRNA-derived species - for preprocessing, discovery, quantification, differential expression, and target prediction. The category leads with the decisions that make or break a small-RNA analysis: ligation bias makes absolute cross-miRNA abundance untrustworthy, the library kit determines which classes are even captured, a miRDeep2 score is a hypothesis not a finding, miRNA normalization is compositionally fragile, and a predicted target is a hypothesis that needs expression evidence.

**Tool type:** mixed | **Primary tools:** cutadapt, miRDeep2, miRge3, DESeq2, miRanda, MINTmap

## Skills

| Skill | Description |
|-------|-------------|
| smrna-preprocessing | Kit-specific adapter, UMI, and 4N handling; size selection; length-histogram QC |
| mirdeep2-analysis | De novo miRNA discovery scored against the Dicer/Drosha biogenesis signature |
| mirge3-analysis | Fast known-miRNA, isomiR, tRF, and A-to-I quantification against curated libraries |
| differential-mirna | DE with compositionally-aware normalization and low-count handling |
| target-prediction | Seed-based prediction and validated databases, filtered by expression evidence |
| trf-pirna-profiling | tRF, piRNA, and small-RNA-class profiling with MINTmap, unitas, and proTRAC |

## Example Prompts

- "Trim the kit adapter and size-select my small RNA-seq, then read the length histogram"
- "Discover novel miRNAs with miRDeep2 and pick a score cutoff from the signal-to-noise table"
- "Quantify known miRNAs and isomiRs with miRge3 using MirGeneDB"
- "Find DE miRNAs, checking whether a few dominant miRNAs distort normalization"
- "Predict targets of my DE miRNAs and keep only those anti-correlated with mRNA-seq"
- "Profile tRNA-derived fragments and test the piRNA ping-pong signature"
- "My TruSeq library shows no tRNA halves - is that real or an assay artifact?"

## Requirements

```bash
# Preprocessing and QC
pip install cutadapt umi_tools
conda install -c bioconda cutadapt fastp seqkit

# Quantification and discovery
conda install -c bioconda mirdeep2 bowtie viennarna
pip install mirge3

# tRF / piRNA profiling
conda install -c bioconda mintmap unitas protrac

# Differential expression (R)
# BiocManager::install(c('DESeq2', 'edgeR', 'apeglm', 'EnhancedVolcano', 'pheatmap'))

# Target prediction
conda install -c bioconda miranda
pip install biopython gseapy
```

## Related Skills

- **read-qc** - General read quality control, adapter trimming, and UMI processing
- **differential-expression** - Bulk RNA-seq DE framework (DESeq2/edgeR mechanics)
- **rna-quantification** - RNA quantification concepts
- **clip-seq** - AGO-CLIP/CLASH direct miRNA-target evidence
- **genome-annotation** - tRNA/rRNA/snoRNA locus annotation underlying tRF/piRNA tools
