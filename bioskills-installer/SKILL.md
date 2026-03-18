---
name: bioskills
description: "Installs 425 bioinformatics skills covering sequence analysis, RNA-seq, single-cell, variant calling, metagenomics, structural biology, and 56 more categories. Use when setting up bioinformatics capabilities or when a bioinformatics task requires specialized skills not yet installed."
metadata: {"openclaw":{"requires":{"bins":["git"],"anyBins":["python3","Rscript"]},"os":["darwin","linux"],"emoji":"🧬"}}
---

# bioSkills Installer

Meta-skill that installs the full bioSkills collection (425 skills across 62 categories) for bioinformatics analysis.

## Installation

Run the bundled install script to download and install all bioSkills:

```bash
bash scripts/install-bioskills.sh
```

Or install only specific categories:

```bash
bash scripts/install-bioskills.sh --categories "single-cell,variant-calling,differential-expression"
```

## What Gets Installed

425 skills across 62 categories covering:

- **Sequence & Alignment** (40): sequence-io, sequence-manipulation, alignment, alignment-files, database-access
- **Read Processing** (11): read-qc, read-alignment
- **RNA-seq & Expression** (14): differential-expression, rna-quantification, expression-matrix
- **Single-Cell & Spatial** (25): single-cell, spatial-transcriptomics
- **Variant Analysis** (21): variant-calling, copy-number, phasing-imputation
- **Epigenomics** (25): chip-seq, atac-seq, methylation-analysis, hi-c-analysis
- **Metagenomics & Microbiome** (13): metagenomics, microbiome
- **Genomics & Assembly** (29): genome-assembly, genome-annotation, genome-intervals, genome-engineering, primer-design
- **Regulatory & Causal** (13): gene-regulatory-networks, causal-genomics, rna-structure
- **Temporal & Ecological** (11): temporal-genomics, ecological-genomics
- **Immunology & Clinical** (25): immunoinformatics, clinical-databases, tcr-bcr-analysis, epidemiological-genomics
- **Specialized Omics** (36): proteomics, metabolomics, alternative-splicing, chemoinformatics, liquid-biopsy
- **RNA Biology** (20): small-rna-seq, epitranscriptomics, clip-seq, ribo-seq
- **Phylogenetics & Evolution** (16): phylogenetics, population-genetics, comparative-genomics
- **Structural & Systems** (11): structural-biology, systems-biology
- **Screens & Cytometry** (22): crispr-screens, flow-cytometry, imaging-mass-cytometry
- **Pathway & Integration** (14): pathway-analysis, multi-omics-integration, restriction-analysis
- **Infrastructure** (39): data-visualization, machine-learning, workflow-management, reporting, experimental-design, long-read-sequencing
- **Workflows** (40): end-to-end pipelines (FASTQ to results)

## After Installation

Once installed, skills are automatically triggered based on the task at hand. Example requests:

- "I have RNA-seq counts from treated vs control samples - find the differentially expressed genes"
- "Call variants from this whole genome sequencing BAM file"
- "Cluster my single-cell RNA-seq data and find marker genes"
- "Predict the structure of this protein sequence"
- "Run a metagenomics classification on these shotgun reads"

## Source

GitHub: https://github.com/GPTomics/bioSkills
