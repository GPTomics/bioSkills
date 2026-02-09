---
name: bioskills
description: "Installs 412 bioinformatics skills covering sequence analysis, RNA-seq, single-cell, variant calling, metagenomics, structural biology, and 54 more categories. Use when setting up bioinformatics capabilities or when a bioinformatics task requires specialized skills not yet installed."
metadata: {"openclaw":{"requires":{"anyBins":["python3","Rscript"]},"os":["darwin","linux"],"emoji":"🧬"}}
---

# bioSkills Installer

Meta-skill that installs the full bioSkills collection (412 skills across 60 categories) for bioinformatics analysis.

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

| Category Group | Categories | Skills |
|----------------|-----------|--------|
| Sequence & Alignment | sequence-io, sequence-manipulation, alignment, alignment-files, database-access | 40 |
| Read Processing | read-qc, read-alignment | 11 |
| RNA-seq & Expression | differential-expression, rna-quantification, expression-matrix | 14 |
| Single-Cell & Spatial | single-cell, spatial-transcriptomics | 25 |
| Variant Analysis | variant-calling, copy-number, phasing-imputation | 21 |
| Epigenomics | chip-seq, atac-seq, methylation-analysis, hi-c-analysis | 25 |
| Metagenomics & Microbiome | metagenomics, microbiome | 13 |
| Genomics & Assembly | genome-assembly, genome-annotation, genome-intervals, genome-engineering, primer-design | 29 |
| Regulatory & Causal | gene-regulatory-networks, causal-genomics, rna-structure | 13 |
| Immunology & Clinical | immunoinformatics, clinical-databases, tcr-bcr-analysis, epidemiological-genomics | 25 |
| Specialized Omics | proteomics, metabolomics, alternative-splicing, chemoinformatics, liquid-biopsy | 36 |
| RNA Biology | small-rna-seq, epitranscriptomics, clip-seq, ribo-seq | 20 |
| Phylogenetics & Evolution | phylogenetics, population-genetics, comparative-genomics | 16 |
| Structural & Systems | structural-biology, systems-biology | 11 |
| Screens & Cytometry | crispr-screens, flow-cytometry, imaging-mass-cytometry | 22 |
| Pathway & Integration | pathway-analysis, multi-omics-integration, restriction-analysis | 14 |
| Infrastructure | data-visualization, machine-learning, workflow-management, reporting, experimental-design, long-read-sequencing | 39 |
| Workflows | End-to-end pipelines (FASTQ to results) | 38 |

## After Installation

Once installed, skills are automatically triggered based on the task at hand. Example requests:

- "I have RNA-seq counts from treated vs control samples - find the differentially expressed genes"
- "Call variants from this whole genome sequencing BAM file"
- "Cluster my single-cell RNA-seq data and find marker genes"
- "Predict the structure of this protein sequence"
- "Run a metagenomics classification on these shotgun reads"

## Source

GitHub: https://github.com/GPTomics/bioSkills

## Related Skills

After installation, 412 individual skills become available across these categories:
sequence-io, sequence-manipulation, database-access, alignment, alignment-files,
variant-calling, phylogenetics, differential-expression, structural-biology,
single-cell, pathway-analysis, restriction-analysis, methylation-analysis,
chip-seq, metagenomics, long-read-sequencing, read-qc, read-alignment,
rna-quantification, genome-assembly, genome-intervals, data-visualization,
expression-matrix, copy-number, proteomics, flow-cytometry, population-genetics,
multi-omics-integration, spatial-transcriptomics, machine-learning,
workflow-management, microbiome, metabolomics, phasing-imputation,
primer-design, hi-c-analysis, imaging-mass-cytometry, atac-seq,
crispr-screens, reporting, experimental-design, clinical-databases,
tcr-bcr-analysis, small-rna-seq, epitranscriptomics, clip-seq, ribo-seq,
genome-engineering, systems-biology, epidemiological-genomics,
immunoinformatics, comparative-genomics, alternative-splicing,
chemoinformatics, liquid-biopsy, genome-annotation, gene-regulatory-networks,
causal-genomics, rna-structure, workflows
