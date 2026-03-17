---
name: bioskills
description: >
  Batch-install 425 bioinformatics agent skills for sequence analysis, RNA-seq
  differential expression, single-cell clustering, variant calling, metagenomics
  classification, and structural biology. Use when setting up bioinformatics
  capabilities, the user says "install bioSkills", or a bioinformatics task
  requires specialized skills not yet available.
metadata: {"openclaw":{"requires":{"bins":["git"],"anyBins":["python3","Rscript"]},"os":["darwin","linux"],"emoji":"🧬"}}
---

# bioSkills Installer

Meta-skill that installs the full bioSkills collection (425 skills across 62 categories) for bioinformatics analysis.

## Installation

```bash
# install all 425 skills
bash scripts/install-bioskills.sh

# install specific categories only
bash scripts/install-bioskills.sh --categories "single-cell,variant-calling,differential-expression"
```

### Verify installation

After running the installer, confirm skills are available:

```bash
ls ~/.claude/skills/ | grep -c "."   # expected: 425 for full install
```

If the install fails partway (network timeout, disk space), re-run the same command — it skips already-installed skills.

## What Gets Installed

425 skills across 62 categories covering:

- **Sequence & Alignment** (40): sequence-io, alignment, alignment-files, database-access
- **RNA-seq & Expression** (14): differential-expression, rna-quantification, expression-matrix
- **Single-Cell & Spatial** (25): single-cell, spatial-transcriptomics
- **Variant Analysis** (21): variant-calling, copy-number, phasing-imputation
- **Epigenomics** (25): chip-seq, atac-seq, methylation-analysis, hi-c-analysis
- **Metagenomics & Microbiome** (13): metagenomics, microbiome
- **Genomics & Assembly** (29): genome-assembly, genome-annotation, genome-intervals
- **Immunology & Clinical** (25): immunoinformatics, clinical-databases, tcr-bcr-analysis
- **Specialized Omics** (36): proteomics, metabolomics, chemoinformatics, liquid-biopsy
- **Infrastructure** (39): data-visualization, machine-learning, workflow-management, reporting
- **Workflows** (40): end-to-end pipelines (FASTQ to results)
- Plus 12 more category groups (RNA biology, phylogenetics, structural biology, screens, etc.)

## After Installation

Skills are automatically triggered based on the task. Example requests:

- "find differentially expressed genes from these RNA-seq counts"
- "call variants from this whole genome sequencing BAM file"
- "cluster my single-cell RNA-seq data and find marker genes"
- "run metagenomics classification on these shotgun reads"

## Source

GitHub: https://github.com/GPTomics/bioSkills
