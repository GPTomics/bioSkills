# ChIP-seq Analysis

ChIP-seq analysis using MACS3 for peak calling, ChIPseeker for annotation, and DiffBind for differential binding analysis.

## Overview

This category covers ChIP-seq data analysis: MACS3 identifies enriched regions (peaks) from aligned reads, ChIPseeker annotates peaks to genomic features and genes, and DiffBind finds differentially bound regions between conditions. ChIP-seq is used to study protein-DNA interactions including transcription factor binding and histone modifications.

**Tool type:** `mixed`
**Primary tools:** MACS3 (CLI), ChIPseeker (R), DiffBind (R)

## Skills

| Skill | Description |
|-------|-------------|
| [peak-calling](peak-calling/) | Call peaks with MACS3 from ChIP-seq BAM files |
| [peak-annotation](peak-annotation/) | Annotate peaks to genes with ChIPseeker |
| [differential-binding](differential-binding/) | Differential binding analysis with DiffBind |
| [chipseq-visualization](chipseq-visualization/) | Visualize ChIP-seq data and peaks |

## Workflow

```
Raw ChIP-seq Reads (FASTQ)
    |
    v
[Alignment: bwa/bowtie2] ----> Align to reference genome
    |
    v
[alignment-files] -----------> Sort, index, deduplicate
    |
    v
[peak-calling] --------------> Call peaks with MACS3
    |
    v
[peak-annotation] -----------> Annotate peaks to genes
    |
    v
[differential-binding] ------> Find differential peaks
    |
    v
[pathway-analysis] ----------> Functional enrichment
```

## ChIP-seq Data Types

| Target | Type | Peak Shape | MACS3 Mode |
|--------|------|------------|------------|
| Transcription factors | Point source | Sharp, narrow | Default |
| H3K4me3, H3K27ac | Active marks | Sharp peaks | Default |
| H3K4me1 | Enhancers | Broad | --broad |
| H3K36me3 | Gene bodies | Broad | --broad |
| H3K9me3, H3K27me3 | Repressive | Very broad | --broad |

## Example Prompts

### Peak Calling
- "Call peaks from my ChIP-seq BAM with MACS3"
- "Run MACS3 with input control"
- "Call broad peaks for H3K27me3"

### Peak Annotation
- "Annotate my peaks with nearby genes"
- "Find peaks in promoters vs enhancers"
- "What genes have peaks in their promoter?"

### Differential Binding
- "Find differential peaks between treatment and control"
- "Compare binding between two conditions with DiffBind"
- "Identify lost and gained peaks"

### Visualization
- "Create a heatmap of peak signal"
- "Plot peak distribution around TSS"
- "Visualize peaks in a genomic region"

## Requirements

### MACS3 (CLI)

```bash
# Conda installation (recommended)
conda install -c bioconda macs3

# Or pip
pip install macs3

# Verify
macs3 --version
```

### R Packages

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('ChIPseeker', 'DiffBind', 'TxDb.Hsapiens.UCSC.hg38.knownGene'))

# For visualization
BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
```

## Key Functions

| Function | Tool | Purpose |
|----------|------|---------|
| macs3 callpeak | MACS3 | Call peaks |
| macs3 bdgcmp | MACS3 | Compare signal tracks |
| annotatePeak | ChIPseeker | Annotate peaks to genes |
| plotAnnoPie | ChIPseeker | Pie chart of annotations |
| plotDistToTSS | ChIPseeker | Distance to TSS plot |
| covplot | ChIPseeker | Coverage plot |
| dba | DiffBind | Create DBA object |
| dba.count | DiffBind | Count reads in peaks |
| dba.contrast | DiffBind | Set up comparison |
| dba.analyze | DiffBind | Differential analysis |

## Notes

- **Input control required** - Always use input/IgG control for accurate peak calling
- **Duplicate removal** - Mark duplicates before peak calling (but don't remove for MACS3)
- **Fragment size** - MACS3 models fragment size, but you can specify with --extsize
- **Broad peaks** - Use --broad for histone marks like H3K27me3
- **Multiple samples** - Use DiffBind for comparing conditions with replicates

## Related Skills

- **alignment-files** - BAM preparation before peak calling
- **pathway-analysis** - Functional enrichment of peak-associated genes
- **differential-expression** - Compare ChIP-seq with RNA-seq
- **methylation-analysis** - Combine with DNA methylation data

## References

- [MACS3 GitHub](https://github.com/macs3-project/MACS)
- [ChIPseeker Bioconductor](https://bioconductor.org/packages/ChIPseeker/)
- [DiffBind Bioconductor](https://bioconductor.org/packages/DiffBind/)
- [ENCODE ChIP-seq Pipeline](https://www.encodeproject.org/chip-seq/)
