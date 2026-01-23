---
name: bio-workflows-merip-pipeline
description: End-to-end MeRIP-seq analysis from FASTQ to m6A peaks and differential methylation. Use when analyzing epitranscriptomic m6A modifications from immunoprecipitation data.
tool_type: mixed
primary_tool: exomePeak2
---

# MeRIP-seq Pipeline

## Pipeline Overview

```
FASTQ → Align IP+Input → Peak calling → Differential → Visualization
```

## Step 1: Alignment

```bash
# Align IP and Input samples
for sample in IP_rep1 IP_rep2 Input_rep1 Input_rep2; do
    STAR --genomeDir index \
        --readFilesIn ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}_
done
```

## Step 2: Peak Calling

```r
library(exomePeak2)

result <- exomePeak2(
    bam_ip = c('IP_rep1.bam', 'IP_rep2.bam'),
    bam_input = c('Input_rep1.bam', 'Input_rep2.bam'),
    gff = 'genes.gtf',
    genome = 'hg38'
)

exportResults(result, format = 'BED')
```

## Step 3: Differential Methylation

```r
design <- data.frame(condition = c('ctrl', 'ctrl', 'treat', 'treat'))

diff_result <- exomePeak2(
    bam_ip = all_ip_bams,
    bam_input = all_input_bams,
    gff = 'genes.gtf',
    experiment_design = design
)
```

## Step 4: Visualization

```r
library(Guitar)
GuitarPlot(peaks, txdb, saveToPDFprefix = 'm6a_metagene')
```

## QC Checkpoints

1. **After alignment**: Check IP/Input correlation
2. **After peaks**: Verify DRACH motif enrichment
3. **After metagene**: Should see stop codon enrichment

## Related Skills

- **epitranscriptomics/** - Individual m6A analysis skills
- **chip-seq/peak-calling** - Similar concepts
