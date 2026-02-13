---
name: bio-epitranscriptomics-m6a-peak-calling
description: Call m6A peaks from MeRIP-seq IP vs input comparisons. Use when identifying m6A modification sites from methylated RNA immunoprecipitation data.
tool_type: mixed
primary_tool: exomePeak2
---

## Version Compatibility

Reference examples tested with: MACS3 3.0+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# m6A Peak Calling

**"Call m6A peaks from my MeRIP-seq data"** → Identify m6A-modified RNA regions by comparing immunoprecipitated (IP) and input samples using statistical enrichment testing.
- R: `exomePeak2::exomePeak2()` for GC-bias aware peak calling
- CLI: `macs3 callpeak` as an alternative broad peak caller

## exomePeak2 (Recommended)

**Goal:** Identify m6A-enriched regions by comparing IP and input samples with GC-bias correction and replicate-aware statistical testing.

**Approach:** Provide IP and input BAM files along with a gene annotation to exomePeak2, which models read counts in sliding windows across the transcriptome and calls significant enrichment peaks.

```r
library(exomePeak2)

# Peak calling with biological replicates
result <- exomePeak2(
    bam_ip = c('IP_rep1.bam', 'IP_rep2.bam'),
    bam_input = c('Input_rep1.bam', 'Input_rep2.bam'),
    gff = 'genes.gtf',
    genome = 'hg38',
    paired_end = TRUE
)

# Export peaks
exportResults(result, format = 'BED')
```

## MACS3 Alternative

```bash
# Call peaks treating input as control
macs3 callpeak \
    -t IP_rep1.bam IP_rep2.bam \
    -c Input_rep1.bam Input_rep2.bam \
    -f BAMPE \
    -g hs \
    -n m6a_peaks \
    --nomodel \
    --extsize 150 \
    -q 0.05
```

## MeTPeak

```r
library(MeTPeak)

# GTF-aware peak calling
metpeak(
    IP_BAM = c('IP_rep1.bam', 'IP_rep2.bam'),
    INPUT_BAM = c('Input_rep1.bam', 'Input_rep2.bam'),
    GENE_ANNO_GTF = 'genes.gtf',
    OUTPUT_DIR = 'metpeak_output'
)
```

## Peak Filtering

```bash
# Filter by fold enrichment and q-value
# FC > 2, q < 0.05 typical thresholds
awk '$7 > 2 && $9 < 0.05' peaks.xls > filtered_peaks.bed
```

## Related Skills

- merip-preprocessing - Prepare data for peak calling
- m6a-differential - Compare peaks between conditions
- chip-seq/peak-calling - Similar concepts
