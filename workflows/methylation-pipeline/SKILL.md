---
name: bio-workflows-methylation-pipeline
description: Orchestrates the end-to-end bisulfite/EM-seq methylation pipeline from FASTQ to differentially methylated regions, chaining Trim Galore/fastp QC, Bismark alignment + deduplication, methylation calling, methylKit coverage-filtering/normalization, and selection-aware DMR detection (dmrseq/DSS). Use when gating the run on bisulfite conversion (lambda + pUC19 controls) BEFORE any beta value, committing the genome build + library directionality once, keeping mate-overlap deduplicated (--no_overlap), M-bias-trimming from the plot, filtering coverage before testing, choosing a count model (beta-binomial/DSS) over a bare-beta t-test, or using a region-selection-aware FDR (dmrseq/DSS) rather than raw methylKit tiles. Hands mechanism to the methylation-analysis component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: Bismark
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - methylation-analysis/bismark-alignment
  - methylation-analysis/methylation-calling
  - methylation-analysis/methylkit-analysis
  - methylation-analysis/differential-cpg-testing
  - methylation-analysis/dmr-detection
qc_checkpoints:
  - after_qc: "Q30 >80%, adapter content removed"
  - after_alignment: "Mapping efficiency >50%, bisulfite conversion >99%"
  - after_calling: "Coverage distribution reasonable, no biased positions"
---

## Version Compatibility

Reference examples tested with: Bismark 0.24+, Bowtie2 2.5.3+, FastQC 0.12+, Trim Galore 0.6.10+, fastp 0.23+, methylKit 1.28+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Methylation Pipeline

**"Analyze my bisulfite sequencing data from FASTQ to DMRs"** -> Chain QC/trim, Bismark alignment + dedup, methylation calling, coverage-filtered per-CpG testing, and selection-aware DMR detection.
- CLI + R: Trim Galore/fastp -> bismark -> deduplicate_bismark -> bismark_methylation_extractor -> methylKit (filter/normalize/unite) -> DSS/dmrseq

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A methylation callset is decided at four seams, not inside the caller.

1. **Bisulfite conversion is verified BEFORE any beta value is trusted — this is the gate that most pipelines skip.** An unmethylated lambda (or spike-in) bounds UNDER-conversion (residual C read as methylated -> false hyper); a methylated pUC19 control bounds OVER-conversion (5mC deaminated -> false hypo). Require conversion >99% from the lambda control before computing a single methylation level; a 98% library silently shifts every call.
2. **Mate overlap must be deduplicated once.** In paired-end WGBS the R1/R2 insert overlaps, and counting a CpG in both mates double-weights it. `bismark_methylation_extractor --paired-end` applies `--no_overlap` by default; running the extractor in single-end mode on paired data (or losing `--no_overlap`) inflates coverage and distorts levels.
3. **M-bias is trimmed from the plot, and coverage is filtered BEFORE testing.** End-repair fill-in biases the first/last few bases (read the M-bias plot, trim positionally — not a fixed number). Then filter low- and extreme-coverage CpGs before any test: variance depends on coverage, so unfiltered low-coverage sites dominate the FDR.
4. **The statistic must respect counts, and region FDR must respect selection.** A bare-beta t-test discards coverage (the precision unique to sequencing); use a beta-binomial/overdispersion model (DSS, methylKit `overdispersion='MN'`) for counts, or limma-on-M-values for arrays. For REGIONS, methylKit fixed tiles do not correct for the region-selection step — use dmrseq (permutation null over selection) or DSS `callDMR` for a rigorous region-level FDR.

## Made-once commitments

| Commitment | Choice | Consequence inherited downstream |
|------------|--------|----------------------------------|
| Genome build + library model | One build; directional (WGBS/EM-seq) vs non-directional/PBAT | Wrong strand model tanks mapping; build fixes all coordinates |
| Conversion controls | Lambda (unmethylated) + pUC19 (methylated) spike-ins | Without them, under/over-conversion is undetectable and biases every call |
| Assay entry | WGBS/EM-seq (this pipeline) vs Infinium array (array-preprocessing) | Array data enters at beta/M matrix, not Bismark |
| Context | CpG (default) vs CHG/CHH (plants/non-CpG) | Non-CpG contexts need conversion-aware calling and separate testing |

## Workflow Overview

```
FASTQ files
    |
    v
[1. QC & Trimming] -----> fastp/Trim Galore
    |
    v
[2. Alignment] ---------> Bismark
    |
    v
[3. Deduplication] -----> deduplicate_bismark
    |
    v
[4. Methylation Calling] -> bismark_methylation_extractor
    |
    v
[5. Per-CpG Analysis] ---> methylKit (R) or scipy (Python)
    |
    v
[6. DMR Detection] ------> methylKit/DSS
    |
    v
Differentially methylated regions
```

## Primary Path: Bismark + methylKit

### Step 1: Quality Control

```bash
# Trim Galore recommended for bisulfite data (handles adapter bias)
trim_galore --paired --fastqc \
    -o trimmed/ \
    sample_R1.fastq.gz sample_R2.fastq.gz

# Or fastp with conservative settings
fastp -i sample_R1.fastq.gz -I sample_R2.fastq.gz \
    -o trimmed/sample_R1.fq.gz -O trimmed/sample_R2.fq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 35 \
    --html qc/sample_fastp.html
```

### Step 2: Bismark Alignment

```bash
# Prepare genome (once)
bismark_genome_preparation --bowtie2 genome/

# Align
bismark --genome genome/ \
    -1 trimmed/sample_R1_val_1.fq.gz \
    -2 trimmed/sample_R2_val_2.fq.gz \
    -o aligned/ \
    --parallel 4 \
    --temp_dir tmp/

# Output: sample_R1_val_1_bismark_bt2_pe.bam
```

**QC Checkpoint:** Check Bismark report
- Mapping efficiency >50% (the 3-letter alphabet lowers uniqueness; 50-70% is normal for WGBS)
- Bisulfite conversion >99% from the unmethylated lambda spike-in (bounds under-conversion -> false hyper); also check a methylated pUC19 control for over-conversion (-> false hypo). With NO spike-in, use the genome-wide non-CpG (CHH) methylation rate as a conversion proxy in mammals (expected near 0)

### Step 3: Deduplication (WGBS / EM-seq ONLY)

Deduplicate WGBS and EM-seq. Do NOT deduplicate RRBS, amplicon, or other target-enrichment libraries: their reads legitimately stack at the MspI cut sites, so positional dedup destroys real coverage (Bismark's own docs say so). Skip this step entirely for RRBS.

```bash
# WGBS / EM-seq only -- skip for RRBS/amplicon
deduplicate_bismark \
    --bam \
    -p \
    --output_dir deduplicated/ \
    aligned/sample_R1_val_1_bismark_bt2_pe.bam
```

### Step 4: Methylation Calling

```bash
# --paired-end enables --no_overlap by DEFAULT (deduplicates the R1/R2 insert overlap so a CpG in
# the overlap is not double-counted). Do NOT run the extractor in single-end mode on paired data.
bismark_methylation_extractor \
    --paired-end \
    --comprehensive \
    --bedGraph \
    --cytosine_report \
    --genome_folder genome/ \
    -o methylation/ \
    deduplicated/sample_R1_val_1_bismark_bt2_pe.deduplicated.bam

# Generate summary report
bismark2report
bismark2summary
```

### Step 5: Analysis with methylKit

**Goal:** Turn per-sample coverage/cytosine reports into a coverage-filtered, normalized, united methylation object ready for testing.

**Approach:** Read each sample with the matching pipeline, drop low-coverage and extreme-coverage CpGs, normalize coverage across libraries, then unite to the sites covered in every sample.

```r
library(methylKit)

# Read methylation calls
files <- list(
    'methylation/control_1.CpG_report.txt',
    'methylation/control_2.CpG_report.txt',
    'methylation/treated_1.CpG_report.txt',
    'methylation/treated_2.CpG_report.txt'
)

sample_ids <- c('control_1', 'control_2', 'treated_1', 'treated_2')
treatment <- c(0, 0, 1, 1)

# Read cytosine reports
meth_obj <- methRead(
    location = as.list(files),
    sample.id = as.list(sample_ids),
    assembly = 'hg38',
    treatment = treatment,
    context = 'CpG',
    pipeline = 'bismarkCytosineReport'
)

# Filter by coverage
meth_filtered <- filterByCoverage(meth_obj, lo.count = 10, hi.perc = 99.9)

# Normalize coverage
meth_norm <- normalizeCoverage(meth_filtered)

# Merge samples (keep sites covered in all)
meth_merged <- unite(meth_norm, destrand = TRUE)

# Sample statistics
getMethylationStats(meth_obj[[1]], plot = TRUE)
getCoverageStats(meth_obj[[1]], plot = TRUE)
```

### Step 5b: Python Alternative for Per-CpG Testing

When methylKit is unavailable or a Python-only workflow is preferred, per-CpG testing can be performed with scipy and statsmodels on beta values computed from the coverage files.

```python
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import numpy as np

# Read Bismark coverage files and compute beta values
# beta = count_methylated / (count_methylated + count_unmethylated)
# Filter CpGs with < 10x coverage in any sample
# Run per-CpG Welch's t-test between groups
# Apply BH FDR correction: multipletests(pvals, method='fdr_bh')
# See methylation-analysis/differential-cpg-testing for full pipeline
```

A bare-beta t-test discards coverage (the precision information unique to sequencing) and is only a quick look. For sequencing counts, route to a beta-binomial / overdispersion-corrected count model (DSS, or methylKit with overdispersion='MN'); for array or continuous data, use limma on M-values. The count-vs-continuous decision is owned by methylation-analysis/differential-cpg-testing.

### Step 6: DMR Detection

methylKit fixed tiles are a fast screen, but their region q-value is not corrected for the region-selection step. For a rigorous region-level FDR use dmrseq (a permutation null over the selection) or DSS callDMR, and confirm with cross-tool overlap - see methylation-analysis/dmr-detection.

```r
# Calculate differential methylation (per CpG). overdispersion='MN' + test='Chisq' applies the
# overdispersion correction seam #4 requires; the default 'none' gives underdispersed p-values.
diff_meth <- calculateDiffMeth(meth_merged, overdispersion = 'MN', test = 'Chisq')

# Get significant DMCs
dmc <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01)

# Tile into regions (DMRs)
tiles <- tileMethylCounts(meth_merged, win.size = 1000, step.size = 1000)
diff_tiles <- calculateDiffMeth(tiles, overdispersion = 'MN', test = 'Chisq')   # same overdispersion correction as per-CpG (seam #4)
dmr <- getMethylDiff(diff_tiles, difference = 25, qvalue = 0.01)

# Export
write.csv(as.data.frame(dmc), 'dmc_results.csv')
write.csv(as.data.frame(dmr), 'dmr_results.csv')

# Annotate with genomic features
library(genomation)
gene_obj <- readTranscriptFeatures('genes.bed')
annotateWithGeneParts(as(dmr, 'GRanges'), gene_obj)
```

## Parameter Recommendations

| Step | Parameter | Value |
|------|-----------|-------|
| Trim Galore | default | Recommended for BS-seq |
| Bismark | --parallel | 4 (per sample parallelization) |
| methylKit | lo.count | 10 (minimum coverage) |
| methylKit | difference | 25 (% methylation difference) |
| methylKit | qvalue | 0.01 |
| DMR tiles | win.size | 500-1000 bp |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Genome-wide hyper- or hypo-methylation shift | Under/over-conversion never checked | Gate on lambda (>99%) + pUC19 controls BEFORE trusting any beta value |
| Coverage inflated, levels off in mate-overlap regions | Extractor run single-end on paired data / lost `--no_overlap` | Use `--paired-end` (applies `--no_overlap`); do not single-end paired data |
| Systematic bias at read ends | M-bias from end-repair fill-in | Trim positionally from the M-bias plot, not a fixed number |
| Low-coverage CpGs dominate the DMC list | No coverage filter before testing | `filterByCoverage(lo.count=10, hi.perc=99.9)` before `calculateDiffMeth` |
| Spurious DMCs / underdispersed p-values | Bare-beta t-test ignores counts/overdispersion | Beta-binomial/DSS or methylKit `overdispersion='MN'`; limma-M for arrays |
| Region q-values too optimistic | methylKit fixed tiles ignore the region-selection step | Use dmrseq (permutation null) or DSS `callDMR` for region-level FDR |
| Very low mapping efficiency | Wrong library directionality (PBAT/non-directional aligned as directional) | Set the correct Bismark strand model (methylation-analysis/bismark-alignment) |

The full per-step chain is shown above; the runnable methylKit analysis is in this skill's examples/ (`methylkit_analysis.R`).

## References

- Krueger F, Andrews SR (2011) Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. *Bioinformatics* 27:1571-1572. DOI 10.1093/bioinformatics/btr167.
- Akalin A, Kormaksson M, Li S, et al (2012) methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. *Genome Biology* 13:R87. DOI 10.1186/gb-2012-13-10-r87.
- Feng H, Conneely KN, Wu H (2014) A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. *Nucleic Acids Research* 42:e69. DOI 10.1093/nar/gku154. (DSS.)
- Korthauer K, Chakraborty S, Benjamini Y, Irizarry RA (2019) Detection and accurate false discovery rate control of differentially methylated regions from whole genome bisulfite sequencing. *Biostatistics* 20:367-383. DOI 10.1093/biostatistics/kxy007. (dmrseq; region-selection-aware FDR.)

## Related Skills

- methylation-analysis/bismark-alignment - Bisulfite/EM-seq alignment, library/strand model, conversion QC
- methylation-analysis/methylation-calling - Calling from BAM (Bismark/MethylDackel), contexts, variant-aware
- methylation-analysis/methylkit-analysis - methylKit object model and overdispersion gotchas
- methylation-analysis/differential-cpg-testing - Per-CpG testing (count-vs-continuous fork)
- methylation-analysis/dmr-detection - Selection-aware region callers (dmrseq/DSS) and PMD segmentation
- methylation-analysis/array-preprocessing - Alternate entry: Infinium IDAT to beta/M matrix
- methylation-analysis/cell-type-deconvolution - Cell-fraction covariates for bulk-tissue EWAS
- methylation-analysis/epigenetic-clocks - DNAm age and age acceleration
- methylation-analysis/ewas-design - EWAS confounding, batch, inflation, and replication
