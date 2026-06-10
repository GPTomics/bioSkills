# methylation-analysis

## Overview

DNA methylation analysis across both major platforms: short-read bisulfite/EM-seq sequencing (WGBS, RRBS, PBAT) and Illumina Infinium arrays (450K, EPIC, EPICv2). Covers alignment of C->T-converted reads, methylation calling, per-CpG differential testing, DMR detection, array preprocessing and QC, cell-type deconvolution, epigenetic clocks, and EWAS study design. Decision-grade framing throughout: bisulfite conversion efficiency masquerades as methylation, an array beta is a two-chemistry ratio that is meaningless until corrected, a per-CpG test is dictated by the data object, a DMR is a region-definition choice, and an EWAS hit is a cell-composition difference until proven otherwise. Native long-read MM/ML modification calling lives in long-read-sequencing/nanopore-methylation.

**Tool type:** mixed | **Primary tools:** Bismark (CLI), methylKit/dmrseq/DSS (R), sesame/minfi (R, arrays), EpiDISH (R, deconvolution), methylclock (R, clocks), meffil/limma (R, EWAS), scipy (Python)

## Skills

| Skill | Description |
|-------|-------------|
| bismark-alignment | Bisulfite/EM-seq read alignment with Bismark, library/strand model, conversion QC |
| methylation-calling | Per-CpG calling from BAM (Bismark/MethylDackel), conversion QC, contexts, variant-aware |
| methylkit-analysis | methylKit object model: import, filter, normalize, unite, calculateDiffMeth |
| differential-cpg-testing | Per-CpG/per-site testing: count-vs-continuous fork, M-values, differential variability |
| dmr-detection | Region callers (dmrseq/DSS/DMRcate/comb-p), post-selection inference, PMD segmentation |
| array-preprocessing | Infinium IDAT to corrected beta/M matrix: Type I/II, normalization, detection masking |
| array-qc-filtering | Probe/sample QC: cross-reactive/SNP/sex filtering, EPICv2 collapse, sample identity |
| cell-type-deconvolution | Reference-based/free cell-fraction estimation; cell-type-resolved EWAS |
| epigenetic-clocks | DNAm age and age acceleration (Horvath/GrimAge/DunedinPACE); reliability, EPICv2 dropout |
| ewas-design | EWAS confounding hierarchy, batch/SVA, genomic inflation/BACON, thresholds, replication |

## Example Prompts

- "Align my bisulfite sequencing reads with Bismark and check conversion efficiency"
- "Extract CpG methylation calls from a bwa-meth BAM with MethylDackel"
- "Load my coverage files into methylKit and correct for overdispersion"
- "Test individual CpG sites for differential methylation on M-values"
- "Find DMRs between conditions with a selection-aware region FDR"
- "Process EPICv2 IDATs into a normalized beta matrix with sesame"
- "Filter cross-reactive and SNP-overlapping probes and check for sample swaps"
- "Estimate blood cell-type proportions to adjust my EWAS"
- "Compute Horvath DNAm age and age acceleration from my array data"
- "Design an EWAS: which covariates, what genome-wide threshold, how to handle batch"

## Requirements

```bash
# Bisulfite alignment + calling (CLI)
conda install -c bioconda bismark bowtie2 hisat2 samtools trim-galore methyldackel

# Python (per-CpG testing, EWAS helpers)
pip install numpy pandas scipy statsmodels
```

```r
# Sequencing downstream
BiocManager::install(c('methylKit', 'bsseq', 'dmrseq', 'DSS', 'DMRcate', 'limma', 'missMethyl', 'GenomicRanges'))

# Arrays, deconvolution, clocks, EWAS
BiocManager::install(c('sesame', 'minfi', 'maxprobes', 'EpiDISH', 'FlowSorted.Blood.EPIC', 'methylclock', 'sva', 'bacon'))
# meffil, dnaMethyAge, methylCIPHER are installed from GitHub
```

## Related Skills

- **long-read-sequencing** - Native long-read MM/ML modification calling (different platform)
- **read-qc** - Adapter trimming and QC before bisulfite alignment
- **alignment-files** - BAM file manipulation after alignment
- **sequence-io** - FASTQ handling before alignment
- **pathway-analysis** - CpG-bias-aware functional enrichment of genes near DMRs
- **cell-type-deconvolution / single-cell** - scRNA atlases for deconvolution references
- **causal-genomics** - mQTL-based Mendelian randomization for EWAS causal orientation
- **machine-learning** - Predictor training behind clocks and methylation risk scores
- **clinical-biostatistics** - Survival/mortality modeling of epigenetic age acceleration
