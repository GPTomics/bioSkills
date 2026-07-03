---
name: bio-workflows-smrna-pipeline
description: End-to-end small RNA-seq analysis from FASTQ to differential miRNA expression and target prediction. Use when analyzing miRNA, isomiR, tRF, or piRNA sequencing data through preprocessing, quantification, discovery, DE, and targets.
tool_type: mixed
primary_tool: miRDeep2
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: cutadapt 4.4+, miRDeep2 2.0.1.3+, miRge3.0 0.1.4+, DESeq2 1.42+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Small RNA-seq Pipeline

**"Analyze my small RNA-seq data from FASTQ to differential miRNAs"** -> Orchestrate adapter trimming (cutadapt), quantification (miRge3) or discovery (miRDeep2), differential expression (DESeq2), and target prediction (miRanda), with the small-RNA-specific decisions surfaced at each step.

## Pipeline Overview

```
FASTQ -> cutadapt trim -> miRge3 (known) / miRDeep2 (novel) -> DESeq2 -> target prediction
                                                |
                                                +-> tRF / piRNA profiling (MINTmap / unitas)
```

## Step 1: Preprocessing

```bash
# The 3' adapter is on EVERY real read (insert ~22 nt << read length), so --discard-untrimmed
# drops adapter-dimer / no-insert junk - the inverse of genomic DNA. For NEXTflex 4N libraries,
# strip the random bases AFTER adapter removal (cutadapt -u 4 -u -4); for QIAseq UMIs, extract
# and dedup with umi_tools. Never position-dedup a non-UMI small-RNA library.
cutadapt -a TGGAATTCTCGGGTGCCAAGG \
    --minimum-length 18 --maximum-length 30 --discard-untrimmed \
    -o trimmed.fastq.gz reads.fastq.gz
```

## Step 2: Quantification or discovery

```bash
# Known-miRNA + isomiR quantification (the common case) - fast, curated libraries:
miRge3.0 annotate -s trimmed.fastq.gz -lib /path/to/miRge3_Lib -on human -db miRBase \
    -a illumina -gff -ai -cpu 8 -o mirge_out

# OR novel discovery (high false-positive; needs a genome and bowtie 1):
mapper.pl trimmed.fastq.gz -e -h -i -j -l 18 -m -p genome_index \
    -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf
miRDeep2.pl reads_collapsed.fa genome.fa reads_collapsed_vs_genome.arf \
    mature_ref.fa none hairpin_ref.fa -t Human
# Choose the miRDeep2 score cutoff from survey.pl signal-to-noise; it is not a fixed rule.
```

## Step 3: Differential expression

```r
library(DESeq2)
# Feed RAW counts (miR.Counts.csv), NOT RPM. A few miRNAs can dominate the library, so
# inspect sizeFactors for compositional distortion before trusting the calls.
counts <- read.csv('mirge_out/miR.Counts.csv', row.names = 1)
dds <- DESeqDataSetFromMatrix(round(counts), colData, ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # lower prefilter than mRNA
dds <- DESeq(dds)
print(sizeFactors(dds))
res <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
```

## Step 4: Target prediction

```bash
# Prediction is a hypothesis: intersect with anti-correlated mRNA DE before trusting targets.
miranda mature_mirnas.fa target_3utrs.fa -sc 140 -en -20 -strict -out targets.txt
```

## QC Checkpoints

0. **Before quantification**: run miRTrace for length/complexity, RNA-class composition (miRNA vs rRNA/tRNA/artifact), and cross-clade contamination (a good mapping rate hides sample swaps and reagent contamination). For a fully reproducible pipeline, nf-core/smrnaseq chains these steps.
1. **After trimming**: read-length distribution should peak at 21-23 nt; a 30+ nt smear means degradation or tRNA/rRNA contamination, a 26-32 nt peak means piRNA (for PANDORA/phospho preps the broad tRF/rRF distribution is expected)
2. **Ligation bias**: absolute cross-miRNA abundance within a sample is untrustworthy; compare the same miRNA across samples, never across kits
3. **After alignment**: check the RNA-class composition (miRNA vs tRF/rRF/piRNA); abundant non-miRNA classes mean this is a tRF/piRNA story (see trf-pirna-profiling)
4. **Before DE**: confirm RAW counts (not RPM) and inspect size factors for compositional distortion
5. **After DE**: report base mean with every call; a significant fold-change on a ~5-count miRNA is noise
6. **Targets**: filter predictions by anti-correlated mRNA DE and validated databases, never enrichment on raw predictions

## Related Skills

- small-rna-seq/smrna-preprocessing - Kit-specific adapter, UMI, and 4N handling
- small-rna-seq/mirdeep2-analysis - Novel miRNA discovery
- small-rna-seq/mirge3-analysis - Known-miRNA and isomiR quantification
- small-rna-seq/differential-mirna - Compositionally-aware DE
- small-rna-seq/target-prediction - Seed prediction filtered by expression
- small-rna-seq/trf-pirna-profiling - tRF and piRNA profiling
- workflow-management/nf-core-pipelines - Run nf-core/smrnaseq as a curated, reproducible pipeline instead of chaining the steps by hand
