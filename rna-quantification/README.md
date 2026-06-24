# rna-quantification

## Overview

Quantify gene and transcript expression from RNA-seq data. Covers BAM-based counting with featureCounts and alignment-free quantification with Salmon/kallisto, plus import to R with the correct length offset and QC of the count matrix before differential expression. Two paths exist: align-then-count (deterministic) and quantify-then-import (probabilistic, resolves multi-mapping). Both feed gene-level DE through tximport; choose counting for specified features (exons, custom intervals) or when an EM model is unwanted, and alignment-free for accuracy on multi-isoform genes and any transcript-level work.

**Tool type:** mixed | **Primary tools:** featureCounts, Salmon, kallisto, tximport

## Skills

| Skill | Description |
|-------|-------------|
| featurecounts-counting | Count reads per gene from BAM files; strandedness, fragments, multimappers, summary QC |
| alignment-free-quant | Salmon/kallisto quantification; decoy index, library type, bias, inferential replicates |
| tximport-workflow | Import transcript estimates with the length offset; countsFromAbundance, version IDs |
| count-matrix-qc | Normalization, VST/rlog, PCA, Cook's outliers, batch and sample-swap checks |

## Example Prompts

- "Count reads per gene from my BAM files"
- "Determine the strandedness of my RNA-seq before counting"
- "Decide whether my paralog-heavy genes need EM quantification instead of counting"
- "Quantify transcripts without aligning to the genome"
- "Build a decoy-aware index and explain why it matters"
- "Generate inferential replicates for transcript-level testing"
- "Import Salmon results into R with the correct length offset"
- "Choose a countsFromAbundance mode for my 3'-tag library"
- "Resolve a transcript-ID version mismatch on import"
- "Check my count matrix for outliers and batch effects before DE"
- "Decide whether to drop a suspected outlier sample"
- "Diagnose whether PC1 is biology or sequencing depth"

## Requirements

```bash
# featureCounts (part of Subread)
conda install -c bioconda subread

# Salmon
conda install -c bioconda salmon

# kallisto
conda install -c bioconda kallisto
```

```r
BiocManager::install(c('tximport', 'tximeta'))
```

## Related Skills

- **read-qc** - Upstream quality control
- **alignment-files** - BAM file processing
- **differential-expression** - Downstream gene-level DE (and catchSalmon transcript DTE)
- **alternative-splicing** - DTU and swish transcript-level testing from quantification
- **expression-matrix** - Count ingestion, sleuth, and ID mapping
- **genome-intervals** - GTF/GFF annotation handling
