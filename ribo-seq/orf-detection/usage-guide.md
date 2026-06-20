# ORF Detection - Usage Guide

## Overview

Detect actively translated open reading frames from Ribo-seq data using 3-nucleotide periodicity (not coverage) as the evidence of translation: canonical CDS, uORFs, internal ORFs, dORFs, and novel ORFs. Classify them by the 2022 community standard, quantify ORF-level translation, and validate calls before trusting them.

## Prerequisites

```bash
pip install RiboCode
```

```r
# Isoform-aware quantification and the general toolkit (distinct packages)
BiocManager::install(c("ORFquant", "ORFik", "DESeq2"))
```

## Quick Start

Tell your AI agent what you want to do:
- "Find translated ORFs in my Ribo-seq data with RiboCode"
- "Detect uORFs including near-cognate (CUG/GUG) starts"
- "Quantify ORF-level translation isoform-aware with ORFquant"
- "Validate my called ORFs (in-frame fraction, FLOSS)"

## Example Prompts

### ORF Detection

> "Run RiboCode on my transcriptome BAM and report ORFs with adjusted p < 0.05"

> "Include CUG and GTG near-cognate start codons in the uORF search"

> "How many novel ORFs and uORFs are detected?"

### Classification and Nomenclature

> "Classify my ORFs by the Mudge 2022 categories (uORF, uoORF, intORF, dORF, doORF, lncRNA-ORF)"

> "Find translated ORFs on long non-coding RNAs"

### Quantification and Validation

> "Quantify ORF-level translation isoform-aware with ORFquant"

> "Compute FLOSS and in-frame fraction to validate my called ORFs"

> "Compare ORF-level translation between conditions with DESeq2"

## What the Agent Will Do

1. Prepare transcript annotation and run metaplots to select periodic read lengths
2. Call ORFs with RiboCode (optionally with near-cognate starts)
3. Classify ORFs by the standard positional categories
4. Quantify isoform-aware with ORFquant or count P-sites per ORF
5. Validate with in-frame fraction, FLOSS, and conservation/peptide evidence

## Tips

- **Periodicity, not coverage** - frame-0 enrichment is the translation signal
- **Near-cognate starts** - uORFs often begin at CUG/GUG; an ATG-only scan misses them
- **RiboCode -l is a toggle** - read lengths come from the metaplots config, not -l
- **ORFik is not ORFquant** - distinct packages; install the one you mean
- **ORF types** - RiboCode emits Overlap_uORF/Overlap_dORF/Internal/novel (no "noncoding")
- **Validate novel ORFs** - FLOSS, in-frame fraction, PhyloCSF, and mass-spec peptides

## Related Skills

- ribosome-periodicity - Calibrate the P-site offsets ORF callers consume
- initiation-site-mapping - Map non-AUG start codons from TI-seq data
- translation-efficiency - Add an RNA-seq denominator for TE
- differential-expression/deseq2-basics - Differential ORF-level translation
