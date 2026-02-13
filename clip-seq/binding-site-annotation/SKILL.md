---
name: bio-clip-seq-binding-site-annotation
description: Annotate CLIP-seq binding sites to genomic features including 3'UTR, 5'UTR, CDS, introns, and ncRNAs. Use when characterizing where an RBP binds in transcripts.
tool_type: mixed
primary_tool: ChIPseeker
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Binding Site Annotation

**"Annotate where my RBP binds in transcripts"** → Map CLIP-seq peaks to genomic features (3'UTR, 5'UTR, CDS, introns, ncRNAs) to characterize RNA-binding protein target regions.
- R: `ChIPseeker::annotatePeak()` with transcript annotation databases
- CLI: `bedtools intersect` with gene model BED files

## Using ChIPseeker (R)

**Goal:** Classify CLIP-seq binding sites by genomic feature (3'UTR, 5'UTR, CDS, intron).

**Approach:** Load peaks and a TxDb transcript database, annotate with annotatePeak, and visualize the feature distribution with a pie chart.

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peaks <- readPeakFile('peaks.bed')
anno <- annotatePeak(peaks, TxDb = txdb)

plotAnnoPie(anno)
```

## Using BEDTools

```bash
# Annotate to UTRs
bedtools intersect -a peaks.bed -b 3utr.bed -wa -wb > peaks_3utr.bed
```

## Python Annotation

```python
import pandas as pd

def annotate_peaks(peaks_bed, annotation_gtf):
    '''Annotate peaks to genomic features'''
    # Load peaks and annotations
    # Intersect and categorize
    pass
```

## Related Skills

- clip-peak-calling - Get peaks
- genome-intervals/interval-arithmetic - Intersect peaks with genomic features
