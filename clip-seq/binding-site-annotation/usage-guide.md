# Binding Site Annotation - Usage Guide

## Overview

Annotate CLIP-seq peaks to RNA features (5' UTR, CDS, 3' UTR, intron, splice site, ncRNA, repeat element). Interpretation is RBP-class-specific: splicing factors at splice junctions, HuR/PUM2 in 3' UTRs, EIF3 in 5' UTRs, FASTKD2 on chrM, repeat-binding RBPs (MATR3, HNRNPK) on Alu/LINE/LTR. ChIPseeker is the fast tool; RBP-Maps (Yeo lab) is mandatory for splicing factors; RepeatMasker overlap is a separate axis. Default ChIPseeker `tssRegion=c(-3000,3000)` over-extends for CLIP - tighten to `c(-100,100)`.

## Prerequisites

```r
BiocManager::install(c('ChIPseeker','RCAS','GenomicFeatures',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene'))
```
```bash
conda install -c bioconda bedtools deeptools rseqc
# RBP-Maps (Yeo)
git clone https://github.com/YeoLab/rbp-maps
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate peaks to 5' UTR / CDS / 3' UTR / intron with ChIPseeker, transcript-level, tight TSS region"
- "Build the RBP-Maps splicing regulatory map for my PTBP1 eCLIP"
- "Add a repeat-element axis to my annotation (Alu / LINE / LTR / SINE)"
- "Metagene of CLIP signal across 3' UTRs with deepTools computeMatrix"
- "Why do 50% of my peaks fall in 'Promoter (2-3kb)'? Tighten TSS region"
- "Sum of region counts > total peaks - enforce hierarchy"
- "FASTKD2 should bind chrM - my TxDb excludes it"

## Example Prompts

### Global Region Distribution

> "ChIPseeker annotation with tssRegion=c(-100,100), level=transcript, then plotAnnoPie and plotDistToTSS"

> "What fraction of my HuR peaks are in 3' UTRs versus introns? Expected > 50% 3' UTR"

### Splicing Factor Maps

> "Run RBP-Maps for PTBP1 - generate the 1400 nt cassette-exon regulatory metagene"

> "Need cassette exon BED from ENCODE shRNA-KD RNA-seq before RBP-Maps will produce signal"

### Repeat-Element Axis

> "MATR3 binds LINE-1; intersect peaks with RepeatMasker BED and report repeat-class fraction"

> "My peaks are 30% Alu - is this real or alignment artifact? Compare to ENCODE Alu binders"

### Metagene

> "deepTools metagene of CLIP signal scaled across 3' UTRs - look for stop-codon-proximal peak"

> "RSeQC geneBody_coverage.py to check for 5' vs 3' positional bias"

### Diagnostics

> "Many peaks labeled 'Distal Intergenic' - is my TxDb missing transcripts?"

> "RBP-Maps metagene is flat - I don't have the cassette exon table"

## What the Agent Will Do

1. Choose ChIPseeker (fast global) + RCAS (RNA-aware details) + RBP-Maps (splicing) based on RBP class
2. Set `tssRegion=c(-100,100)`, `level='transcript'`, customized priority hierarchy
3. Annotate peaks; report per-region pie, distance-to-TSS, distance-to-stop-codon
4. Cross-check expected vs observed dominant region (HuR -> 3' UTR; PTBP1 -> intron; FASTKD2 -> chrM)
5. Add repeat-element axis via RepeatMasker intersect (`-s` for strand)
6. For splicing factors: generate the RBP-Maps cassette-exon regulatory metagene
7. Flag anomalies: too many "Promoter", too many "Intergenic", missing chrM, repeat over-counting

## Tips

- **Default `tssRegion` is for ChIP, not CLIP.** Tighten to `c(-100,100)` or `c(-50,50)`.
- **`level='transcript'` not `level='gene'`.** Isoform context matters for splicing factors.
- **RBP-Maps is mandatory for splicing factors.** PTBP1, U2AF2, RBFOX, SRSF1-9, MBNL, HNRNPC.
- **RepeatMasker is a separate axis.** Don't lose repeat-binding biology in the region axis.
- **chrM peaks can be lost.** Verify TxDb includes mt-transcripts (newer GENCODE does).
- **`bedtools intersect -s` always.** CLIP is strand-specific.
- **Hierarchy must be enforced.** Sum of region counts == total peaks (no double-counting).
- **Cassette exon BED needed for RBP-Maps.** From companion RNA-seq KD or ENCODE shRNA tables.
- **Metagene reveals position-dependent biology.** m6A readers stop-codon-proximal; PUM2 3' UTR proximal.
- **GENCODE v38+ for hg38.** Older GTFs miss lincRNAs and snoRNAs.

## Related Skills

- clip-seq/clip-peak-calling - Upstream peaks
- clip-seq/crosslink-site-detection - Fine-grained metagene
- clip-seq/clip-motif-analysis - Motifs in annotated regions
- clip-seq/differential-clip - Cross-condition regulatory maps
- clip-seq/ago-clip-mirna-targets - miRNA seed sites
- clip-seq/m6a-clip - Stop-codon m6A annotation
- genome-intervals/gtf-gff-handling - GTF prep
- genome-intervals/interval-arithmetic - bedtools patterns
- alternative-splicing/differential-splicing - Cassette tables
