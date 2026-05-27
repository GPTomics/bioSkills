# m6A Peak Calling - Usage Guide

## Overview

Call m6A peaks from MeRIP-seq paired IP-vs-input data with multiple callers (exomePeak2 / MeTPeak / MACS3), confirm DRACH motif enrichment as a sanity check on the peak set, reconcile differing peak sets across tools via intersection, and flag the m6A-vs-m6Am ambiguity at 5'UTR peaks that antibody-based methods cannot resolve. Provides BED12 / narrowPeak output ready for downstream differential analysis, motif scanning, and visualisation.

## Prerequisites

```r
BiocManager::install(c('exomePeak2', 'GenomicFeatures', 'BSgenome.Hsapiens.UCSC.hg38',
                       'rtracklayer', 'GenomicRanges', 'ggseqlogo'))
devtools::install_github('compgenomics/MeTPeak')
```

```bash
conda install -c bioconda macs3 macs2 samtools homer bedtools
```

Reference inputs:

- Coordinate-sorted indexed GENOME BAM files for IP and input libraries (from merip-preprocessing)
- Matching GENCODE / Ensembl GTF
- BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38) for exomePeak2 GC correction and motif sequence retrieval

## Quick Start

- "Call m6A peaks with exomePeak2 from paired IP/input BAMs"
- "Run MeTPeak as a second caller and intersect the peak sets"
- "Call broad m6A peaks with MACS3 (viral genome or kilobase-scale peaks)"
- "Confirm DRACH enrichment on my called peaks with HOMER and ggseqlogo"
- "Build a multi-tool consensus peak set with 2-of-3 caller agreement"
- "Flag 5'UTR peaks as m6A-vs-m6Am ambiguous"

## Example Prompts

### Standard Peak Calling

> "Run exomePeak2 with 3 IP and 3 matched input BAMs against GRCh38 + GENCODE v44 GTF; use BSgenome for GC correction; export BED12 + RDS to exomepeak2_output/."

> "Call m6A peaks with MeTPeak using default parameters (window 50, step 50, fragment 100) plus FDR cutoff 0.05; output to metpeak_output/."

> "Call broad MeRIP peaks with MACS3 using --nomodel --extsize 150 --keep-dup all --broad --broad-cutoff 0.1."

### Cross-Caller Reconciliation

> "Intersect exomePeak2, MeTPeak, and MACS3 broadPeak outputs requiring at least 2-of-3 caller agreement; report per-tool counts and consensus count."

> "Compare exomePeak2 peaks against published m6A-Atlas peaks at common transcripts; report concordance percentage."

### DRACH Sanity Check

> "Run HOMER findMotifsGenome.pl on my exomePeak2 peaks in RNA mode; confirm DRACH-like enrichment in top motifs with E-value < 1e-50."

> "Render a ggseqlogo of peak-centre 5-mers to visually confirm DRACH consensus."

### m6A vs m6Am Disambiguation

> "Flag peaks within 50 nt of TSS as m6A-or-m6Am ambiguous; export ambiguous peaks separately from internal peaks."

> "Cross-reference my 5' peaks with published PCIF1-KO MeRIP data to identify m6Am-dominant sites."

### Orthogonal Validation

> "Cross-validate my top 100 exomePeak2 peaks against published GLORI sites in HEK293T; report which are stoichiometry-confirmed."

> "Compare my MeRIP peaks against published miCLIP single-nucleotide m6A sites at common gene loci."

## What the Agent Will Do

1. Build a TxDb from the matched GTF (exomePeak2 / MeTPeak input)
2. Verify BAM and GTF chromosome naming match (chr1 vs 1 reconciliation)
3. Run exomePeak2 with paired IP/input vectors + TxDb + BSgenome (GC-correction)
4. Run MeTPeak with either `GENE_ANNO_GTF=` (file path) or `TXDB=` (uppercase, TxDb object) + IP/Input BAM vectors
5. Run MACS3 with `--nomodel --extsize 150 --keep-dup all --broad --broad-cutoff 0.1`
6. Intersect peak sets via bedtools; build 2-of-3 consensus
7. Run HOMER `findMotifsGenome.pl` on peak set for DRACH sanity check
8. Render ggseqlogo of peak-centre 5-mers
9. Flag peaks within 50 nt of TSS as m6A-or-m6Am ambiguous
10. Annotate peaks against transcript features (5'UTR / CDS / 3'UTR / stop-codon) via ChIPseeker
11. Cross-reference against published m6A-Atlas / REPIC if requested
12. Report per-tool peak counts, consensus count, cross-replicate overlap, DRACH enrichment E-value

## Tips

- exomePeak2 is the field default for transcript-aware GC-corrected calls. MeTPeak HMM smoothing helps at low coverage. MACS3 broad mode is for viral / kb-scale peaks.
- DRACH is a SANITY CHECK on the peak set, NEVER a per-peak filter. Filtering drops 5-10% of real m6A sites.
- Peaks within ~50 nt of TSS are m6A-or-m6Am ambiguous because anti-m6A antibodies cross-react with m6Am (cap-adjacent N6,2'-O-dimethyladenosine installed by PCIF1).
- Always pass `--keep-dup all` to MACS3 / MACS2 for MeRIP. Default `--keep-dup 1` destroys signal.
- MeTPeak accepts `GENE_ANNO_GTF=` (file path) OR `TXDB=` (uppercase, TxDb object); exomePeak2 uses lowercase `txdb=`. Don't confuse argument casing between the two APIs.
- Cross-tool peak overlap is typically ~70% between exomePeak2 and MeTPeak; report multi-tool consensus to bound the false-positive rate.
- For absolute stoichiometry claims (not just enrichment), use GLORI (Liu 2023) or SAC-seq (Hu 2022); MeRIP is qualitative.
- METTL3-KO is the gold-standard control for "is this peak m6A?"; peaks that remain in KO are antibody artifacts or non-METTL3 modifications.
- McIntyre 2020 documented ~45% median between-lab peak overlap. Treat single-study peak sets as one observation among many.

## Related Skills

- merip-preprocessing - Upstream IP/input BAM preparation
- m6a-differential - Compare peak sets between conditions
- m6anet-analysis - Orthogonal validation via ONT direct-RNA
- modification-visualization - Metagene, browser tracks, heatmaps
- clip-seq/peak-calling - miCLIP / m6A-CLIP single-nucleotide validation
- chip-seq/peak-calling - General IP-vs-input peak calling concepts
- chip-seq/peak-annotation - Gene-feature annotation
- rna-quantification/featurecounts-counting - Peak count matrices for downstream differential
- pathway-analysis/go-enrichment - GO enrichment on peak-bearing genes
- data-visualization/multipanel-figures - Figure assembly
