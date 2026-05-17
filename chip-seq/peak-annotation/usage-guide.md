# Peak Annotation - Usage Guide

## Overview

Annotate ChIP-seq peaks to genomic features (promoter, exon, intron, intergenic), nearest or host genes, ENCODE candidate cis-regulatory elements (cCREs: PLS, pELS, dELS, CTCF-only, DNase-H3K4me3), regulatory domains via GREAT, and cell-type-specific enhancer-gene targets via ENCODE-rE2G. Supports ChIPseeker (R), HOMER annotatePeaks.pl (CLI), pyranges (Python), rGREAT, chipenrich, Broad-Enrich, and direct cCRE BED intersection.

## Prerequisites

```r
BiocManager::install(c('ChIPseeker', 'GenomicFeatures', 'rtracklayer',
                       'rGREAT', 'chipenrich',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'))
```

```bash
conda install -c bioconda homer bedtools
pip install pandas pyranges
```

## Quick Start

Tell the agent what to do:
- "Annotate my peaks to nearest genes with the provided GTF, host-gene convention"
- "Classify each peak as promoter / exon / intron / intergenic with 2 kb promoter window"
- "Intersect peaks with ENCODE cCREs and tell me what fraction are PLS / pELS / dELS / CTCF-only"
- "Run GREAT regulatory-domain gene-set enrichment for distal H3K27ac peaks"
- "Run ChIP-Enrich with locus-length adjustment for GO enrichment on TF peaks"
- "Annotate broad H3K27me3 domains using broadenrich method"
- "Plot the peak distribution by feature type and signed distance to TSS"

## Example Prompts

### Standard ChIPseeker annotation
> "Annotate my H3K4me3 narrowPeak file using ChIPseeker with TxDb.Hsapiens.UCSC.hg38.knownGene, 2 kb promoter window, host-gene convention (`overlap='all'`)."

### Custom GTF
> "Annotate peaks using a GENCODE v44 GTF I'm providing. Build TxDb from the GTF, map gene symbols from the GTF directly (annoDb doesn't work for custom TxDb), strip Ensembl version suffixes when joining."

### ENCODE cCRE classification
> "Download the ENCODE cCRE BED from SCREEN and intersect my peaks to classify them as PLS, pELS, dELS, CTCF-only, or DNase-H3K4me3. Report counts per class."

### GREAT regulatory-domain enrichment
> "I have ~5000 distal H3K27ac peaks. Run rGREAT with default basal+extension regulatory domain rules for GO:BP enrichment. Filter out hyper-ChIPable peaks first."

### ChIP-Enrich (locus-length-adjusted)
> "Run chipenrich on my CTCF peaks with `locusdef='nearest_tss'` and `method='chipenrich'` for GO:BP enrichment. The locus-length adjustment is important for fair comparison."

### Broad-Enrich for broad histones
> "Run broadenrich on H3K27me3 broadPeak file; the method accounts for region width which is critical for broad domains."

### Enhancer-gene linking
> "Link H3K27ac enhancer peaks to target genes using ENCODE-rE2G for K562 cells."

### Decoupling diagnostic
> "Show me peaks where the nearest-TSS gene differs from the host gene (peak in gene body)."

## What the Agent Will Do

1. **Determine annotation source**: pre-built TxDb (standard genomes) or makeTxDbFromGFF (custom GTF)
2. **Verify genome assembly match**: TxDb genome must match BAM/peak alignment genome
3. **Choose convention**: nearest-TSS (HOMER, ChIPseeker default) or host-gene (ChIPseeker `overlap='all'`) based on biology
4. **Coordinate-system conversion**: GTF 1-based to BED 0-based for TSS positions
5. **Compute signed distance to TSS**: strand-aware (negative = upstream)
6. **Classify feature category**: promoter (within window) > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic
7. **Optional: ENCODE cCRE intersection**: assign each peak to a cCRE class if it overlaps
8. **Optional: enhancer-gene linking**: ENCODE-rE2G or ABC for distal regulatory peaks
9. **Optional: gene-set enrichment**: GREAT (regulatory domain) or ChIP-Enrich (locus-length-adjusted)
10. **Visualize**: annotation pie chart, distance-to-TSS distribution, per-feature stacked bar
11. **Export**: annotated TSV with peak coords, gene, distance, feature, cCRE class
12. **Document**: annotation source + version, convention used, promoter window, locusdef for enrichment

## Tips

- **Match genome assembly explicitly.** Silent assembly mismatch produces wrong distances and feature calls.
- **Custom GTF needs custom symbol mapping.** ChIPseeker's `annoDb` parameter doesn't work with `makeTxDbFromGFF`; map symbols from the GTF directly.
- **Promoter window matters.** ChIPseeker default is `c(-3000, 3000)`; many studies use `c(-2000, 2000)` or `c(-1000, 500)`. Match to the analysis requirements.
- **HOMER promoter window is hard-coded at -1 kb / +100 bp.** Reclassify by Distance to TSS column for custom windows.
- **Strip Ensembl version suffixes when joining.** GENCODE GTFs have `ENSG00000142192.25`; the TxDb may store version-less IDs.
- **Use host-gene convention for gene-body marks (H3K36me3, H3K27me3) and exonic peaks.** Use nearest-TSS for distal TF binding.
- **For enhancer-gene linking, use ENCODE-rE2G or ABC, not nearest-TSS.** Enhancers regulate non-nearest genes in 30-50% of cases.
- **Filter hyper-ChIPable peaks before gene-set enrichment.** rRNA / mtDNA / housekeeping artifacts inflate "ribosomal" and "translation" GO terms.
- **GREAT for distal regulatory; ChIP-Enrich for promoter-focused; broadenrich for broad histone marks.** Different statistical models address different bias structures.
- **Cross-reference against ENCODE cCRE atlas.** Useful sanity check: a TF that should bind enhancers should have peaks dominated by dELS / pELS, not PLS or CTCF-only.

## Troubleshooting

### `seqlevels` mismatch in ChIPseeker

Chromosome naming convention (chr1 vs 1) differs between peaks and TxDb. Fix with:
```r
seqlevelsStyle(peaks) <- 'UCSC'  # or 'Ensembl'
```

### Gene symbols all NA

Using `annoDb='org.Hs.eg.db'` with a custom `makeTxDbFromGFF` TxDb. Map symbols from the original GTF post-hoc; strip Ensembl version suffixes before joining.

### HOMER "no annotation"

Genome not configured for HOMER. Run `perl configureHomer.pl -install hg38` once.

### GREAT timeout / network errors

rGREAT v2.x runs locally without biomart; rGREAT v1.x used biomart over the network. For large peak sets, use `rGREAT::great()` with pre-loaded TxDb and gene sets.

### Unexpected "ribosomal" / "translation" GO terms dominate

Hyper-ChIPable artifacts at ribosomal protein genes and rRNA processing factors. Filter peaks against blacklist + custom hyper-ChIPable BED (top-1% input signal) before enrichment.

### ENCODE-rE2G cell type not in registry

ENCODE-rE2G has cell-type-specific weights for a fixed set of cell types. For other cell types, use ABC with cell-type-matched ATAC + H3K27ac + Hi-C contact freq.

### Differential peaks annotated to different genes than original peak set

When peak boundaries shift (e.g., differential peaks are summit-recentered ±250 bp), nearest-TSS can change. Run annotation consistently on the original consensus peakset, then map differential calls.

## Related Skills

- chip-seq/peak-calling - Generate peaks for annotation
- chip-seq/chipseq-qc - Filter hyper-ChIPable peaks before enrichment
- chip-seq/super-enhancers - Annotate SE-associated genes
- chip-seq/differential-binding - Annotate differential peaks
- chip-seq/cut-and-run-tag - CUT&RUN/CUT&Tag peak annotation (same tools, different upstream)
- atac-seq/enhancer-gene-linking - ENCODE-rE2G workflow for distal regulatory peaks
- pathway-analysis/go-enrichment - Downstream GO enrichment
- pathway-analysis/reactome-pathways - Reactome pathway enrichment
- genome-intervals/gtf-gff-handling - GTF parsing utilities
- genome-intervals/proximity-operations - bedtools closest / window
