---
name: bio-chipseq-peak-annotation
description: Annotate ChIP-seq peaks to genomic features and nearest genes. Classify peaks as promoter, exon, intron, or intergenic using ChIPseeker (R), HOMER annotatePeaks.pl (CLI), or Python (pandas/pyranges). Supports pre-built annotation databases and custom GTF files. Handles promoter definition, feature priority, category collapsing, and signed distance-to-TSS. Use when assigning genomic context to ChIP-seq peaks or linking peaks to target genes.
tool_type: mixed
primary_tool: ChIPseeker
---

## Version Compatibility

Reference examples tested with: ChIPseeker 1.38+, GenomicFeatures 1.54+, rtracklayer 1.62+, HOMER 4.11+, pyranges 0.0.129+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `annotatePeaks.pl` (HOMER prints version on run)

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Peak Annotation

**"Annotate my peaks to genes and genomic features"** -> Assign each peak to a genomic feature category (promoter, exon, intron, intergenic), find the nearest gene, and calculate signed distance to TSS.
- R: `ChIPseeker::annotatePeak(peaks, TxDb=txdb)`
- CLI: `annotatePeaks.pl peaks.bed hg38 -gtf annotation.gtf`
- Python: Parse GTF, compute intervals, classify by overlap

## Choosing an Annotation Approach

| Context | Recommended | Why |
|---------|-------------|-----|
| Standard genome, pre-built annotations available | ChIPseeker with TxDb package | Simplest setup; automatic gene symbol mapping via annoDb |
| Custom or project-specific GTF provided | ChIPseeker + makeTxDbFromGFF, HOMER -gtf, or Python | All three handle custom annotations; choose based on pipeline |
| HOMER already in pipeline | HOMER annotatePeaks.pl | Reuses tag directory; combined annotation + motif workflow |
| Fine-grained control over classification logic | Python | Full control over priority rules, distance calculation, output format |
| Quick annotation with standard categories | HOMER annotatePeaks.pl | Single command; no R environment required |

**Critical:** Always use the same annotation source for peak annotation as was used for the analysis. Mixing databases (e.g., UCSC knownGene TxDb with a GENCODE GTF alignment) produces gene name mismatches and incorrect feature assignments. When a specific GTF is provided, use it directly rather than a pre-built TxDb package.

## Coordinate Systems and TSS

BED files use 0-based half-open coordinates `[start, end)`. GTF files use 1-based closed coordinates `[start, end]`. Mixing these without conversion shifts annotations by one base.

**Peak center** (from BED): `(start + end) // 2`

**TSS from GTF:**
- Plus-strand genes: TSS = `start` (1-based) -> 0-based: `start - 1`
- Minus-strand genes: TSS = `end` (1-based) -> 0-based: `end`

**Signed distance:** Negative = upstream of TSS, positive = downstream.
- Plus-strand: `distance = peak_center - tss`
- Minus-strand: `distance = -(peak_center - tss)`

## Annotation with ChIPseeker (R)

### Pre-built TxDb (Standard Genomes)

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile('peaks.narrowPeak')
peak_anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = 'org.Hs.eg.db')
anno_df <- as.data.frame(peak_anno)
```

### Custom GTF Annotations

**Goal:** Annotate peaks using a project-specific GTF rather than a pre-built annotation package.

**Approach:** Build a TxDb from the GTF with `makeTxDbFromGFF()`, annotate peaks against it, then map gene symbols from the original GTF since custom TxDb objects lack the ID mappings that `annoDb` requires.

```r
library(ChIPseeker)
library(GenomicFeatures)
library(rtracklayer)

txdb <- makeTxDbFromGFF('genes.gtf.gz', format = 'gtf')
peaks <- readPeakFile('peaks.bed')
peak_anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-2000, 2000))
anno_df <- as.data.frame(peak_anno)

# Map gene symbols from GTF (annoDb does not work with custom TxDb)
gtf <- import('genes.gtf.gz')
gene_map <- unique(data.frame(
    gene_id = sub('\\..*', '', gtf$gene_id),
    symbol = gtf$gene_name, stringsAsFactors = FALSE))
gene_map <- gene_map[!is.na(gene_map$symbol), ]
anno_df$geneId_base <- sub('\\..*', '', anno_df$geneId)
anno_df$SYMBOL <- gene_map$symbol[match(anno_df$geneId_base, gene_map$gene_id)]
```

The version-suffix stripping (`sub('\\..*', '', ...)`) handles GENCODE gene IDs like `ENSG00000142192.25` where the TxDb may store the full ID but the GTF attribute has it without the version.

### Custom Promoter Definition

The `tssRegion` parameter defines the promoter window around each TSS:

| Window | Use case |
|--------|----------|
| c(-1000, 1000) | Strict core promoter |
| c(-2000, 2000) | Common custom definition |
| c(-3000, 3000) | ChIPseeker default; broader capture |
| c(-2000, 500) | Asymmetric; emphasizes upstream regulatory elements |

Match this to analysis requirements. Many studies define specific windows (e.g., 2kb symmetric) -- always check.

### Annotation Priority

ChIPseeker resolves overlapping features using `genomicAnnotationPriority`:

Default: `Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic`

A peak overlapping both a promoter of gene A and an intron of gene B receives "Promoter". To customize:

```r
peak_anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-2000, 2000),
    genomicAnnotationPriority = c('Promoter', '5UTR', '3UTR', 'Exon', 'Intron',
                                   'Downstream', 'Intergenic'))
```

### Collapse Annotation Categories

ChIPseeker returns detailed subcategories. To collapse to four standard categories:

| ChIPseeker Output | Collapsed |
|-------------------|-----------|
| Promoter (<=1kb), Promoter (1-2kb), Promoter (2-3kb) | promoter |
| 5' UTR, 3' UTR, 1st Exon, Other Exon | exon |
| 1st Intron, Other Intron | intron |
| Downstream (<=300), Downstream (<=1kb), Distal Intergenic | intergenic |

```r
collapse_annotation <- function(ann) {
    ifelse(grepl('Promoter', ann), 'promoter',
    ifelse(grepl("5' UTR|3' UTR|Exon", ann), 'exon',
    ifelse(grepl('Intron', ann), 'intron', 'intergenic')))
}
anno_df$feature <- collapse_annotation(anno_df$annotation)
```

To suppress subcategories at the source:

```r
options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_promoter_subcategory = TRUE)
```

### Export Results

```r
output <- data.frame(chr = anno_df$seqnames, start = anno_df$start, end = anno_df$end,
    nearest_gene = anno_df$SYMBOL, distance_to_tss = anno_df$distanceToTSS,
    feature = anno_df$feature)
write.table(output, 'annotations.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
```

## Annotation with HOMER (CLI)

### Standard Genome

```bash
annotatePeaks.pl peaks.bed hg38 > annotated.txt
```

### Custom GTF

```bash
# With installed genome
annotatePeaks.pl peaks.bed hg38 -gtf genes.gtf > annotated.txt

# Without installed genome (annotation from GTF only)
annotatePeaks.pl peaks.bed none -gtf genes.gtf > annotated.txt
```

For gzipped GTFs, decompress first: `gunzip -k genes.gtf.gz`

### Annotation Statistics

```bash
annotatePeaks.pl peaks.bed hg38 -gtf genes.gtf -annStats ann_stats.txt > annotated.txt
```

### Parse HOMER Output

HOMER produces 19 tab-delimited columns. Key columns for annotation:

| Column | Name | Content |
|--------|------|---------|
| 8 | Annotation | Feature category (promoter-TSS, exon, intron, Intergenic) |
| 10 | Distance to TSS | Signed distance (negative = upstream) |
| 16 | Gene Name | Gene symbol |

```bash
awk -F'\t' 'NR>1 {print $2, $3, $4, $16, $10, $8}' OFS='\t' annotated.txt > summary.tsv
```

### HOMER Category Mapping

| HOMER Category | Collapsed |
|---------------|-----------|
| promoter-TSS | promoter |
| 5' UTR, 3' UTR, exon, non-coding | exon |
| intron | intron |
| Intergenic, TTS | intergenic |

**Limitation:** HOMER defines promoter as -1kb to +100bp from TSS; this window is not configurable via flags. For custom promoter windows, reclassify using the Distance to TSS column:

```bash
awk -F'\t' 'NR>1 {
    dist = ($10 < 0) ? -$10 : $10
    feat = (dist <= 2000) ? "promoter" : $8
    print $2, $3, $4, $16, $10, feat
}' OFS='\t' annotated.txt > reclassified.tsv
```

## Annotation with Python

### Parse GTF and Extract Gene Models

**Goal:** Build gene, exon, and TSS tables from a GTF file for peak annotation.

**Approach:** Parse GTF line by line, extract gene and exon features, compute TSS positions from strand and coordinates, converting from 1-based GTF to 0-based BED coordinates.

```python
import gzip, pandas as pd

def parse_gtf(gtf_path):
    records = []
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            attrs = {}
            for item in fields[8].strip().rstrip(';').split(';'):
                item = item.strip()
                if ' ' in item:
                    key, val = item.split(' ', 1)
                    attrs[key] = val.strip('"')
            records.append({'chrom': fields[0], 'feature': fields[2],
                            'start': int(fields[3]) - 1, 'end': int(fields[4]),
                            'strand': fields[6], **attrs})
    return pd.DataFrame(records)

gtf = parse_gtf('genes.gtf.gz')
genes = gtf[gtf['feature'] == 'gene'].copy()
genes['tss'] = genes.apply(lambda r: r['start'] if r['strand'] == '+' else r['end'], axis=1)
exons = gtf[gtf['feature'] == 'exon']
```

### Find Nearest Gene and Classify Features

**Goal:** For each peak, find the nearest gene by TSS distance, compute signed distance, and assign a feature category with proper priority.

**Approach:** Calculate absolute TSS distance for all genes on the same chromosome, select the nearest, compute strand-aware signed distance, then classify using promoter window -> exon overlap -> gene body overlap -> intergenic, in priority order.

```python
peaks = pd.read_csv('peaks.bed', sep='\t', header=None,
                     names=['chr', 'start', 'end', 'name', 'score'])
peaks['center'] = (peaks['start'] + peaks['end']) // 2
promoter_window = 2000  # bp from TSS; match to analysis requirements

results = []
for _, peak in peaks.iterrows():
    chrom_genes = genes[genes['chrom'] == peak['chr']]
    abs_dists = (chrom_genes['tss'] - peak['center']).abs()
    nearest = chrom_genes.loc[abs_dists.idxmin()]

    raw_dist = peak['center'] - nearest['tss']
    signed_dist = -raw_dist if nearest['strand'] == '-' else raw_dist

    # Feature classification: promoter > exon > intron > intergenic
    if abs(raw_dist) <= promoter_window:
        feature = 'promoter'
    else:
        chrom_exons = exons[exons['chrom'] == peak['chr']]
        in_exon = ((chrom_exons['start'] <= peak['center']) & (peak['center'] < chrom_exons['end'])).any()
        in_gene = ((chrom_genes['start'] <= peak['center']) & (peak['center'] < chrom_genes['end'])).any()
        feature = 'exon' if in_exon else ('intron' if in_gene else 'intergenic')

    results.append({'chr': peak['chr'], 'start': peak['start'], 'end': peak['end'],
                    'nearest_gene': nearest['gene_name'], 'distance_to_tss': int(signed_dist),
                    'feature': feature})

result_df = pd.DataFrame(results)
result_df.to_csv('annotations.tsv', sep='\t', index=False)
```

### Alternative: pyranges

For projects with pyranges installed, GTF parsing is simpler:

```python
import pyranges as pr

gtf = pr.read_gtf('genes.gtf.gz')  # auto-converts to 0-based half-open
genes = gtf[gtf.Feature == 'gene']
peaks = pr.read_bed('peaks.bed')
nearest = peaks.nearest(genes)  # adds Distance column
```

Feature classification and signed distance still require manual logic on the result DataFrame.

## Visualization

```r
plotAnnoPie(peak_anno)
plotAnnoBar(peak_anno)
plotDistToTSS(peak_anno, title = 'Distribution of peaks relative to TSS')
```

### Compare Multiple Peak Sets

```r
peak_files <- list(H3K4me3 = 'h3k4me3_peaks.bed', H3K27ac = 'h3k27ac_peaks.bed')
anno_list <- lapply(peak_files, function(f) annotatePeak(readPeakFile(f), TxDb = txdb))
plotAnnoBar(anno_list)
plotDistToTSS(anno_list)
```

## Key Parameters

| Parameter (ChIPseeker) | Default | Description |
|------------------------|---------|-------------|
| tssRegion | c(-3000, 3000) | Promoter window around TSS |
| level | "transcript" | "transcript" or "gene"; gene-level merges all isoforms |
| genomicAnnotationPriority | Promoter > ... > Intergenic | Feature priority for overlapping annotations |
| overlap | "TSS" | "TSS" (nearest by TSS) or "all" (any gene body overlap) |
| sameStrand | FALSE | Restrict to same-strand genes only |
| addFlankGeneInfo | FALSE | Include neighboring gene information |

## Related Skills

- peak-calling - Generate peak files with MACS3 or HOMER
- motif-analysis - De novo and known motif enrichment in peak regions
- differential-binding - Compare peaks between conditions
- chipseq-visualization - Signal tracks, heatmaps, profile plots
- genome-intervals/gtf-gff-handling - Parse and convert GTF/GFF annotation files
- genome-intervals/proximity-operations - bedtools closest and window operations
- pathway-analysis/go-enrichment - Functional enrichment of peak-associated genes
