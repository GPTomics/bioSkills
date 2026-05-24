---
name: bio-chipseq-peak-annotation
description: Annotates ChIP-seq peaks to genomic features, nearest genes, ENCODE candidate cis-regulatory elements (cCREs), and regulatory domains. Uses ChIPseeker (R), HOMER annotatePeaks.pl (CLI), pyranges (Python), GREAT/rGREAT (regulatory domain gene-set enrichment), ChIP-Enrich (locus-length-adjusted), ENCODE SCREEN cCRE classification (PLS/pELS/dELS/CTCF-only/DNase-H3K4me3), and ENCODE-rE2G for cell-type-specific enhancer-gene linking. Handles nearest-TSS vs host-gene ambiguity, promoter window definition, and feature priority. Use when assigning genomic context to peaks, linking enhancer peaks to target genes, classifying peaks against ENCODE cCRE registry, or running gene-set enrichment on peak-associated genes.
tool_type: mixed
primary_tool: ChIPseeker
---

## Version Compatibility

Reference examples tested with: ChIPseeker 1.38+, GenomicFeatures 1.54+, rtracklayer 1.62+, HOMER 4.11+, rGREAT 2.4+, chipenrich 2.26+, pyranges 0.0.129+, pandas 2.2+.

ENCODE cCRE registry expanded to 2.35M human and 927k mouse elements (Nature 2025; Vu Ernst expansion). SCREEN web app at screen.encodeproject.org provides browser access; ENCODE provides bed files for batch annotation.

# Peak Annotation

**"What genes and regulatory elements do my peaks correspond to?"** -> Assign each peak to a genomic feature (promoter, exon, intron, intergenic), its target gene (via nearest-TSS or host-gene), and where applicable an ENCODE cCRE class (PLS/pELS/dELS/CTCF-only/DNase-H3K4me3).

- R (gene-feature): `ChIPseeker::annotatePeak(peaks, TxDb=txdb)`
- CLI (gene-feature): `annotatePeaks.pl peaks.bed hg38 -gtf annotation.gtf`
- Python (custom): pyranges + pandas
- R (cCRE classification): intersect peaks with ENCODE cCRE BED from SCREEN
- R (gene-set enrichment): `rGREAT::great()` or `chipenrich::chipenrich()`

The single biggest source of misinterpretation is the **nearest-TSS vs host-gene** distinction (see below). For enhancer-driven biology, ENCODE-rE2G or ABC (in atac-seq/enhancer-gene-linking) is more accurate than nearest-TSS.

## Choosing an Annotation Approach

| Context | Recommended | Why |
|---------|-------------|-----|
| Standard genome, pre-built annotations available | ChIPseeker with TxDb package | Simplest; automatic gene symbol mapping via annoDb |
| Custom or project-specific GTF | ChIPseeker + makeTxDbFromGFF, HOMER -gtf, or pyranges | All three handle custom annotations |
| HOMER already in pipeline | HOMER annotatePeaks.pl | Reuses tag directory; combined with motif workflow |
| Fine-grained control | pyranges (Python) | Full control over priority rules, distance calculation |
| Enhancer peaks (distal regulatory) | GREAT / rGREAT | Regulatory domain assignment (basal + extension), not just nearest |
| Cell-type-specific enhancer-gene linking | ENCODE-rE2G | Modern (2024); ABC-trained logistic regression with chromatin context |
| Gene-set enrichment with locus-length adjustment | chipenrich / Broad-Enrich | Corrects for systematic gene-length bias in peak assignment |
| Compare against ENCODE cCRE atlas | SCREEN cCRE BED intersect | Cross-reference standard regulatory registry |
| Promoter-coverage decomposition | bedtools intersect with TSS windows | Quick stats per peak set |

**Critical:** Use the same annotation source as the alignment (UCSC knownGene TxDb with GENCODE GTF alignment causes mismatches). When a specific GTF is provided, use it directly via `makeTxDbFromGFF` rather than a mismatched pre-built TxDb package.

## Nearest-TSS vs Host-Gene Convention

Peak annotation involves two decisions that should be coupled but often aren't:
1. Which gene to assign (target gene)
2. What feature the peak overlaps (promoter / exon / intron / intergenic)

Default tools decouple these, producing internally inconsistent annotations.

| Convention | Gene from | Feature from | Tools |
|------------|-----------|---------------|-------|
| Nearest-TSS (default) | Gene with closest TSS | Physical overlap at peak center | ChIPseeker `overlap='TSS'` (default), HOMER |
| Host-gene priority | Gene whose body contains the peak | Same gene's features | ChIPseeker `overlap='all'` |

**Example failure:** Peak inside gene A's intron, near gene B's TSS. Default tools report `nearest_gene=B, feature=intron` — but the intron belongs to gene A, not gene B. The annotation is internally inconsistent.

### Choosing per Biology

| Context | Convention | Rationale |
|---------|-----------|-----------|
| Distal TF binding (enhancers) | Nearest-TSS, but prefer ENCODE-rE2G / ABC | Enhancers can regulate gene A despite sitting in gene B's intron |
| Histone marks in gene bodies (H3K36me3, H3K27me3) | Host-gene | Mark reflects host transcriptional state |
| Promoter-associated marks (H3K4me3, H3K27ac at promoters) | Either | Most peaks at promoters where conventions agree |
| Custom annotation against project GTF | Host-gene | Internal consistency |
| Reproducing published HOMER results | Nearest-TSS | Matches HOMER default |

When a task says "nearest gene," clarify which definition. For most annotation purposes where gene + feature should be consistent, use host-gene; for distal enhancer biology, use a proper enhancer-gene linker (ENCODE-rE2G, ABC).

## Coordinate Systems and TSS

BED uses 0-based half-open `[start, end)`. GTF uses 1-based closed `[start, end]`. Mixing without conversion shifts annotations by one base.

**Peak center (BED):** `(start + end) // 2`

**TSS from GTF (1-based to 0-based):**
- Plus-strand: `tss_0based = start - 1`
- Minus-strand: `tss_0based = end`

**Signed distance** (negative = upstream of TSS):
- Plus-strand: `distance = peak_center - tss`
- Minus-strand: `distance = -(peak_center - tss)`

## ChIPseeker (R)

**Goal:** Assign each ChIP-seq peak to a gene and a feature category using a transcript database.

**Approach:** Load the TxDb (pre-built or custom-built from GTF), pass peaks to `annotatePeak()` with the desired `tssRegion` window and `overlap` convention (host-gene vs nearest-TSS), then export the annotated data frame with gene symbols mapped from `annoDb` or the original GTF.

**Standard genome:**

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

peaks <- readPeakFile('peaks.narrowPeak')
peak_anno <- annotatePeak(peaks,
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           tssRegion = c(-2000, 2000),
                           annoDb = 'org.Hs.eg.db',
                           overlap = 'all')   # host-gene convention
anno_df <- as.data.frame(peak_anno)
```

**Custom GTF** (use makeTxDbFromGFF; map symbols from original GTF since custom TxDb objects lack annoDb mappings):

```r
library(GenomicFeatures)
library(rtracklayer)

txdb <- makeTxDbFromGFF('genes.gtf.gz', format = 'gtf')
peaks <- readPeakFile('peaks.bed')
peak_anno <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-2000, 2000),
                           overlap = 'all')

gtf <- import('genes.gtf.gz')
gene_map <- unique(data.frame(
    gene_id = sub('\\..*', '', gtf$gene_id),
    symbol = gtf$gene_name, stringsAsFactors = FALSE))
gene_map <- gene_map[!is.na(gene_map$symbol), ]
anno_df <- as.data.frame(peak_anno)
anno_df$gene_id_base <- sub('\\..*', '', anno_df$geneId)
anno_df$SYMBOL <- gene_map$symbol[match(anno_df$gene_id_base, gene_map$gene_id)]
```

GENCODE gene IDs have version suffixes (`ENSG00000142192.25`); strip before joining.

**Promoter window:** `tssRegion = c(-2000, 2000)` is common; `c(-3000, 3000)` is ChIPseeker default. Match to analysis requirements.

**Feature priority:** Default `Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic`. A peak in both a promoter (gene A) and an intron (gene B) receives "Promoter (gene A)" by default.

## HOMER annotatePeaks.pl (CLI)

```bash
# Standard genome (HOMER's installed annotation)
annotatePeaks.pl peaks.bed hg38 > annotated.txt

# Custom GTF (overrides HOMER's default)
annotatePeaks.pl peaks.bed hg38 -gtf genes.gtf > annotated.txt

# Without installed genome, GTF only
annotatePeaks.pl peaks.bed none -gtf genes.gtf > annotated.txt

# Generate annotation statistics
annotatePeaks.pl peaks.bed hg38 -gtf genes.gtf -annStats stats.txt > annotated.txt
```

HOMER's 19-column output: columns 8 (Annotation), 10 (Distance to TSS), 16 (Gene Name) are the primary annotation columns.

**HOMER promoter window is fixed at -1kb / +100bp** — not configurable via flags. For custom windows, reclassify using the Distance to TSS column post-hoc.

## ENCODE cCRE Classification

The ENCODE Registry of candidate cis-Regulatory Elements (cCREs) provides 2.35M human + 927k mouse elements classified into 5 categories:

| Class | Definition | Marker pattern |
|-------|------------|-----------------|
| **PLS** (Promoter-Like Signature) | ≤ 200 bp of annotated TSS; high DNase + high H3K4me3 | DNase + H3K4me3 |
| **pELS** (Proximal Enhancer-Like Signature) | ≤ 2 kb of TSS; enhancer-like (DNase + H3K27ac, low H3K4me3) | DNase + H3K27ac |
| **dELS** (Distal Enhancer-Like Signature) | > 2 kb of TSS; enhancer-like | DNase + H3K27ac |
| **DNase-H3K4me3** | Promoter signature without annotated TSS | DNase + H3K4me3 (no TSS) |
| **CTCF-only** | Potential boundary; DNase + CTCF | DNase + CTCF |

```bash
# Download ENCODE cCRE BED from SCREEN (hg38)
wget https://api.wenglab.org/screen_v13/screen_human_ccres_simple.bed.gz
gunzip screen_human_ccres_simple.bed.gz

# Intersect peaks with cCRE; -wa preserves peak coords, -wb adds cCRE class
bedtools intersect -a peaks.narrowPeak -b screen_human_ccres_simple.bed -wa -wb \
    > peaks_ccre.tsv
```

Cross-referencing peaks against cCREs:
- Indicates whether peaks overlap canonical regulatory elements
- Provides the cCRE class (PLS / pELS / dELS / CTCF-only / DNase-H3K4me3)
- Cell-type-specific activity profiles available via SCREEN web app

## GREAT / rGREAT (Regulatory Domain Gene-Set Enrichment)

GREAT (McLean 2010) addresses two problems with standard gene-set enrichment on peaks:
1. Peak-to-gene assignment via regulatory domains (not nearest TSS)
2. Statistical correction for region-locus length bias

**Regulatory domain rules** (default):
- Basal domain: -5 kb / +1 kb of TSS
- Extension: up to 1 Mb in each direction, OR until reaching neighbor's basal domain
- Each peak is assigned to ALL genes whose regulatory domain it overlaps (not just nearest)

```r
library(rGREAT)

# Submit peaks for regulatory-domain gene-set enrichment
res <- great(gr = peaks, gene_sets = 'GO:BP', tss_source = 'TxDb.Hsapiens.UCSC.hg38.knownGene',
              biomart_dataset = 'hsapiens_gene_ensembl')

# Top enriched gene sets
table_results <- getEnrichmentTable(res)
head(table_results)

# Visualization
plotVolcano(res)
plotRegionGeneAssociationGraphs(res)
```

GREAT is most appropriate for distal regulatory elements (enhancer ChIP, ATAC). For promoter-focused marks (H3K4me3), ChIP-Enrich is more standard.

## ChIP-Enrich (Locus-Length-Adjusted Gene-Set Enrichment)

Welch 2014: standard gene-set enrichment on peak-associated genes systematically over-counts long genes. ChIP-Enrich models locus length as a covariate.

```r
library(chipenrich)

res <- chipenrich(peaks = 'peaks.bed', genome = 'hg38',
                   genesets = 'GOBP', locusdef = 'nearest_tss',
                   out_name = 'chipenrich_out', n_cores = 4)
# Locus definitions: nearest_tss, nearest_gene, 1kb, 5kb, 10kb, gene_body, exon
# Methods: chipenrich (default), polyenrich (for high peak counts), broadenrich (broad marks)
```

For broad marks (H3K27me3, H3K9me3): use `method = 'broadenrich'` which accounts for region width.

## ENCODE-rE2G (Modern Enhancer-Gene Linking)

ENCODE-rE2G (2024) replaces ABC for cell types with ENCODE data. Cell-type-specific logistic-regression weights map distal enhancer peaks to target genes with higher accuracy than nearest-TSS or basal+extension.

See atac-seq/enhancer-gene-linking for full workflow; the same model applies to ChIP-seq enhancer marks (H3K27ac, H3K4me1, H3K4me2).

## Per-Tool Failure Modes

### ChIPseeker -- TxDb / annoDb genome mismatch

**Trigger:** Using hg19 TxDb on hg38-aligned BAMs / peaks.

**Mechanism:** Silent; ChIPseeker doesn't verify genome assembly.

**Symptom:** Annotated gene symbols look reasonable but distance-to-TSS is wrong; promoter / intron classifications drift.

**Fix:** Match TxDb to BAM alignment genome explicitly; verify with `seqlevels(peaks) == seqlevels(txdb)`.

### ChIPseeker -- Default `overlap='TSS'` decouples gene from feature

**Trigger:** Default annotation call on peaks in gene bodies.

**Mechanism:** `overlap='TSS'` assigns nearest gene by TSS; feature classification is independent of that gene.

**Symptom:** Annotation reports `nearest_gene=X, feature=intron` where the intron belongs to a different gene.

**Fix:** Pass `overlap='all'` for host-gene-consistent annotation; or accept TSS-only convention and clarify in methods.

### ChIPseeker -- Custom TxDb has no annoDb

**Trigger:** Building TxDb from GTF and passing `annoDb='org.Hs.eg.db'`.

**Mechanism:** Custom TxDb lacks the gene_id-to-symbol mapping that org.Hs.eg.db provides; ChIPseeker silently returns NA for symbols.

**Fix:** Map symbols separately from the original GTF after annotation; strip Ensembl version suffixes before joining.

### HOMER -- Hard-coded promoter window

**Trigger:** Needing a 2 kb or 5 kb promoter window with HOMER.

**Mechanism:** HOMER's promoter classification is hard-coded to -1 kb / +100 bp; not configurable.

**Fix:** Post-hoc reclassify using `Distance to TSS` column:
```bash
awk -F'\t' 'NR>1 { dist = ($10 < 0) ? -$10 : $10; \
    feat = (dist <= 2000) ? "promoter_custom" : $8; \
    print $2, $3, $4, $16, $10, feat }' OFS='\t' annotated.txt
```

### GREAT -- Default regulatory domain inappropriate for some species / cell types

**Trigger:** Using default basal+extension on insect or compact-genome data.

**Mechanism:** 1 Mb maximum extension assumes vertebrate-scale enhancer-target distances; not appropriate for organisms with shorter regulatory ranges.

**Fix:** Adjust `extension` parameter; for non-default species, configure regulatory domain explicitly.

### GREAT -- Hyper-ChIPable peaks inflate enrichment

**Trigger:** Including unfiltered peaks at rRNA / housekeeping / mtDNA in GREAT analysis.

**Mechanism:** Hyper-ChIPable artifacts are enriched at highly-transcribed loci; GREAT assigns them to associated genes, inflating GO terms for "translation" and "ribosomal" categories.

**Symptom:** Top enriched GO terms always include "ribosomal", "translation", "mitochondrion" regardless of biology.

**Fix:** Blacklist filter + custom hyper-ChIPable filter (top-1% input signal) before GREAT.

### ENCODE cCRE -- Cell-type-agnostic vs specific

**Trigger:** Using the master cCRE BED (cell-type-agnostic) to claim cell-type-specific regulatory activity.

**Mechanism:** Master cCRE BED is the union across all cell types. Specific activity profile per cell type is a separate dataset.

**Fix:** Use SCREEN web app or per-cell-type activity profiles for cell-type-specific claims.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ChIPseeker nearest-TSS gene ≠ HOMER nearest gene | Different TSS reference; HOMER uses RefSeq | Verify both use same TxDb / RefSeq + UCSC knownGene |
| GREAT enrichment ≠ ChIP-Enrich enrichment | GREAT uses regulatory domain; ChIP-Enrich uses locus length adjustment | Both are valid; use GREAT for distal regulatory, ChIP-Enrich for promoter-focused |
| Peak overlaps cCRE but classified differently than expected | Cell-type-specific activity profile not used | Check SCREEN per-cell-type profile |
| Enhancer peak's nearest gene differs from ENCODE-rE2G target | ENCODE-rE2G uses cell-type chromatin context | Use ENCODE-rE2G for cell-type-specific enhancer-gene claims |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `seqlevels` mismatch in ChIPseeker | chr vs no-chr naming | `seqlevelsStyle(peaks) <- 'UCSC'` |
| Gene symbols all NA in ChIPseeker | Custom TxDb without annoDb | Map symbols from original GTF |
| HOMER reports "no annotation" | Genome not installed | `perl configureHomer.pl -install hg38` |
| rGREAT timeout | Large peak set + slow biomart | Use pre-computed gene sets; lower peak count |
| chipenrich slow | Default locusdef computed on-the-fly | Use built-in locusdef shortcuts (`nearest_tss`, `1kb`) |
| pyranges feature-overlap result missing strand | pyranges 0.x conversion drops strand by default | Pass `strandedness='same'` to overlap operations |

## References

- Yu G et al 2015 Bioinformatics 31:2382 (ChIPseeker)
- Heinz S et al 2010 Mol Cell 38:576 (HOMER annotatePeaks)
- McLean CY et al 2010 Nat Biotechnol 28:495 (GREAT)
- Gu Z 2023 Bioinformatics 39:btac745 (rGREAT)
- Welch RP et al 2014 Nucleic Acids Res 42:e105 (ChIP-Enrich)
- Cavalcante RG, Lee C, Welch RP, ... Sartor MA 2014 Bioinformatics 30:i393-i400 (Broad-Enrich)
- ENCODE Project Consortium 2020 Nature 583:699 (cCRE registry v1)
- ENCODE Project Consortium 2025 Nature (expanded cCRE registry, 2.35M elements)
- Fulco CP et al 2019 Nat Genet 51:1664 (ABC model precursor to ENCODE-rE2G)
- Kundaje lab / ENCODE 2024 (ENCODE-rE2G)
- SCREEN: screen.encodeproject.org

## Related Skills

- chip-seq/peak-calling - Generate peaks for annotation
- chip-seq/chipseq-qc - Filter hyper-ChIPable peaks before GREAT / chipenrich
- chip-seq/super-enhancers - Annotate SE-associated genes
- chip-seq/differential-binding - Annotate differential peaks
- atac-seq/enhancer-gene-linking - ENCODE-rE2G workflow for cell-type-specific E-G linking
- pathway-analysis/go-enrichment - Standard GO enrichment on peak-associated genes
- pathway-analysis/reactome-pathways - Reactome pathway enrichment
- genome-intervals/gtf-gff-handling - Parse and convert GTF/GFF
- genome-intervals/proximity-operations - bedtools closest and window operations
