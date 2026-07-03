---
name: bio-epitranscriptomics-modification-visualization
description: Visualises RNA-modification data with transcript-feature metagene plots (Guitar GuitarPlot; MetaPlotR; deepTools computeMatrix scale-regions), peak-centred heatmaps (ComplexHeatmap; deepTools plotHeatmap), IP-vs-input paired browser tracks (log2 IP/input bigWig via deepTools bamCompare; pyGenomeTracks; Gviz; IGV/UCSC track hubs), DRACH sequence-logo plots (ggseqlogo; MEME), feature-distribution stacked bars, and volcano/MA plots for differential modification. Establishes stop-codon enrichment in the metagene plot as the biological QC anchor for any MeRIP dataset (Dominissini 2012; Meyer 2012). Use when producing the canonical metagene plot with stop-codon enrichment as a QC anchor, building paired IP/input genome-browser tracks at single-locus resolution, plotting peak-centred heatmaps clustered by condition, summarising peak distribution across transcript features, generating DRACH motif logos as sanity checks, rendering volcano plots of differential m6A, or reproducing the stop-codon enrichment plot.
tool_type: mixed
primary_tool: Guitar
---

## Version Compatibility

Reference examples tested with: Guitar 2.18+ (Bioconductor), MetaPlotR (GitHub, unversioned), deepTools 3.5+, ggcoverage 1.4+, pyGenomeTracks 3.8+, Gviz 1.46+, ComplexHeatmap 2.18+, ggseqlogo 0.2+, ggplot2 3.5+, rtracklayer 1.62+, GenomicFeatures 1.54+, BSgenome.Hsapiens.UCSC.hg38 1.4+, TxDb.Hsapiens.UCSC.hg38.knownGene 3.18+, ChIPseeker 1.38+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('Guitar')` then `?GuitarPlot` to verify parameters
- CLI: `<tool> --help`; `deepTools <tool> --help`

If R throws `unused argument` or `argument is missing`, the API moved between Bioconductor minor releases. Guitar's older API used `txdb=` while newer uses `txTxdb=`; verify with `?GuitarPlot`. deepTools `bamCompare --operation log2` is the modern syntax (older `--ratio log2` is being phased out).

# RNA Modification Visualisation

**"Make the canonical m6A metagene plot for my paper"** -> Render the 5'UTR / CDS / 3'UTR transcript-feature metagene with stop-codon enrichment for visual confirmation of the canonical m6A topology (THE smoke test that the antibody-IP captured real m6A), the peak-feature-distribution stacked bar ("where do my peaks land?" Figure 1), peak-centred heatmaps clustered by condition, IP-over-Input paired browser tracks at specific loci, the DRACH sequence logo as a sanity check on antibody specificity, and the volcano / MA plots for differential modification. CRITICAL: the stop-codon enrichment plot is a biological-QC anchor — a MeRIP dataset that does NOT show enrichment at and around the stop codon indicates IP failure, wrong antibody, wrong protocol, or the assay captured a different modification (e.g., m1A is centred at TSS, not stop). Build this plot FIRST as a smoke test BEFORE any downstream visualisation.

- R: `Guitar::GuitarPlot(peakBedFiles, txTxdb=...)` -- canonical transcript-feature metagene (THE field-standard plot)
- CLI: `deeptools computeMatrix scale-regions ...` + `plotProfile` -- genome-coordinate metagene
- R: `ComplexHeatmap::Heatmap()` on peak-centred signal matrix -- peak-centred heatmap clustered by condition
- CLI: `pyGenomeTracks --tracks tracks.ini --region chr:start-end` -- multi-track browser plot
- R: `ggcoverage::ggcoverage(data=track, mark.region=peaks)` + `geom_gene()` -- ggplot2-based browser tracks
- R: `ggseqlogo::ggseqlogo(seqs)` -- DRACH sequence logo from peak-centre 5-mers
- R: `ChIPseeker::annotatePeak()` + `ggplot2` stacked bar -- 5'UTR / CDS / 3'UTR feature distribution

## The Single Most Important Modern Insight -- Stop-codon enrichment in the metagene plot is the biological QC anchor — a MeRIP dataset without it is suspect

Dominissini 2012 *Nature* 485:201 and Meyer 2012 *Cell* 149:1635 — concurrent papers from different labs using different protocols (MeRIP-seq and m6A-seq, respectively) — both independently showed m6A enrichment at and around the stop codon (3'UTR-proximal end of the CDS). This is the most robust biological signal in MeRIP-seq, reproduced across cell types, conditions, and decades. A metagene plot from a MeRIP library that does NOT show stop-codon enrichment indicates: (1) IP failure (antibody did not bind), (2) wrong antibody, (3) wrong protocol (e.g., the assay actually captured a different modification — m1A is centred at TSS, not stop), or (4) sample-RNA degradation. Build the Guitar metagene plot FIRST after peak calling, BEFORE any downstream visualisation, as a smoke test. The corollary: the 5'UTR / CDS / 3'UTR feature-distribution stacked bar (the "where do m6A peaks land?" Figure 1) should always be PAIRED with the metagene plot, because the stacked bar can look right (peaks land in 3'UTR / stop area) while the metagene is off (peak DENSITY not concentrated at the codon itself), or vice versa. The two plots are complementary biology-QC anchors, not redundant.

## Algorithmic Taxonomy

| Tool / plot | Mechanism | Inputs | Output | Strength | Fails when |
|-------------|-----------|--------|--------|----------|------------|
| Guitar GuitarPlot (Cui 2016 *Biomed Res Int* 2016:8367534) | Per-peak distance computation relative to TSS / start / stop / TES; feature-scaled rendering | BED + TxDb | PDF + per-feature density | THE field-standard m6A metagene; transcript-feature-aware | Older / newer API uses `txdb=` vs `txTxdb=` argument |
| MetaPlotR (Olarerin-George 2017 *Bioinformatics* 33:1563) | Alternative metagene approach with different segment-rescaling | BED + GTF | PDF | Alternative segment rescaling philosophy | GitHub-only; less actively maintained |
| deepTools computeMatrix + plotProfile (Ramírez 2016 *NAR* 44:W160) | Genome-coordinate signal aggregation over scaled gene regions | bigWig + BED | PDF + matrix | Fast; flexible region/anchor choice; well-tested ChIP-seq lineage | Not transcript-feature-aware; misses 5'UTR / CDS / 3'UTR distinction |
| deepTools plotHeatmap | Same matrix as plotProfile; rendered as heatmap with row clustering | bigWig matrix | PDF | Standard peak-centred heatmap | Single-condition view (use ComplexHeatmap for multi-condition clustering) |
| ComplexHeatmap (Gu 2016 *Bioinformatics* 32:2847) | General-purpose heatmap with multi-dimensional clustering | numeric matrix | PDF / interactive | Multi-condition cluster + annotation; publication-quality | Requires pre-computed signal matrix |
| pyGenomeTracks (Lopez-Delisle 2021 *Bioinformatics* 37:422) | Config-file (INI) driven multi-track browser plot | bigWig / bed / GTF + INI config | PDF / PNG | Reproducible browser figures via config; multi-track stacking | INI config syntax error-prone for new users |
| ggcoverage (Song & Wang 2023 *BMC Bioinformatics* 24:309) | ggplot2-native track plotting | bigWig / BAM + GTF | ggplot2 object | ggplot2 syntax; combinable with annotation layers | Newer tool; smaller user base than Gviz |
| Gviz (Hahne & Ivanek 2016 *Methods Mol Biol* 1418:335) | R/Bioconductor general genome track | bigWig / BAM / GRanges + TxDb | PDF / R plot | Most flexible R-native; long Bioconductor history | Heavier than ggcoverage; steeper learning curve |
| IGV / IGV.js (Robinson 2011 *Nat Biotechnol* 29:24) | Interactive browser via Java / JS | bigWig / BAM / bed | interactive | Standard for ad-hoc inspection | Not reproducible for figure generation |
| UCSC Track Hubs | UCSC genome browser display | bigWig + hub.txt + genomes.txt | URL | Public-display standard | Setup overhead for short projects |
| ggseqlogo (Wagih 2017 *Bioinformatics* 33:3645) | Sequence-logo rendering in ggplot2 | character vector of equal-length sequences | ggplot2 object | Native ggplot2; method='bits' or 'probability' | Requires equal-length sequences |
| MEME-ChIP (Machanick & Bailey 2011 *Bioinformatics* 27:1696) | Motif discovery + logo | BED + genome | HTML + PWM | Comprehensive motif suite; gold standard | Heavier than ggseqlogo for simple visualisation |
| HOMER findMotifsGenome.pl (Heinz 2010 *Mol Cell* 38:576) | Motif discovery via cumulative hypergeometric on shuffled background | BED + genome | HTML + motif files | Battle-tested; RNA mode via `-rna` | Less flexible output formatting than MEME |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| Canonical Figure 1 metagene plot for m6A paper | Guitar GuitarPlot with TxDb -- transcript-feature scaled | deepTools computeMatrix is genome-coordinate only; misses 5'UTR / CDS / 3'UTR distinction |
| QC: does my MeRIP show stop-codon enrichment? | Guitar GuitarPlot FIRST -- the canonical biological-QC anchor | Skipping this and rushing to peak counting is the most common QC failure |
| Browser-track Figure for a specific locus | pyGenomeTracks (config-file driven; reproducible) OR ggcoverage (ggplot2-native) | IGV is interactive but not reproducible for figures |
| Multi-condition peak-centred heatmap with clustering | ComplexHeatmap; row-cluster k-means; annotate by condition | deepTools plotHeatmap is single-condition view |
| DRACH sequence logo from peak centres | ggseqlogo on resized peak ranges (width=5, fix='center') | MEME-ChIP heavier than needed for visualisation only |
| Differential m6A volcano plot | ggplot2 native (see m6a-differential) -- modification-visualization should NOT duplicate | Duplicating volcano recipe across skills inflates content |
| Genome-coordinate metagene over custom features | deepTools computeMatrix scale-regions; plotProfile | Guitar restricted to transcript features; doesn't scale to arbitrary BED |
| Cross-sample paired IP/input track display | bamCompare -> log2 bigWig; pyGenomeTracks multi-track | Single bigWig hides the IP/input pairing structure |
| 5'UTR / CDS / 3'UTR stacked bar | ChIPseeker annotatePeak -> ggplot2 stacked bar | Manual annotation error-prone; ChIPseeker handles edge cases |
| Reproducing Dominissini 2012 / Meyer 2012 stop-codon plot | Guitar GuitarPlot; the canonical recipe | Other tools produce similar but not identical plots |

Methodology evolves; before high-stakes figure generation, web-search "Guitar Bioconductor release notes" and "pyGenomeTracks vs ggcoverage" for current best practice.

## Guitar Transcript-Feature Metagene (THE Canonical Plot)

**Goal:** Render the canonical m6A metagene plot showing peak density along scaled transcript features (5'UTR, CDS, 3'UTR), with the stop-codon-proximal enrichment that is the field's biological QC anchor.

**Approach:** Load called peaks from BED; pass to GuitarPlot with the matched TxDb; export per-feature density and the PDF.

```r
library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Older Guitar versions: argument is `txdb=`; newer: `txTxdb=`. Verify with ?GuitarPlot.
GuitarPlot(
    txTxdb           = txdb,
    stBedFiles       = list('exomepeak2_output/m6a_run1/peaks.bed'),
    miscOutFilePrefix = 'figures/m6a_metagene',
    enableCI         = FALSE,
    saveToPDFprefix  = 'figures/m6a_metagene'
)
```

If `txTxdb=` is rejected, fall back to `txdb=` and consult `?GuitarPlot`. The expected output: stop-codon-proximal enrichment, with peak density rising toward and peaking near the stop codon in the 3'UTR. If the plot does NOT show this pattern, do not proceed to downstream visualisation — investigate IP / antibody / protocol failure in merip-preprocessing.

## deepTools Genome-Coordinate Metagene

**Goal:** Render a genome-coordinate metagene over scaled gene regions using deepTools, useful when transcript-feature scaling is not needed or when the regions of interest are not standard genes (custom BED).

**Approach:** Build a signal matrix with `computeMatrix scale-regions` from a bigWig (typically IP-over-Input log2 from merip-preprocessing); render as profile plot + heatmap.

```bash
mkdir -p figures

computeMatrix scale-regions \
    --regionsFileName refs/protein_coding.bed \
    --scoreFileName tracks/IP_rep1_over_Input_rep1.bw \
    --regionBodyLength 2000 \
    --upstream 500 \
    --downstream 500 \
    --skipZeros \
    --numberOfProcessors 8 \
    --outFileName figures/m6a_matrix.gz

plotProfile \
    --matrixFile figures/m6a_matrix.gz \
    --outFileName figures/m6a_profile.pdf \
    --plotTitle 'm6A IP over Input metagene' \
    --plotType lines

plotHeatmap \
    --matrixFile figures/m6a_matrix.gz \
    --outFileName figures/m6a_heatmap.pdf \
    --colorMap RdBu_r \
    --plotTitle 'm6A IP over Input heatmap'
```

`scale-regions` is the right mode for gene-body metagene (5'-end-to-3'-end scaled to common length); `reference-point` is for peak-centred plots. For peak-centred, see the next section.

## Peak-Centred Heatmap

**Goal:** Render a peak-centred heatmap of IP-over-Input signal at +/-window around each peak, clustered by condition or by signal pattern, for cross-condition comparison.

**Approach:** Use deepTools `computeMatrix reference-point --referencePoint center` for the matrix; ComplexHeatmap for the clustered render.

```bash
computeMatrix reference-point \
    --regionsFileName exomepeak2_output/m6a_run1/peaks.bed \
    --scoreFileName tracks/IP_rep1_over_Input_rep1.bw tracks/IP_rep2_over_Input_rep2.bw tracks/IP_rep3_over_Input_rep3.bw \
    --referencePoint center \
    --upstream 500 \
    --downstream 500 \
    --binSize 25 \
    --skipZeros \
    --numberOfProcessors 8 \
    --outFileName figures/peak_centred_matrix.gz \
    --outFileNameMatrix figures/peak_centred_matrix.tab

plotHeatmap \
    --matrixFile figures/peak_centred_matrix.gz \
    --outFileName figures/peak_centred_heatmap.pdf \
    --kmeans 3 \
    --colorMap viridis \
    --plotTitle 'm6A signal centred at peaks'
```

For multi-condition heatmap with annotations, parse the matrix into R and use ComplexHeatmap:

```r
library(ComplexHeatmap)
library(circlize)

# deepTools --outFileNameMatrix has a 3-line JSON-style header before per-bin numeric columns.
raw <- read.delim('figures/peak_centred_matrix.tab', skip=3, header=FALSE)
mat <- as.matrix(raw[, -(1:6)])

col_fun <- colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red'))

Heatmap(
    mat,
    name              = 'log2 IP / Input',
    col               = col_fun,
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    row_km            = 3,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    column_title      = 'Peak-centred (+/- 500 bp)'
)
```

## IP-vs-Input Paired Browser Tracks via pyGenomeTracks

**Goal:** Render publication-quality genome-browser tracks at a specific locus showing paired IP / Input bigWig tracks, the IP/input log2 ratio, peak calls, and gene annotation; reproducible via INI config.

**Approach:** Build a `tracks.ini` config file listing each track type (bigwig, bed, gtf); invoke `pyGenomeTracks --tracks tracks.ini --region chr:start-end`.

```ini
[x-axis]
where = top
fontsize = 12

[IP rep1]
file = tracks/IP_rep1.bw
title = IP rep1
color = #d62728
height = 2
min_value = 0

[Input rep1]
file = tracks/Input_rep1.bw
title = Input rep1
color = #1f77b4
height = 2
min_value = 0

[IP / Input log2]
file = tracks/IP_rep1_over_Input_rep1.bw
title = log2 IP / Input
color = #2ca02c
height = 2

[spacer]

[m6A peaks]
file = exomepeak2_output/m6a_run1/peaks.bed
title = m6A peaks (exomePeak2)
color = #ff7f0e
height = 1
labels = false

[genes]
file = refs/annotation.gtf
title = GENCODE genes
height = 2
prefered_name = gene_name
merge_transcripts = true
```

```bash
pyGenomeTracks \
    --tracks tracks.ini \
    --region chr19:54,792,000-54,799,000 \
    --outFileName figures/browser_metti3_locus.pdf \
    --width 14 \
    --plotWidth 12
```

## ggcoverage ggplot2-Native Browser Track

**Goal:** Render genome-browser tracks in ggplot2 syntax for combining with other ggplot2 layers (annotations, peak highlights, custom theming).

**Approach:** `LoadTrackFile()` parses a bigWig / bigBed / BAM input into the dataframe ggcoverage expects; chain `ggcoverage()` with `geom_gene()` for transcript annotation; `mark.region` requires columns `start`, `end`, and `label`.

```r
library(ggcoverage)
library(rtracklayer)

peaks <- as.data.frame(import('exomepeak2_output/m6a_run1/peaks.bed'))

track.df <- LoadTrackFile(
    track.file = 'tracks/IP_rep1_over_Input_rep1.bw',
    format     = 'bw',
    region     = 'chr19:54792000-54799000'
)

mark.df <- data.frame(
    start = peaks$start,
    end   = peaks$end,
    label = peaks$name
)

ggcoverage(
    data        = track.df,
    color       = 'auto',
    mark.region = mark.df
) +
    geom_gene(gtf.file = 'refs/annotation.gtf') +
    ggplot2::theme_classic()
```

## DRACH Sequence Logo

**Goal:** Render a sequence logo of peak-centre 5-mers as a sanity check that the called peak set is enriched for the DRACH consensus motif.

**Approach:** Resize peak GRanges to fixed 5-nt centred windows; extract genomic sequences; pass to ggseqlogo with `method='probability'`.

```r
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggseqlogo)
library(rtracklayer)

peaks <- import('exomepeak2_output/m6a_run1/peaks.bed')

peak_centres <- resize(peaks, width=5, fix='center')

genome <- BSgenome.Hsapiens.UCSC.hg38
seqs <- as.character(getSeq(genome, peak_centres))

ggseqlogo(seqs, method='probability') +
    ggplot2::labs(title='Peak-centre 5-mer (DRACH consensus expected)',
                  subtitle='Sanity check on antibody specificity — NOT a per-peak filter')
```

The expected output: a logo showing approximately D-R-A-C-H consensus (D=A/G/U, R=A/G, A=methylated, C, H=A/C/U) with A clearly dominating position 3. If the logo does NOT show DRACH-like enrichment, the IP failed OR the wrong antibody was used.

## 5'UTR / CDS / 3'UTR Stacked Bar

**Goal:** Render the Figure 1 "where do my peaks land?" stacked bar showing the fraction of peaks in 5'UTR vs CDS vs 3'UTR vs intron, paired with the metagene to confirm canonical m6A topology.

**Approach:** Use ChIPseeker `annotatePeak()` with a matched TxDb; aggregate to per-feature counts; render as ggplot2 stacked bar.

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(ggplot2)
library(dplyr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- import('exomepeak2_output/m6a_run1/peaks.bed')

peak_anno <- annotatePeak(peaks, TxDb=txdb, level='transcript', verbose=FALSE)
anno_df <- as.data.frame(peak_anno@anno)

anno_df$feature <- gsub(' \\(.*\\)', '', anno_df$annotation)

feature_counts <- anno_df %>%
    group_by(feature) %>%
    summarise(n=n()) %>%
    mutate(fraction=n/sum(n)) %>%
    arrange(desc(fraction))

ggplot(feature_counts, aes(x='m6A peaks', y=fraction, fill=feature)) +
    geom_col(width=0.5) +
    scale_y_continuous(labels=scales::percent_format()) +
    labs(x=NULL, y='Fraction of peaks', title='m6A peak distribution across transcript features') +
    theme_minimal()
```

## Per-Method Failure Modes

### Metagene without stop-codon enrichment

**Trigger:** Guitar GuitarPlot of called peaks does NOT show enrichment at and around the stop codon.

**Mechanism:** Canonical m6A topology (Dominissini 2012 / Meyer 2012) shows stop-codon-proximal enrichment. Absence indicates (1) IP failure, (2) wrong antibody, (3) wrong protocol, (4) sample-RNA degradation, OR (5) the assay captured a different modification (m1A is TSS-centred; m5C distribution differs).

**Symptom:** Metagene is flat OR peaks at TSS (m1A signature) OR peaks in introns (unusual).

**Fix:** Do NOT proceed to downstream analysis. Diagnose at the merip-preprocessing layer: plotFingerprint, per-transcript IP/input distribution, antibody-lot QC. Re-do the IP if necessary. The metagene plot is the biological QC anchor; without it, all downstream interpretation is suspect.

### DRACH logo not enriched

**Trigger:** ggseqlogo of peak-centre 5-mers shows no consensus, OR shows a non-DRACH-like motif.

**Mechanism:** Antibody specificity failure OR wrong protocol. Anti-m6A antibodies have ~70% DRACH-context enrichment on real m6A peaks; a failed IP captures random sequences with no consensus.

**Fix:** Re-inspect IP enrichment in merip-preprocessing. If IP is clean but DRACH logo is absent, the assay may have captured a different modification — investigate before claiming m6A.

### Guitar `txdb=` vs `txTxdb=` argument confusion

**Trigger:** `GuitarPlot(stBedFiles=..., txdb=txdb)` rejected with "unused argument" error.

**Mechanism:** Guitar changed the argument name between Bioconductor releases — `txdb=` (older) vs `txTxdb=` (newer).

**Fix:** Try alternative; consult `?GuitarPlot` for installed version. Pin Guitar version in reproducible analyses.

### deepTools genome-coordinate metagene used where transcript-feature is needed

**Trigger:** `computeMatrix scale-regions` over a BED of genes used to show "5'UTR / CDS / 3'UTR enrichment".

**Mechanism:** deepTools scales by genomic length, NOT by transcript-feature length. A gene with long 5'UTR and short CDS will scale 5'UTR more than CDS; the metagene loses 5'UTR / CDS / 3'UTR semantics.

**Fix:** Use Guitar (transcript-feature-aware) for 5'UTR / CDS / 3'UTR semantics. Use deepTools for genome-coordinate metagene over arbitrary BED.

### IGV-only browser figure in published paper

**Trigger:** Figure 4 of paper shows an IGV screenshot of one locus; not reproducible from code.

**Mechanism:** IGV is interactive; the figure cannot be re-generated from a config file. Reviewers cannot reproduce the figure if track files change.

**Fix:** Use pyGenomeTracks (INI config) or ggcoverage (ggplot2 script) for figures intended for publication. IGV is for ad-hoc inspection only.

### bamCompare `--ratio log2` deprecation

**Trigger:** Older deepTools syntax `--ratio log2` used; newer requires `--operation log2`.

**Fix:** Switch to `--operation log2`. Both currently work but `--ratio` is being phased out.

### ggseqlogo "sequences must be equal length"

**Trigger:** Mixed-length peak sequences passed to ggseqlogo.

**Mechanism:** ggseqlogo requires equal-length input strings.

**Fix:** Resize peak ranges to fixed width before extraction: `resize(peaks, width=5, fix='center')`.

### pyGenomeTracks INI typo

**Trigger:** Track section header missing brackets, OR `file =` instead of `file=` (whitespace inconsistency).

**Fix:** Use the documented INI syntax precisely; bracketed section headers; consistent whitespace around `=`. Validate with `pyGenomeTracks --tracks tracks.ini --region 'chr:start-end' --outFileName out.pdf` and iterate on errors.

## Reconciliation: When Plots Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Metagene shows stop-codon enrichment but stacked bar has many 5'UTR peaks | 5'UTR peaks are present but density is concentrated at stop codon | Both plots are correct; report together |
| Guitar metagene and deepTools metagene look very different | Guitar transcript-feature-scaled vs deepTools genome-scaled | Both correct; Guitar for 5'UTR / CDS / 3'UTR semantics; deepTools for genome coordinates |
| Peak-centred heatmap shows two clusters; condition annotation crosses cluster boundary | Biological signal not condition-aligned; OR clustering driven by per-peak coverage variance | Inspect per-cluster fold-change; consider removing low-coverage peaks before re-clustering |
| pyGenomeTracks shows IP track higher than Input but bamCompare bigWig shows log2 near zero | Track signal magnitudes are CPM-normalised; absolute counts and log2 ratios are different summaries | Report log2 ratio as the primary; per-track CPM as supplement |
| DRACH logo enriched but stacked bar shows mostly intronic peaks | Intronic m6A (Louloupi 2018 nascent transcripts) | Genuine biology; report intronic vs exonic separately; consider library prep (poly-A vs ribo-depleted) |
| Volcano plot symmetric but most differential peaks in 3'UTR | Differential is feature-restricted | Cross-check with per-feature differential testing; biological interpretation |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| Guitar metagene -- expected pattern | Stop-codon-proximal enrichment, 3'UTR > CDS > 5'UTR density | Dominissini 2012 *Nature* 485:201; Meyer 2012 *Cell* 149:1635 |
| deepTools metagene scaled-region length | 2000 bp body + 500 bp flanks | Convention; covers typical mammalian gene span |
| Peak-centred heatmap window | +/-500 bp around peak centre | Standard convention; covers MeRIP fragment width |
| Peak-centred heatmap k-means clusters | 3-5 typical | Cluster count informed by condition count + signal heterogeneity |
| ggseqlogo expected DRACH pattern | A dominant at position 3 (the methylated A); C dominant at position 4 | DRACH consensus |
| pyGenomeTracks figure width | 10-14 inches | Standard publication width |
| Browser-track region width | 5-50 kb | Locus-context dependent |
| 5'UTR / CDS / 3'UTR / intron+other expected distribution | ~10% / ~30% / ~50% / ~10% for canonical m6A | Per published m6A atlases |
| DRACH logo "method" parameter | 'probability' for visual; 'bits' for information-theoretic | ggseqlogo convention |
| Per-feature stacked bar minimum peaks | >=100 | Below this, fractions are noisy |
| Cross-replicate metagene divergence | Should be near-identical within condition | If replicates' metagenes diverge, failed-IP suspect in merip-preprocessing |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Guitar `txTxdb` rejected | Older version uses `txdb` | Switch argument name; consult `?GuitarPlot` |
| Guitar PDF blank | BED file empty OR chromosome mismatch with TxDb | Verify peak BED has peaks; reconcile chromosome naming |
| deepTools `computeMatrix` errors with "no regions" | BED file empty OR file path wrong | Verify BED is non-empty and path correct |
| pyGenomeTracks INI parse error | Bracket / whitespace inconsistency | Match documented INI syntax exactly |
| ggcoverage region not displayed | Region beyond bigWig coverage; OR chromosome naming mismatch | Verify bigWig contains region; reconcile chr1 vs 1 |
| ggseqlogo throws "equal length" error | Mixed-length sequences | `resize(peaks, width=5, fix='center')` before extracting |
| ComplexHeatmap clustering hangs | Very large peak set (>100k) | Subset to top N peaks; or use deepTools plotHeatmap k-means |
| bamCompare `--ratio log2` deprecation warning | Newer syntax | Switch to `--operation log2` |
| IGV screenshot not reproducible | Interactive tool | Switch to pyGenomeTracks / ggcoverage |
| ggseqlogo logo blank | Empty `seqs` input OR all sequences are gaps / Ns | Verify peak BED resolves to valid genomic sequences via `getSeq()` |
| Stacked bar shows >100% | Peak annotation overlap counted multiple times | Use `annotatePeak` with explicit hierarchy |
| ComplexHeatmap colour scale wrong | colorRamp2 breaks at outliers | Set scale based on robust quantiles (quantile(x, 0.05), 0, quantile(x, 0.95)) |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Does the MeRIP show the canonical stop-codon enrichment?" | Yes — Guitar metagene plot shows expected stop-codon-proximal enrichment; cited Dominissini 2012 / Meyer 2012 |
| "Why Guitar and not deepTools?" | Guitar is transcript-feature-scaled (5'UTR / CDS / 3'UTR semantics); deepTools is genome-coordinate; both reported when needed |
| "Is the DRACH motif enriched?" | ggseqlogo of peak-centre 5-mers; OR HOMER findMotifsGenome.pl E-value reported |
| "Is the browser figure reproducible?" | Yes — pyGenomeTracks INI config OR ggcoverage R script; not IGV screenshot |
| "How were peaks annotated to features?" | ChIPseeker annotatePeak with matched TxDb; hierarchy explicit |
| "Why these specific clusters in the heatmap?" | k-means with k=3 chosen via elbow / silhouette; clusters reflect signal heterogeneity |
| "Does the metagene differ between conditions?" | Per-condition metagenes plotted alongside; differences quantified at the feature level |
| "Why is the colour scheme red-blue?" | Standard convention for log2 ratios (red = up, blue = down); colour-blind-safe palette via viridis available |
| "Was a cross-check with published m6A-Atlas peaks done?" | Common-core overlap reported; cited m6A-Atlas v2 |
| "Are the browser track signal magnitudes comparable across samples?" | bigWig CPM-normalised in merip-preprocessing; documented |

## References

- Dominissini D, Moshitch-Moshkovitz S, Schwartz S et al (2012) Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. *Nature* 485(7397):201-206. doi:10.1038/nature11112
- Meyer KD, Saletore Y, Zumbo P, Elemento O, Mason CE, Jaffrey SR (2012) Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. *Cell* 149(7):1635-1646. doi:10.1016/j.cell.2012.05.003
- Cui X, Wei Z, Zhang L et al (2016) Guitar: an R/Bioconductor package for gene annotation guided transcriptomic analysis of RNA-related genomic features. *Biomed Res Int* 2016:8367534. doi:10.1155/2016/8367534
- Olarerin-George AO, Jaffrey SR (2017) MetaPlotR: a Perl/R pipeline for plotting metagenes of nucleotide modifications and other transcriptomic sites. *Bioinformatics* 33(10):1563-1564. doi:10.1093/bioinformatics/btx002
- Ramírez F, Ryan DP, Grüning B et al (2016) deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Res* 44(W1):W160-W165. doi:10.1093/nar/gkw257
- Gu Z, Eils R, Schlesner M (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics* 32(18):2847-2849. doi:10.1093/bioinformatics/btw313
- Lopez-Delisle L, Rabbani L, Wolff J et al (2021) pyGenomeTracks: reproducible plots for multivariate genomic datasets. *Bioinformatics* 37(3):422-423. doi:10.1093/bioinformatics/btaa692
- Song Y, Wang J (2023) ggcoverage: an R package to visualize and annotate genome coverage for various NGS data. *BMC Bioinformatics* 24(1):309. doi:10.1186/s12859-023-05438-2
- Robinson JT, Thorvaldsdóttir H, Winckler W et al (2011) Integrative genomics viewer. *Nat Biotechnol* 29(1):24-26. doi:10.1038/nbt.1754
- Wagih O (2017) ggseqlogo: a versatile R package for drawing sequence logos. *Bioinformatics* 33(22):3645-3647. doi:10.1093/bioinformatics/btx469
- Yu G, Wang LG, He QY (2015) ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization. *Bioinformatics* 31(14):2382-2383. doi:10.1093/bioinformatics/btv145
- Hahne F, Ivanek R (2016) Visualizing Genomic Data Using Gviz and Bioconductor. *Methods Mol Biol* 1418:335-351. doi:10.1007/978-1-4939-3578-9_16
- Heinz S, Benner C, Spann N et al (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. *Mol Cell* 38(4):576-589. doi:10.1016/j.molcel.2010.05.004
- Machanick P, Bailey TL (2011) MEME-ChIP: motif analysis of large DNA datasets. *Bioinformatics* 27(12):1696-1697. doi:10.1093/bioinformatics/btr189

## Related Skills

- merip-preprocessing - Generates the bigWig tracks (bamCompare log2 IP/input) used here
- m6a-peak-calling - Generates the peak BED used for metagene + stacked bar + heatmap
- m6a-differential - Differential results visualised via volcano + MA (those plots live in m6a-differential, not duplicated here)
- m6anet-analysis - Per-site DRS modification calls; visualisation analogous via metagene
- data-visualization/ggplot2-fundamentals - General ggplot2 grammar
- data-visualization/multipanel-figures - Combining metagene + heatmap + volcano into figures
- data-visualization/heatmaps-clustering - General heatmap clustering patterns
- data-visualization/volcano-and-ma-plots - General volcano / MA recipes (modification-specific volcano lives in m6a-differential)
- data-visualization/genome-tracks - General genome-track rendering (this skill adds the IP/input pairing specifics)
- data-visualization/sequence-logos - General sequence-logo plotting (this skill adds DRACH-specific context)
- chip-seq/chipseq-visualization - Closest sibling for browser-track + peak-centred heatmap patterns
- pathway-analysis/enrichment-visualization - For visualising downstream pathway results
