# data-visualization

## Overview

Publication-quality data visualization for bioinformatics at PhD/post-doc decision-grade depth. Encodes the perceptual-effectiveness theory (Cleveland-McGill 1984; Munzner 2014), CVD-safe palette selection (Wong 2011; Crameri 2020; Nuñez 2018 cividis), statistical-graphics rigor (Weissgerber 2015 bar critique; Allen 2019 raincloud), and the version-specific API gotchas that distinguish reproducible publication figures from notebook scratch.

**Tool type:** mixed | **Primary tools:** ggplot2, matplotlib, ComplexHeatmap, plotly, patchwork

## Skills

| Skill | Description |
|-------|-------------|
| ggplot2-fundamentals | Grammar of graphics; theme_classic baseline; programmatic tidy-eval aes; cairo_pdf TrueType embedding |
| matplotlib-fundamentals | Object-oriented Figure/Axes API; Type-42 fonts; constrained_layout; seaborn integration |
| color-palettes | Crameri scientific colormaps; cividis; Okabe-Ito; grayscale-monotonicity and CVD-simulation tests |
| multipanel-figures | patchwork ≥1.2.0 axes='collect'; cowplot alignment; GridSpec/subfigures; journal-spec sizing |
| heatmaps-clustering | ComplexHeatmap draw() trap; ward.D vs ward.D2 (Murtagh-Legendre 2014); OLO (Bar-Joseph 2001); robust bounds |
| dimensionality-reduction-plots | PCA loadings; UMAP/t-SNE/PHATE choice; Kobak-Berens 2019 init; Chari-Pachter 2023 critique |
| volcano-and-ma-plots | apeglm/ashr LFC shrinkage; padj-aware threshold lines; combined-rank labeling; raster for large N |
| manhattan-qq-locuszoom | Analysis-appropriate thresholds; λGC + LDSC; y-cap with indicator; LD reference matched to ancestry |
| oncoprint-mutation-matrices | ComplexHeatmap oncoPrint stacking; memoSort; somaticInteractions; DISCOVER for pan-cancer mutex |
| lollipop-protein-maps | maftools / trackViewer; HGVSp + isoform; hotspot validation against TCGA + ICGC |
| sequence-logos | ggseqlogo / Logomaker / WebLogo; bits vs probability; background-composition correction |
| forest-funnel-plots | metafor REML random-effects; I² / τ² / prediction interval; contour-enhanced funnel; Egger / trim-and-fill |
| distribution-plots | Boxplot / violin / beeswarm / raincloud; N-based decision tree; Sheather-Jones bandwidth |
| statistical-annotation | ggpubr / ggsignif / statannotations; test selection; Holm adjustment; effect-size reporting |
| flow-and-transition-plots | Sankey vs alluvial; ggalluvial; CONSORT trial-flow; pyCirclize chord diagrams |
| upset-plots | ComplexUpset (active) vs UpSetR (unmaintained); cardinality vs degree sort; attribute panels |
| network-visualization | NetworkX / PyVis / Cytoscape; layout-as-artifact warning; hive plots; edge bundling; ForceAtlas2 |
| circos-plots | circlize.clear() trap; chromosome.index ordering; when NOT to use circular (Cleveland-McGill 1984) |
| interactive-visualization | plotly / bokeh / gganimate; Kaleido v1 static-export pipeline (orca EOL); always export static alongside |
| genome-tracks | pyGenomeTracks / Gviz / IGV batch; bamCoverage --normalizeUsing UNDOES spike-in; shared y-axis |

## Example Prompts

- "Build a volcano plot from DESeq2 results with apeglm shrinkage, padj threshold, top-N by combined rank"
- "Cluster heatmap of top 500 variable genes with ward.D2 + Optimal Leaf Ordering + Crameri vik palette"
- "PCA + UMAP side-by-side for sample QC, with hyperparameters explicit"
- "Manhattan plot for trans-ancestry GWAS at 5e-9 threshold, y-cap at 50 with indicators"
- "OncoPrint of top 25 mutated genes split by molecular subtype, with log-transformed TMB bar"
- "Lollipop plot for TP53 across cohort with R175/R248/R273 labeled and cohort-comparison overlay"
- "ggseqlogo of TF binding sites using bits encoding with human genome background composition"
- "Raincloud comparing biomarker levels across 3 arms, Holm-adjusted Wilcoxon brackets"
- "Forest plot from REML random-effects meta-analysis with I², τ², and 95% prediction interval"
- "Alluvial of cell-state transitions across 3 timepoints, colored by starting cluster"
- "CONSORT 2010 diagram for trial enrollment flow"
- "ComplexUpset for 6 gene sets with cardinality sort and metadata stacks"
- "Multi-panel publication figure at 183mm (Nature double-column) with patchwork axes/guides collection"
- "pyGenomeTracks INI for ChIP-Rx with spike-in scale factor via --scaleFactor (NOT --normalizeUsing)"

## Requirements

```r
install.packages(c('ggplot2', 'patchwork', 'cowplot', 'scales', 'ggrepel', 'ggtext', 'viridis',
                   'scico', 'khroma', 'colorspace', 'ggsci', 'ggalluvial', 'ggbeeswarm',
                   'ggdist', 'gghalves', 'lvplot', 'ggpubr', 'ggsignif', 'rstatix', 'ggrastr',
                   'metafor', 'forestplot', 'survminer', 'MendelianRandomization',
                   'ComplexUpset', 'pheatmap', 'circlize', 'seriation',
                   'qqman', 'CMplot', 'locuszoomr', 'LDlinkR',
                   'maftools', 'trackViewer', 'g3viz', 'consort', 'networkD3',
                   'igraph', 'tidygraph', 'ggraph', 'plotly', 'htmlwidgets', 'DT',
                   'ggseqlogo', 'survival', 'lme4'))
BiocManager::install(c('ComplexHeatmap', 'Gviz', 'EnhancedVolcano', 'EnsDb.Hsapiens.v86',
                       'TxDb.Hsapiens.UCSC.hg38.knownGene', 'GenomicRanges',
                       'apeglm', 'ashr', 'DESeq2', 'PCAtools'))
```

```bash
pip install matplotlib seaborn plotly bokeh kaleido cmcrameri colorcet colorspacious \
            scanpy umap-learn openTSNE phate scikit-learn \
            statannotations adjustText \
            upsetplot pyCirclize networkx pyvis py4cytoscape \
            pyGenomeTracks deepTools comut logomaker \
            ptitprince ggalluvial datashader
```

## Related Skills

- **differential-expression/de-visualization** - DESeq2 / edgeR built-in plot helpers
- **pathway-analysis/enrichment-visualization** - clusterProfiler / enrichplot functions
- **chip-seq/chipseq-visualization** - ChIP-seq-specific track + heatmap patterns
- **hi-c-analysis/hic-visualization** - Hi-C contact maps
- **copy-number/cnv-visualization** - karyoploteR + cohort CNV heatmaps
- **phylogenetics/tree-visualization** - ggtree / ETE4 phylogenetic trees
- **spatial-transcriptomics/spatial-visualization** - Tissue overlay + spatial scatter
- **single-cell/markers-annotation** - scRNA dotplot / matrix plot / stacked violin
- **clinical-biostatistics/survival-analysis** - KM curves with risk tables
- **reporting/figure-export** - DPI / format / journal-spec compliance
