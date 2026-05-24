# Reference: ComplexHeatmap 2.18+, maftools 2.18+, circlize 0.4.16+ | Verify API if version differs

# PhD-level OncoPrint encoding the four correctness traps:
# (1) ;-separated cells preserve multi-class stacking, (2) default memoSort retained,
# (3) remove_empty_columns=FALSE preserves cohort N, (4) log10 TMB to avoid hypermutator saturation.

library(ComplexHeatmap)
library(circlize)
library(maftools)

# 1. INPUT -- MAF + CNV calls; convert each (gene, sample) to a ;-separated class string
# Result: gene-by-sample matrix `mat`, where cells look like 'Missense;Amp' or '' (empty)

# 2. ALTERATION-CLASS COLOR -- CVD-safe distinguishable palette
col <- c(Missense   = '#56B4E9',
         Truncating = '#000000',
         Splice     = '#CC79A7',
         Amp        = '#D55E00',
         HomDel     = '#0072B2',
         Fusion     = '#009E73')

# 3. alter_fun -- one renderer per class; order = bottom-to-top stacking
alter_fun <- list(
    background = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h - unit(0.5, 'mm'),
                  gp = gpar(fill = '#EEEEEE', col = NA)),
    Amp = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h - unit(0.5, 'mm'),
                  gp = gpar(fill = col['Amp'], col = NA)),
    HomDel = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h - unit(0.5, 'mm'),
                  gp = gpar(fill = col['HomDel'], col = NA)),
    Missense = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h * 0.5,
                  gp = gpar(fill = col['Missense'], col = NA)),
    Truncating = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h * 0.33,
                  gp = gpar(fill = col['Truncating'], col = NA)),
    Splice = function(x, y, w, h)
        grid.rect(x, y, w - unit(0.5, 'mm'), h * 0.25,
                  gp = gpar(fill = col['Splice'], col = NA)),
    Fusion = function(x, y, w, h)
        grid.points(x, y, pch = 17, size = unit(2, 'mm'),
                    gp = gpar(col = col['Fusion'])))

# 4. CLINICAL ANNOTATION -- top track; TMB log-transformed for hypermutator robustness
ha_top <- HeatmapAnnotation(
    TMB     = anno_barplot(log10(clinical$tmb + 1),
                            gp = gpar(fill = '#56B4E9', col = NA),
                            axis_param = list(at = log10(c(1, 10, 100, 1000) + 1),
                                              labels = c('1', '10', '100', '1000'))),
    Subtype = clinical$subtype,
    Stage   = clinical$stage,
    col = list(Subtype = c(Luminal='#0072B2', Basal='#D55E00', HER2='#009E73'),
               Stage   = c(I='#FFFFCC', II='#FED976', III='#FD8D3C', IV='#BD0026')),
    annotation_name_gp = gpar(fontsize = 8),
    show_legend = TRUE)

# 5. RENDER -- remove_empty_columns=FALSE preserves cohort denominators
pdf('oncoprint.pdf', width = 12, height = 7)
draw(oncoPrint(mat,
               alter_fun = alter_fun,
               col = col,
               top_annotation = ha_top,
               column_title = 'TCGA-BRCA mutation landscape (N=600)',
               row_names_gp = gpar(fontsize = 8),
               pct_gp = gpar(fontsize = 7),
               show_pct = TRUE,
               remove_empty_columns = FALSE,                  # CRITICAL: cohort denominator
               remove_empty_rows = FALSE,
               heatmap_legend_param = list(
                   title = 'Alteration',
                   at = names(col),
                   labels = names(col))),
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()

# 6. MUTUAL EXCLUSIVITY / CO-OCCURRENCE -- maftools (Fisher) or DISCOVER (rate-aware)
# Use DISCOVER for pan-cancer; Fisher is fine for homogeneous cohort
maf <- read.maf(maf = 'cohort.maf', clinicalData = clinical)
si <- somaticInteractions(maf = maf, top = 20,
                          pvalue = c(0.05, 0.01),
                          fontSize = 0.7)
# si returns a matrix of -log10(p) with sign = direction (+ co-occur, - mutex)

# 7. SUBTYPE-SPLIT layout
# pass column_split = clinical$subtype to oncoPrint() for per-subtype panels

# 8. SAMPLE-ORDER OVERRIDE warning
# Default memoSort produces the canonical staircase; do NOT override column_order
# unless intentional. Document the sort criterion in the caption.
