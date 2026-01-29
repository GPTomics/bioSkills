library(UpSetR)

# --- ALTERNATIVE: Use real gene sets ---
# For realistic examples, try GO term gene sets or DE results:
#
# library(clusterProfiler)
# library(org.Hs.eg.db)
# ego <- enrichGO(gene_list, OrgDb = org.Hs.eg.db, ont = 'BP')
# gene_sets <- geneInCategory(ego)[1:6]
#
# Or use MSigDB gene sets:
# library(msigdbr)
# hallmark <- msigdbr(species = 'Homo sapiens', category = 'H')

# Simulated gene sets representing different analyses
# Typical scenario: comparing DE genes across conditions/timepoints
set.seed(42)

all_genes <- paste0('Gene', 1:500)

# Simulate overlapping gene sets (e.g., DE genes from different comparisons)
# Real gene sets often have ~10-30% overlap between related conditions
gene_sets <- list(
    Treatment_vs_Control = sample(all_genes, 150),
    Timepoint_6h = sample(all_genes, 120),
    Timepoint_24h = sample(all_genes, 180),
    Drug_A = sample(all_genes, 100),
    Drug_B = sample(all_genes, 90),
    Combined_Treatment = sample(all_genes, 200)
)

# Add some forced overlaps for biological realism
# Core response genes present in multiple conditions
core_genes <- sample(all_genes, 30)
for (i in 1:4) {
    gene_sets[[i]] <- unique(c(gene_sets[[i]], core_genes))
}

# Basic UpSet plot sorted by frequency
# Shows largest intersections first - most useful default
pdf('upset_basic.pdf', width = 12, height = 7)
upset(fromList(gene_sets),
      nsets = 6,
      order.by = 'freq',
      mainbar.y.label = 'Genes in Intersection',
      sets.x.label = 'Total Genes per Set')
dev.off()

# Customized with colors and adjusted ratios
# mb.ratio controls matrix-to-bar height ratio (default 0.7, 0.3)
pdf('upset_customized.pdf', width = 12, height = 8)
upset(fromList(gene_sets),
      nsets = 6,
      nintersects = 30,
      order.by = 'freq',
      decreasing = TRUE,
      mb.ratio = c(0.55, 0.45),
      point.size = 3.5,
      line.size = 1.2,
      sets.bar.color = c('#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4'),
      main.bar.color = '#7E6148',
      matrix.color = '#7E6148',
      text.scale = c(1.5, 1.2, 1.2, 1, 1.5, 1.2),
      set_size.show = TRUE)
dev.off()

# With query highlights
# Highlight specific intersections of interest
pdf('upset_queries.pdf', width = 12, height = 8)
upset(fromList(gene_sets),
      nsets = 6,
      order.by = 'freq',
      queries = list(
          # Highlight genes unique to combined treatment
          list(query = intersects,
               params = list('Combined_Treatment'),
               color = '#E64B35',
               active = TRUE,
               query.name = 'Combined only'),
          # Highlight core genes in multiple timepoints
          list(query = intersects,
               params = list('Timepoint_6h', 'Timepoint_24h'),
               color = '#4DBBD5',
               active = TRUE,
               query.name = 'Both timepoints')
      ),
      query.legend = 'bottom')
dev.off()

message('UpSet plots saved: upset_basic.pdf, upset_customized.pdf, upset_queries.pdf')

# Print intersection statistics
cat('\nSet sizes:\n')
for (name in names(gene_sets)) {
    cat(sprintf('  %s: %d genes\n', name, length(gene_sets[[name]])))
}
