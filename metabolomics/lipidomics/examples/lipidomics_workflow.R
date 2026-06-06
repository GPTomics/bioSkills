# Reference: lipidr 2.16+ | Verify API if version differs
# End-to-end lipidr workflow on the bundled data_normalized dataset:
# class-based differential analysis, lipid set enrichment, and honest
# resolution-level reporting. Demonstrates the skill's central claim that
# a lipid name carries a structural-resolution level the tool may overstate.
library(lipidr)

out_dir <- tempdir()

# data_normalized ships with lipidr (PQN-normalized, log2). For real data:
#   d <- read_skyline(list.files(datadir, 'data.csv', full.names = TRUE))
#   d <- add_sample_annotation(d, 'clinical.csv')
#   d <- normalize_pqn(d, measure = 'Area', exclude = 'blank', log = TRUE)
# Class-based quantification instead of PQN uses normalize_istd (one IS per class):
#   d <- normalize_istd(d, measure = 'Area', exclude = 'blank', log = TRUE)
data(data_normalized)

cat('Lipid classes present:', paste(unique(rowData(data_normalized)$Class), collapse = ', '), '\n')

# Contrast references sample-group labels directly; group_col defaults to the first annotation column
de_results <- de_analysis(data_normalized, HighFat_water - NormalDiet_water, measure = 'Area')

# logFC.cutoff is on the limma log2 scale; 1 = two-fold, 0.05 is the conventional FDR gate
sig <- significant_molecules(de_results, p.cutoff = 0.05, logFC.cutoff = 1)
cat('Significant lipids (|log2FC| > 1, adj.P < 0.05):', length(sig), '\n')

volcano <- plot_results_volcano(de_results, show.labels = FALSE)
ggplot2::ggsave(file.path(out_dir, 'lipid_volcano.png'), volcano, width = 8, height = 6)

# Lipid set enrichment tests class, total chain length, and unsaturation sets automatically
enrich <- lsea(de_results, rank.by = 'logFC')
sig_sets <- significant_lipidsets(enrich, p.cutoff = 0.05, size.cutoff = 2)

enrich_plot <- plot_enrichment(de_results, sig_sets, annotation = 'class', measure = 'logFC')
ggplot2::ggsave(file.path(out_dir, 'lipid_class_enrichment.png'), enrich_plot, width = 8, height = 6)

# Honest resolution-level reporting: keep limma columns, never invent an sn claim the data lacks.
# de_analysis returns a tidy data.frame: Molecule, Class, total_cl, total_cs, logFC, P.Value, adj.P.Val.
results_table <- de_results[, c('Molecule', 'Class', 'total_cl', 'total_cs', 'logFC', 'P.Value', 'adj.P.Val')]
write.csv(results_table, file.path(out_dir, 'lipidomics_de_results.csv'), row.names = FALSE)
cat('Results written to', out_dir, '\n')
