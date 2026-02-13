# Reference: CEMiTool 1.26+ | Verify API if version differs
# Automated co-expression analysis with CEMiTool

library(CEMiTool)

expr_matrix <- read.csv('normalized_counts.csv', row.names = 1)
sample_annot <- data.frame(
    SampleName = colnames(expr_matrix),
    Class = c(rep('control', 15), rep('treated', 15))
)

cem <- cemitool(expr_matrix, sample_annot, filter = TRUE, plot = TRUE, verbose = TRUE)
cat('Found', nmodules(cem), 'modules\n')

mod_genes <- module_genes(cem)
for (mod in names(mod_genes)) {
    cat(mod, ':', nrow(mod_genes[[mod]]), 'genes\n')
}

# Optional: ORA enrichment with gene sets
# gene_sets <- read_gmt('pathways.gmt')
# cem <- mod_ora(cem, gene_sets)
# cem <- plot_ora(cem)

generate_report(cem, directory = 'cemitool_report')
cat('Report saved to cemitool_report/\n')

save_plots(cem, 'cemitool_plots', force = TRUE)
cat('Plots saved to cemitool_plots/\n')
