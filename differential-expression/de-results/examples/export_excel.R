# Export DE results to Excel with multiple sheets

library(openxlsx)
library(dplyr)

# Simulate results
set.seed(42)
n_genes <- 1000
res_df <- data.frame(
    gene = paste0('gene', 1:n_genes),
    baseMean = 10^runif(n_genes, 1, 4),
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    pvalue = 10^(-runif(n_genes, 0, 5))
)
res_df$padj <- p.adjust(res_df$pvalue, method = 'BH')

# Filter results
all_results <- res_df %>% arrange(padj)
sig_genes <- res_df %>% filter(padj < 0.05) %>% arrange(padj)
up_genes <- res_df %>% filter(padj < 0.05, log2FoldChange > 0) %>% arrange(padj)
down_genes <- res_df %>% filter(padj < 0.05, log2FoldChange < 0) %>% arrange(padj)

# Create workbook
wb <- createWorkbook()

# Add sheets
addWorksheet(wb, 'All Results')
writeData(wb, 'All Results', all_results)

addWorksheet(wb, 'Significant (padj<0.05)')
writeData(wb, 'Significant (padj<0.05)', sig_genes)

addWorksheet(wb, 'Up-regulated')
writeData(wb, 'Up-regulated', up_genes)

addWorksheet(wb, 'Down-regulated')
writeData(wb, 'Down-regulated', down_genes)

# Add summary sheet
summary_df <- data.frame(
    Category = c('Total genes', 'Significant (padj < 0.05)', 'Up-regulated', 'Down-regulated'),
    Count = c(nrow(all_results), nrow(sig_genes), nrow(up_genes), nrow(down_genes))
)
addWorksheet(wb, 'Summary')
writeData(wb, 'Summary', summary_df)

# Save
saveWorkbook(wb, 'de_results.xlsx', overwrite = TRUE)
cat('Results exported to de_results.xlsx\n')
cat(sprintf('  - All results: %d genes\n', nrow(all_results)))
cat(sprintf('  - Significant: %d genes\n', nrow(sig_genes)))
cat(sprintf('  - Up-regulated: %d genes\n', nrow(up_genes)))
cat(sprintf('  - Down-regulated: %d genes\n', nrow(down_genes)))
