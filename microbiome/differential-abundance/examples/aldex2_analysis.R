# Differential abundance testing with ALDEx2
library(ALDEx2)
library(phyloseq)
library(ggplot2)

# Public microbiome datasets:
# - Bioconductor: microbiomeDataSets package (curated phyloseq objects)
# - curatedMetagenomicData package (1000s of samples)
# - HMP: https://hmpdacc.org (Human Microbiome Project)
# - Example: data('GlobalPatterns', package = 'phyloseq')

ps <- readRDS('phyloseq_object.rds')

# Filter low-abundance taxa
# 0.1 (10%): Require taxa present in at least 10% of samples. Standard prevalence filter.
# mean(x) > 10: Require mean abundance >10 reads. Removes rare taxa with unreliable estimates.
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) > 0.1 * nsamples(ps), TRUE)
ps_filt <- filter_taxa(ps_filt, function(x) mean(x) > 10, TRUE)
cat('Taxa after filtering:', ntaxa(ps_filt), '\n')

# Prepare data for ALDEx2
otu <- as.data.frame(otu_table(ps_filt))
if (!taxa_are_rows(ps_filt)) otu <- t(otu)

groups <- as.character(sample_data(ps_filt)$Group)
cat('Groups:', paste(unique(groups), collapse = ', '), '\n')

# Run ALDEx2
# mc.samples=128: Monte Carlo samples for posterior. 128 is standard; use 256+ for publication.
cat('\nRunning ALDEx2 (this may take a few minutes)...\n')
aldex_out <- aldex(otu, groups, mc.samples = 128, test = 'welch',
                   effect = TRUE, include.sample.summary = FALSE, denom = 'all')

# Summarize results
aldex_out$taxon <- rownames(aldex_out)
# q<0.05: FDR threshold. Use 0.1 for exploratory, 0.01 for stringent.
# |effect|>1: Minimum effect size. 1.0 is moderate; 0.5 for small effects, 2.0 for large.
aldex_out$significant <- aldex_out$we.eBH < 0.05 & abs(aldex_out$effect) > 1

cat('\nResults summary:\n')
cat('  Total taxa tested:', nrow(aldex_out), '\n')
cat('  Significant (q<0.05, |effect|>1):', sum(aldex_out$significant), '\n')
cat('  Enriched in', unique(groups)[1], ':', sum(aldex_out$significant & aldex_out$effect > 0), '\n')
cat('  Enriched in', unique(groups)[2], ':', sum(aldex_out$significant & aldex_out$effect < 0), '\n')

# Add taxonomy
if (!is.null(tax_table(ps_filt))) {
    tax <- as.data.frame(tax_table(ps_filt))
    aldex_out <- merge(aldex_out, tax, by.x = 'taxon', by.y = 'row.names', all.x = TRUE)
}

# Sort by effect size
aldex_out <- aldex_out[order(-abs(aldex_out$effect)), ]
write.csv(aldex_out, 'aldex2_results.csv', row.names = FALSE)

# Effect plot
p <- ggplot(aldex_out, aes(x = effect, y = -log10(we.eBH))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'grey50') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed', color = 'grey50') +
    scale_color_manual(values = c('grey60', 'firebrick'), name = 'Significant') +
    theme_minimal() +
    labs(x = 'Effect Size', y = '-log10(Adjusted P-value)', title = 'ALDEx2 Differential Abundance')

ggsave('aldex2_effect_plot.pdf', p, width = 7, height = 6)
cat('\nSaved: aldex2_results.csv, aldex2_effect_plot.pdf\n')
