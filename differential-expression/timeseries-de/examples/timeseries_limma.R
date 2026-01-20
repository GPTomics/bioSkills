#!/usr/bin/env Rscript
# Time-series differential expression with limma splines

library(limma)
library(edgeR)
library(splines)
library(ggplot2)

run_timeseries_de <- function(counts_file, metadata_file, output_prefix = 'timeseries') {
    cat('Loading data...\n')
    counts <- read.csv(counts_file, row.names = 1)
    meta <- read.csv(metadata_file)

    # Design with natural splines
    cat('Creating design matrix...\n')
    time <- meta$time
    group <- factor(meta$group)

    # Use 3 degrees of freedom for spline (adjust based on timepoints)
    design <- model.matrix(~ group * ns(time, df = 3))

    # Voom transformation
    cat('Running voom...\n')
    dge <- DGEList(counts = counts)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)

    # Fit model
    cat('Fitting linear model...\n')
    fit <- lmFit(v, design)
    fit <- eBayes(fit)

    # Test for time effect
    time_coefs <- grep('ns\\(time', colnames(design))
    results <- topTable(fit, coef = time_coefs, number = Inf, sort.by = 'F')

    # Test for group:time interaction
    interaction_coefs <- grep('group.*ns\\(time', colnames(design))
    if (length(interaction_coefs) > 0) {
        interaction_results <- topTable(fit, coef = interaction_coefs, number = Inf, sort.by = 'F')
        write.csv(interaction_results, paste0(output_prefix, '_interaction.csv'))
    }

    # Write results
    write.csv(results, paste0(output_prefix, '_time_effect.csv'))

    # Significant genes
    sig_genes <- rownames(results)[results$adj.P.Val < 0.05]
    cat(sprintf('Significant time-varying genes: %d\n', length(sig_genes)))

    # Plot top genes
    cat('Plotting top genes...\n')
    pdf(paste0(output_prefix, '_top_genes.pdf'), width = 10, height = 8)
    top_genes <- head(rownames(results), 9)

    par(mfrow = c(3, 3))
    for (gene in top_genes) {
        expr <- as.numeric(v$E[gene, ])
        plot(time, expr, col = as.numeric(group), pch = 16,
             main = gene, xlab = 'Time', ylab = 'Expression')
        for (g in levels(group)) {
            idx <- group == g
            lines(sort(time[idx]), expr[idx][order(time[idx])], col = which(levels(group) == g))
        }
    }
    dev.off()

    cat('Results saved to:', paste0(output_prefix, '_time_effect.csv\n'))
    return(results)
}

# Run if executed directly
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    counts_file <- args[1]
    metadata_file <- args[2]
    output_prefix <- if (length(args) > 2) args[3] else 'timeseries'
    run_timeseries_de(counts_file, metadata_file, output_prefix)
} else {
    cat('Usage: Rscript timeseries_limma.R counts.csv metadata.csv [output_prefix]\n')
    cat('\nMetadata should have columns: sample, time, group\n')
}
