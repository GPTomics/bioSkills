#!/usr/bin/env Rscript
# Isoform switching analysis with IsoformSwitchAnalyzeR
# Identifies switches and predicts functional consequences

library(IsoformSwitchAnalyzeR)

# Step 1: Import transcript quantification from Salmon
import_salmon_data <- function(salmon_dir, design_matrix) {
    salmonQuant <- importIsoformExpression(
        parentDir = salmon_dir,
        addIsofomIdAsColumn = TRUE
    )

    # Create switch analysis object
    switchAnalyzeRlist <- importRdata(
        isoformCountMatrix = salmonQuant$counts,
        isoformRepExpression = salmonQuant$abundance,
        designMatrix = design_matrix,
        isoformExonAnnoation = 'annotation.gtf',
        isoformNtFasta = 'transcripts.fa',
        showProgress = TRUE
    )

    return(switchAnalyzeRlist)
}

# Step 2: Filter and test for switches
run_switch_analysis <- function(switchList) {
    # Filter lowly expressed isoforms
    # geneExpressionCutoff: minimum gene-level TPM
    switchList <- preFilter(
        switchList,
        geneExpressionCutoff = 1,
        isoformExpressionCutoff = 0,
        removeSingleIsoformGenes = TRUE
    )

    cat(sprintf('Genes after filtering: %d\n', length(unique(switchList$isoformFeatures$gene_id))))

    # Test for isoform switches using DEXSeq
    # reduceToSwitchingGenes keeps only genes with significant switches
    switchList <- isoformSwitchTestDEXSeq(
        switchList,
        reduceToSwitchingGenes = TRUE
    )

    return(switchList)
}

# Step 3: Extract sequences for external annotation
extract_sequences <- function(switchList, output_dir) {
    switchList <- extractSequence(
        switchList,
        pathToOutput = output_dir,
        writeToFile = TRUE
    )

    cat(sprintf('Sequences written to: %s\n', output_dir))
    cat('Run external tools:\n')
    cat('  - CPC2: http://cpc2.gao-lab.org/\n')
    cat('  - Pfam: https://www.ebi.ac.uk/Tools/hmmer/\n')
    cat('  - SignalP: https://services.healthtech.dtu.dk/service.php?SignalP\n')
    cat('  - IUPred2A: https://iupred2a.elte.hu/\n')

    return(switchList)
}

# Step 4: Import external annotations and analyze consequences
analyze_consequences <- function(switchList, annotation_dir) {
    # Import CPC2 coding potential predictions
    if (file.exists(file.path(annotation_dir, 'cpc2_results.txt'))) {
        switchList <- analyzeCPC2(
            switchList,
            pathToCPC2resultFile = file.path(annotation_dir, 'cpc2_results.txt'),
            removeNoncodinORFs = TRUE
        )
    }

    # Import Pfam domain predictions
    if (file.exists(file.path(annotation_dir, 'pfam_results.txt'))) {
        switchList <- analyzePFAM(
            switchList,
            pathToPFAMresultFile = file.path(annotation_dir, 'pfam_results.txt')
        )
    }

    # Analyze functional consequences
    # dIFcutoff: minimum isoform fraction change (default 0.1)
    switchList <- analyzeSwitchConsequences(
        switchList,
        consequencesToAnalyze = c(
            'intron_retention',
            'coding_potential',
            'ORF_seq_similarity',
            'NMD_status',
            'domains_identified',
            'signal_peptide_identified'
        ),
        dIFcutoff = 0.1,
        showProgress = TRUE
    )

    return(switchList)
}

# Step 5: Summarize and visualize results
summarize_results <- function(switchList) {
    # Extract summary of significant switches
    summary <- extractSwitchSummary(
        switchList,
        filterForConsequences = TRUE
    )

    cat('\nSwitch Summary:\n')
    print(summary)

    # Consequence enrichment analysis
    enrichment <- extractConsequenceEnrichment(
        switchList,
        consequencesToAnalyze = 'all'
    )

    return(list(summary = summary, enrichment = enrichment))
}

# Step 6: Plot individual gene switches
plot_gene_switch <- function(switchList, gene_name, condition1, condition2, output_file) {
    pdf(output_file, width = 10, height = 8)
    switchPlot(
        switchList,
        gene = gene_name,
        condition1 = condition1,
        condition2 = condition2
    )
    dev.off()
    cat(sprintf('Plot saved to: %s\n', output_file))
}

# Example usage
cat('IsoformSwitchAnalyzeR workflow\n')
cat('==============================\n\n')

# Example design matrix
# design <- data.frame(
#     sampleID = c('ctrl1', 'ctrl2', 'ctrl3', 'treat1', 'treat2', 'treat3'),
#     condition = c('control', 'control', 'control', 'treatment', 'treatment', 'treatment')
# )
#
# switchList <- import_salmon_data('salmon_quant/', design)
# switchList <- run_switch_analysis(switchList)
# switchList <- extract_sequences(switchList, 'sequences/')
# # ... run external tools ...
# switchList <- analyze_consequences(switchList, 'annotations/')
# results <- summarize_results(switchList)
# plot_gene_switch(switchList, 'TP53', 'control', 'treatment', 'tp53_switch.pdf')

cat('Provide paths to Salmon quantification and annotation files to run analysis\n')
