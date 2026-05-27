# Reference: exomePeak2 1.14+ (Bioconductor 3.18+), GenomicFeatures 1.54+, BSgenome.Hsapiens.UCSC.hg38 1.4+, ggplot2 3.5+ | Verify with packageVersion('exomePeak2'); ?exomePeak2 if installed releases differ.
# Differential m6A analysis with exomePeak2: populate bam_ip + bam_input (control) AND bam_treated_ip + bam_treated_input (treatment).
# Effect-size filter (|log2FC| >= 0.5) AND FDR < 0.05 are non-negotiable per McIntyre 2020 *Sci Rep* 10:6590.
# For batch / antibody-lot covariates, fall through to featureCounts-on-peaks -> DESeq2 (exomePeak2 top-level API does not accept covariates).

library(exomePeak2)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

ctrl_ip <- c('aligned/ctrl_IP1.bam',
             'aligned/ctrl_IP2.bam',
             'aligned/ctrl_IP3.bam')

ctrl_input <- c('aligned/ctrl_Input1.bam',
                'aligned/ctrl_Input2.bam',
                'aligned/ctrl_Input3.bam')

treat_ip <- c('aligned/treat_IP1.bam',
              'aligned/treat_IP2.bam',
              'aligned/treat_IP3.bam')

treat_input <- c('aligned/treat_Input1.bam',
                 'aligned/treat_Input2.bam',
                 'aligned/treat_Input3.bam')

txdb <- makeTxDbFromGFF('refs/annotation.gtf', format='gtf')

result <- exomePeak2(
    bam_ip            = ctrl_ip,
    bam_input         = ctrl_input,
    bam_treated_ip    = treat_ip,
    bam_treated_input = treat_input,
    txdb              = txdb,
    genome            = BSgenome.Hsapiens.UCSC.hg38,
    paired_end        = TRUE,
    library_type      = 'unstranded',
    peak_calling_mode = 'exon',
    save_dir          = 'exomepeak2_diff_output',
    experiment_name   = 'ctrl_vs_treat'
)

diff_table <- as.data.frame(result)
cat('Total tested peaks:', nrow(diff_table), '\n')

LOG2FC_THRESHOLD <- 0.5
FDR_THRESHOLD <- 0.05

diff_table$direction <- ifelse(diff_table$log2FC > 0, 'up_in_treat', 'down_in_treat')
diff_table$significant <- with(diff_table, padj < FDR_THRESHOLD & abs(log2FC) > LOG2FC_THRESHOLD)

sig <- diff_table[diff_table$significant, ]
cat('Differential peaks (padj <', FDR_THRESHOLD, ', |log2FC| >', LOG2FC_THRESHOLD, '):', nrow(sig), '\n')
cat('  up in treat:', sum(sig$direction == 'up_in_treat'), '\n')
cat('  down in treat:', sum(sig$direction == 'down_in_treat'), '\n')

write.csv(sig, 'exomepeak2_diff_output/ctrl_vs_treat/differential_peaks_significant.csv', row.names=FALSE)
write.csv(diff_table, 'exomepeak2_diff_output/ctrl_vs_treat/differential_peaks_all.csv', row.names=FALSE)

volcano <- ggplot(diff_table, aes(x=log2FC, y=-log10(padj), colour=significant)) +
    geom_point(alpha=0.5, size=0.8) +
    geom_vline(xintercept=c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype='dashed') +
    geom_hline(yintercept=-log10(FDR_THRESHOLD), linetype='dashed') +
    scale_colour_manual(values=c(`TRUE`='red', `FALSE`='grey60'), name='differential') +
    labs(x='log2 (treat / ctrl) MeRIP enrichment ratio',
         y='-log10 (FDR)',
         title='Differential m6A: ctrl vs treat',
         caption='MeRIP IP fold-change is RELATIVE enrichment, not absolute stoichiometry.\nValidate top hits with GLORI / m6Anet mod_ratio / SAC-seq.') +
    theme_minimal()

ggsave('exomepeak2_diff_output/ctrl_vs_treat/volcano.pdf', volcano, width=8, height=6)

cat('Wrote differential results and volcano to exomepeak2_diff_output/ctrl_vs_treat/\n')
cat('For batch / antibody-lot covariates: featureCounts-on-peaks -> DESeq2 with lot in design.\n')
