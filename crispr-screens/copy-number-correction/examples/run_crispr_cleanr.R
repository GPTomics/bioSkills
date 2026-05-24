# Reference: CRISPRcleanR 3.0+, BiocManager | Verify API if version differs
#
# CRISPRcleanR pre-hoc copy-number bias correction for a single cancer-line screen.
# Outputs corrected counts compatible with MAGeCK / BAGEL2 / drugZ downstream.

library(CRISPRcleanR)

# === INPUTS ===
counts_file <- 'counts.txt'                       # tab-separated: sgRNA, gene, sample columns
plasmid_col <- 'Plasmid'
sample_cols <- c('Veh_r1', 'Veh_r2', 'Drug_r1', 'Drug_r2')
data(KY_Library_v1.0)                              # built-in KY library; replace with custom

# === LOAD AND NORMALIZE ===
raw_counts <- read.table(counts_file, header=TRUE, sep='\t', stringsAsFactors=FALSE)
norm <- ccr.NormfoldChanges(raw_counts,
                              min_reads = 30,                # min reads/sgRNA: drop low-coverage
                              EXPname   = 'cancer_screen',
                              libraryAnnotation = KY_Library_v1.0)

# === GENOME-SORTED LFC ===
gw_lfc <- ccr.logFCs2chromPos(norm$norm_fold_changes,
                                KY_Library_v1.0)

# === CN-AWARE SEGMENTATION + CORRECTION ===
# ccr.GWclean detects spatial enrichment / depletion patterns and shifts segments
# toward the global mean. Unsupervised: no CN profile required.
cleaned <- ccr.GWclean(gw_lfc,
                         display = TRUE,
                         label   = 'cancer_screen')

# === DERIVE CORRECTED COUNTS FOR DOWNSTREAM TOOLS ===
corrected_counts <- ccr.correctCounts(my_screen           = norm,
                                        correction          = cleaned,
                                        outprefix           = 'cancer_screen_cleanr',
                                        libraryAnnotation   = KY_Library_v1.0)
# Output: cancer_screen_cleanr_corrected_counts.txt (MAGeCK-compatible)

# === DIAGNOSTIC: PRE vs POST CORRECTION ===
# (assuming CN profile is available)
# Compare Spearman correlation of gene LFC with CN profile pre and post

cat('Pre-correction logFC mean (amplified):',
    mean(gw_lfc$gw_lfc[gw_lfc$CN > 4]), '\n')
cat('Post-correction logFC mean (amplified):',
    mean(cleaned$corrected_logFCs$correctedFC[cleaned$corrected_logFCs$CN > 4]), '\n')
cat('Correction effectiveness: difference =',
    mean(cleaned$corrected_logFCs$correctedFC[cleaned$corrected_logFCs$CN > 4]) -
    mean(gw_lfc$gw_lfc[gw_lfc$CN > 4]), '\n')

# === HIT CALLING ON CORRECTED COUNTS ===
# Feed cancer_screen_cleanr_corrected_counts.txt to MAGeCK / BAGEL2 / drugZ
# system('mageck test -k cancer_screen_cleanr_corrected_counts.txt ...')
