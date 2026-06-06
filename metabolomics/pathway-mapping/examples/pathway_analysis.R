# Reference: R 4.3+ (base stats only) | Verify API if version differs
#
# Two teaching points, both runnable with no database downloads:
#   (1) ORA on a small identified-compound list via the hypergeometric test,
#       showing that the BACKGROUND choice (assay-coverage vs 'all of KEGG')
#       flips the p-value.
#   (2) The mummichog background pitfall: building the null from the FULL
#       feature table (R_all) gives an honest test, but using the significant
#       features alone as the background reports spurious significance.
# The real MetaboAnalystR/FELLA chains are in SKILL.md; this script demonstrates
# the STATISTICS that those tools implement, on synthetic data.

set.seed(1)

# ---- (1) ORA: the background silently controls the p-value ----------------

# A toy pathway of 10 compounds; our identified hit list overlaps it by 4.
pathway_size <- 10
hits_in_pathway <- 4
hit_list_size <- 18

# A correct ORA background is the set of compounds the assay could detect and
# identify (assay coverage), NOT the whole database. Wieder 2021 measured this:
# a KEGG-human background held ~3,373 compounds vs 286-1,110 actually measurable.
assay_coverage_background <- 300
all_of_kegg_background <- 3373

ora_pvalue <- function(hits, pathway, hit_list, universe) {
    # P(>= hits successes) under the hypergeometric null
    phyper(hits - 1, m = pathway, n = universe - pathway, k = hit_list, lower.tail = FALSE)
}

p_assay <- ora_pvalue(hits_in_pathway, pathway_size, hit_list_size, assay_coverage_background)
p_kegg <- ora_pvalue(hits_in_pathway, pathway_size, hit_list_size, all_of_kegg_background)

cat('--- ORA: same hits, different background ---\n')
cat(sprintf('assay-coverage background (n=%d):  p = %.4g\n', assay_coverage_background, p_assay))
cat(sprintf("'all of KEGG' background (n=%d):    p = %.4g  (inflated)\n", all_of_kegg_background, p_kegg))
cat(sprintf('inflation factor: %.1fx more significant with the oversized background\n\n', p_assay / p_kegg))

# ---- (2) Mummichog: the permutation-null background pitfall ----------------

# Synthetic untargeted feature table. Some features mass-match a compound in the
# pathway of interest. There is NO genuine perturbation here: pathway membership
# is INDEPENDENT of significance, so an honest test must return non-significant.
n_features <- 2000
feature_table <- data.frame(
    feature = seq_len(n_features),
    in_pathway = rbinom(n_features, 1, 0.08),
    significant = rbinom(n_features, 1, 0.10)
)

# Mummichog scores how over-represented pathway-mapped compounds are among the
# SIGNIFICANT features, against a background of the full feature table. We use a
# hypergeometric tail as the enrichment test; the only thing that changes
# between the two runs is the BACKGROUND POOL.
n_sig <- sum(feature_table$significant)
hits <- sum(feature_table$significant == 1 & feature_table$in_pathway == 1)

enrich_p <- function(pool) {
    pathway_in_pool <- sum(pool$in_pathway)
    sig_in_pool <- sum(pool$significant)
    phyper(hits - 1, m = pathway_in_pool, n = nrow(pool) - pathway_in_pool, k = sig_in_pool, lower.tail = FALSE)
}

# Correct background: the ENTIRE feature table (R_all). The significant features
# are a random slice of it w.r.t. pathway membership, so enrichment is null.
p_full <- enrich_p(feature_table)

# Cardinal error: collapse the background to the significant features alone. The
# background now equals the query, the test has no proper null to compare
# against, and the result is meaningless (here it degenerates to p = 1).
p_sig <- enrich_p(feature_table[feature_table$significant == 1, ])

cat('--- Mummichog: which pool is the background ---\n')
cat(sprintf('background = FULL feature table (R_all):  p = %.4g  (correct: no real signal)\n', p_full))
cat(sprintf('background = significant features only:    p = %.4g  (broken null)\n', p_sig))
cat('When the background equals the significant set, the test loses its null;\n')
cat('the correct background is the entire feature table (R_all).\n')

out <- file.path(tempdir(), 'pathway_analysis_demo.txt')
writeLines(c(
    sprintf('ORA p (assay background):    %.4g', p_assay),
    sprintf('ORA p (all-of-KEGG):         %.4g', p_kegg),
    sprintf('mummichog p (R_all bg):      %.4g', p_full),
    sprintf('mummichog p (sig-only bg):   %.4g', p_sig)
), out)
cat(sprintf('\nSummary written to %s\n', out))
