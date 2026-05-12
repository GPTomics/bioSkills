# Reference: MAGMA 1.10+, Open Targets Genetics API, PoPS, cS2G, ABC | Verify API if version differs
#
# Multi-evidence concordance scoring for effector-gene prioritization at GWAS loci.
# Combines six orthogonal evidence streams into a per-(locus, gene) concordance score:
#   1. Fine-mapping (SuSiE PIP + credible set purity)
#   2. Colocalization (coloc.abf or coloc.susie PP.H4)
#   3. Distance (variant-to-TSS within regulatory window)
#   4. PoPS polygenic priority (top-decile per locus)
#   5. Open Targets L2G score (gradient-boosting classifier)
#   6. ABC / ENCODE-rE2G enhancer-gene linking
# Operational rule (Mountjoy 2021; Open Targets Genetics):
#   >= 3 of 6 = 'high-confidence'; >= 4 = 'strong'; >= 5 = 'near-certain'

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
})

candidates <- read_tsv('locus_candidates.tsv', show_col_types = FALSE)

# Expected columns per row (one row per (locus, gene)):
#   locus, gene, pip_top_variant, credible_set_purity, coloc_pph4,
#   distance_to_tss, pops_score, pops_decile_rank, l2g_score,
#   abc_score, encode_re2g_score

thresh <- list(
    pip = 0.5,
    purity = 0.5,
    pph4 = 0.7,
    distance = 100000,
    pops_decile = 1,
    l2g = 0.5,
    abc = 0.02,
    re2g = 0.5)

scored <- candidates %>%
    mutate(
        pass_finemap = pip_top_variant > thresh$pip & credible_set_purity > thresh$purity,
        pass_coloc = coloc_pph4 >= thresh$pph4,
        pass_distance = distance_to_tss <= thresh$distance,
        pass_pops = pops_decile_rank == thresh$pops_decile,
        pass_l2g = l2g_score >= thresh$l2g,
        pass_abc_re2g = abc_score >= thresh$abc | encode_re2g_score >= thresh$re2g,
        concordance = as.integer(pass_finemap) + as.integer(pass_coloc) +
                      as.integer(pass_distance) + as.integer(pass_pops) +
                      as.integer(pass_l2g) + as.integer(pass_abc_re2g),
        confidence_tier = case_when(
            concordance >= 5 ~ 'near_certain',
            concordance == 4 ~ 'strong',
            concordance == 3 ~ 'high',
            concordance == 2 ~ 'suggestive',
            TRUE ~ 'associational_only'))

high_conf <- scored %>%
    filter(concordance >= 3) %>%
    arrange(locus, desc(concordance), desc(l2g_score)) %>%
    select(locus, gene, concordance, confidence_tier, l2g_score,
           pops_decile_rank, coloc_pph4, pip_top_variant,
           distance_to_tss, abc_score, encode_re2g_score)

write_tsv(scored, 'effector_gene_concordance_all.tsv')
write_tsv(high_conf, 'effector_gene_high_confidence.tsv')

# Per-locus winner: top gene by concordance (tiebreak L2G)
per_locus_winner <- scored %>%
    group_by(locus) %>%
    slice_max(order_by = concordance + l2g_score / 100, n = 1, with_ties = FALSE) %>%
    ungroup()

write_tsv(per_locus_winner, 'effector_gene_per_locus_winner.tsv')

# L2G vs PoPS concordance flag at each locus
l2g_pops_check <- scored %>%
    group_by(locus) %>%
    mutate(
        is_l2g_top = l2g_score == max(l2g_score, na.rm = TRUE),
        is_pops_top = pops_score == max(pops_score, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(is_l2g_top | is_pops_top) %>%
    group_by(locus) %>%
    summarise(
        l2g_top_gene = gene[is_l2g_top][1],
        pops_top_gene = gene[is_pops_top][1],
        methods_concordant = l2g_top_gene == pops_top_gene,
        .groups = 'drop')

write_tsv(l2g_pops_check, 'l2g_pops_concordance.tsv')

cat('Effector-gene prioritization complete.\n')
cat('High-confidence candidates (concordance >= 3):', nrow(high_conf), '\n')
cat('L2G + PoPS concordant loci:', sum(l2g_pops_check$methods_concordant, na.rm = TRUE),
    'of', nrow(l2g_pops_check), '\n')
