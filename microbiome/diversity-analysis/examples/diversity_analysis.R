# Reference: phyloseq 1.46+, vegan 2.6+, picante 1.8+, GUniFrac 1.8+ | Verify API if version differs
# Alpha and beta diversity of an amplicon ASV table.
# The three knobs are made explicit: the rarefaction DEPTH (declared + dropped samples reported),
# the TREE (must be attached; prefer SEPP/Greengenes2 over de novo), and the METRIC (weighted AND unweighted).
library(phyloseq)
library(vegan)
library(ggplot2)

# Public amplicon datasets / built-in demo:
# - data('GlobalPatterns', package = 'phyloseq')  # ships a tree, used below
# - Qiita https://qiita.ucsd.edu ; Earth Microbiome Project https://earthmicrobiome.org

data('GlobalPatterns')
ps <- GlobalPatterns
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
cat('Samples:', nsamples(ps), '| ASVs:', ntaxa(ps), '\n')

# --- Knob 1: the sampling depth. Read the distribution; do NOT default to min(). ---
depths <- sort(sample_sums(ps))
cat('Per-sample depth range:', min(depths), '-', max(depths), '| median:', median(depths), '\n')
# chosen_depth sits below the median to retain most samples while saturating richness;
# min(sample_sums) is the worst choice (one tiny library drags everyone to noise).
chosen_depth <- as.integer(quantile(depths, 0.10))
dropped <- names(sample_sums(ps))[sample_sums(ps) < chosen_depth]
cat('Chosen depth:', chosen_depth, '| samples dropped:', length(dropped), '->',
    if (length(dropped)) paste(dropped, collapse = ', ') else 'none', '\n')

ps_rare <- rarefy_even_depth(ps, sample.size = chosen_depth, rngseed = 42, replace = FALSE)

# --- Alpha diversity: span the Hill spectrum, report effective species (base-invariant) ---
alpha <- estimate_richness(ps_rare, measures = c('Observed', 'Shannon', 'InvSimpson'))
alpha$Shannon_eff <- exp(alpha$Shannon)   # Hill q=1; QIIME2 log2 vs R ln differ, exp(H') dodges the base
alpha$Group <- sample_data(ps_rare)$SampleType

cat('\nMean Shannon (nats) by group:\n')
print(aggregate(Shannon ~ Group, data = alpha, FUN = mean))

kw <- kruskal.test(Shannon ~ Group, data = alpha)   # non-parametric; escalate to a mixed model for covariates
cat(sprintf('\nKruskal-Wallis (Shannon): chi-sq = %.2f, p = %.4g\n', kw$statistic, kw$p.value))

p_alpha <- ggplot(alpha, aes(x = Group, y = Shannon, fill = Group)) +
    geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'Shannon diversity', y = 'Shannon (nats)')
ggsave('alpha_diversity.pdf', p_alpha, width = 6, height = 5)

# --- Knob 2: the tree. UniFrac/Faith PD inherit it; GlobalPatterns ships one. ---
# Faith PD (picante) if installed; phylo metrics need the phy_tree slot.
if (requireNamespace('picante', quietly = TRUE) && !is.null(phy_tree(ps_rare, errorIfNULL = FALSE))) {
    faith <- picante::pd(as(t(otu_table(ps_rare)), 'matrix'), phy_tree(ps_rare), include.root = TRUE)
    alpha$Faith_PD <- faith$PD
    cat('\nFaith PD computed for', nrow(faith), 'samples\n')
}

# --- Knob 3: the metric. Report weighted AND unweighted UniFrac (they can flip the story). ---
wu  <- UniFrac(ps_rare, weighted = TRUE)    # abundant-lineage view
uwu <- UniFrac(ps_rare, weighted = FALSE)   # rare-lineage + topology view
bray <- phyloseq::distance(ps_rare, method = 'bray')

meta <- data.frame(sample_data(ps_rare))
# permutations=999: resolution floor for p ~ 0.001; use 9999 for publication.
perm_wu  <- adonis2(wu  ~ Group, data = meta, permutations = 999)
perm_uwu <- adonis2(uwu ~ Group, data = meta, permutations = 999)
cat('\nPERMANOVA weighted UniFrac:   R2 =', round(perm_wu$R2[1], 3),  '| p =', perm_wu$`Pr(>F)`[1], '\n')
cat('PERMANOVA unweighted UniFrac: R2 =', round(perm_uwu$R2[1], 3), '| p =', perm_uwu$`Pr(>F)`[1], '\n')

# MANDATORY dispersion check: a significant PERMANOVA can be spread, not location.
bd <- betadisper(wu, meta$Group)
pt <- permutest(bd)
cat('betadisper (weighted UniFrac) p =', pt$tab$`Pr(>F)`[1],
    '-> if significant, the PERMANOVA is location-vs-dispersion ambiguous\n')

# PCoA on weighted UniFrac
pcoa <- ordinate(ps_rare, method = 'PCoA', distance = wu)
p_beta <- plot_ordination(ps_rare, pcoa, color = 'SampleType') +
    stat_ellipse(level = 0.95) + theme_minimal() +
    labs(title = sprintf('PCoA (weighted UniFrac)  R2=%.2f p=%.3f', perm_wu$R2[1], perm_wu$`Pr(>F)`[1]))
ggsave('beta_diversity_wunifrac.pdf', p_beta, width = 7, height = 6)

cat('\nDone. Reported: chosen depth + dropped samples; weighted AND unweighted UniFrac; betadisper.\n')
cat('Plots: alpha_diversity.pdf, beta_diversity_wunifrac.pdf\n')
