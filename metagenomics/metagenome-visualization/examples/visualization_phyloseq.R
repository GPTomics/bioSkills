# Reference: phyloseq 1.46+, vegan 2.6+, microViz 0.12+, ALDEx2 1.34+ | Verify API if version differs
# Honest community stats from a MetaPhlAn table: CLR ordination, PERMANOVA PAIRED WITH betadisper,
# Hill-number diversity, and a two-tool differential-abundance consensus. Bray-Curtis is incoherent
# on relative abundance; pair every adonis2 with a dispersion test.
library(phyloseq)
library(vegan)

abundance <- read.table('merged_abundance.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
species <- abundance[grepl('s__', rownames(abundance)) & !grepl('t__', rownames(abundance)), ]
rownames(species) <- gsub('s__', '', sapply(strsplit(rownames(species), '\\|'), tail, 1))

otu <- otu_table(as.matrix(species), taxa_are_rows = TRUE)
meta <- data.frame(Group = rep(c('Control', 'Treatment'), length.out = ncol(species)),
                   row.names = colnames(species))
ps <- phyloseq(otu, sample_data(meta))

# Beta diversity: declare the metric. Bray-Curtis is the field default but compositionally incoherent;
# for a compositional ordination use microViz CLR-PCA. Show the conclusion survives both.
mat <- t(as(otu_table(ps), 'matrix'))               # samples as rows
dist_bc <- vegdist(mat, method = 'bray')
pm <- adonis2(dist_bc ~ Group, data = meta, permutations = 999, by = 'terms')

# ALWAYS pair PERMANOVA with a dispersion test - a significant adonis2 can be a spread difference,
# not a location shift (Anderson & Walsh 2013).
bd <- betadisper(dist_bc, meta$Group)
disp <- permutest(bd, permutations = 999)

# Alpha diversity as Hill numbers (effective species). estimate_richness Observed/Chao1 assume INTEGER
# counts - valid for Bracken, meaningless on MetaPhlAn percentages, so prefer evenness-weighted q=1/q=2.
alpha <- estimate_richness(ps, measures = c('Shannon', 'InvSimpson'))
alpha$hill_q1 <- exp(alpha$Shannon)
alpha$hill_q2 <- alpha$InvSimpson

cat('PERMANOVA p =', pm[['Pr(>F)']][1], '| betadisper p =', disp$tab[['Pr(>F)']][1], '\n')
cat('If betadisper is significant, the PERMANOVA is ambiguous - report both.\n')

# Differential abundance: never one tool, never uncorrected Wilcoxon on relative abundance.
# Run >=2 compositional tools (e.g. ALDEx2 + ANCOM-BC), prevalence-filter, BH-correct, intersect.
# library(ALDEx2)
# aldex_res <- aldex(round(as(otu_table(ps), 'matrix')), meta$Group, test = 't', effect = TRUE)
# Report taxa with we.eBH < 0.05 AND |effect| > 1, intersected with a second tool's hits.
