# Reference: bcftools 1.19+ | Verify API if version differs
library(hierfstat)
library(adegenet)
library(ggplot2)

# --- Load genotype data ---
# genepop format: populations separated by 'Pop' keyword
data_genind <- read.genepop('populations.gen')

# Alternatively, from VCF via vcfR:
# library(vcfR); vcf <- read.vcfR('variants.vcf')
# data_genind <- vcfR2genind(vcf)
# pop(data_genind) <- pop_assignments

cat('Populations:', nlevels(pop(data_genind)), '\n')
cat('Individuals:', nInd(data_genind), '\n')
cat('Loci:', nLoc(data_genind), '\n')

# --- Convert to hierfstat format ---
data_hf <- genind2hierfstat(data_genind)

# --- Basic F-statistics (Weir & Cockerham estimators) ---
bstats <- basic.stats(data_hf)

cat('\n--- Overall Statistics ---\n')
cat('Fis:', round(bstats$overall['Fis'], 4), '\n')
cat('Fst:', round(bstats$overall['Fst'], 4), '\n')
cat('Fit:', round(bstats$overall['Fit'], 4), '\n')

# --- Per-population diversity ---
pop_names <- levels(pop(data_genind))
pop_diversity <- data.frame(
    population = pop_names,
    Ho = colMeans(bstats$Ho, na.rm = TRUE),
    He = colMeans(bstats$Hs, na.rm = TRUE),
    Fis = colMeans(bstats$Fis, na.rm = TRUE),
    n = as.numeric(table(pop(data_genind)))
)

cat('\n--- Per-Population Diversity ---\n')
print(pop_diversity)

# --- Rarefied allelic richness ---
# Rarefaction corrects for unequal sample sizes across populations
ar <- allelic.richness(data_hf)
pop_diversity$AR <- colMeans(ar$Ar, na.rm = TRUE)

cat('\nRarefied allelic richness:\n')
print(setNames(pop_diversity$AR, pop_diversity$population))

# --- Private alleles ---
library(poppr)
pa <- private_alleles(data_genind, count.alleles = TRUE)
pop_diversity$private_alleles <- rowSums(pa > 0)

# --- Pairwise Fst ---
pw_fst <- pairwise.WCfst(data_hf)
cat('\n--- Pairwise Fst ---\n')
print(round(pw_fst, 4))

# Bootstrap confidence intervals
# nboot=1000: standard for publication-quality intervals
boot_result <- boot.ppfst(data_hf, nboot = 1000)

cat('\n--- Pairwise Fst 95% CI (lower) ---\n')
print(round(boot_result$ll, 4))
cat('\n--- Pairwise Fst 95% CI (upper) ---\n')
print(round(boot_result$ul, 4))

# --- Visualization: diversity comparison ---
div_long <- data.frame(
    population = rep(pop_diversity$population, 2),
    metric = rep(c('Ho', 'He'), each = nrow(pop_diversity)),
    value = c(pop_diversity$Ho, pop_diversity$He)
)

p1 <- ggplot(div_long, aes(x = population, y = value, fill = metric)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.7) +
    labs(x = 'Population', y = 'Heterozygosity',
         title = 'Observed vs Expected Heterozygosity', fill = '') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('heterozygosity_comparison.pdf', p1, width = 8, height = 5)

# Allelic richness barplot
p2 <- ggplot(pop_diversity, aes(x = population, y = AR)) +
    geom_bar(stat = 'identity', fill = 'steelblue', width = 0.6) +
    labs(x = 'Population', y = 'Rarefied Allelic Richness',
         title = 'Allelic Richness by Population') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('allelic_richness.pdf', p2, width = 8, height = 5)

# Fst heatmap
fst_df <- as.data.frame(as.table(pw_fst))
colnames(fst_df) <- c('Pop1', 'Pop2', 'Fst')
fst_df <- fst_df[!is.na(fst_df$Fst), ]

p3 <- ggplot(fst_df, aes(x = Pop1, y = Pop2, fill = Fst)) +
    geom_tile() +
    geom_text(aes(label = round(Fst, 3)), size = 3) +
    scale_fill_gradient(low = 'white', high = 'darkred') +
    labs(title = 'Pairwise Fst') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('pairwise_fst_heatmap.pdf', p3, width = 7, height = 6)

# --- Export summary ---
write.csv(pop_diversity, 'population_diversity_summary.csv', row.names = FALSE)
cat('\nSummary written to population_diversity_summary.csv\n')
