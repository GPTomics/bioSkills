# Reference: ggpubr 0.6+, ggsignif 0.6+, rstatix 0.7+ | Verify API if version differs

# PhD-level statistical annotation encoding the four correctness traps:
# (1) Wilcoxon for non-normal small N, (2) Holm adjustment for multiple pairs,
# (3) paired vs unpaired distinction, (4) effect size in caption.

library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

# 1. NORMALITY CHECK -- pick test based on data, not by default
shapiro_results <- df %>% group_by(group) %>%
    summarise(shapiro_p = shapiro.test(value)$p.value)
# If any group p < 0.05, prefer Wilcoxon over t-test for that comparison

# 2. PAIRWISE TESTS with Holm adjustment
pairs <- list(c('Control', 'Treatment'),
              c('Control', 'Vehicle'),
              c('Treatment', 'Vehicle'))

# Wilcoxon for non-normal; t.test for normal
stat_test <- df %>%
    pairwise_wilcox_test(value ~ group,
                          comparisons = pairs,
                          p.adjust.method = 'holm') %>%
    add_xy_position(x = 'group', step.increase = 0.1)

# 3. EFFECT SIZE -- Cliff's delta for non-parametric (Cohen's d for parametric)
effect_sizes <- df %>% wilcox_effsize(value ~ group, ci = TRUE)
print(effect_sizes)                              # report in caption / supplementary

# 4. PLOT with brackets and overall test
p <- ggboxplot(df, x = 'group', y = 'value', color = 'group',
               add = 'jitter', palette = c('#0072B2', '#D55E00', '#009E73'),
               outlier.shape = NA) +
    stat_pvalue_manual(stat_test, label = 'p.adj.signif',     # asterisks per adjusted p
                       tip.length = 0.01, bracket.size = 0.3) +
    stat_compare_means(method = 'kruskal.test',                # overall non-parametric
                       label.y = 1.15 * max(df$value),
                       label = 'p.format', size = 3) +
    labs(x = NULL, y = 'Value',
         caption = sprintf('Wilcoxon pairwise, Holm-adjusted. Cliff d ranges: %.2f-%.2f',
                           min(effect_sizes$effsize), max(effect_sizes$effsize))) +
    theme_classic(base_size = 10) + theme(legend.position = 'none')

ggsave('stat_annotated.pdf', p, width = 89, height = 75, units = 'mm', device = cairo_pdf)

# 5. PAIRED EXAMPLE -- before/after measurements per subject
paired_test <- df_paired %>%
    pairwise_wilcox_test(value ~ time, paired = TRUE, p.adjust.method = 'holm')

ggpaired(df_paired, x = 'time', y = 'value', id = 'subject_id',
         color = 'time', line.color = '#888888', line.size = 0.2,
         palette = c('#0072B2', '#D55E00')) +
    stat_pvalue_manual(paired_test, label = 'p.adj') +
    labs(caption = 'Wilcoxon signed-rank, paired by subject_id')

# 6. NESTED DATA -- LMM, NOT pairwise t-test on cells
# Wrong: pairwise tests treating each cell as independent
# Right: aggregate to per-subject value OR use mixed model
library(lme4)
nested_lmm <- lmer(value ~ group + (1 | subject_id), data = df_nested)
summary(nested_lmm)
# For pairwise contrasts: emmeans::emmeans(lmm, pairwise ~ group, adjust = 'holm')

# 7. ALTERNATIVE: ggsignif (lighter API)
library(ggsignif)
ggplot(df, aes(group, value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_signif(comparisons = pairs,
                test = 'wilcox.test',
                map_signif_level = TRUE,                    # asterisks
                step_increase = 0.1) +
    scale_fill_manual(values = c('#0072B2', '#D55E00', '#009E73'))
