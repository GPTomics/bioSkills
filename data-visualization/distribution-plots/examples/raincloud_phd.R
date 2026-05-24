# Reference: ggplot2 3.5+, ggbeeswarm 0.7+, ggdist 3.3+, gghalves 0.1.4+ | Verify API if version differs

# PhD-level distribution plot showing the N-based decision tree:
# (1) raw points for N<30, (2) raincloud for N 30-200, (3) letter-value for N>200,
# (4) Sheather-Jones bandwidth, (5) N annotated on axis.

library(ggplot2)
library(ggbeeswarm)
library(gghalves)
library(dplyr)

# Compute per-group N for axis annotation
n_per_group <- df %>% count(group) %>%
    mutate(label = paste0(group, '\n(n=', n, ')'))
df <- df %>% left_join(n_per_group, by = 'group')

okabe_subset <- c('#0072B2', '#D55E00', '#009E73')

# 1. SMALL N (<30) -- beeswarm + median crossbar; show every point
p_small <- ggplot(df_small, aes(group, value, color = group)) +
    geom_quasirandom(method = 'quasirandom', width = 0.3, alpha = 0.7, size = 1.5) +
    stat_summary(fun = median, geom = 'crossbar', width = 0.5,
                 color = 'black', linewidth = 0.4) +
    scale_color_manual(values = okabe_subset) +
    scale_x_discrete(labels = function(x) n_per_group$label[match(x, n_per_group$group)]) +
    labs(x = NULL, y = 'Expression') +
    theme_classic(base_size = 10) + theme(legend.position = 'none')

# 2. MEDIUM N (30-200) -- horizontal raincloud (Allen 2019)
p_raincloud <- ggplot(df_med, aes(group, value, fill = group, color = group)) +
    geom_half_violin(side = 'r', alpha = 0.7,
                     trim = FALSE,                          # don't cut at data range
                     adjust = 1, bw = 'SJ',                 # Sheather-Jones bandwidth
                     position = position_nudge(x = 0.15)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7,
                 position = position_nudge(x = -0.05)) +
    geom_half_point(side = 'l', alpha = 0.5, size = 1.2,
                    range_scale = 0.4,
                    position = position_nudge(x = -0.2)) +
    scale_fill_manual(values = okabe_subset) +
    scale_color_manual(values = okabe_subset) +
    coord_flip() +
    labs(x = NULL, y = 'Biomarker level') +
    theme_classic(base_size = 10) + theme(legend.position = 'none')

# 3. LARGE N (>200) -- letter-value plot preserves tails better than boxplot
library(lvplot)
p_lv <- ggplot(df_large, aes(group, value, fill = group)) +
    geom_lv(k = 5, alpha = 0.7) +                          # 5 letter-value levels
    scale_fill_manual(values = okabe_subset) +
    labs(x = NULL, y = 'Expression', caption = sprintf('N >= 200 per group')) +
    theme_classic(base_size = 10) + theme(legend.position = 'none')

# 4. SPLIT VIOLIN for 2-condition comparison within each group (paired data)
library(introdataviz)                                       # remotes::install_github('PsyTeachR/introdataviz')
p_split <- ggplot(df_paired, aes(cluster, expression, fill = condition)) +
    geom_split_violin(alpha = 0.7, trim = FALSE, bw = 'SJ') +
    geom_boxplot(width = 0.15, position = position_dodge(0.5), outlier.shape = NA) +
    scale_fill_manual(values = c(Control = '#56B4E9', Treatment = '#D55E00')) +
    labs(x = 'Cluster', y = 'Expression', fill = NULL) +
    theme_classic(base_size = 10) + theme(legend.position = 'top')

# 5. EXPORT -- 89mm = Nature single column width
ggsave('raincloud.pdf', p_raincloud, width = 89, height = 70, units = 'mm',
       device = cairo_pdf)
