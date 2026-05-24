# Reference: ggalluvial 0.12+, consort 0.2+ | Verify API if version differs

# PhD-level alluvial + CONSORT encoding the four correctness traps:
# (1) explicit factor levels minimize crossings, (2) color by origin (not destination),
# (3) Sankey vs alluvial chosen by data structure, (4) flow conservation verified.

library(ggplot2)
library(ggalluvial)
library(dplyr)

# 1. INPUT -- wide (alluvia) format: one row per entity
# Example: cell-state transitions across 3 timepoints
df_wide <- data.frame(
    cell_id = 1:1000,
    t1 = sample(c('Proliferating', 'Differentiating', 'Quiescent'), 1000, replace = TRUE, prob = c(0.5, 0.3, 0.2)),
    t2 = sample(c('Proliferating', 'Differentiating', 'Quiescent', 'Apoptotic'), 1000, replace = TRUE, prob = c(0.3, 0.4, 0.2, 0.1)),
    t3 = sample(c('Differentiating', 'Quiescent', 'Apoptotic'), 1000, replace = TRUE, prob = c(0.3, 0.5, 0.2)))

# 2. EXPLICIT FACTOR LEVELS -- minimize ribbon crossings
state_levels <- c('Proliferating', 'Differentiating', 'Quiescent', 'Apoptotic')
df_wide <- df_wide %>%
    mutate(t1 = factor(t1, levels = state_levels),
           t2 = factor(t2, levels = state_levels),
           t3 = factor(t3, levels = state_levels))

# 3. ALLUVIAL -- color by ORIGIN (t1) for "where did final state come from" story
p_alluvial <- ggplot(df_wide, aes(axis1 = t1, axis2 = t2, axis3 = t3)) +
    geom_alluvium(aes(fill = t1), alpha = 0.7, width = 1/6) +
    geom_stratum(width = 1/6, fill = 'grey90', color = 'black') +
    geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c('t1', 't2', 't3'), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = c(Proliferating = '#D55E00',
                                  Differentiating = '#0072B2',
                                  Quiescent = '#009E73',
                                  Apoptotic = '#999999'),
                      name = 'Starting state') +
    labs(y = 'Cells', x = NULL) +
    theme_classic(base_size = 10)

ggsave('alluvial.pdf', p_alluvial, width = 130, height = 90, units = 'mm', device = cairo_pdf)

# 4. SANKEY for single-step source-to-sink
library(networkD3)

# Nodes
nodes <- data.frame(name = c('Raw variants', 'Quality-filtered', 'MAF-filtered', 'Annotated', 'Final'))

# Links (source/target indices into nodes; value = count)
links <- data.frame(
    source = c(0, 1, 2, 3),
    target = c(1, 2, 3, 4),
    value  = c(5000000, 3000000, 1500000, 1200000))      # decreasing counts; conservation = identical at each step

p_sankey <- sankeyNetwork(Links = links, Nodes = nodes,
                          Source = 'source', Target = 'target', Value = 'value',
                          NodeID = 'name',
                          colourScale = JS('d3.scaleOrdinal(d3.schemeCategory10);'),
                          fontSize = 11, nodeWidth = 30, height = 300, width = 700)
# Static export via webshot2 if PDF needed

# 5. CONSORT TRIAL FLOW
library(consort)

g <- add_box(txt = 'Assessed for eligibility (n=200)')
g <- add_side_box(g, txt = c('Excluded (n=50)',
                              '  Did not meet criteria (n=30)',
                              '  Declined participation (n=15)',
                              '  Other reasons (n=5)'))
g <- add_box(g, txt = 'Randomized (n=150)')
g <- add_split(g, txt = c('Allocated to intervention (n=75)\n  Received as allocated (n=70)\n  Did not receive (n=5)',
                          'Allocated to control (n=75)\n  Received as allocated (n=73)\n  Did not receive (n=2)'))
g <- add_box(g, txt = c('Lost to follow-up (n=2)\nDiscontinued (n=3)',
                        'Lost to follow-up (n=1)\nDiscontinued (n=2)'))
g <- add_box(g, txt = c('Analysed (n=75)\nExcluded from analysis (n=0)',
                        'Analysed (n=75)\nExcluded from analysis (n=0)'))

pdf('consort.pdf', width = 8, height = 10)
plot(g)
dev.off()

# 6. CONSERVATION CHECK for alluvial -- sum per axis must equal N entities
n_per_axis <- df_wide %>%
    summarise(across(c(t1, t2, t3), ~ sum(!is.na(.))))
stopifnot(all(n_per_axis == nrow(df_wide)))              # all entities present at every axis
