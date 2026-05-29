# Reference: flowCore 2.14+, ggplot2 3.5+ | Verify API if version differs
library(flowCore)
library(ggplot2)

# Singlet discrimination uses the FSC-A vs FSC-H DIAGONAL (Area-Height non-proportionality),
# applied as a polygon gate - not an arbitrary percentile cutoff. This script simulates the
# geometry, then gates with flowCore::polygonGate exactly as the SKILL recommends.

# === 1. SIMULATE FSC-A vs FSC-H (in practice: ff <- read.FCS(...)) ===
set.seed(42)
n <- 50000

# singlets: Area proportional to Height (tight diagonal)
fsc_h_s <- rnorm(n * 0.95, 100000, 20000)
fsc_a_s <- fsc_h_s * 1.1 + rnorm(n * 0.95, 0, 5000)

# doublets: ~same Height, elevated Area (deflect ABOVE the diagonal)
n_d <- n * 0.05
fsc_h_d <- rnorm(n_d, 100000, 20000)
fsc_a_d <- fsc_h_d * 1.5 + rnorm(n_d, 20000, 5000)

mat <- cbind('FSC-A' = pmax(c(fsc_a_s, fsc_a_d), 0),
             'FSC-H' = pmax(c(fsc_h_s, fsc_h_d), 0))
ff <- flowFrame(mat)
cat('Total events:', nrow(ff), '\n')

# === 2. DIAGONAL POLYGON SINGLET GATE ===
# matrix dimnames preserve 'FSC-A'/'FSC-H'; data.frame() would mangle them to FSC.A
singlet_gate <- polygonGate(filterId = 'singlets', .gate = matrix(
    c(20000, 15000, 260000, 245000, 260000, 270000, 20000, 35000),
    ncol = 2, byrow = TRUE, dimnames = list(NULL, c('FSC-A', 'FSC-H'))))

singlets <- Subset(ff, singlet_gate)
pct_singlet <- nrow(singlets) / nrow(ff) * 100
cat('Singlets:', nrow(singlets), '(', round(pct_singlet, 1), '%)\n')
cat('Doublet rate:', round(100 - pct_singlet, 1), '%\n')

# === 3. VISUALIZE THE GATE ON THE DIAGONAL ===
df <- as.data.frame(exprs(ff))
verts <- as.data.frame(singlet_gate@boundaries)
p <- ggplot(df, aes(`FSC-H`, `FSC-A`)) +
    geom_point(alpha = 0.15, size = 0.4, color = 'gray40') +
    geom_polygon(data = verts, aes(`FSC-H`, `FSC-A`),
                 fill = NA, color = 'red', linewidth = 0.8) +
    theme_bw() +
    labs(title = 'Singlet gate on the FSC-A vs FSC-H diagonal',
         subtitle = 'doublets deflect above the diagonal (high Area, normal Height)')
ggsave('singlet_gate.png', p, width = 7, height = 6)
cat('Saved singlet_gate.png\n')
