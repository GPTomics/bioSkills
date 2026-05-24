# Reference: scico 1.5+, viridis 0.6+, RColorBrewer 1.1+, colorspace 2.1+ | Verify API if version differs

# PhD-level palette selection encoding the four criteria:
# (1) perceptual uniformity, (2) CVD safety, (3) grayscale monotonicity, (4) symmetric bounds for diverging.

library(ggplot2)
library(scico)
library(viridis)
library(colorspace)

# 1. SEQUENTIAL -- Crameri batlow OR viridis cividis
ggplot(df, aes(x, y, fill = expression)) + geom_tile() +
    scale_fill_scico(palette = 'batlow')                          # Crameri 2020
# alt: scale_fill_viridis_c(option = 'cividis')                   # Nuñez 2018, CVD-optimal

# 2. DIVERGING -- Crameri vik or roma; ALWAYS symmetric bounds
vmax <- quantile(abs(df$lfc), 0.99, na.rm = TRUE)                  # robust 1-99% bounds
ggplot(df, aes(x, y, fill = lfc)) + geom_tile() +
    scale_fill_scico(palette = 'vik', midpoint = 0,
                     limits = c(-vmax, vmax), oob = scales::squish)

# 3. CYCLIC (phase, time-of-day, angle)
ggplot(df, aes(x, y, color = phase)) + geom_point() +
    scale_color_scico(palette = 'romaO')                          # cyclic variant of roma

# 4. CATEGORICAL <=8 -- Okabe-Ito (Wong 2011 Nat Methods)
okabe_ito <- c('#E69F00', '#56B4E9', '#009E73', '#F0E442',
               '#0072B2', '#D55E00', '#CC79A7', '#000000')
# Built into R 4.0+: palette.colors(8, 'Okabe-Ito')
ggplot(df, aes(x, y, color = cell_type)) + geom_point() +
    scale_color_manual(values = okabe_ito)

# 5. CATEGORICAL DE convention -- Up/Down/NS
de_colors <- c(Up = '#D55E00', Down = '#0072B2', NS = '#999999')
ggplot(de_df, aes(log2FC, neg_log10_p, color = significance)) + geom_point() +
    scale_color_manual(values = de_colors)

# 6. CVD SIMULATION -- mandatory check
demoplot(scico(8, palette = 'batlow'), type = 'heatmap')           # normal vision
demoplot(deutan(scico(8, palette = 'batlow')), type = 'heatmap')   # deuteranopia simulation
demoplot(protan(scico(8, palette = 'batlow')), type = 'heatmap')   # protanopia simulation
# also: demoplot(tritan(palette), 'heatmap') for the rare tritan deficiency

# 7. GRAYSCALE MONOTONICITY TEST
library(scales)
show_col(scico(10, palette = 'batlow'))                            # full color
show_col(desaturate(scico(10, palette = 'batlow')))                # grayscale equivalent
# If the grayscale gradient is still monotonic, the colormap is luminance-monotonic.

# 8. COMPARE -- the jet / rainbow trap visualized
show_col(rainbow(10))                                              # bad: non-monotonic luminance
show_col(desaturate(rainbow(10)))                                  # confirms banding artifact
show_col(viridis(10))                                              # good: monotonic
show_col(desaturate(viridis(10)))                                  # confirms luminance ramp

# 9. JOURNAL-BRAND PALETTES (stylistic, not accessible)
library(ggsci)
ggplot(df, aes(x, y, color = group)) + geom_point() +
    scale_color_npg()                                              # Nature Publishing Group
# also: scale_color_aaas(), _lancet(), _jama(), _jco(), _nejm()
# NB: not CVD-validated. Use journal brand for house style, Okabe-Ito for accessibility.
