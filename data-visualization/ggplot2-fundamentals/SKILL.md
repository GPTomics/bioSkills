---
name: bio-data-visualization-ggplot2-fundamentals
description: Build publication-quality figures in R with ggplot2 using the grammar of graphics (data + aesthetics + geometries + scales + facets + themes) with CVD-safe palettes, cairo_pdf TrueType embedding, programmatic aes via tidy evaluation, and the theme_classic publication baseline. Use when producing static figures in R for papers, presentations, or reports.
tool_type: r
primary_tool: ggplot2
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: ggplot2 3.5+, scales 1.3+, ggrepel 0.9.5+, ggtext 0.1.2+, viridis 0.6+, scico 1.5+, patchwork 1.2+ (axes='collect' requires 1.2.0+).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# ggplot2 Fundamentals

**"Build a publication figure in R"** -> Express the figure as **data + aesthetic mappings + one or more geometries + scales + facets + theme**. The grammar of graphics (Wilkinson 2005; Wickham 2010 *J Comput Graph Stat* 19:3) makes each visual element separately addressable — change scales without rewriting geoms; swap geom_point for geom_violin without touching aesthetics.

- R: `ggplot(data, aes(x, y)) + geom_point() + scale_color_manual(...) + theme_classic()`
- Programmatic: `aes(x = .data[[var]])` for tidy-eval; `!!sym(var)` for older base R style

## The Three Modern Defaults

1. **theme_classic() + remove panel grid + Okabe-Ito palette** as the publication baseline. `theme_minimal` adds light gridlines; `theme_bw` adds a panel border; both work but `theme_classic` is the cleanest for journals.

2. **cairo_pdf for export** — `ggsave('out.pdf', device = cairo_pdf)` embeds TrueType fonts (searchable PDFs); default `ggsave('.pdf')` uses pdf() which produces journal-incompatible fonts on some systems.

3. **Tidy evaluation for programmatic aes** — `aes(x = .data[[var]])` is the modern idiom (ggplot2 3.0+); the older `aes_string(x = var)` is deprecated. For dplyr-style symbol evaluation, use `!!sym(var)` with `aes(x = !!sym(var))`.

## Grammar in Layers

```r
library(ggplot2)

# data + aes + geom is the minimum
ggplot(df, aes(x = condition, y = expression)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    # scales
    scale_y_continuous(trans = 'log10', labels = scales::label_log()) +
    scale_color_manual(values = c('#0072B2', '#D55E00')) +
    # labels
    labs(x = NULL, y = 'Expression (log10)',
         title = 'Gene X across conditions',
         caption = 'Source: ...') +
    # facets
    facet_wrap(~ tissue, ncol = 3, scales = 'free_y') +
    # theme
    theme_classic(base_size = 10) +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = 'bold'))
```

## Common Geoms

```r
geom_point(alpha = 0.7, size = 1, rasterize = TRUE)   # rasterize: ggplot2 3.5+ inline OR ggrastr::rasterize()
geom_line(linewidth = 0.5)                             # linewidth replaces size for lines (ggplot2 3.4+)
geom_col()                                              # bar with y values (use this; geom_bar(stat='identity') is older)
geom_bar()                                              # bar with counts
geom_boxplot(outlier.shape = NA)                       # always suppress when overlaying jitter
geom_violin(bw = 'SJ', trim = FALSE)                   # Sheather-Jones bandwidth; show full tails
geom_histogram(bins = 30)                              # bins NOT binwidth for control
geom_density(alpha = 0.5)
geom_tile(aes(fill = z))                               # heatmap building block
geom_text(aes(label = label), check_overlap = TRUE)
geom_text_repel(aes(label = label), max.overlaps = Inf)   # ggrepel; max.overlaps = Inf prevents silent label drops
```

## Aesthetic Mappings

```r
aes(x, y, color, fill, shape, size, alpha, linetype, linewidth, group)

# Color vs fill: color = stroke; fill = interior (boxplot, bar, area, polygon)
# Use both when needed: geom_point(aes(color = group, fill = group), shape = 21)
```

**Constant inside vs mapping inside aes** is a common confusion:
```r
geom_point(color = 'red')             # constant: every point red
geom_point(aes(color = group))        # mapping: color varies with group
```

## Scales

```r
# Continuous
scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2),
                    labels = scales::label_number(scale = 1e-6, suffix = 'M'))
scale_y_log10()
scale_y_continuous(trans = 'sqrt')

# Discrete
scale_x_discrete(limits = c('Control', 'Treatment', 'Vehicle'))   # explicit order
scale_color_manual(values = c(Control = '#0072B2', Treatment = '#D55E00'))

# Colormap (sequential, diverging, cyclic) -- see color-palettes
scale_color_viridis_c(option = 'viridis')
scale_color_scico(palette = 'batlow')                              # Crameri
scale_fill_gradient2(low = '#0072B2', mid = 'white', high = '#D55E00', midpoint = 0)

# Date / time
scale_x_date(date_breaks = '1 year', date_labels = '%Y')
```

## Facets

```r
facet_wrap(~ var, ncol = 3, scales = 'free_y')
facet_grid(rows = vars(condition), cols = vars(timepoint), scales = 'free_x')
facet_grid(condition ~ timepoint)                                  # formula syntax
```

`scales = 'free_y'` lets each panel have its own y-range — appropriate when biological scales differ across facets. `scales = 'fixed'` (default) is the right choice when comparing across panels.

## Theme

```r
# Publication baseline
theme_pub <- theme_classic(base_size = 10) +
    theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black', linewidth = 0.3),
        axis.line = element_line(color = 'black', linewidth = 0.3),
        legend.position = 'right',
        legend.key.size = unit(0.4, 'cm'),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold', size = 9),
        plot.title = element_text(face = 'bold', size = 11),
        plot.tag = element_text(face = 'bold', size = 11))

# Save as a function for re-use across project
```

## Programmatic Plots (Tidy Evaluation)

```r
# Pass variable name as a string
plot_var <- function(df, x_var, y_var) {
    ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]])) +
        geom_point()
}
plot_var(df, 'PC1', 'PC2')

# Alternative: bare names via embracing
plot_var2 <- function(df, x_var, y_var) {
    ggplot(df, aes(x = {{ x_var }}, y = {{ y_var }})) +
        geom_point()
}
plot_var2(df, PC1, PC2)
```

`aes_string` is deprecated as of ggplot2 3.0. `.data[[var]]` is the modern programmatic idiom.

## Labels with ggtext (rich-text)

```r
library(ggtext)
ggplot(df, aes(x, y)) + geom_point() +
    labs(x = 'log<sub>2</sub> fold change',
         y = '\\u2212log<sub>10</sub>(*p*)') +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
```

ggtext renders inline HTML / Markdown in titles, captions, axis labels — much better than `expression(...)` for italics + subscripts + special characters.

## Saving — TrueType Embedding

```r
# cairo_pdf for TrueType embedded; portable across systems
ggsave('figure.pdf', plot = p,
       width = 89, height = 70, units = 'mm',
       device = cairo_pdf)

# Vector + raster mix via ggrastr (for large scatter)
library(ggrastr)
ggplot(df, aes(x, y)) +
    rasterise(geom_point(alpha = 0.5), dpi = 300) +
    theme_pub
ggsave('out.pdf', device = cairo_pdf)

# PNG for raster
ggsave('figure.png', p, width = 89, height = 70, units = 'mm', dpi = 300)

# TIFF for some journals
ggsave('figure.tiff', p, width = 89, height = 70, units = 'mm', dpi = 300,
       compression = 'lzw')
```

## Common Failure Modes

### Default ggsave fonts not embedded

**Trigger:** `ggsave('out.pdf', p)` without `device = cairo_pdf`.

**Mechanism:** Default pdf() device on some systems produces non-embedded fonts.

**Symptom:** Reviewer or coauthor opens PDF; text renders in wrong font; journal rejects.

**Fix:** Always `device = cairo_pdf` for PDF saves.

### Mapping vs constant aesthetic confusion

**Trigger:** `geom_point(aes(color = 'red'))` — string 'red' becomes a categorical mapping.

**Mechanism:** `aes()` interprets its arguments as variables; 'red' becomes a 1-level factor and gets mapped to the FIRST default color.

**Symptom:** Points appear blue (or whatever default) with a legend showing "red" as a category.

**Fix:** Move outside aes: `geom_point(color = 'red')` for a constant; keep inside for a mapping.

### linewidth vs size for lines

**Trigger:** `geom_line(size = 0.5)` in ggplot2 3.4+.

**Mechanism:** ggplot2 3.4+ renamed line-width control from `size` to `linewidth`; `size` still works for points.

**Symptom:** Warning "Using `size` aesthetic for lines was deprecated"; lines render but warning.

**Fix:** `geom_line(linewidth = 0.5)`. `geom_point(size = 1)` is correct.

### facet_wrap scales = 'free' confuses cross-panel comparison

**Trigger:** `facet_wrap(~ var, scales = 'free')` for figures intended to compare across panels.

**Mechanism:** Each panel has its own scale; visual comparison invalid.

**Symptom:** Reviewer asks "why are these heights different?"

**Fix:** Use `scales = 'fixed'` (default) when cross-panel comparison matters; use `'free_y'` only when panels are inherently different scales.

### aes_string deprecated

**Trigger:** `aes_string(x = 'PC1', y = 'PC2')` for programmatic plotting.

**Mechanism:** Deprecated since ggplot2 3.0; emits warning.

**Symptom:** Deprecation warning in script log.

**Fix:** `aes(x = .data[['PC1']], y = .data[['PC2']])` OR `aes(x = !!sym(x_var))`.

### ggrepel max.overlaps default drops labels

**Trigger:** `geom_text_repel(aes(label = label))` with N > 10 labels.

**Mechanism:** Default `max.overlaps = 10`; labels exceeding this are silently dropped with a warning.

**Symptom:** Some labeled genes are silently missing; warning buried in log.

**Fix:** `geom_text_repel(aes(label = label), max.overlaps = Inf)` OR `options(ggrepel.max.overlaps = Inf)` at script top.

### Saving with size in inches but intended mm

**Trigger:** `ggsave('out.pdf', p, width = 89, height = 70)` thinking mm.

**Mechanism:** Default `units = 'in'`.

**Symptom:** Figure is 89 inches wide — too large to open in Illustrator.

**Fix:** `units = 'mm'` explicit. Nature single column = 89mm; double column = 183mm.

## References

- Wickham H. 2016. *ggplot2: Elegant Graphics for Data Analysis* (2nd ed). Springer.
- Wickham H. 2010. A layered grammar of graphics. *J Comput Graph Stat* 19(1):3-28.
- Wilkinson L. 2005. *The Grammar of Graphics* (2nd ed). Springer.

## Related Skills

- data-visualization/color-palettes - Scale_color/_fill palette selection
- data-visualization/multipanel-figures - patchwork composition
- data-visualization/distribution-plots - Box / violin / raincloud geoms
- data-visualization/volcano-and-ma-plots - ggplot2 volcano with ggrepel
- data-visualization/heatmaps-clustering - ComplexHeatmap and ggplot2 geom_tile
