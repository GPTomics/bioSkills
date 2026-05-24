---
name: bio-data-visualization-multipanel-figures
description: Compose multi-panel publication figures with patchwork, cowplot, gridExtra (R), or matplotlib GridSpec/subfigures (Python) including shared axes/legends/guides collection, panel labels in Nature/Cell convention, and journal-spec sizing. Covers patchwork ≥1.2.0 axes='collect' feature, Type-42 font embedding, and the cairo_pdf save path. Use when composing 2+ subpanels into a single figure for journal submission.
tool_type: mixed
primary_tool: patchwork
---

## Version Compatibility

Reference examples tested with: patchwork 1.2+ (axes='collect' requires this version, released 2024-01-05), cowplot 1.1+, ggplot2 3.5+, matplotlib 3.8+ (subfigures stable since 3.4).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Multi-Panel Figures

**"Combine plots into a multi-panel figure"** -> Arrange individual plots into a single composed figure with consistent sizing, shared legends/axes, and panel labels (a, b, c) in the Nature/Cell convention. The decision space: which composition library (patchwork most modern in R; matplotlib subfigures in Python), how to share legends and axes, and how to size at journal specifications.

- R: `patchwork` (modern; supports axes/guides collection since 1.2), `cowplot` (older; align_plots), `gridExtra` (basic grid arrange)
- Python: `matplotlib.gridspec.GridSpec`, `fig.subfigures()` (matplotlib 3.4+)

## The Single Most Important Modern Insight -- Axes Collection Requires patchwork ≥ 1.2.0

patchwork 1.2.0 (released 2024-01-05) added `axes = 'collect'` and `axis_titles = 'collect'` to `plot_layout()`. These collect repeated axes / titles across subplots into a single shared axis label — the same way `guides = 'collect'` (available since patchwork 1.0) collects legends.

Without this, multi-panel figures with shared axes show redundant labels on every subplot (visually cluttered AND non-Nature compliant). Verify patchwork version is ≥ 1.2.0; older versions silently ignore the `axes` argument.

## patchwork -- Modern R Composition

**Goal:** Compose 4 ggplot objects into a 2×2 panel figure with shared legend, collected axes, and bold panel labels (a, b, c, d) in upper-left of each subplot.

**Approach:** Combine plots with `+`, `/`, `|` operators; apply `plot_layout(guides='collect', axes='collect')` for shared elements; add `plot_annotation(tag_levels='a')` for Nature-style panel labels.

```r
library(patchwork)
library(ggplot2)

p1 <- ggplot(df, aes(x, y)) + geom_point() + theme_classic()
p2 <- ggplot(df, aes(group, value)) + geom_boxplot() + theme_classic()
p3 <- ggplot(df, aes(x)) + geom_histogram() + theme_classic()
p4 <- ggplot(df, aes(x, y, color = group)) + geom_point() + theme_classic()

# 2x2 grid
fig <- (p1 + p2) / (p3 + p4) +
    plot_annotation(tag_levels = 'a',
                    theme = theme(plot.tag = element_text(face = 'bold', size = 10))) +
    plot_layout(guides = 'collect',         # share legends
                axes = 'collect',           # share axes (patchwork >= 1.2.0)
                axis_titles = 'collect')

ggsave('figure1.pdf', fig, width = 180, height = 140, units = 'mm', device = cairo_pdf)
```

## patchwork Operators

```r
p1 + p2                                     # side-by-side
p1 / p2                                     # vertical stack
(p1 | p2) / p3                              # mixed: top row two, bottom one
p1 + p2 + p3 + plot_layout(ncol = 3)
p1 + p2 + plot_layout(widths = c(2, 1))     # 2:1 width ratio

# Complex grid via design string
design <- "
AAB
AAB
CCC
"
p1 + p2 + p3 + plot_layout(design = design)

# Inset
p1 + inset_element(p2, left = 0.6, bottom = 0.6, right = 1, top = 1)
```

## cowplot -- Alternative with Alignment Focus

```r
library(cowplot)

# plot_grid is the workhorse
combined <- plot_grid(p1, p2, p3, p4,
                       ncol = 2, labels = 'AUTO',         # 'AUTO' = A, B, C, D
                       label_size = 12, label_fontface = 'bold',
                       align = 'hv',                       # align horizontally + vertically
                       rel_widths = c(1, 1), rel_heights = c(1, 1))

# Nested grids
top_row <- plot_grid(p1, p2, ncol = 2, labels = c('A', 'B'))
bottom <- plot_grid(p3, p4, ncol = 2, labels = c('C', 'D'))
combined <- plot_grid(top_row, bottom, nrow = 2, rel_heights = c(1, 1.2))

ggsave('figure.pdf', combined, width = 180, height = 140, units = 'mm', device = cairo_pdf)
```

cowplot is older but its alignment behavior is sometimes more reliable than patchwork on edge cases (axes-with-titles of different lengths).

## matplotlib GridSpec (Python)

```python
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(180/25.4, 120/25.4), constrained_layout=True)
gs = GridSpec(2, 3, figure=fig)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1:])         # top right, spans columns 1-2
ax3 = fig.add_subplot(gs[1, :])           # bottom row, spans all columns

ax1.scatter(x, y, s=4, rasterized=True)
ax2.plot(x, y)
ax3.bar(cats, vals)

# Panel labels at (-0.15, 1.05) of each axes
for ax, lbl in zip([ax1, ax2, ax3], 'abc'):
    ax.text(-0.15, 1.05, lbl, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top')

fig.savefig('figure.pdf', dpi=300, bbox_inches='tight')
```

## matplotlib Subfigures

```python
fig = plt.figure(figsize=(180/25.4, 120/25.4), constrained_layout=True)
subfigs = fig.subfigures(1, 2, width_ratios=[2, 1])

# Left subfigure has 2 stacked panels
axs_left = subfigs[0].subplots(2, 1)
axs_left[0].plot(x, y)
axs_left[1].scatter(x, y, rasterized=True)

# Right subfigure has one panel
ax_right = subfigs[1].subplots(1, 1)
ax_right.imshow(matrix)
subfigs[1].colorbar(ax_right.images[0], ax=ax_right, shrink=0.5)
```

Subfigures are stronger than GridSpec for complex compositions because each subfigure has its own constrained_layout.

## Journal Sizing

| Journal | Single col | Double col | Max height |
|---------|------------|------------|------------|
| Nature | 89 mm | 183 mm | 247 mm |
| Cell | 85 mm | 174 mm | 235 mm |
| Science | 55 mm | 120 mm | 220 mm |
| PNAS | 87 mm | 178 mm | 225 mm |
| eLife | 86 mm | 175 mm | ~240 mm |

Always set explicit units in mm; default inches is the most common source of "figure too large" errors.

## Panel Labels — Nature/Cell Convention

- **Nature**: lowercase bold serif (a, b, c) in upper-left corner of each panel; 8 pt
- **Cell**: uppercase bold sans-serif (A, B, C); placed flush left at panel top
- **Science**: capital bold (A, B, C)

```r
# patchwork tag_levels for lowercase (Nature)
plot_annotation(tag_levels = 'a',
                theme = theme(plot.tag = element_text(face = 'bold', size = 9)))
# 'A' for uppercase (Cell)
plot_annotation(tag_levels = 'A')
# 'i' for roman numerals (sometimes for sub-panels)
```

```r
# cowplot
plot_grid(..., labels = 'AUTO')   # auto uppercase A, B, C
plot_grid(..., labels = 'auto')   # auto lowercase a, b, c
```

## Per-Method Failure Modes

### patchwork axes='collect' silently ignored

**Trigger:** Using `plot_layout(axes='collect')` with patchwork < 1.2.0.

**Mechanism:** Older versions silently accept the argument but don't act on it.

**Symptom:** Redundant axes on each subplot; no warning or error.

**Fix:** `packageVersion('patchwork')` must be ≥ 1.2.0. Update with `install.packages('patchwork')`.

### Default ggsave produces non-portable PDF

**Trigger:** `ggsave('out.pdf', fig)` without `device = cairo_pdf`.

**Mechanism:** Default pdf() device produces fonts that journals reject on some systems.

**Symptom:** Submission rejected at automated check; "non-embedded fonts."

**Fix:** Always `device = cairo_pdf`.

### Figure dimensions in inches when mm intended

**Trigger:** `ggsave('out.pdf', fig, width = 180, height = 140)`.

**Mechanism:** Default `units = 'in'`.

**Symptom:** Figure file rejected for being 180 × 140 inches.

**Fix:** Explicit `units = 'mm'`.

### Panel labels not aligned to panel content

**Trigger:** patchwork `plot_annotation(tag_levels)` with subplots of different y-axis label widths.

**Mechanism:** Tag is positioned relative to the plot canvas, including the y-axis label area.

**Symptom:** Labels are at different horizontal positions in each panel.

**Fix:** Either standardize y-label widths (pad with whitespace) OR move tags inside the plotting area: `theme(plot.tag.position = c(0.02, 0.98))`.

### cowplot align='v' fails on plots of different widths

**Trigger:** `plot_grid(p_wide, p_narrow, align = 'v')`.

**Mechanism:** Vertical alignment requires same x-axis widths.

**Symptom:** Plots align at the y-axis but x-axis labels are offset.

**Fix:** Use `align = 'hv'` if both alignments needed; otherwise patchwork's `axes='collect'` handles this more gracefully.

### Shared legend lost in patchwork

**Trigger:** `(p1 + p2) + plot_layout(guides = 'collect')` but p1 and p2 use different scales.

**Mechanism:** `guides='collect'` merges identical guides; different scales produce duplicate (not merged) legends.

**Symptom:** Two legends still appear.

**Fix:** Standardize the scales across subplots (same `scale_color_manual(values=...)`); OR drop one legend via `& theme(legend.position = 'none')` on the redundant plot.

### matplotlib GridSpec with constrained_layout=False

**Trigger:** Older code with `plt.subplots` no constrained_layout; tight_layout fails on colorbars.

**Mechanism:** tight_layout doesn't know about post-hoc colorbars.

**Symptom:** Colorbar overlaps adjacent subplot.

**Fix:** `plt.figure(constrained_layout=True)` and use `fig.add_subplot(gs[...])`. constrained_layout is the default-on choice in matplotlib 3.6+.

## Reconciliation

| Pattern | Cause | Action |
|---------|-------|--------|
| patchwork and cowplot align differently | Different alignment algorithms | Try both; cowplot's `align='hv'` and patchwork's `axes='collect'` rarely produce identical results |
| Panel labels position differs between sessions | Different y-axis label widths | Standardize across panels |
| Shared legend duplicated | Scales differ across subplots | Use identical scales OR drop legend from N-1 panels |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Nature single column | 89 mm | Nature figure guidelines |
| Nature double column | 183 mm | Nature figure guidelines |
| Body text size | 5-7 pt | Nature rejects outside range |
| Panel label size | 8 pt bold | Nature convention |
| patchwork axes='collect' minimum version | 1.2.0 (2024-01-05) | patchwork release notes |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Redundant axis labels per panel | patchwork < 1.2.0 OR axes='collect' not set | Update + add to plot_layout |
| Non-embedded font rejection | Default ggsave device | `device = cairo_pdf` |
| Figure 180 in × 140 in | Default units = 'in' | `units = 'mm'` |
| Panel tags misaligned | Different y-label widths | Standardize or move tag inside |
| Cowplot vertical alignment fails | Different x-axis widths | Use 'hv' OR switch to patchwork |
| Two legends instead of shared | Scales differ across subplots | Unify scales |
| matplotlib colorbar overlaps subplot | No constrained_layout | constrained_layout=True |

## References

- Pedersen TL. 2024. patchwork: the composer of plots. *CRAN package* (v1.2.0 release notes).
- Wilke CO. 2017. cowplot: streamlined plot theme and plot annotations for ggplot2. *CRAN package*.
- Hunter JD. 2007. Matplotlib: A 2D graphics environment. *Comput Sci Eng* 9(3):90-95.

## Related Skills

- data-visualization/ggplot2-fundamentals - Individual ggplot objects
- data-visualization/matplotlib-fundamentals - Python equivalent
- reporting/figure-export - DPI / format / journal-spec compliance
- data-visualization/color-palettes - Consistent palette across subpanels
