---
name: bio-reporting-figure-export
description: Exports publication-ready figures in various formats with proper resolution, sizing, and typography. Use when preparing figures for journal submission, creating vector graphics for presentations, or ensuring consistent figure styling across analyses.
tool_type: mixed
primary_tool: matplotlib
---

## Version Compatibility

Reference examples tested with: ggplot2 3.5+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Publication-Ready Figure Export

**"Export figures for publication"** → Save plots as high-resolution PDF/SVG/TIFF with journal-required DPI and dimensions.
- Python: `fig.savefig('fig.pdf', dpi=300, bbox_inches='tight')`
- R: `ggsave('fig.pdf', width=7, height=5, units='in')`

## Python (matplotlib)

```python
import matplotlib.pyplot as plt

# Set publication defaults
plt.rcParams.update({
    'font.size': 8,
    'font.family': 'Arial',
    'axes.linewidth': 0.5,
    'lines.linewidth': 1,
    'figure.dpi': 300
})

fig, ax = plt.subplots(figsize=(3.5, 3))  # Single column width
# ... create plot ...

# Save in multiple formats
fig.savefig('figure1.pdf', bbox_inches='tight', dpi=300)
fig.savefig('figure1.png', bbox_inches='tight', dpi=300)
fig.savefig('figure1.svg', bbox_inches='tight')
```

## R (ggplot2)

```r
library(ggplot2)

p <- ggplot(data, aes(x, y)) + geom_point() +
  theme_classic(base_size = 8) +
  theme(text = element_text(family = 'Arial'))

# PDF for vector graphics
ggsave('figure1.pdf', p, width = 3.5, height = 3, units = 'in')

# High-res PNG
ggsave('figure1.png', p, width = 3.5, height = 3, units = 'in', dpi = 300)

# TIFF (some journals require)
ggsave('figure1.tiff', p, width = 3.5, height = 3, units = 'in',
       dpi = 300, compression = 'lzw')
```

## Journal Requirements

| Journal Type | Format | Resolution | Width |
|--------------|--------|------------|-------|
| Most journals | PDF/EPS | Vector | 3.5" (1-col), 7" (2-col) |
| Online-only | PNG | 300 DPI | Variable |
| Print | TIFF | 300-600 DPI | Column width |

## Multi-panel Figures

```python
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(7, 5))  # Two-column width
gs = GridSpec(2, 3, figure=fig)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1:])
ax3 = fig.add_subplot(gs[1, :])

# Add panel labels
for ax, label in zip([ax1, ax2, ax3], ['A', 'B', 'C']):
    ax.text(-0.1, 1.1, label, transform=ax.transAxes,
            fontsize=10, fontweight='bold')

fig.savefig('figure_multipanel.pdf', bbox_inches='tight')
```

## Color Considerations

- Use colorblind-friendly palettes (viridis, cividis)
- Ensure sufficient contrast for grayscale printing
- Maintain consistency across all figures

## Related Skills

- data-visualization/ggplot2-fundamentals - Creating plots in R
- data-visualization/heatmaps-clustering - Complex visualizations
- data-visualization/multipanel-figures - Figure composition
