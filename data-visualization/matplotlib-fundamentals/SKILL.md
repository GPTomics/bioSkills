---
name: bio-data-visualization-matplotlib-fundamentals
description: Build publication-quality figures with matplotlib using the object-oriented Figure/Axes API, constrained_layout, rcParams customization, TrueType (Type-42) font embedding for journal submission, and CVD-safe palettes. Covers seaborn integration, common chart types, axis formatting, and the small gotchas that distinguish reproducible matplotlib from notebook scratch. Use when producing publication figures in Python — RNA-seq scatter, single-cell embeddings, generic biological plotting.
tool_type: python
primary_tool: matplotlib
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, seaborn 0.13+, numpy 1.26+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# matplotlib Fundamentals

**"Make a publication figure in Python"** -> Build via the **object-oriented Figure/Axes API** (not pyplot state-machine), with `constrained_layout` for axes alignment, `pdf.fonttype=42` for journal-compliant TrueType fonts, CVD-safe palettes, and rasterized point layers for large scatter. The pyplot interface is for notebook scratch; the Figure/Axes API is for reproducible figures.

- Python: `fig, ax = plt.subplots()` -> `ax.scatter` / `ax.plot` / `ax.bar`; `seaborn.objects` (new grammar API) for ggplot-like

## The Three Modern Defaults

1. **Object-oriented API** — `fig, ax = plt.subplots(figsize=(4, 3))` then `ax.scatter(x, y)`, `ax.set_xlabel(...)`. The pyplot state-machine (`plt.scatter`, `plt.xlabel`) hides which axes are being modified and breaks in multi-subplot figures.

2. **constrained_layout** — `plt.subplots(constrained_layout=True)` automatically prevents axis-label clipping and tight-packs subplots. Replaces the older `tight_layout()` and is the default in matplotlib 3.6+.

3. **Type-42 (TrueType) font embedding** — `plt.rcParams['pdf.fonttype']=42` produces searchable/editable PDF text. Default Type-3 PostScript glyphs are not searchable and **rejected by Nature, IEEE, ACM, and many other publishers**.

## Standard Setup for Publication

```python
import matplotlib.pyplot as plt
import matplotlib as mpl

# rcParams for publication compliance
mpl.rcParams.update({
    'pdf.fonttype': 42,                 # TrueType -- searchable PDFs
    'ps.fonttype': 42,                  # TrueType in EPS
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,                     # Nature requires 5-7 pt body text
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'figure.dpi': 100,                  # display
    'savefig.dpi': 300,                 # save
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'lines.linewidth': 1.0,
    'patch.linewidth': 0.5,
})
```

## Figure / Axes API

```python
import matplotlib.pyplot as plt

# Single axes
fig, ax = plt.subplots(figsize=(89/25.4, 70/25.4),       # 89mm x 70mm in inches; Nature single col
                       constrained_layout=True)
ax.scatter(x, y, c='#0072B2', s=10, alpha=0.7, edgecolors='none', rasterized=True)
ax.set_xlabel('PC1 (45%)')
ax.set_ylabel('PC2 (12%)')
ax.spines[['top', 'right']].set_visible(False)
fig.savefig('scatter.pdf')

# Grid of axes
fig, axes = plt.subplots(2, 3, figsize=(180/25.4, 100/25.4),  # 180mm double col
                          constrained_layout=True)
for ax, (label, panel_data) in zip(axes.flat, data.items()):
    ax.plot(panel_data['x'], panel_data['y'])
    ax.set_title(label, fontsize=8)
```

## Common Chart Types

```python
# Scatter -- always rasterized for >1000 points
ax.scatter(x, y, c=values, cmap='viridis', s=8, alpha=0.6,
           edgecolors='none', rasterized=True)
plt.colorbar(ax.collections[0], ax=ax, label='Expression', shrink=0.8)

# Line
ax.plot(x, y1, color='#0072B2', label='Control', linewidth=1)
ax.plot(x, y2, color='#D55E00', label='Treatment', linewidth=1)
ax.fill_between(x, y_low, y_high, color='#0072B2', alpha=0.2)
ax.legend(frameon=False, fontsize=6)

# Bar
ax.bar(categories, values, color='#0072B2', edgecolor='black', linewidth=0.5)

# Box / violin (prefer seaborn for these -- see distribution-plots)
ax.boxplot([group_a, group_b, group_c], labels=['A', 'B', 'C'],
           patch_artist=True, boxprops=dict(facecolor='#0072B2', alpha=0.7))

# Histogram
ax.hist(values, bins=30, color='#0072B2', edgecolor='white', linewidth=0.5)

# Heatmap (prefer seaborn for clustered; see heatmaps-clustering)
im = ax.imshow(matrix, cmap='RdBu_r', aspect='auto', vmin=-vmax, vmax=vmax)
plt.colorbar(im, ax=ax, label='Z-score')
```

## seaborn Integration

```python
import seaborn as sns

# seaborn shares the matplotlib Figure/Axes -- pass ax= argument
fig, ax = plt.subplots(figsize=(4, 3), constrained_layout=True)
sns.scatterplot(data=df, x='log_fold_change', y='neg_log_p',
                hue='significance', palette=['#999999', '#0072B2', '#D55E00'],
                s=10, alpha=0.7, ax=ax, rasterized=True)

# seaborn 0.13+ has the `objects` grammar interface (ggplot-like)
import seaborn.objects as so
(so.Plot(df, x='log_fold_change', y='neg_log_p')
   .add(so.Dots(pointsize=2), color='significance')
   .scale(color=['#999999', '#0072B2', '#D55E00']))
```

**Return-type gotcha:** seaborn axes-level functions (`scatterplot`, `boxplot`, `barplot`) return Axes. Figure-level (`displot`, `relplot`, `catplot`) return FacetGrid — needs `.set_axis_labels(x, y)` not `.set_xlabel(x)`.

## Axis Formatting

```python
# Log scale
ax.set_yscale('log')

# Scientific notation
from matplotlib.ticker import ScalarFormatter
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))

# Date axis
import matplotlib.dates as mdates
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))

# Tick frequency
ax.set_xticks(np.arange(0, 10, 2))
ax.set_xticklabels(['A', 'B', 'C'], rotation=45, ha='right')

# Grid
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
```

## Color and Palette

```python
# CVD-safe categorical
okabe_ito = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']

# Perceptually-uniform sequential (Crameri batlow / viridis cividis)
from cmcrameri import cm as cmc
plt.imshow(data, cmap=cmc.batlow)
plt.imshow(data, cmap='viridis')                          # built-in

# Diverging symmetric for LFC / z-score
vmax = np.quantile(np.abs(data), 0.99)
plt.imshow(data, cmap='RdBu_r', vmin=-vmax, vmax=vmax)   # symmetric
```

See `data-visualization/color-palettes` for full palette decision tree.

## Saving

```python
# PDF for vector text + raster scatter (best of both)
fig.savefig('figure.pdf', dpi=300, bbox_inches='tight')

# PNG for raster (web, presentations)
fig.savefig('figure.png', dpi=300, bbox_inches='tight')

# TIFF for some journals
fig.savefig('figure.tiff', dpi=300, pil_kwargs={'compression': 'tiff_lzw'})

# SVG for editable vector
fig.savefig('figure.svg', bbox_inches='tight')
```

## Common Failure Modes

### Default Type-3 fonts rejected by journals

**Trigger:** Default `pdf.fonttype=3` (PostScript Type 3 glyphs as drawing operators).

**Mechanism:** Type-3 glyphs are not searchable or selectable; many journals reject.

**Symptom:** Submission rejected at automated check; "Type 3 fonts not permitted."

**Fix:** `mpl.rcParams['pdf.fonttype']=42` AND `ps.fonttype=42`. Verify with `pdffonts figure.pdf` showing `TrueType`.

### tight_layout fails on complex grids

**Trigger:** `plt.tight_layout()` on a figure with colorbars or shared axes.

**Mechanism:** tight_layout doesn't account for axes added after-the-fact (colorbars).

**Symptom:** Labels clipped; subplots overlap colorbar.

**Fix:** Use `constrained_layout=True` in `plt.subplots()` instead; or `fig.set_constrained_layout(True)` after creation.

### pyplot state-machine in multi-subplot

**Trigger:** `plt.xlabel(...)` after `plt.subplots(2, 3)`.

**Mechanism:** pyplot calls modify the *current* axes — usually the last created. Multi-subplot code becomes order-dependent.

**Symptom:** Wrong subplot gets the label.

**Fix:** Use `ax.set_xlabel(...)` with explicit axes reference.

### Scatter of 100000 points crashes PDF viewer

**Trigger:** Vector scatter at large N; one PDF page becomes 50 MB.

**Mechanism:** Each scatter point is a vector circle.

**Symptom:** PDF takes 30 seconds to open; Illustrator crashes; reviewer files complaint.

**Fix:** `rasterized=True` on the scatter call. Keep axes and text vector.

### seaborn FacetGrid vs Axes return-type confusion

**Trigger:** `g = sns.displot(...)`; calling `g.set_xlabel('x')` fails.

**Mechanism:** displot returns FacetGrid; needs `.set_axis_labels(x, y)` or per-axes iteration.

**Symptom:** AttributeError on .set_xlabel.

**Fix:** Use `set_axis_labels` for FacetGrid; `set_xlabel` for Axes. Switch to axes-level `sns.histplot(ax=ax)` to get Axes-API behavior.

### figsize in inches when mm was intended

**Trigger:** `figsize=(89, 70)` thinking mm; matplotlib expects inches.

**Mechanism:** Default figure unit is inches.

**Symptom:** Figure is 89 inches wide.

**Fix:** Convert: `figsize=(89/25.4, 70/25.4)` for mm input.

### Colorbar over-fills the axes

**Trigger:** Default `plt.colorbar(im, ax=ax)`.

**Mechanism:** Colorbar takes the same height as the axes; on small subplots dominates.

**Symptom:** Subplot looks squished.

**Fix:** `plt.colorbar(im, ax=ax, shrink=0.6, aspect=20)`; or use `make_axes_locatable` for fine control.

### Vector grid + rasterized scatter mixed properly

**Trigger:** Want vector axes + raster scatter; save as PDF.

**Mechanism:** Default rasterization can include axes if not controlled.

**Symptom:** Whole plot rasterized; axis text blurry on zoom.

**Fix:** Per-element `rasterized=True` on scatter only; axes and text stay vector. Set `fig.set_rasterization_zorder(0)` to globally control.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| PDF rejected by journal | Type-3 fonts | `pdf.fonttype=42` |
| Subplots overlap | No constrained_layout | `plt.subplots(constrained_layout=True)` |
| Wrong subplot labeled | pyplot state-machine | Use ax.set_xlabel explicitly |
| 50 MB PDF | Vector scatter at large N | `rasterized=True` on scatter |
| Figure too big | mm interpreted as inches | Divide by 25.4 |
| Colorbar dominates | Default size | `shrink=0.6, aspect=20` |
| seaborn .set_xlabel fails | FacetGrid not Axes | `g.set_axis_labels(x, y)` |
| Axes spine missing | Wrong API | `ax.spines[['top','right']].set_visible(False)` |

## References

- Hunter JD. 2007. Matplotlib: A 2D graphics environment. *Comput Sci Eng* 9(3):90-95.
- Rougier NP, Droettboom M, Bourne PE. 2014. Ten simple rules for better figures. *PLOS Comp Biol* 10(9):e1003833.
- Waskom ML. 2021. seaborn: statistical data visualization. *J Open Source Softw* 6(60):3021.

## Related Skills

- data-visualization/color-palettes - Palette selection
- data-visualization/multipanel-figures - GridSpec and patchwork-equivalent layouts
- data-visualization/distribution-plots - seaborn boxplot/violin/raincloud
- data-visualization/heatmaps-clustering - seaborn.clustermap
- data-visualization/volcano-and-ma-plots - matplotlib scatter for volcano
- reporting/figure-export - DPI / format / journal-spec details
