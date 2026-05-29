---
name: bio-data-visualization-color-palettes
description: Select colormaps and qualitative palettes for scientific figures using perceptual-uniformity, color-vision-deficiency safety, and luminance-monotonicity criteria. Covers Crameri scientific colormaps, viridis/cividis/magma, Okabe-Ito categorical, ColorBrewer, and the rainbow/jet critique. Use when choosing palettes for heatmaps, scatter, networks, or any encoding where color carries quantitative or categorical meaning.
tool_type: mixed
primary_tool: viridis
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: viridis 0.6+, RColorBrewer 1.1+, scico 1.5+ (Crameri colormaps in R), khroma 1.12+ (Tol/Crameri palettes in R), matplotlib 3.8+, colorcet 3.0+, ggsci 3.0+, colorspace 2.1+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Color Palettes for Scientific Visualization

**"Pick a color palette"** -> Choose a colormap that (a) is perceptually uniform along the relevant data axis, (b) remains interpretable under common color-vision deficiencies, (c) prints correctly to grayscale, and (d) matches the data type — sequential, diverging, cyclic, or qualitative.

- R: `viridis::viridis`, `scico::scale_color_scico`, `khroma::color`, `RColorBrewer::brewer.pal`
- Python: `matplotlib.colormaps`, `colorcet`, `seaborn.color_palette`, `cmcrameri.cm`

## The Three Modern Standards

1. **Perceptual uniformity** -- equal data steps produce equal perceived color steps. viridis (van der Walt 2015), cividis (Nuñez 2018), and the Crameri family (batlow, roma, vik) are designed for this. Jet, rainbow, and red->green are not.

2. **Color vision deficiency safety** -- ~6% of males have deuteranopia / protanopia (red-green deficiency). cividis was explicitly designed to be near-identical under normal and CVD viewing (Nuñez 2018 *PLOS ONE* 13:e0199239). The Okabe-Ito 8-color qualitative palette (popularized in Wong 2011 *Nat Methods* 8:441) is the CVD-safe categorical default.

3. **Grayscale monotonicity** -- a perceptually-uniform sequential colormap has monotonically increasing luminance. Convert the figure to grayscale; if the order is still readable, the colormap is luminance-monotonic. This is the single most actionable test.

## Palette Type by Data Type

| Data type | Use | Avoid |
|-----------|-----|-------|
| Sequential (expression, coverage, density) | viridis, magma, cividis, batlow, lipari | jet, rainbow, hsv |
| Diverging (log fold change, z-score, signed correlation) | vik, roma, RdBu, BrBG, PiYG | jet, rainbow |
| Cyclic (phase, time-of-day, angle) | romaO, vikO, twilight | linear sequential (wrap creates artifactual jump) |
| Categorical (≤8 groups) | Okabe-Ito (Wong 2011), Tol bright, Dark2 | rainbow with N=20, Set1 if CVD matters |
| Categorical (9-20 groups) | tab20, Paired, Polychrome | too-many categorical hues -- consider faceting |
| Categorical (>20) | None -- reconsider design | More colors will not help |

## The Crameri Scientific Colormaps

Crameri 2020 *Nat Commun* 11:5444 documented the prevalence of misleading palettes (rainbow, red-green) across published science and released a family of perceptually-uniform CVD-safe colormaps via Zenodo (doi:10.5281/zenodo.8409685). Key entries:

| Crameri name | Type | Use case |
|--------------|------|----------|
| `batlow` | sequential | Default jet replacement; runs through dark-blue -> ochre -> light-yellow |
| `lipari` | sequential | Higher-saturation alternative; better for projection |
| `vik` | diverging | Blue -> white -> red equivalent, perceptually uniform |
| `roma` | diverging | Slightly warmer than vik |
| `bam` | diverging | Brown -> white -> green |
| `romaO` | cyclic | Phase, time-of-day, angle data |
| `vikO` | cyclic | Diverging cyclic |

```r
library(scico)
# Sequential
ggplot(df, aes(x, y, fill = value)) + geom_tile() +
    scale_fill_scico(palette = 'batlow')
# Diverging
ggplot(df, aes(x, y, fill = lfc)) + geom_tile() +
    scale_fill_scico(palette = 'vik', midpoint = 0)
```

```python
from cmcrameri import cm
import matplotlib.pyplot as plt
plt.imshow(data, cmap=cm.batlow)         # sequential
plt.imshow(data, cmap=cm.vik, vmin=-vmax, vmax=vmax)   # diverging, symmetric
```

## viridis Family (matplotlib default since 3.0)

```r
library(viridis)
scale_color_viridis_c(option = 'viridis')   # default: dark blue -> yellow
scale_color_viridis_c(option = 'magma')     # black -> red -> yellow
scale_color_viridis_c(option = 'inferno')   # black -> purple -> yellow
scale_color_viridis_c(option = 'plasma')    # purple -> pink -> yellow
scale_color_viridis_c(option = 'cividis')   # CVD-optimized
scale_color_viridis_c(option = 'turbo')     # jet-like but perceptually uniform
```

```python
plt.imshow(data, cmap='viridis')   # 'magma', 'inferno', 'plasma', 'cividis', 'turbo'
```

**cividis is the only viridis-family colormap optimized for CVD.** Use it for any figure intended to remain interpretable under deuteranopia/protanopia.

## Okabe-Ito Categorical Palette (Wong 2011)

The 8-color CVD-safe categorical palette. Memorize the hexes:

```r
okabe_ito <- c(
    '#E69F00',  # orange
    '#56B4E9',  # sky blue
    '#009E73',  # bluish green
    '#F0E442',  # yellow
    '#0072B2',  # blue
    '#D55E00',  # vermilion
    '#CC79A7',  # reddish purple
    '#000000'   # black
)
scale_color_manual(values = okabe_ito)
```

Available as `palette.colors(8, 'Okabe-Ito')` in R 4.0+, `scale_color_manual(values = palette.colors(8, 'Okabe-Ito'))`. In matplotlib, `colorblind` style or manual hex list.

For DE plots, the canonical assignment is Up = `#D55E00` (vermilion), Down = `#0072B2` (blue), NS = `#999999` (grey).

## ColorBrewer (Harrower & Brewer 2003)

```r
library(RColorBrewer)
display.brewer.all()                    # interactive palette browser
display.brewer.all(colorblindFriendly = TRUE)   # CVD-safe subset only
brewer.pal(n = 8, name = 'Dark2')       # qualitative
brewer.pal(n = 9, name = 'YlOrRd')      # sequential
brewer.pal(n = 11, name = 'RdBu')       # diverging
```

ColorBrewer's CVD-safe sequential and diverging palettes are publication-defaults. For qualitative beyond 8 colors, switch to Tol/Polychrome — ColorBrewer qualitative tops out at 12 (Set3).

## Scientific Journal Brand Palettes

```r
library(ggsci)
scale_color_npg()       # Nature Publishing Group
scale_color_aaas()      # Science (AAAS)
scale_color_lancet()    # Lancet
scale_color_jama()      # JAMA
scale_color_jco()       # JCO
scale_color_nejm()      # NEJM
```

These are CVD-imperfect — use journal palettes for stylistic compliance, not for accessibility. Verify by colorblindness simulation (below).

## CVD Simulation -- The Mandatory Check

```r
library(colorspace)
# Simulate deuteranopia / protanopia on a palette
cvd_emulator(palette, type = 'deutan')
cvd_emulator(palette, type = 'protan')
cvd_emulator(palette, type = 'tritan')

# Visual side-by-side
demoplot(palette, type = 'heatmap')
```

```python
# colorspacious provides CVD simulation
from colorspacious import cspace_converter
# or use a CVD-safe palette by construction (cividis, Okabe-Ito, Crameri)
```

If a palette is unreadable under deutan simulation, do not use it for accessible figures. Period.

## Grayscale Monotonicity Test

```r
library(scales)
show_col(viridis(10))           # full color
show_col(grey(seq(0, 1, length = 10)))   # equivalent grayscale gradient
```

In practice: save the figure as PNG, open in an image editor, desaturate. If the data order is still readable, the colormap is luminance-monotonic. If it shows arbitrary "rings" or "bands," the colormap is non-monotonic — fix before submitting.

Rainbow / jet fails this test catastrophically. viridis and cividis pass.

## Diverging Palette Setup (LFC, z-score)

```r
library(circlize)
col_fun <- colorRamp2(c(-2, 0, 2), c('#0072B2', 'white', '#D55E00'))
# Symmetric around 0; ALWAYS use symmetric bounds for signed data
```

```python
import matplotlib.pyplot as plt
plt.imshow(data, cmap='RdBu_r', vmin=-2, vmax=2)   # symmetric
# do NOT use vmin=data.min(), vmax=data.max() for diverging data
```

The most common diverging-palette error is asymmetric bounds (`vmin=min, vmax=max`) which mis-aligns zero with the white midpoint.

## Custom Palette Construction

```r
# Discrete categorical
my_palette <- c('Control' = '#0072B2', 'Treatment' = '#D55E00', 'Vehicle' = '#009E73')
scale_color_manual(values = my_palette)

# Continuous gradient between custom colors
colorRampPalette(c('#0072B2', 'white', '#D55E00'))(100)
```

```python
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('cvd_div', ['#0072B2', '#FFFFFF', '#D55E00'])
```

When building a custom diverging palette: pick endpoints with similar luminance (so neither side dominates), pass through pure white at the midpoint (NOT light gray), and verify with the grayscale test.

## Common Failure Modes

### Asymmetric bounds on diverging data

**Trigger:** `vmin=data.min()`, `vmax=data.max()` on signed data with skewed distribution.

**Mechanism:** Zero no longer maps to the midpoint (white) of the diverging palette.

**Symptom:** Half the cells visually look "below zero" but are actually positive; reviewer confusion.

**Fix:** `vmax = max(abs(data.min()), abs(data.max()))`; then `vmin = -vmax`. Or pre-clip data to a fixed range.

### Categorical palette with too many colors

**Trigger:** 15+ groups all on one colormap.

**Mechanism:** Human color discrimination saturates around 8-10 distinct hues.

**Symptom:** Groups look identical; legend has no information value.

**Fix:** Facet by category, or aggregate small groups into "Other," or use a categorical+marker-shape combination.

### Rainbow / jet still in use

**Trigger:** Default colormaps in older matplotlib (<2.0), MATLAB-derived code, or `colorRampPalette(rainbow(...))`.

**Mechanism:** Rainbow has non-monotonic luminance and includes a perceptual "yellow band" that creates artifactual boundaries.

**Symptom:** Figures show banding that doesn't exist in the data; CVD viewers cannot interpret.

**Fix:** Migrate to viridis (sequential) or vik/roma (diverging). For nostalgic jet-like appearance with perceptual properties, use `turbo` (matplotlib 3.3+).

### Light gray midpoint instead of pure white

**Trigger:** `colorRamp2(c(-2, 0, 2), c('blue', '#EEEEEE', 'red'))`.

**Mechanism:** Light gray reads as "weakly significant" rather than zero — the visual "where is zero" anchor is lost.

**Symptom:** Zero values appear muted, drawing the eye away from the actual midpoint.

**Fix:** Use pure white `'#FFFFFF'` or `'white'` at the midpoint.

### Wrong cmap for non-numeric data

**Trigger:** Continuous colormap applied to a categorical variable (e.g., cluster ID as a continuous gradient).

**Mechanism:** Cluster IDs are nominal — ordering is meaningless; gradient implies false ordering.

**Symptom:** Cluster 2 "looks closer to" cluster 1 than cluster 8, but the cluster numbering is arbitrary.

**Fix:** Use a qualitative categorical palette (Okabe-Ito ≤8; tab20 for more).

### CVD-unsafe palette in a CVD-sensitive figure

**Trigger:** Red-green Set1 in a clinical figure intended for broad audience.

**Mechanism:** ~6% of male readers cannot distinguish red from green.

**Symptom:** Reviewer or colleague reports the figure is unreadable.

**Fix:** Pre-flight with `colorspace::cvd_emulator`; switch to Okabe-Ito for categorical, cividis for sequential.

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| Max distinguishable categorical hues | 8-10 | Wong 2011; perceptual research |
| CVD prevalence (males of European descent) | ~6% deutan/protan, ~0.5% tritan | Various; Nuñez 2018 cites figures |
| Diverging midpoint | pure white (#FFFFFF) | Not light gray; preserves zero anchor |
| Crameri batlow / lipari -- general-purpose sequential | – | Crameri 2020 |
| cividis -- CVD-optimal sequential | – | Nuñez 2018 |
| Okabe-Ito -- 8-color qualitative CVD-safe | – | Wong 2011 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Diverging plot with zero not at white | Asymmetric bounds | Use symmetric `vmin = -vmax` |
| Rainbow "bands" visible in heatmap | Non-monotonic luminance of rainbow | Replace with viridis or turbo |
| Categorical plot with indistinguishable groups | Too many hues | Facet, aggregate, or shape+color |
| CVD viewer reports unreadable figure | Red-green palette | Switch to Okabe-Ito or cividis |
| Light gray at diverging midpoint | Wrong center color | Use pure white |
| Grayscale conversion shows banding | Non-luminance-monotonic colormap | Use viridis family or Crameri |
| Heatmap with one cell saturating the scale | No quantile clipping | See data-visualization/heatmaps-clustering for robust bounds |

## References

- Crameri F, Shephard GE, Heron PJ. 2020. The misuse of colour in science communication. *Nat Commun* 11:5444. doi:10.1038/s41467-020-19160-7
- Harrower M, Brewer CA. 2003. ColorBrewer.org: an online tool for selecting colour schemes for maps. *Cartogr J* 40(1):27-37.
- Nuñez JR, Anderton CR, Renslow RS. 2018. Optimizing colormaps with consideration for color vision deficiency to enable accurate interpretation of scientific data. *PLOS ONE* 13(7):e0199239.
- Borland D, Taylor RM II. 2007. Rainbow color map (still) considered harmful. *IEEE Comput Graph Appl* 27(2):14-17.
- Light A, Bartlein PJ. 2004. The end of the rainbow? Color schemes for improved data graphics. *Eos Trans AGU* 85(40):385,391.
- Wong B. 2010. Points of view: Color coding. *Nat Methods* 7(8):573.
- Wong B. 2011. Points of view: Color blindness. *Nat Methods* 8(6):441.
- Gehlenborg N, Wong B. 2012. Points of view: Mapping quantitative data to color. *Nat Methods* 9(8):769.

## Related Skills

- data-visualization/heatmaps-clustering - Robust diverging bounds for heatmaps
- data-visualization/volcano-and-ma-plots - Okabe-Ito Up/Down/NS conventions
- data-visualization/ggplot2-fundamentals - Applying palettes in ggplot2 scales
- data-visualization/dimensionality-reduction-plots - Categorical palette for cluster labels
