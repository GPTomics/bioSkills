---
name: bio-reporting-figure-export
description: Exports publication-ready figures with the correct vector/raster split, embedded editable fonts, color-space-robust palettes, and journal-correct sizing and resolution in matplotlib and ggplot2. Use when preparing figures for journal submission, exporting a dense single-cell or GWAS plot without producing an unopenable vector file, or fixing fonts and colors that break in print.
tool_type: mixed
primary_tool: matplotlib
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, ggplot2 3.5+, ggrastr 1.0+, ragg 1.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show matplotlib` then `help(matplotlib.figure.Figure.savefig)`; introspect `matplotlib.rcParams` if a key is renamed
- R: `packageVersion('ggplot2')` then `?ggsave`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Publication-Ready Figure Export

**"Export this figure for the journal"** -> Save the plot so each part of it survives production: vector structure stays crisp, the dense data layer stays a reasonable file size, text stays editable, and color survives the print conversion.
- Python: `fig.savefig('fig.pdf', dpi=600)` with publication rcParams set
- R: `ggsave('fig.pdf', width=89, units='mm', device=cairo_pdf)`

## The Load-Bearing Idea: A Figure Is Four Layers

A publication figure is not one object; it is four superimposed layers, and export means giving each the representation it needs:

1. **Vector structure** - axes, ticks, spines, gridlines, fit lines, error bars, annotations. Resolution-independent; must stay vector so the typesetter can scale it to 89 mm without pixelation.
2. **Raster data layer** - the dense part: a scatter with 10^5-10^7 points, a heatmap, a micrograph. Drawing a million points as a million vector circles makes a hundreds-of-MB PDF that crashes Illustrator. This layer wants to be pixels.
3. **Type** - every glyph. Editors need it to stay selectable text, not flattened paths or pixels baked into the raster.
4. **Color encoding** - the data-to-color mapping. A scientific choice (perceptual uniformity, color-vision-deficiency safety, grayscale survival) that also interacts with the RGB->CMYK conversion the journal performs without asking.

The expert move is the **hybrid figure**: rasterize only layer 2, keep layers 1 and 3 vector, embed editable fonts, and pick a colormap that survives CMYK and grayscale. Everything below serves that.

A reproducibility framing: a figure is a pure function of (data, code, theme, font availability). If any of those is unpinned, the figure is not reproducible.

## The Hybrid Figure (rasterize only the dense layer)

The single most important export skill for single-cell (UMAP/tSNE) and GWAS (Manhattan) figures. A fully-vector million-point scatter is unopenable and gets rejected by the typesetter's RIP; rasterizing just the data layer keeps file size sane while axes and text stay crisp at any zoom.

```python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(3.5, 3.0))     # physical inches = the real size control
ax.scatter(x, y, s=2, rasterized=True, zorder=0)   # 10^6 points -> embedded raster
ax.plot(xfit, yfit, color='black', zorder=2)        # stays vector
ax.set_xlabel('UMAP1')                              # stays vector text
fig.savefig('fig.pdf', dpi=600)                     # dpi governs ONLY the rasterized layer
```

In a vector container, `savefig(dpi=...)` sets the resolution of the embedded raster patch and does nothing to the vector parts. `ax.set_rasterization_zorder(0)` rasterizes every artist below a z-order cleanly. In R use `ggrastr::rasterise(geom_point(size=0.3), dpi=600)` (wraps any geom since 0.2.0), keeping `theme_*` vector.

## Fonts: Keep the Text Editable

The most common typesetter complaint about matplotlib PDFs is un-editable text. matplotlib's default `pdf.fonttype` is **3** (Type 3), which embeds glyphs as PostScript procedures that import into Illustrator as ungrouped paths and cannot be re-selected as text. Set **42** (TrueType, wrapped) so text stays selectable and editable. Fix it once, globally:

```python
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42      # TrueType, editable text in PDF
mpl.rcParams['ps.fonttype']  = 42      # same for EPS/PS
mpl.rcParams['svg.fonttype'] = 'none'  # SVG: emit real <text>, reference the font by name
```

Type 42/TrueType/Type 3 fonts are subsetted (only used glyphs embedded); Type 1 are not. With large glyph sets (CJK) Type 42 can bloat the PDF - a real tradeoff. `svg.fonttype='none'` keeps words editable in Inkscape/Illustrator but the viewer must have the font (else it substitutes); the default `'path'` outlines every glyph (portable, uneditable). Default to editable text; only outline if a specific production desk asks.

In R, the base `pdf()` device has weak font handling and inconsistent cross-OS rendering; use `cairo_pdf` (embeds fonts and supports alpha). For raster, `ragg::agg_png()`/`agg_tiff()` render anti-aliased text better than cairo/base devices. Rule of thumb: `showtext` for vector devices, `ragg` for raster.

## Color Space: The Author Works in RGB, Print Is CMYK

Screens are additive RGB; offset print is subtractive CMYK. matplotlib and ggplot2 author in RGB only - there is no honest path to a true CMYK figure from them. The journal's pipeline converts RGB->CMYK, and because the CMYK gamut is smaller than sRGB, **saturated out-of-gamut colors shift**: pure RGB blue (#0000FF) and vivid green/cyan come back muddier and darker on paper. The neon scatter that pops on screen can print gray.

What to do: pull colors slightly off full saturation (they survive conversion better); soft-proof downstream in Illustrator/Photoshop with a CMYK profile if it matters; and submit RGB - Nature, Science, Cell, and PLOS all explicitly want RGB, not CMYK, because their pipeline does the conversion and online is RGB anyway. If a legacy desk demands CMYK, convert downstream with an explicit profile and re-check that nothing shifted; do not fake it with a colorspace flag.

## Transparency: EPS Has No Alpha

EPS/PostScript do not support alpha - matplotlib's PS backend renders partially-transparent artists as opaque, so an alpha-blended overplotted scatter loses its density information on EPS export. PDF and SVG support alpha natively; prefer them when transparency carries meaning. If a journal forces EPS and transparency is needed, rasterize that layer (`rasterized=True`) or re-encode density as hexbin/2D-KDE. `savefig(transparent=True)` makes the background transparent for slide overlays, not for print.

## DPI Is Meaningless for Vector

A vector PDF has no inherent resolution; it renders sharp at any zoom. DPI governs only raster formats (PNG/TIFF) and the rasterized data layer inside a vector file. The real size control is the **physical figure size** in inches/mm - design at the journal's exact column width from the start; rescaling a 300-dpi raster to 200% halves its effective resolution. Font sizes are in points (1 pt = 1/72 in) independent of DPI.

The DPI tiers follow the IMAGE CLASS, because print reproduces tone via halftone dots (follow the target journal's own numbers; these are the common production convention):

| Image class | DPI | Why |
|-------------|-----|-----|
| Halftone / grayscale / color photo | 300 | continuous tone matches typical screen rulings |
| Combination (halftone + line/text) | 500-600 | thin lines and small type must not jag against toned background |
| Line art (pure black/white) | 1000-1200 | hard edges alias badly at low DPI - or keep it vector and DPI is moot |

`savefig.dpi` is the file resolution; `figure.dpi` is the on-screen resolution. Independently of DPI, very thin strokes (below ~0.25 pt / 0.1 mm) can drop out or thicken unpredictably at the printer's RIP even in a vector file - keep hairlines at or above the journal's minimum line weight.

## Format Decision

| Format | Type | Use for | Avoid for |
|--------|------|---------|-----------|
| PDF | vector(+raster) | default for most journals; hybrid figures; alpha works | - |
| EPS | vector(+raster) | legacy journal requirement | anything with alpha (flattened opaque) |
| SVG | vector | web; handoff to Illustrator/Inkscape for editing | final print at some desks (support varies) |
| TIFF (LZW) | raster, lossless | print production when a journal demands raster (Cell, PLOS) | large vector-friendly line figures |
| PNG | raster, lossless | online, previews, slides, README figures | print where vector is accepted |
| JPEG | raster, LOSSY | photographs only | any figure with text/lines/edges (DCT ringing) |

Never JPEG for line/text figures - block compression rings along high-contrast edges (gray halos on text, fringing on thin lines). For TIFF, use LZW (near-universal reader support) for 8-bit figures; use ZIP for 16-bit (LZW can inflate 16-bit files).

## Colormaps Are a Scientific Choice, Not Taste

Perceptually-uniform maps (viridis, cividis, magma, inferno, plasma) are constructed in CAM02-UCS so equal data steps map to equal perceived steps with monotonically increasing lightness - which is exactly why they survive grayscale and avoid inventing false gradients. **cividis** is additionally optimized so viewers with and without red-green color-vision deficiency see nearly the same image. By contrast, jet/rainbow has non-monotonic luminance that invents bright/dark bands the data does not have (false edges at yellow/cyan) and collapses to mush in grayscale - a correctness failure, not an aesthetic one.

- Sequential map for ordered data; diverging map (with a meaningful midpoint) for signed data; a categorical CVD-safe palette for discrete classes - never a continuous rainbow for categories.
- Use the **Okabe-Ito** 8-color Color Universal Design palette for categories. Red-green CVD affects up to ~8% of males (population-dependent), so add **redundant encoding** (shape + color, linetype + color, direct labels) so color is never the sole channel.
- Run the grayscale-photocopy test: convert to grayscale and confirm the figure still reads.

## Reproducible Export

- **Byte-stable PDFs:** matplotlib stamps a `CreationDate` into every PDF, so two identical runs differ byte-for-byte (noisy git diffs). Pass `metadata={'CreationDate': None}` to `savefig`, or set the `SOURCE_DATE_EPOCH` env var, for deterministic output.
- **`bbox_inches='tight'` breaks exact widths.** It recomputes the bounding box from drawn content, so output dimensions depend on tick-label lengths and the renderer's font metrics - the same script on two machines (different fonts) yields different-sized PDFs, and it can clip annotations. For camera-ready figures at an exact 89 mm, design to size with `constrained_layout=True` and save without `bbox_inches='tight'`; if cropping is unavoidable, pair it with explicit `pad_inches`.
- **Font-availability nondeterminism:** Helvetica on a Mac vs DejaVu Sans on CI gives different glyph widths, line breaks, and (with tight bbox) different sizes. Pin the font or accept the default and don't crop-to-content.
- **Headless rendering:** call `matplotlib.use('Agg')` before importing pyplot on a cluster/CI box, or just use the file backends, so no display is required.

## Journal Specs (verify against the target journal at submission)

Specs change and vary by sub-journal; re-pull the target's author-guideline page. Snapshot, June 2026:

| Journal | Widths | Min DPI | Formats | Color |
|---------|--------|---------|---------|-------|
| Nature | 89 mm single / 183 mm double; <=170 mm tall | 300 photo, 600+ line | vector AI/EPS/PDF preferred; TIFF raster | RGB |
| Science | 5.7 / 12.1 / 18.4 cm | >=300 at final size; vector preferred | Illustrator-openable vector; no PowerPoint | RGB (not CMYK) |
| Cell Press | 85 / 114 / 174 mm | 300 color, 500 grayscale, 1000 line | TIFF (LZW) or vector | RGB |
| PLOS | 789-2250 px wide; <=2625 px tall | 300-600 (do not exceed 600) | TIFF or EPS only; flattened, LZW, no alpha/layers | RGB or grayscale, 8-bit; no CMYK |

PLOS is strictest (8-bit RGB/grayscale TIFF, no alpha channel, no layers). Cell's grayscale (500) and line (1000) tiers exceed the generic numbers - follow the journal's own.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Typesetter: "supply editable text" | default Type 3 fonts | `pdf.fonttype=42`, `ps.fonttype=42`, `svg.fonttype='none'` |
| PDF won't open / hundreds of MB | fully-vector dense scatter | rasterize the data layer (`rasterized=True` / `ggrastr::rasterise`) |
| Colors muddy in print | saturated RGB out of CMYK gamut | desaturate slightly; soft-proof; submit RGB |
| Transparency gone on EPS | EPS has no alpha | rasterize that layer, or use PDF/SVG, or hexbin |
| Figure not exactly 89 mm | `bbox_inches='tight'` non-deterministic size | design to size + `constrained_layout`, drop tight bbox |
| Heatmap shows false bands | jet/rainbow non-monotonic luminance | viridis/cividis (sequential), diverging map for signed data |
| Noisy git diff on identical figure | PDF `CreationDate` timestamp | `metadata={'CreationDate': None}` or `SOURCE_DATE_EPOCH` |

## Related Skills

- data-visualization/ggplot2-fundamentals - Building the plots in R
- data-visualization/matplotlib-fundamentals - Building the plots in Python
- data-visualization/multipanel-figures - Composing multi-panel layouts
- data-visualization/color-palettes - Choosing perceptual and CVD-safe palettes
- reporting/publication-tables - The table counterpart to figure export

## References

- Borland D, Taylor RM 2nd. Rainbow Color Map (Still) Considered Harmful. IEEE Comput Graph Appl. 2007;27(2):14-17. doi:10.1109/MCG.2007.323435
- Nuñez JR, Anderton CR, Renslow RS. Optimizing colormaps with consideration for color vision deficiency to enable accurate interpretation of scientific data. PLoS ONE. 2018;13(7):e0199239. doi:10.1371/journal.pone.0199239
- Okabe M, Ito K. Color Universal Design (CUD): how to make figures and presentations friendly to colorblind people. jfly.uni-koeln.de/color/ (8-color CVD-safe palette)
- van der Walt S, Smith N. A Better Default Colormap for Matplotlib. SciPy 2015 (conference talk; viridis/magma/inferno/plasma constructed in CAM02-UCS). bids.github.io/colormap
