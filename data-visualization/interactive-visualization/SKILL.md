---
name: bio-data-visualization-interactive-visualization
description: Build interactive HTML/web visualizations with plotly (Python/R), bokeh (Python), and gganimate/plotly frames for animation, with awareness of current Kaleido static-export model (post-orca-EOL), HTML file-size bloat, and the limits of interactive-only output for journal submission. Use when producing zoomable/hoverable plots for notebook EDA, supplementary HTML, dashboards, or animated time-course / iteration visualizations.
tool_type: mixed
primary_tool: plotly
---

## Version Compatibility

Reference examples tested with: plotly 5.24+, plotly R 4.10+, bokeh 3.4+, kaleido 1.0+ (note: v1 dropped bundled Chrome), gganimate 1.0.9+, altair 5.4+, htmlwidgets 1.6+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)`
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Interactive Visualization

**"Build an interactive plot"** -> Render a zoomable, hoverable, panable HTML/web visualization, knowing that interactive output is a SUPPLEMENT to (not replacement for) the static figure needed for journal submission. Choose plotly for fastest onboarding and ggplot2 conversion (`ggplotly`); bokeh for streaming/server-side; altair for grammar-of-graphics; D3.js for full custom.

- Python: `plotly.graph_objects`, `plotly.express`, `bokeh`, `altair`
- R: `plotly` (via `ggplotly`), `htmlwidgets` ecosystem (leaflet, networkD3, DT)

## The Single Most Important Modern Insight -- Kaleido v1 and the Static-Export Pipeline

Interactive plots produce HTML, but journals need static PDF/PNG. The plotly static-export pipeline changed materially in 2025:

- **Orca is end-of-life** (deprecated 2021, removed pipeline 2025)
- **`fig.write_image(..., engine='orca')` removed in plotly 6.2** (post-Sept 2025)
- **Kaleido v1+ is the current standard** — pass no `engine=` argument
- **Kaleido v1 dropped bundled Chrome** — requires installed Chrome / Chromium
- **EPS export removed in Kaleido v1** (was supported via orca's bundled Chromium)

For static export of plotly figures in 2026: `pip install kaleido`; verify Chrome installed; `fig.write_image('out.pdf')`. Test by writing to a known path and inspecting file size; silent failure on missing Chrome was a 2024-2025 pain point that v1 partially addresses with clearer errors.

## Interactive vs Static — The Reproducibility Cost

Interactive HTML has hidden trade-offs:

- **File size**: a 5000-point plotly HTML is 3-5 MB (embedded JS bundle). 50000 points crashes browsers without WebGL acceleration.
- **Non-citable**: a paper figure must be static. Always export static alongside.
- **Browser version drift**: HTML from 2020 plotly may not render in 2026 browsers.
- **Cannot be alt-text described**: accessibility weaker than static.

Use interactive for notebooks (exploration), supplementary HTML (online journal supplement), dashboards (Streamlit/Dash/Shiny). For the journal figure, always also produce static.

## plotly (Python) — Standard Interactive

**Goal:** Build an interactive HTML plot with zoom, pan, and hover-tooltip behavior; export both interactive HTML for supplements and static PDF for the journal figure.

**Approach:** Use `plotly.express` for declarative high-level plots OR `graph_objects` for fine control; enable WebGL via `render_mode='webgl'` or `Scattergl` for >5000 points; export HTML with `write_html()` and static with `write_image()` after installing Kaleido v1+ and Chrome.

```python
import plotly.express as px
import plotly.graph_objects as go

# Express: high-level, declarative
fig = px.scatter(df, x='PC1', y='PC2', color='cluster',
                  hover_data=['gene_count', 'sample_id'],
                  color_discrete_sequence=['#0072B2', '#D55E00', '#009E73'],
                  title='PCA')
fig.update_layout(template='plotly_white', width=600, height=500)

# WebGL acceleration for >5000 points
fig = px.scatter(df, x='PC1', y='PC2', color='cluster', render_mode='webgl')

# Save
fig.write_html('pca.html')
fig.write_image('pca.pdf')                   # requires kaleido + Chrome

# Graph_objects: low-level
fig = go.Figure(go.Scattergl(                 # Scattergl == WebGL scatter
    x=df['PC1'], y=df['PC2'],
    mode='markers',
    marker=dict(color=df['cluster_code'], colorscale='Tab10', size=4),
    text=df['sample_id'], hoverinfo='text'))
```

## plotly (R) — ggplotly Conversion

```r
library(plotly)
library(ggplot2)

p <- ggplot(df, aes(x = PC1, y = PC2, color = cluster, text = sample_id)) +
    geom_point() + theme_classic()

# Convert ggplot to interactive plotly
p_int <- ggplotly(p, tooltip = c('text', 'x', 'y', 'colour'))

# Save
htmlwidgets::saveWidget(p_int, 'pca.html', selfcontained = TRUE)
```

ggplotly is the lowest-friction R interactive path — write ggplot, get plotly.

## bokeh (Python) — Server-Side / Streaming

```python
from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, HoverTool

output_file('pca_bokeh.html')

source = ColumnDataSource(df)
p = figure(title='PCA', x_axis_label='PC1', y_axis_label='PC2',
           tools='pan,wheel_zoom,box_zoom,reset,hover,save')
p.scatter('PC1', 'PC2', source=source, size=8, alpha=0.7,
          color={'field': 'cluster', 'transform': cluster_cmap})
p.add_tools(HoverTool(tooltips=[('Sample', '@sample_id'), ('Cluster', '@cluster')]))
save(p)
```

bokeh is stronger than plotly for streaming dashboards and server-side aggregation. Static export via `bokeh.io.export_png` requires selenium + Chrome.

## Animation — gganimate (R) and plotly frames (Python)

```r
library(gganimate)
p <- ggplot(df, aes(x, y, color = condition)) +
    geom_point(size = 3) +
    theme_classic() +
    transition_time(time) +                  # animate over time
    labs(title = 'Time: {frame_time}')

anim <- animate(p, nframes = 100, fps = 20, width = 600, height = 400,
                 renderer = gifski_renderer())
anim_save('time_course.gif', anim)
```

```python
import plotly.express as px
fig = px.scatter(df, x='x', y='y', color='condition',
                  animation_frame='time',
                  animation_group='entity_id',
                  range_x=[xmin, xmax], range_y=[ymin, ymax])
fig.write_html('time_course.html')
```

Animation suits time-course data, iterative algorithm visualization, before-after comparisons. Limit to ≤100 frames; longer animations bloat file size and tax viewer attention.

## htmlwidgets Ecosystem (R)

```r
library(DT)                                 # interactive tables
datatable(df, filter = 'top', extensions = 'Buttons',
          options = list(dom = 'Bfrtip', buttons = c('csv', 'excel')))

library(leaflet)                            # interactive maps
leaflet(spatial_df) %>% addTiles() %>% addCircles()

library(networkD3)                          # interactive networks
sankeyNetwork(...) %>% saveWidget('sankey.html')
```

htmlwidgets is the R answer to plotly's JavaScript wrapping — many specialized packages for tables, maps, networks, all producing standalone HTML.

## Per-Method Failure Modes

### plotly static export silently fails

**Trigger:** `fig.write_image('out.pdf')` without kaleido installed.

**Mechanism:** plotly previously fell back to orca (now removed); current versions raise ValueError but older versions silently skipped.

**Symptom:** No file written; OR file written with default settings.

**Fix:** `pip install kaleido`; verify Chrome is installed (kaleido v1+ requires it); test with `fig.write_image('test.pdf')` after install.

### orca dependency in older code

**Trigger:** Following 2020-2022 plotly tutorials with `engine='orca'`.

**Mechanism:** orca is EOL; `engine=` parameter deprecated in plotly 6.2 (post-Sep 2025).

**Symptom:** ValueError or DeprecationWarning.

**Fix:** Remove `engine=` argument; use Kaleido v1 (default).

### EPS export needed but Kaleido v1 dropped it

**Trigger:** Journal requires EPS; Kaleido v1 only supports PDF/PNG/SVG/JPG/WebP.

**Mechanism:** Bundled Chromium in v0 supported EPS; v1 unbundled and dropped it.

**Symptom:** kaleido error on EPS export.

**Fix:** Export PDF, then convert via `pdf2ps` (ghostscript). For complex figures may produce raster EPS — verify acceptability with journal.

### HTML file > 10 MB

**Trigger:** Plotly scatter of 50000 points exported as HTML.

**Mechanism:** Each point + hover data embedded; JS bundle ~3 MB; data scales linearly.

**Symptom:** Browser hangs opening; reviewer's network throttles upload.

**Fix:** Use Scattergl (WebGL); OR Datashader pre-aggregation; OR ship static + small HTML supplement.

### gganimate slow on large frames

**Trigger:** `transition_time` with 100+ frames and 10000+ points per frame.

**Mechanism:** Each frame rendered independently.

**Symptom:** Animation takes hours.

**Fix:** Downsample frames; pre-aggregate per-frame data; OR use plotly animation (in-browser interpolation faster).

### Interactive plot shown as figure in paper

**Trigger:** Manuscript references interactive HTML as Figure 2.

**Mechanism:** Journals require static; interactive HTML is supplement.

**Symptom:** Submission requires figure resubmission as static.

**Fix:** Always produce both static (figure) + interactive (supplement) versions.

## Reconciliation

| Pattern | Cause | Action |
|---------|-------|--------|
| Kaleido / orca confusion in plotly | Pipeline changed 2024-2025 | Use Kaleido v1+; no `engine=` |
| ggplotly drops some custom theme | Conversion loses non-translatable ggplot elements | Manually re-add via `plotly::layout()` |
| bokeh static export fails | selenium not installed | `pip install selenium`; Chrome required |
| htmlwidgets self-contained doesn't work offline | CDN-linked resources by default | `saveWidget(..., selfcontained = TRUE)` |

## Quantitative Thresholds

| Threshold | Value | Source |
|-----------|-------|--------|
| HTML file size warning | >10 MB | Practical |
| Scattergl trigger | >5000 points | plotly performance |
| Animation max frames | ~100 | Viewer attention + file size |
| Selfcontained HTML on | always for portability | htmlwidgets best practice |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Static export silent failure | kaleido / Chrome missing | Install both |
| HTML bloated | Large N points | Scattergl or Datashader |
| orca DeprecationWarning | Following old tutorial | Remove engine=, use Kaleido v1 |
| EPS export fails | Kaleido v1 dropped EPS | PDF + pdf2ps |
| ggplotly tooltips show wrong fields | Default `tooltip` argument | Specify `tooltip = c(...)` |
| Animation file too large | Too many frames | Downsample / pre-aggregate |
| Interactive cited as paper figure | Journal requires static | Produce both |

## References

- Sievert C. 2020. *Interactive Web-Based Data Visualization with R, plotly, and shiny.* Chapman and Hall/CRC.
- Plotly Python — Static Image Generation Changes (2024-2025). https://plotly.com/python/static-image-generation-changes/
- Bostock M, Ogievetsky V, Heer J. 2011. D³ Data-Driven Documents. *IEEE TVCG* 17(12):2301-2309.
- Wickham H, Pedersen TL, Seidel D. 2022. gganimate (CRAN). https://gganimate.com

## Related Skills

- reporting/quarto-reports - Embed interactive HTML in scientific reports
- reporting/rmarkdown-reports - htmlwidgets in Rmd
- data-visualization/ggplot2-fundamentals - ggplot input for ggplotly
- data-visualization/dimensionality-reduction-plots - Interactive UMAP/PCA exploration
- data-visualization/network-visualization - PyVis interactive networks
