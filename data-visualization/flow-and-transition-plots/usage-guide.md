# Flow and Transition Plots - Usage Guide

## Overview

Sankey and alluvial plots show how entities flow or transition between categories. Sankey shows aggregate flows from source to sink (single step); alluvial tracks individual entity trajectories across multiple ordered axes (multi-step). CONSORT diagrams are the formal vertical-box trial-flow standard. ggalluvial (R) is the modern default for alluvial; networkD3 / plotly for Sankey; consort for CONSORT.

## Prerequisites

```r
install.packages(c('ggalluvial', 'networkD3', 'consort'))
```

```bash
pip install plotly kaleido pySankey
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot cell-state transitions across 3 timepoints as an alluvial diagram"
- "Build a CONSORT diagram for the trial enrollment flow"
- "Sankey of variants flowing through the QC filter stages"
- "Color ribbons by starting cluster to show 'where this came from' story"
- "Interactive Sankey of metabolic pathway flux"

## Example Prompts

### Cell-state alluvial

> "Alluvial of N=1000 cells across 3 timepoints (t1, t2, t3). Color by t1 cluster identity. Use ggalluvial. Order categories alphabetically within each axis."

### CONSORT trial diagram

> "CONSORT 2010-compliant trial flow: 200 assessed, 50 excluded with reasons, 150 randomized 1:1, then per-arm allocation/loss/analyzed boxes."

### Variant QC pipeline

> "Sankey showing how N variants flow through QC stages (raw -> quality-filtered -> MAF-filtered -> annotated -> final). Counts at each transition."

### Interactive multi-omics flow

> "Interactive Sankey (networkD3 or plotly) showing pathway -> gene -> variant relationships."

### Drug-response trajectory

> "Alluvial of patient response class (CR/PR/SD/PD) across 3 assessments. Highlight responders -> non-responders."

## What the Agent Will Do

1. Determine flow type: single-step source-to-sink (Sankey) vs multi-step entity trajectories (alluvial) vs trial-flow (CONSORT).
2. Reshape data to required format: long (lodes) or wide (alluvia) for ggalluvial; nodes/links for networkD3/plotly Sankey.
3. Set explicit factor levels for category ordering within each axis (minimizes ribbon crossings).
4. Choose color encoding: by origin (where ribbons started) OR destination (where they ended).
5. Apply ggalluvial::geom_alluvium + geom_stratum for static; networkD3::sankeyNetwork for interactive HTML.
6. For CONSORT: use consort package which scaffolds required boxes.
7. Verify conservation: in-flows = out-flows at each internal node.
8. Export PDF via cairo_pdf (ggalluvial) or kaleido (plotly).

## Tips

- **Sankey vs alluvial are different** - Sankey collapses to aggregate flows; alluvial tracks individual entity paths. Pick based on the scientific story.

- **Set explicit factor levels** for categories within each axis (`factor(x, levels = c('A','B','C'))`). Default ordering shuffles position across axes, creating spaghetti ribbons.

- **Color by origin OR destination, not both.** Common pattern: `aes(fill = first_axis)` for "where did these end-state cells come from?"

- **Max 5-7 categories per column for legibility.** Above this, ribbons become unreadable; aggregate small categories into "Other."

- **Max 4-5 axes for alluvial.** Above this, even with optimal ordering, crossings dominate.

- **CONSORT 2010 is required** for randomized clinical trials. Use the `consort` R package to ensure all required boxes are present.

- **plotly Sankey static export requires kaleido.** `pip install kaleido`; orca is EOL.

- **Flow conservation at internal nodes** - in-flows must equal out-flows. Verify the input before plotting; misalignment renders silently with wrong node sizes.

- **Add count annotations to ribbons** when the number of entities per transition is the key information.

- **For pipeline filtering** (variants through QC, reads through trimming/alignment/dedup), CONSORT-style vertical flow is more readable than horizontal Sankey.

## Related Skills

- data-visualization/upset-plots - Set-intersection alternative
- clinical-biostatistics/trial-reporting - CONSORT in trial publication
- single-cell/trajectory-inference - Cell-state transitions to visualize
- workflows/biomarker-pipeline - Pipeline filtering flows
