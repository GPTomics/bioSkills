# Lollipop / Needle Protein Maps - Usage Guide

## Overview

Lollipop plots visualize per-gene mutation distributions on a protein domain map. Each mutation appears as a vertical stem at its amino-acid position, capped with a circle whose size reflects recurrence count and color encodes variant class. The biological story is hotspot identification - a tall stack at a single residue (KRAS G12, TP53 R175) is the canonical driver-mutation signature. maftools::lollipopPlot is the most-used implementation; trackViewer::lolliplot is more customizable; g3viz produces interactive HTML.

## Prerequisites

```r
install.packages(c('g3viz'))
BiocManager::install(c('maftools', 'trackViewer'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Lollipop plot for TP53 mutations across the cohort, with hotspots R175H, R248Q, R273H labeled"
- "Compare KRAS mutations between cohort A and cohort B as paired lollipops"
- "Use Pfam domain map for TP53 with kinase = blue, DNA-binding = green"
- "Show actual mutation counts on each lollipop, not just size-encoded"
- "Verify that the specified protein isoform matches our variant annotation"

## Example Prompts

### Single-gene cohort lollipop

> "Build a lollipop plot for TP53 from the cohort MAF. Mark canonical hotspots R175H, R248Q, R273H. Color stems by variant class (missense, truncating, splice). Show actual count on each lollipop."

### Cohort comparison

> "Paired lollipop comparing TP53 mutations between BRCA-Luminal and BRCA-Basal cohorts using maftools::lollipopPlot2."

### Custom domain map

> "Use UniProt domain coordinates instead of maftools' cached Pfam. Render with trackViewer::lolliplot."

### Interactive supplement

> "Build an interactive HTML lollipop with g3viz for online supplement; nature theme."

### Isoform-aware annotation

> "Mutations called against ENST00000269305. Annotate the figure with isoform ID; verify residue numbering matches."

## What the Agent Will Do

1. Load MAF / variant table; verify `HGVSp_Short` column for per-mutation AA position.
2. Confirm protein isoform between variant calls and reference; pass explicit `proteinID` if needed.
3. Decide tool: maftools for fast TCGA-style; trackViewer for custom layouts; g3viz for interactive HTML.
4. Set color palette by variant class (CVD-safe; missense one color, truncating another).
5. Mark canonical hotspots if known (TP53 R175/R248/R273; KRAS G12/G13; PIK3CA E545/H1047).
6. Annotate counts (`printCount = TRUE`) so the reader can read actual recurrence.
7. Verify hotspots are recurrent in independent cohort if making novel claims.

## Tips

- **HGVSp_Short column required.** maftools defaults to `HGVSp_Short` (e.g., `p.R175H`). If the MAF has `HGVSp` only, reformat or pass `AACol = 'HGVSp'`.

- **Isoform consistency is critical.** Mutations annotated against one isoform plotted on another are off-by-residue. Document the isoform in the caption; pass `proteinID` explicitly.

- **maftools' Pfam cache may lag.** For current Pfam: pull from UniProt API and pass via `trackViewer::lolliplot` features.

- **`printCount = TRUE`** annotates each lollipop with its mutation count. Size encoding saturates beyond ~10 recurrences; numbers make the difference between 30 and 300 visible.

- **Domain colors should encode functional class** (kinase, binding, regulatory), not random hue. Map manually.

- **lollipopPlot2 for cohort comparison** - one cohort above, one below, sharing the same domain map. Useful for subtype-specific hotspot identification.

- **3D hotspots require 3D-aware tests** (HotMAPS, 3D Hotspots) - linear lollipop misses spatial proximity in folded proteins.

- **Novel hotspot claims require validation** in independent cohorts (TCGA Pan-Cancer + ICGC) plus formal hotspot test (MutSig - Lawrence 2014).

- **g3viz outputs HTML.** Interactive for supplements; static export (`htmlwidgets::saveWidget`) for figures.

## Related Skills

- data-visualization/oncoprint-mutation-matrices - Cohort mutation matrix
- variant-calling/variant-annotation - HGVSp annotation upstream
- clinical-databases/variant-prioritization - Filter to driver candidates
- data-visualization/color-palettes - Variant-class palettes
- structural-biology/structure-navigation - 3D protein structure for hotspot interpretation
