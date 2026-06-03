# Hi-C Visualization - Usage Guide

## Overview

This skill renders Hi-C contact data into honest, reproducible figures. The central decision is the colorscale: the same matrix under raw, ICE-balanced, or observed-over-expected (O/E) tells three different biological stories, and a reviewer-grade figure states its normalization, its color scale (LogNorm vs symmetric-diverging) with limits or clip percentile, and its bin resolution. It covers square and rotated-triangle matrix maps, virtual 4C profiles, APA/saddle/on-diagonal pileups, two-condition comparisons, NaN/white-stripe handling, and when to use interactive exploration (HiGlass) versus scripted-static publication tools (matplotlib/cooltools, pyGenomeTracks, FAN-C, CoolBox, plotgardener).

## Prerequisites

```bash
pip install cooler cooltools bioframe matplotlib
conda install -c bioconda hicexplorer pygenometracks
# Optional: pip install fanc coolbox coolpuppy   # alternative plotting tools
```

A cooler must be balanced (`cooler balance matrix.cool`) before `matrix(balance=True)` returns anything but NaN. For interactive HiGlass, build a multiresolution tileset first (`cooler zoomify --balance base.cool`).

## Quick Start

Tell your AI agent what you want to do:
- "Plot my balanced Hi-C matrix for chr1:50-60Mb"
- "Make an observed/expected map so compartments and loops are visible"
- "Pile up my loop calls and report the APA score"
- "Compare treatment vs control as a log2-ratio map"
- "Build a virtual 4C profile from a viewpoint"

## Example Prompts

### Single matrix, right transform

> "I have a balanced .mcool. Plot chr2:30-40Mb at 10kb resolution with a fall colormap and LogNorm, set vmax to the 99th off-diagonal percentile, and render masked bins in light gray instead of white."

> "Compartments don't show up on my balanced map. Compute the cis expected, plot log2(observed/expected) with coolwarm and symmetric limits, and tell me what resolution to use for the checkerboard."

### Track-stacked figures

> "Make a rotated-triangle Hi-C panel with genes, a CTCF ChIP bigWig, and the insulation score stacked beneath it on a shared x-axis for chr1:50-60Mb. Pick a depth that covers a 2 Mb TAD."

### Pileups

> "Pile up my 8,000 loop calls as an observed/expected APA, average over the loops, plot log2 with symmetric limits, and compute the center-to-corner APA score."

> "Make a saddle plot to show compartment strength; remember to phase the eigenvector first."

### Comparisons and 4C

> "Compare KO vs WT Hi-C as a single log2-ratio map. Balance and depth-match both libraries first, use RdBu_r with white at no-change, and grey out the very-distal noisy bins."

> "Extract a virtual 4C profile from chr1:55Mb on a log-y axis so the near-cis spike doesn't swamp the distal signal."

## What the Agent Will Do

1. Open the cooler at a single-resolution URI and confirm it is balanced.
2. Choose the transform (raw / balanced / O/E / log2-ratio) for the biological claim.
3. Pick the resolution to match the feature (compartments 100-500kb, TADs 10-40kb, loops 5-10kb).
4. Set the norm and colormap: LogNorm + sequential for counts/balanced; symmetric diverging centered at 0 for O/E and ratios; vmax from a stated off-diagonal percentile.
5. Render NaN bins explicitly with `set_bad`; never quantify on interpolated/imputed matrices.
6. For track stacks, delegate the 45deg shear and depth crop to pyGenomeTracks/HiCExplorer/FAN-C/plotgardener.
7. For pileups, build the stack with `cooltools.pileup`, average over `axis=0`, and report the APA/saddle score.
8. For comparisons, balance and depth-match both libraries before ratioing.
9. Save a static figure whose legend records normalization, scale, and resolution.

## Tips

- Balanced is not the science figure for compartments/loops -- switch to log2(O/E) before any negative claim.
- Symmetric limits (`vmin=-vmax`) are mandatory on every divergent/ratio map, or zero drifts off white.
- Report the vmax percentile; the apparent strength of TADs/loops moves visibly as it slides.
- White stripes at a known blacklist/centromere are expected; an unexpected stripe is a coverage problem worth investigating before trusting nearby features.
- The triangle y-axis is genomic separation, not a second coordinate -- and `depth` silently truncates features larger than it.
- HiGlass is for finding the region/resolution; reproduce that exact view in a scripted tool for the manuscript.
- Use the `file.mcool::/resolutions/<bp>` URI -- a bare `.mcool` gives empty or wrong-resolution output.

## Related Skills

- hic-data-io - Load the cooler/.mcool files and zoomify for HiGlass tilesets
- matrix-operations - Balancing and O/E that the divergent map depends on
- compartment-analysis - Eigenvector phasing behind the saddle plot
- tad-detection - Insulation/boundary tracks stacked under the triangle
- loop-calling - Loop calls and peak-anchored pileup conventions visualized here
- hic-differential - Replicate-aware testing behind the two-condition comparison
- data-visualization/genome-tracks - Config-driven multi-track stacks (pyGenomeTracks hic_matrix)
- genome-intervals/bigwig-tracks - Export eigenvector/insulation as bigWig for the track stack
