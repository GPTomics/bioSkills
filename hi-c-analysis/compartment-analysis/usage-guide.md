# A/B Compartment Analysis - Usage Guide

## Overview

This skill calls A/B chromatin compartments from balanced Hi-C contact matrices using cooltools eigenvector decomposition. The active (A) and inactive (B) compartments are read from an eigenvector of the distance-normalized, Pearson-correlated cis matrix -- but two decisions make or break the result: which eigenvector is actually the compartment track (E1 often captures a centromere arm gradient instead), and the eigenvector's sign, which is arbitrary until phased against a GC or gene-density track. The skill covers per-arm eigendecomposition, GC phasing with bioframe, resolution choice, saddle plots and saddle_strength for compartmentalization strength, the counterintuitive cohesin-loss-strengthens-compartments result, subcompartment tooling (SNIPER/Calder/dcHiC), and cross-condition switching.

## Prerequisites

```bash
pip install cooler cooltools bioframe matplotlib
```

The cooler must already be balanced (a stored `weight` column) before compartment analysis -- balance it in matrix-operations first. A genome FASTA with a `.fai` index is needed for the GC phasing track. For `.mcool` files, use a single-resolution URI such as `matrix.mcool::/resolutions/100000`.

## Quick Start

Tell your AI agent what you want to do:
- "Call A/B compartments from my balanced Hi-C cooler at 100kb"
- "Phase the compartment eigenvector by GC content from hg38"
- "Run the eigenvector per chromosome arm so the centromere gradient does not dominate"
- "Build a saddle plot and report compartment strength"
- "Compare compartmentalization between treatment and control"

## Example Prompts

### Compartment calling with phasing
> "I have a balanced .mcool. Call A/B compartments at 100kb, run eigs_cis per chromosome arm, phase E1 by GC content from hg38, and save the eigenvector and A/B calls as a BED file."

> "My compartment track is monotonic across chr1 with no sign flips -- I think I captured the arm gradient. Split the view at centromeres, request 3 eigenvectors, and pick the one with the strongest correlation to GC."

### Compartment strength
> "Compute a saddle plot from my balanced cooler and report the compartment strength as the (AA+BB)/(AB+BA) corner ratio at a fixed extent."

> "I knocked out Nipbl and the compartments look stronger -- is that expected? Quantify strength in both samples with matched saddle settings."

### Cross-condition comparison
> "Compare A/B compartments between two conditions and flag bins that switch -- use dcHiC so the signs are coherent across samples rather than hand-diffing two eigenvectors."

### Subcompartments
> "I want A1/A2/B1/B2/B3 subcompartments, not just A/B -- which tool should I use given my coverage?"

## What the Agent Will Do

1. Confirm the cooler is balanced and select a single-resolution URI at 100kb-1Mb.
2. Build a per-arm `view_df` by fetching centromeres and calling `bioframe.make_chromarms`, subset to the cooler's chromosomes.
3. Compute a GC-content phasing track with `bioframe.frac_gc` on the cooler's bins at the matched resolution.
4. Run `cooltools.eigs_cis` with the GC track, `n_eigs=3`, and `sort_metric='pearsonr'`, then confirm E1 correlates with GC (else pick E2/E3).
5. Assign A (positive E1) / B (negative E1) and export the track.
6. Optionally compute `expected_cis`, build the saddle with `cooltools.saddle`, and read `saddle_strength` at a fixed extent.
7. For multiple conditions, route to dcHiC (hic-differential) for sign-coherent differential compartments.

## Tips

- Balance the cooler first -- an unbalanced matrix makes O/E all-NaN and the eigenvector meaningless.
- Run per chromosome arm, not per whole chromosome: the centromere/arm gradient otherwise hijacks E1.
- Always phase the eigenvector with a GC (or gene-density / H3K27ac) track at the cooler's exact binning; wrong phasing flips A<->B silently.
- Set `sort_metric='pearsonr'` so eigenvectors are ordered by GC correlation, not eigenvalue.
- Call compartments at 100kb-1Mb; finer resolutions mix in TAD/loop structure.
- A preserved-or-stronger saddle after cohesin/Nipbl/RAD21 loss is the expected result -- compartment strength and TAD strength are antagonistic.
- Subcompartments cannot be obtained by raising `n_eigs`; use SNIPER, Calder, or dcHiC.
- Keep `n_bins`, `qrange`, resolution, and the corner extent identical across all compared samples for strength.

## Related Skills

- matrix-operations - Balancing and distance-normalized expected that compartment calling depends on
- hic-data-io - Load and access the cooler files this skill operates on
- hic-differential - dcHiC differential compartments and cross-condition switching
- tad-detection - The loop-extrusion partner of the two-mechanism framework
- hic-visualization - Render the eigenvector track and saddle plot
- chip-seq/chromatin-state-segmentation - Overlay ChromHMM/histone states on A/B compartments
- chip-seq/peak-annotation - Annotate switched bins with TF/histone peaks
- genome-intervals/bigwig-tracks - Export the eigenvector as a bigWig track
- single-cell/scatac-analysis - Single-cell chromatin context for scHi-C compartment work
