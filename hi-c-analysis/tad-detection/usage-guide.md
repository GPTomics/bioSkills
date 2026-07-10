# TAD Detection - Usage Guide

## Overview

This skill calls topologically associating domain (TAD) boundaries from balanced Hi-C contact matrices using the diamond-window insulation score (cooltools `insulation`) and the multi-window TAD-separation score (HiCExplorer `hicFindTADs`). The key idea is that the boundary, not the domain, is the reproducible unit: a population TAD is the ensemble average over stochastic, cell-specific domains (Bintu 2018), and TAD number/size vary 2-5x across caller, resolution, and normalization with no ground truth (Forcato 2017; Zufferey 2018). The skill therefore favors the continuous boundary-strength track over hard domain partitions, runs a list of window sizes (the window is the scale dial), and compares conditions on the differential SCORE rather than by intersecting unstable domain calls.

## Prerequisites

```bash
pip install cooler cooltools bioframe
# For HiCExplorer:
conda install -c bioconda hicexplorer
```

The cooler MUST be balanced (a stored `weight` column) before insulation; analysis takes a single-resolution URI from an `.mcool` (`file.mcool::/resolutions/10000`).

## Quick Start

Tell your AI agent what you want to do:

- "Compute the insulation score on my balanced cooler and find TAD boundaries"
- "Run insulation at several window sizes and report multi-scale boundaries"
- "Rank boundaries by strength and export the strong ones as BED"
- "Compare boundary strength between WT and KO without intersecting domain calls"
- "Call FDR-controlled domains with hicFindTADs"

## Example Prompts

### Multi-scale boundary calling
> "I have a balanced .mcool at 10kb. Compute the insulation score at windows of 3x, 5x, 10x, and 25x the bin size, flag boundaries at the 100kb window, and save the strong boundaries (ranked by boundary_strength) as a BED file."

### Ranking and inspecting boundaries
> "Rank my TAD boundaries by boundary_strength at the 100kb window, drop any that fall in low-coverage bins by checking n_valid_pixels, and tell me how many strong boundaries survive."

### FDR-controlled domains via HiCExplorer
> "Run hicFindTADs on my corrected matrix, sweeping diamond depths from 30kb to 100kb in 10kb steps with FDR correction, and return the domains and boundaries BED files."

### Cross-condition comparison
> "I have insulation tracks for two conditions at matched resolution and depth. Compute the per-bin delta of the log2 insulation score at the 100kb window and tell me which boundaries strengthen or weaken - do NOT intersect domain partitions."

### Boundary annotation
> "Overlap my strong TAD boundaries with a CTCF ChIP-seq peak set and report what fraction of boundaries carry CTCF."

## What the Agent Will Do

1. Load the cooler at a TAD-scale resolution (10-40kb) from a single-resolution `.mcool` URI and confirm it is balanced.
2. Run `cooltools.insulation` with a LIST of window sizes (3-25x the bin), producing `log2_insulation_score_{W}`, `boundary_strength_{W}`, and `is_boundary_{W}` per window.
3. Rank boundaries by the continuous `boundary_strength` (not the per-dataset `is_boundary` flag) and inspect `n_valid_pixels` in sparse loci.
4. Optionally run `hicFindTADs` on the corrected matrix with a depth sweep and FDR correction for an alternative domain/boundary set.
5. For two conditions, compute the bin-matched delta of the continuous insulation score against a permutation/replicate null rather than set-differencing domain BEDs.
6. Export boundaries and the insulation track as BED/bedGraph for visualization and overlap analysis.

## Tips

- Balance the matrix first; insulation on a raw matrix produces coverage-driven garbage valleys.
- The window is the scale dial: small (~3x bin) = sub-TADs, large (~25x bin) = compartment-domains; sweep a list and report multi-scale.
- Use `boundary_strength` (valley prominence) for ranking and cross-sample comparison; `is_boundary` uses a per-dataset Li threshold and is NOT comparable across samples.
- Never set the window below ~3x the bin size; the diamond is then pure noise.
- Compare conditions on the differential SCORE; intersecting domain partitions manufactures spurious gain/loss.
- A TAD boundary is not an A/B compartment switch - insulation and compartmentalization are orthogonal mechanisms.
- Strong boundaries are typically convergent-CTCF + cohesin; use that as a sanity anchor, not a per-boundary filter.

## Related Skills

- matrix-operations - Balancing and O/E that insulation scoring depends on
- hic-data-io - Load and access the cooler files this skill operates on
- compartment-analysis - The orthogonal Mb-scale mechanism; a boundary is not a compartment switch
- loop-calling - Convergent-CTCF loops anchor the strongest boundaries
- hic-differential - Replicate-aware cross-condition contact comparison
- hic-visualization - Render domains/boundaries on the contact matrix
- chip-seq/peak-annotation - Annotate boundaries with CTCF/cohesin peaks
- genome-intervals/interval-arithmetic - Overlap boundary BEDs with features
- genome-intervals/overlap-significance - Test boundary/CTCF co-localization against a matched null
