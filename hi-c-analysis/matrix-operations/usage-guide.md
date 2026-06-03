# Matrix Operations - Usage Guide

## Overview

This skill covers balancing, normalizing, and transforming Hi-C contact matrices with cooler and cooltools. Balancing (ICE/KR/SCALE) removes per-bin coverage bias so a single map is internally self-consistent under the equal-visibility assumption; computing distance-decay expected (`expected_cis` as a per-diagonal P(s) curve, `expected_trans` as a scalar) and dividing observed by it produces an observed/expected (O/E) matrix that exposes loops and compartments above the polymer background. The skill emphasizes two reframes most tutorials skip: balancing is a within-matrix operation and is NOT a cross-sample normalizer (subtracting two balanced maps reads depth as biology), and the equal-visibility model is silently violated by copy-number variation, so raw counts -- not ICE -- are correct for CNV/SV work in aneuploid genomes. It also owns the P(s) log-derivative as a polymer-state diagnostic and the resolution-vs-depth budget.

## Prerequisites

```bash
pip install cooler cooltools bioframe numpy
```

A `.cool`/`.mcool` matrix (see hic-data-io). Analysis functions take a single-resolution URI (`file.mcool::/resolutions/10000`), never the bare `.mcool`. Most operations require a balanced matrix (a stored `weight` column); balancing is a prerequisite for O/E, compartments, insulation, and dots.

## Quick Start

Tell your AI agent what you want to do:
- "Balance my Hi-C cooler with ICE, cis-only"
- "Compute the P(s) contact-decay curve and its log-derivative"
- "Make an observed/expected matrix for chr1"
- "Why does my balanced matrix look wrong on a tumor sample?"
- "What resolution can I afford with this many valid pairs?"

## Example Prompts

### Balancing
> "Balance my 10kb .mcool with ICE, cis-only, masking low-coverage and blacklisted bins, and store the weights. Tell me whether it converged."

> "KR balancing won't converge on my sparse high-resolution map. What should I switch to?"

### Expected and O/E
> "Compute cis expected with smoothing for hg38 chromosome arms, then build a log2(O/E) matrix for chr2 by dividing out the distance-matched expected."

> "Compute the trans expected scalar per chromosome pair so I can make a genome-wide O/E."

### P(s) diagnostics
> "Plot the contact-probability-vs-distance curve and its log-log derivative; I want to see whether cohesin depletion flattened the loop-extrusion bump."

### Copy-number and cross-sample
> "My sample is an aneuploid tumor. Should I balance before calling CNV from the Hi-C coverage, or use raw counts?"

> "I have two conditions at different depth and want to compare them. Is balancing both then subtracting correct?"

## What the Agent Will Do

1. Open the cooler at the correct single-resolution URI and check whether a `weight` column exists.
2. Mask low-coverage and blacklisted bins (mad_max + blacklist), then run cis-only ICE with `cooler.balance_cooler`, reporting convergence; fall back from KR to ICE/SCALE if a juicer map failed to converge.
3. For aneuploid genomes, route CNV/SV work to raw counts and flag that plain ICE erases copy-number (LOIC/CAIC for 3D structure).
4. Compute `expected_cis` (per-diagonal P(s), smoothed) and/or `expected_trans` (scalar), then build O/E by mapping expected onto the dense matrix by diagonal offset.
5. Derive the P(s) log-derivative for polymer-state/loop-extrusion readout, smoothing in logspace first.
6. Advise a resolution given the valid-pair depth (~1000 contacts/bin), and route cross-sample comparison to hic-differential.

## Tips

- Balance cis-only for compartment/TAD/loop work; trans is weak, noisy ambient ligation that distorts the marginals.
- Mask BEFORE balancing -- mad_max filters on log-marginals and blacklisting removes centromere/rDNA/unmappable bins; balancing amplifies a no-signal bin rather than fixing it.
- `ignore_diags=2` removes ligation chemistry (self-ligation, dangling ends, religated fragments), not 3D contact -- that is why P(s)/expected starts at diagonal 2.
- Never subtract two balanced matrices to compare conditions; balanced magnitude scales with depth and is arbitrary-scaled. Downsample, use O/E, and a replicate-aware tool.
- Use raw counts for CNV/SV from Hi-C; ICE forces equal marginals and erases real copy-number in aneuploids.
- cooler weights are multiplicative (raw*w*w); juicer KR/VC weights are divisive (raw/w/w). Check `divisive_weights` before applying an imported vector.
- KR and ICE share the same fixed point -- pick by convergence, not by "which is better". ICE/SCALE are the robust solvers on sparse maps.
- Smooth P(s) in logspace before taking the derivative, or it is pure noise.

## Related Skills

- hic-data-io - Load the cooler files this skill balances
- compartment-analysis - Consumes the O/E this skill produces
- tad-detection - Insulation needs a cis-balanced matrix
- loop-calling - Dots need balanced + expected
- hic-differential - Cross-sample comparison; subtracting/ratioing conditions belongs there
- hic-visualization - Render balanced/O/E/log matrices
- copy-number/cnv-visualization - Raw-count CNV from Hi-C
- genome-intervals/bigwig-tracks - Export expected/eigenvector tracks
