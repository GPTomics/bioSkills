# Isotope Tracing / SIRM Usage Guide

## Overview
Stable-isotope-resolved metabolomics (SIRM) feeds a 13C/15N/2H-labeled tracer and measures how label propagates into downstream metabolites, reporting metabolic ACTIVITY (flux) rather than pool size. This skill guards against the two errors that silently fabricate flux: reading a pool concentration as activity (pool and labeling can move in opposite directions), and interpreting raw isotopologue intensities without the mandatory natural-abundance / tracer-purity correction and steady-state check.

## Prerequisites
```bash
pip install isocor numpy
# high-resolution R alternative:
# R: install.packages('accucor')
```
Conceptual prerequisites: the chosen tracer and its purity, the molecular (and any derivatization) formula per metabolite, the instrument resolution (low-res QqQ vs high-res Orbitrap), and whether the design is single-timepoint (needs verified steady state) or a labeling time course.

## Quick Start
Tell your AI agent what you want to do:
- "Correct these raw isotopologue areas for natural abundance and tracer purity."
- "Compute the MID and fractional enrichment for this metabolite."
- "Help me choose a tracer for measuring glutamine carbon flux."
- "Decide whether my single-timepoint labeling can be read as flux."
- "Explain why my unlabeled control still shows an M+1 peak."

## Example Prompts
### Correction
> "I have M+0..M+5 areas for lactate from a U-13C6-glucose experiment on a triple quad. Apply natural-abundance and 99% tracer-purity correction and give me the corrected MID."
> "My GC-MS data is TBDMS-derivatized; correct the isotopologues accounting for the derivative formula."

### Design and interpretation
> "I want to know whether the pentose phosphate pathway or glycolysis dominates glucose handling in these cells - which tracer should I use and what readout?"
> "My intermediate pool went up but I think flux went down - how do I tell these apart with a tracer?"
> "I sampled labeling at 5, 15, and 30 minutes - is this at isotopic steady state, and can I run classical 13C-MFA?"

## What the Agent Will Do
1. Establish the question: amount (-> targeted-analysis) vs activity/route (-> tracing).
2. Recommend a tracer (13C/15N/2H, uniform vs positional) and a steady-state vs non-stationary design.
3. Apply natural-abundance + tracer-purity correction with IsoCor (any resolution) or AccuCor (high-res), using the correct formula and purity.
4. Compute the corrected MID and fractional enrichment, and visualize as a stacked-bar MID per condition.
5. Check whether labeling has plateaued before any flux interpretation; flag transients for INST-MFA.
6. Hand off flux fitting (13C-MFA / INST-MFA) to a modeling tool (INCA) and pool quantification to targeted-analysis.

## Tips
- Never plot or model raw isotopologue areas; correct first or the MID is wrong by construction.
- Fractional enrichment is concentration-independent - robust to recovery/matrix effects but silent about amount.
- A rising intermediate pool can mean LESS downstream flux; report pool and labeling separately.
- Quench fast and cold; high-turnover metabolites scramble labeling between harvest and extraction.
- Mass spectra resolve isotopologues (count of heavy atoms), not isotopomers (position) - use positional tracers or NMR for position.

## Related Skills
- metabolomics/targeted-analysis - Absolute pool quantification and MRM/SRM mechanics
- metabolomics/xcms-preprocessing - Upstream LC-MS feature detection
- metabolomics/pathway-mapping - Pathway enrichment that interprets pools, not flux
- systems-biology/flux-balance-analysis - Constraint-based predicted flux, distinct from empirical tracing
