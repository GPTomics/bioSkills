# Gating Analysis - Usage Guide

## Overview
Gating defines cell populations by drawing boundaries in marker space, organized as a hierarchy. This skill covers manual gates (rectangle, polygon, quadrant, boolean) and reproducible automated gating (openCyto templates, flowDensity, flowClust), the canonical gate order that operates as a contamination funnel (time -> debris -> singlets -> live -> lineage), and the expert essentials: FMO (not isotype) controls set positive/negative boundaries because spreading error dominates them, hierarchical gates must be recomputed and live on the transformed scale, and rare-event/MRD work stays supervised because clustering fails for ultra-rare populations.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowWorkspace', 'openCyto', 'flowDensity', 'flowCore', 'CytoML', 'ggcyto'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Gate lymphocytes, singlets, live cells, then CD4/CD8 T cells"
- "Automate my FlowJo gating strategy across all samples with openCyto"
- "Set the CD25 positive boundary using my FMO control"
- "Gate a rare MRD population and tell me how many cells I need to acquire"

## Example Prompts
### Manual hierarchy
> "Build a standard T-cell hierarchy on my GatingSet: debris on FSC/SSC, singlets on FSC-A vs FSC-H, live on the viability channel, then a CD4/CD8 quadrant on CD3+ cells - and remember to recompute before pulling counts."
> "Set up a boolean gate for CD4+ AND NOT CD8+ and extract its frequency per sample."

### Automated and reproducible
> "Turn this manual gating scheme into an openCyto CSV template and apply it to all 40 samples; use mindensity for the bimodal markers and flowClust for the CD4/CD8 split."
> "Use flowDensity to reproduce my predefined sequential gating strategy automatically."

### Boundaries and rare events
> "Use my FMO controls to set the positive boundaries for CD25 and CD69 - explain why isotype controls would be wrong here."
> "I need MRD sensitivity of 1e-5 - how many cells must I acquire, and why shouldn't I use clustering for this?"

## What the Agent Will Do
1. Build a GatingSet and add gates in the canonical funnel order.
2. Set positive/negative boundaries from FMO controls where available.
3. Recompute the hierarchy and confirm gates are on the data's (transformed) scale.
4. For automation, write an openCyto template or apply flowDensity, verifying method names against the installed version.
5. Extract per-population statistics or export gated populations / FlowJo workspaces.

## Tips
- Gate order is a funnel: time -> debris -> singlets -> live -> lineage; reordering bakes in artifacts.
- Use FMO controls (full panel minus one) to set boundaries; isotype controls are for qualitative reagent checks only.
- After `gs_pop_add`, call `recompute(gs)` or child populations stay empty.
- Gate coordinates must match the data's scale - on a transformed GatingSet, use transformed values.
- Prefer diagonal polygon gates for singlets over rectangles.
- openCyto method names drift across versions; check `gt_list_methods()`.
- Rare events (<1e-4) are invisible to clustering - use supervised/template gating and size acquisition for the ~50-event Poisson floor.
- FlowJo interop is via CytoML and `.wsp` only (not legacy `.jo`).

## Related Skills
- compensation-transformation - Preprocess before gating; gate on the transformed scale
- doublet-detection - The singlet step of the gating funnel
- clustering-phenotyping - Unsupervised alternative for high-dimensional discovery
- differential-analysis - Compare gated population frequencies between conditions
- fcs-handling - Load FCS and import FlowJo workspaces via CytoML
