---
name: bio-flow-cytometry-gating-analysis
description: Defines cell populations in flow and spectral cytometry through manual gates (rectangle, polygon, quadrant, boolean) and reproducible automated gating (openCyto gating templates, flowDensity data-driven thresholds, flowClust model-based gates), organized as a hierarchical GatingSet (flowWorkspace) and round-tripped with FlowJo via CytoML. Covers the canonical gate order (time -> debris -> singlets -> live -> lineage), FMO-vs-isotype boundary setting, gate-order dependence and recompute semantics, rare-event/MRD gating, and per-population statistics. Use when building a gating strategy, automating a manual FlowJo scheme across samples, choosing manual vs data-driven gates, or extracting population frequencies.
tool_type: r
primary_tool: flowWorkspace
---

## Version Compatibility

Reference examples tested with: flowWorkspace 4.14+, openCyto 2.14+, flowDensity 1.36+, flowCore 2.14+, CytoML 2.14+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

openCyto gating-method names drift across versions - confirm with `gt_list_methods()` on the installed package (e.g. `gate_flowclust_2d` vs `flowClust.2d`). Adapt rather than retrying.

# Gating Analysis

**"Gate my data to identify cell populations"** -> Define populations by drawing boundaries in marker space, organized as a hierarchy, manually or with reproducible data-driven methods.
- R (manual + hierarchy): `flowCore` gates -> `flowWorkspace::GatingSet` -> `gs_pop_add` -> `recompute`
- R (automated): `openCyto` gating template (CSV) or `flowDensity::deGate`

## The Single Most Important Modern Insight -- FMO, Not Isotype, Sets the Boundary; and Gate Order Is a Funnel

The position of a positive/negative boundary is governed by SPREADING ERROR - the variance that every other bright fluorophore spills into the channel of interest - NOT by nonspecific antibody binding (Roederer 2001 *Cytometry* 45:194). An FMO control (full panel minus the one channel) reproduces exactly that spreading and is the correct way to set the gate; an isotype control addresses only nonspecific binding, has a different total fluorochrome load, and sits in the wrong place. Isotypes are deprecated for boundary-setting (still fine for a qualitative new-reagent check). Equally load-bearing is gate ORDER: time -> debris (FSC/SSC) -> singlets (FSC-A vs FSC-H) -> live/dead -> lineage. This is a funnel that removes the broadest, least-specific contaminants first (time instability corrupts ALL channels; doublets are scatter-normal AND viable AND double-positive; dead cells bind antibody nonspecifically) so each narrower downstream gate operates on clean input. Reorder it - gate lineage before singlets - and artifacts are baked into the result that no later gate can remove.

## Automated-Gating Taxonomy

| Method | Citation | Mechanism | When to use |
|--------|----------|-----------|-------------|
| openCyto | Finak 2014 *PLoS Comput Biol* 10:e1003806 | CSV gatingTemplate + per-gate algorithms | reproduce a manual SOP across many samples; human-readable + automated |
| `mindensity` (openCyto) | - | KDE valley between two peaks | clear bimodal marker, 1D cut |
| `tailgate` (openCyto) | - | KDE-derivative tail onset | rare positive tail, no clean second peak |
| `quantileGate` (openCyto) | - | cut at a fixed event quantile | threshold should track a fraction |
| flowDensity | Malek 2015 *Bioinformatics* 31:606 | sequential bivariate density cutoffs | reproduce an entire predefined manual strategy |
| flowClust / `gate_flowclust_2d` | Lo 2009 *BMC Bioinformatics* 10:145 | t-mixture + Box-Cox, K by BIC | overlapping elliptical populations |
| DAFi | Lee 2018 *Cytometry A* 93:597 | recursive filter + clustering on a hierarchy | discovery WITH interpretability |

Rule of thumb: 1D bimodal -> `mindensity`; rare tail -> `tailgate`; overlapping ellipses -> `flowClust.2d`; replicate a full manual SOP -> flowDensity; discovery-with-interpretability -> DAFi.

## Build a Gating Hierarchy

**Goal:** Apply gates in the canonical order and extract population statistics.

**Approach:** Build a GatingSet, add gates parent-by-parent, then `recompute()` - WITHOUT it, child populations are empty. Gates apply on the TRANSFORMED scale if the GatingSet is transformed.

```r
library(flowWorkspace); library(flowCore)

gs <- GatingSet(fs)
# matrix dimnames preserve 'FSC-A'/'FSC-H'; data.frame() would mangle them to FSC.A
singlet <- polygonGate('singlets', .gate = matrix(
  c(2e4, 1e4, 25e4, 2e5, 25e4, 26e4, 2e4, 4e4), ncol = 2, byrow = TRUE,
  dimnames = list(NULL, c('FSC-A', 'FSC-H'))))
gs_pop_add(gs, singlet, parent = 'root')
gs_pop_add(gs, rectangleGate('CD3+', CD3 = c(1.5, Inf)), parent = 'singlets')  # transformed scale
recompute(gs)                                   # REQUIRED - else children are empty
gs_pop_get_stats(gs, type = 'count')
```

## Automated Gating with an openCyto Template

**Goal:** Apply a reproducible, declarative gating strategy across all samples.

**Approach:** A CSV template (alias/pop/parent/dims/gating_method/gating_args) defines the hierarchy; `gt_gating` applies it. Confirm method names with `gt_list_methods()`.

```r
library(openCyto); library(data.table)

tmpl <- fread('
alias,pop,parent,dims,gating_method,gating_args
nonDebris,+,root,FSC-A,mindensity,
singlets,+,nonDebris,"FSC-A,FSC-H",singletGate,
live,-,singlets,"Live_Dead",mindensity,
CD3,+,live,CD3,mindensity,
CD4CD8,+,CD3,"CD4,CD8",gate_flowclust_2d,K=2
')
gt <- gatingTemplate(tmpl)
gs <- GatingSet(fs)
gt_gating(gt, gs)
```

## Rare-Event / MRD Gating

**Goal:** Detect a rare population (e.g. MRD at 1e-4 to 1e-5).

**Approach:** Unsupervised clustering FAILS here (a 1e-5 population is ~10 events, invisible to density/SOM); MRD stays supervised/template-gated. Compute the acquisition depth needed from the target sensitivity and the ~50-event Poisson rule BEFORE acquiring; never downsample.

```r
# Need ~50-60 target events for CV < ~15%; sensitivity 1e-5 => acquire ~1e6 cells.
target_sensitivity <- 1e-5
events_needed <- ceiling(50 / target_sensitivity)   # cells to acquire
# Gate the rare population with a prespecified template; report observed LOD from cells acquired.
```

## Per-Method Failure Modes

### Empty child populations
**Trigger:** querying stats right after `gs_pop_add`. **Mechanism:** membership not computed. **Symptom:** zero counts. **Fix:** `recompute(gs)`.

### Gate coordinates on the wrong scale
**Trigger:** raw-scale gate values on a transformed GatingSet (or vice versa). **Mechanism:** scale mismatch. **Symptom:** gate in the wrong place / empty. **Fix:** set gate values on the same (transformed) scale the GS uses.

### Isotype-defined boundary
**Trigger:** isotype control to set positivity. **Mechanism:** spreading error, not nonspecific binding, sets the edge. **Symptom:** wrong negative boundary. **Fix:** use FMO.

### Clustering used for rare events
**Trigger:** FlowSOM for a 1e-5 population. **Mechanism:** too few events. **Symptom:** rare pop absorbed into a neighbor. **Fix:** supervised/template gating; size acquisition for the Poisson floor.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| ~50-60 events for CV < 15% | Poisson statistics | rare-event detection floor |
| sensitivity 1e-5 needs ~1e6 cells | Poisson floor | to collect ~50 events at that frequency |
| FMO for boundary, not isotype | Roederer 2001; Maecker & Trotter 2006 | spreading error dominates the boundary |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| zero counts in children | no `recompute()` | call it after adding gates |
| `gt_gating` method not found | version-renamed method | check `gt_list_methods()` |
| `filter()` vs `Subset()` confusion | `filter` returns a mask, `Subset` the data | use `Subset(ff, gate)` for the population |
| FlowJo `.jo` won't import | only `.wsp` supported | re-save as wsp; use CytoML |

## References

- Roederer 2001 *Cytometry* 45(3):194-205 — spreading error sets the gate boundary.
- Maecker & Trotter 2006 *Cytometry A* 69(9):1037-1042 — FMO doctrine, controls, positivity.
- Finak 2014 *PLoS Comput Biol* 10(8):e1003806 — openCyto automated gating templates.
- Malek 2015 *Bioinformatics* 31(4):606-607 — flowDensity data-driven gating.
- Lo 2009 *BMC Bioinformatics* 10:145 — flowClust model-based gating.
- Lee 2018 *Cytometry A* 93(6):597-610 — DAFi directed filtering + clustering.
- Spidlen 2015 *Cytometry A* 87(7):683-687 — Gating-ML 2.0 portable gate standard.

## Related Skills

- compensation-transformation - Preprocess before gating; gate on the transformed scale
- doublet-detection - The singlet step of the gating funnel
- clustering-phenotyping - Unsupervised alternative for high-dim discovery
- differential-analysis - Compare gated population frequencies between conditions
- fcs-handling - Load FCS and import FlowJo workspaces via CytoML
