---
name: bio-flow-cytometry-doublet-detection
description: Detects and removes doublets/aggregates from flow, spectral, and mass cytometry before clustering or quantification. Covers FSC-A vs FSC-H singlet discrimination (the Area-Height non-proportionality, not a 1D area gate), FSC-W/SSC width gating, CyTOF Gaussian discrimination parameters (Center/Offset/Width/Residual/Event_length) and DNA intercalator gating, and the residual heterotypic conjugates that survive scatter gating and masquerade as double-positive populations. Use when filtering aggregates before phenotyping, choosing a doublet method for flow vs CyTOF, or diagnosing a suspicious double-positive cluster.
tool_type: r
primary_tool: flowCore
---

## Version Compatibility

Reference examples tested with: flowCore 2.14+, CATALYST 1.26+, ggplot2 3.5+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package and adapt rather than retrying.

# Doublet Detection

**"Remove doublets from my cytometry data"** -> Discriminate single cells from aggregates using pulse geometry (flow) or ion-cloud parameters (CyTOF), before any clustering or quantification.
- R (flow/spectral): `flowCore` gate on the FSC-A vs FSC-H diagonal (+ FSC-W/SSC-W)
- R (mass/CyTOF): gate on DNA intercalator + Gaussian/Event_length parameters

## The Single Most Important Modern Insight -- Doublets Are Caught by Area-vs-Height Non-Proportionality, and Scatter Gating Is Necessary but Not Sufficient

A doublet has roughly double the pulse AREA of a singlet but NOT double the Height, and a longer Width/transit time - so singlets fall on a tight FSC-A vs FSC-H diagonal and doublets deflect above it. A 1D area histogram therefore does NOT remove doublets; the discriminating signal is the Area-Height relationship (plus Width). This matters because an unremoved doublet of a CD3+ and a CD19+ cell reads as an artifactual CD3+CD19+ "double-positive," and clustering will faithfully (and wrongly) carve it out as a real population. Crucially, scatter gating is necessary but NOT sufficient: heterotypic conjugates (e.g. a CD3+CD14+ T:monocyte) survive standard FSC-A/H gates and present as genuine double-positives whose lineage-marker levels look COMPARABLE to true single-positives - the tell is an ELEVATED shared marker (e.g. CD45) and a high bright-field aspect ratio, so the definitive resolver is imaging flow cytometry, not a lineage-intensity check (Stadinski 2020 *Cytometry A* 97:1102). On CyTOF there is no scatter at all - doublets are removed by ion-cloud Gaussian parameters and DNA intercalator content (Bagwell 2020 *Cytometry A* 97:184).

## Method Taxonomy

| Method | Instrument | Principle | Caveat |
|--------|-----------|-----------|--------|
| FSC-A vs FSC-H | flow/spectral | singlets on the A-H diagonal | the standard; the discriminator is non-proportionality, not area |
| FSC-W / SSC-W | flow/spectral | doublets have longer pulse Width | complementary to A-vs-H |
| DNA intercalator (Ir191/193) | CyTOF | doublets show ~2N+ DNA | also separates cells from beads/debris |
| Gaussian params + Event_length | CyTOF | ion-cloud fit residual/length flags fusions | catches fusions DNA alone misses (Bagwell 2020) |
| imaging cytometry | imaging flow | bright-field aspect ratio | the only clean resolver of heterotypic conjugates |

Note: cytometry doublet removal is GATING-based. DoubletFinder/Scrublet/scDblFinder are scRNA-seq DROPLET methods (they simulate artificial doublets) - limited transfer, because cytometry has direct physical doublet signals.

## FSC-A vs FSC-H Singlet Gating (flow/spectral)

**Goal:** Keep events on the singlet diagonal.

**Approach:** A polygon along the A=H diagonal (preferred over a rectangle, which keeps off-diagonal doublets); visualize with the gate overlaid.

```r
library(flowCore); library(ggcyto)

# matrix dimnames preserve 'FSC-A'/'FSC-H'; data.frame() would mangle them to FSC.A
singlet <- polygonGate(filterId = 'singlets', .gate = matrix(
  c(20000, 10000, 250000, 200000, 250000, 260000, 20000, 40000), ncol = 2, byrow = TRUE,
  dimnames = list(NULL, c('FSC-A', 'FSC-H'))))
singlets <- Subset(fs, singlet)
autoplot(fs[[1]], 'FSC-A', 'FSC-H') + ggcyto::geom_gate(singlet)
```

## CyTOF Doublet Removal

**Goal:** Keep intercalator-positive single ion clouds.

**Approach:** Gate DNA intercalator (nucleated, ~2N) and Event_length/Gaussian residual; CATALYST exposes these as channels in the SCE.

```r
library(CATALYST)
# prepData moves Time/Event_length to int_colData by default - keep them in the assay with FACS=TRUE
sce <- prepData(fs, panel, md, transform = TRUE, cofactor = 5, FACS = TRUE)
e <- assay(sce, 'exprs')

dna <- e['DNA1', ]                                   # intercalator-positive = nucleated single cells
keep <- dna > quantile(dna, 0.05) & dna < quantile(dna, 0.95)
if ('Event_length' %in% rownames(sce))               # retained by FACS=TRUE (now on the arcsinh scale)
  keep <- keep & e['Event_length', ] <= quantile(e['Event_length', ], 0.99)   # quantile-relative, so scale is fine
sce_singlets <- sce[, keep]
```

## Per-Method Failure Modes

### 1D area gate leaves doublets
**Trigger:** gating only FSC-A. **Mechanism:** doublets overlap singlets in area. **Symptom:** double-positive clusters persist. **Fix:** gate the FSC-A vs FSC-H diagonal (+ Width).

### Heterotypic conjugate survives scatter gating
**Trigger:** a surprising double-positive between two single-positive clusters. **Mechanism:** T:monocyte conjugate is scatter-normal, lineage markers comparable to singlets. **Symptom:** "novel" DP population with an elevated shared marker (e.g. CD45). **Fix:** treat as suspected doublet; check the shared-marker signal; confirm/resolve by imaging flow (bright-field aspect ratio) when load-bearing.

### CyTOF "doublet gate" using scatter
**Trigger:** porting flow logic to CyTOF. **Mechanism:** no FSC/SSC exists. **Symptom:** no scatter channels. **Fix:** use DNA + Gaussian/Event_length.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| expected doublet rate ~1-5% (PBMC), higher in tissue | community | flag samples far above as prep issues - not a removal cutoff |
| Gaussian + DNA gating improves CV (3.45 -> ~2.04) | Bagwell 2020 *Cytometry A* 97:184 | combined DNA + Gaussian over baseline (Gaussian alone ~2.41) |

Note: a fixed "95th-percentile residual" cutoff is arbitrary; prefer a visual diagonal gate or the instrument's Gaussian parameters over an unjustified quantile.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| double-positive cluster that "shouldn't" exist | residual heterotypic doublets | check for an elevated shared marker (CD45); confirm by imaging flow |
| no FSC/SSC channels (CyTOF) | mass data has no scatter | use DNA/Gaussian/Event_length |
| over-removal of large cells | rectangle gate clips real large singlets | use a diagonal polygon, not a box |

## References

- Stadinski 2020 *Cytometry A* 97(11):1102-1104 — heterotypic doublets survive scatter gating.
- Bagwell 2020 *Cytometry A* 97(2):184-198 — automated CyTOF cleanup via Gaussian/Event_length.
- Finck 2013 *Cytometry A* 83(5):483-494 — CyTOF DNA/event parameters in normalization context.

## Related Skills

Workflow order (CyTOF): EQ-bead drift normalization (raw, FIRST) -> cytometry-qc -> doublet-detection -> clustering -> CytoNorm cross-batch (LAST)

- cytometry-qc - Run first: flow-rate/signal/margin cleaning
- bead-normalization - CyTOF drift correction after doublet removal
- fcs-handling - Load FCS files
- gating-analysis - Where singlet discrimination sits in the hierarchy
- clustering-phenotyping - Downstream analysis after doublet removal
- single-cell/doublet-detection - Droplet scRNA-seq doublet methods (different principle)
