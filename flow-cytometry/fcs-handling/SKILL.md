---
name: bio-flow-cytometry-fcs-handling
description: Reads, inspects, and writes Flow Cytometry Standard (FCS) files from conventional, spectral, and mass cytometry (CyTOF), and parses FlowJo/Cytobank/Diva workspaces. Covers FCS 2.0/3.0/3.1/3.2 internals ($PnE linear-vs-log, $DATATYPE, $SPILLOVER vs SPILL vs $COMP, $TIMESTEP), channel/parameter metadata, the silent linearize/truncate defaults, and R (flowCore, flowWorkspace, CytoML) plus Python (FlowKit, readfcs) readers. Use when loading flow or mass cytometry data, mapping detector channels to antibodies, extracting the event matrix, choosing a reader, or bridging FCS to the scanpy/AnnData ecosystem before preprocessing.
tool_type: mixed
primary_tool: flowCore
---

## Version Compatibility

Reference examples tested with: flowCore 2.14+, flowWorkspace 4.14+, CytoML 2.14+; Python flowkit 1.1+, readfcs 1.1+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# FCS File Handling

**"Load my FCS files and inspect the channels"** -> Parse FCS format into event matrix + parameter metadata, map detector channels to antibodies, and choose a reader appropriate to the instrument and downstream ecosystem.
- R: `flowCore::read.FCS()` / `read.flowSet()` -> `flowFrame`/`flowSet`; `CytoML::flowjo_to_gatingset()` for FlowJo workspaces
- Python: `flowkit.Sample()` (full workflow) or `readfcs.read()` -> AnnData (scanpy/scverse bridge)

## The Single Most Important Modern Insight -- read.FCS Silently Transforms by Default

`flowCore::read.FCS()` defaults to `transformation = "linearize"`, which APPLIES the `$PnE` log-amplification scaling on read. Two pipelines reading "the same raw FCS" (flowCore default vs `fcsparser`/`transformation=FALSE`) therefore return different numbers, and a compensation matrix computed on one will silently mismatch the other. For any preprocessing pipeline that will compensate and transform downstream, read with `transformation = FALSE` (or `NULL`) to get the genuinely raw values, and set `truncate_max_range = FALSE` so out-of-`$PnR` events (common on CyTOF and some digital instruments) are not silently clipped. Decide the read settings deliberately; they are not nuisance defaults.

## FCS Standard Internals (what the keywords mean)

| Keyword | Meaning | Decision-relevant nuance |
|---------|---------|--------------------------|
| `$PnE` | amplification type `"decades,offset"` | `"0,0"` = linear; FCS 3.1 FORBIDS log-stored floats (a float param must be `"0,0"`); log `$PnE` survives only on legacy integer analog-log data |
| `$DATATYPE` | I (uint) / F (float) / D (double) / A (ASCII, deprecated 3.1) | FCS 3.2 allows MIXED types per parameter via `$PnDATATYPE` (integer Time + float fluorescence) |
| `$PnR` | parameter range | for integers defines the bit mask via next power of two (`$PnR=1024` -> 10-bit), NOT a value clamp |
| `$SPILLOVER` | standardized compensation matrix (3.1+) | digital BD instruments wrote non-standard `SPILL` (no `$`); 3.0 `$COMP` stored a matrix WITHOUT naming parameters (ambiguous -> why `$SPILLOVER` exists) |
| `$TIMESTEP` | seconds per Time-channel unit | the master axis for all time-based QC; missing/wrong `$TIMESTEP` silently breaks flow-rate/drift checks |

FCS standards: 3.0 (Seamer 1997 *Cytometry* 28:118), 3.1 (Spidlen 2010 *Cytometry A* 77:97), 3.2 (Spidlen 2021 *Cytometry A* 99:100). Area/Height/Width = pulse integral/peak/duration; FSC-A vs FSC-H is the doublet axis. CyTOF channels are `<Metal><Mass>Di` (e.g. `Yb176Di`) and report dual counts (pulse-counting at low signal, intensity at high).

## Reader Taxonomy

| Reader | Language | What it does | When to use |
|--------|----------|--------------|-------------|
| `flowCore::read.FCS`/`read.flowSet` | R | core FCS -> flowFrame/flowSet | the default for any R/Bioconductor pipeline |
| `flowWorkspace` GatingSet | R | gated hierarchy container | when carrying gates/populations |
| `CytoML` | R | FlowJo (wsp) / Cytobank / Diva import-export | round-tripping a manual analysis (Finak 2018 *Cytometry A* 93:1189) |
| `flowkit` (Session/Sample) | Python | FCS + GatingML 2.0 + FlowJo wsp + compensation/transforms | Python pipelines, FlowJo interop (White 2021 *Front Immunol* 12:768541) |
| `readfcs` | Python | FCS -> AnnData | bridge to scanpy/scverse and the single-cell categories |
| `fcsparser` / `FlowCal` | Python | low-level reader / reader + MEF calibration | quick parse; FlowCal for MESF/MEF work |

## Load and Inspect FCS (R)

**Goal:** Read one file (or a directory) raw, inspect parameters, and map channels to antibodies.

**Approach:** Read with `transformation=FALSE, truncate_max_range=FALSE`; the channel->antibody map lives in `pData(parameters(fcs))` (`name` = detector, `desc` = antibody).

```r
library(flowCore)

fcs <- read.FCS('sample.fcs', transformation = FALSE, truncate_max_range = FALSE)
params <- pData(parameters(fcs))          # name (detector), desc (antibody), range, minRange
channel_map <- setNames(params$desc, params$name)

fs <- read.flowSet(list.files('data', pattern = '\\.fcs$', full.names = TRUE),
                   transformation = FALSE, truncate_max_range = FALSE)
expr <- exprs(fcs)                         # cells x channels
```

## Access the Compensation Matrix from Keywords

**Goal:** Retrieve the acquisition-recorded spillover matrix, handling the three keyword conventions.

**Approach:** Try `$SPILLOVER`, then the legacy `SPILL`, then `$COMP`; `flowCore::spillover()` resolves the standard slots.

```r
kw <- keyword(fcs)
spill <- kw$`$SPILLOVER`
if (is.null(spill)) spill <- kw$SPILL          # digital BD convention
if (is.null(spill)) spill <- kw$`$COMP`        # legacy FCS 3.0 (unnamed columns)
```

## Load FCS in Python (FlowKit / readfcs)

**Goal:** Read FCS in a Python pipeline, either for FlowKit's compensation/gating or as an AnnData for scanpy.

**Approach:** `flowkit.Sample` exposes raw/compensated/transformed events as DataFrames; `readfcs.read` returns AnnData with channels in `var`.

```python
import flowkit as fk
import readfcs

sample = fk.Sample('sample.fcs')
events = sample.as_dataframe(source='raw')   # source in {'raw','comp','xform'}

adata = readfcs.read('sample.fcs')           # AnnData; adata.var has channel + antibody names
```

## Rename Channels, Subset, Write, Annotate Samples

**Goal:** Standardize channel names to antibodies and attach sample-level metadata for downstream tools.

**Approach:** Replace blank `desc` with `name`; attach a `pData` table keyed by `sampleNames(fs)` (CATALYST/diffcyt require this).

```r
new <- ifelse(is.na(params$desc) | params$desc == '', params$name, params$desc)
colnames(fcs) <- new

fcs_markers <- fcs[, c('CD4', 'CD8', 'CD3')]          # subset channels
write.FCS(fcs, 'out.fcs')

pData(fs) <- data.frame(name = sampleNames(fs),
                        condition = c('Control','Control','Treatment','Treatment'),
                        patient = c('P1','P2','P1','P2'),
                        row.names = sampleNames(fs))
```

## Per-Method Failure Modes

### Silent log-linearization on read
**Trigger:** `read.FCS('x.fcs')` with default args. **Mechanism:** `transformation="linearize"` applies `$PnE` scaling. **Symptom:** values differ from `fcsparser`; compensation matrix mismatch. **Fix:** `transformation = FALSE`.

### Out-of-range clipping
**Trigger:** instrument wrote values above `$PnR` (common CyTOF). **Mechanism:** `truncate_max_range=TRUE` (default) clamps them. **Symptom:** a ceiling artifact at the channel max. **Fix:** `truncate_max_range = FALSE`.

### Channel names break formulas
**Trigger:** channels like `FSC-A`, `Pacific Blue-A`. **Mechanism:** hyphens/spaces are not syntactic R names. **Symptom:** formula/gating errors. **Fix:** `alter.names = TRUE` on read.

### FlowJo parsing in the wrong package
**Trigger:** looking for FlowJo import in flowWorkspace. **Mechanism:** parsing lives in CytoML. **Symptom:** function-not-found. **Fix:** `CytoML::open_flowjo_xml()` -> `flowjo_to_gatingset()`; only `.wsp` (FlowJo 10+), not legacy `.jo`.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `exprs()` numbers differ across tools | default linearize | read with `transformation=FALSE` everywhere |
| spillover keyword is `NULL` | instrument used `SPILL`/`$COMP` | try all three keyword names |
| editing `exprs(ff)` corrupts ranges | direct reassignment skips `parameters()` update | use transform/`Subset` workflows |
| readfcs compensation not applied | matrix names don't match `var_names` | align channel names before relying on it |

## References

- Seamer 1997 *Cytometry* 28(2):118-122 — FCS 3.0 standard.
- Spidlen 2010 *Cytometry A* 77(1):97-100 — FCS 3.1 standard.
- Spidlen 2021 *Cytometry A* 99(1):100-102 — FCS 3.2 standard.
- Finak 2018 *Cytometry A* 93(12):1189-1196 — CytoML cross-platform gating import/export.
- White 2021 *Front Immunol* 12:768541 — FlowKit Python toolkit.
- Lee 2008 *Cytometry A* 73(10):926-930 — MIFlowCyt minimum reporting standard.

## Related Skills

- compensation-transformation - Compensate and transform after loading
- cytometry-qc - Assess acquisition quality on the loaded data
- gating-analysis - Define populations from the loaded GatingSet
- clustering-phenotyping - Unsupervised analysis of the event matrix
- single-cell/data-io - readfcs bridges FCS to the AnnData/scanpy ecosystem
- imaging-mass-cytometry/data-preprocessing - Shared metal-channel and FCS conventions
