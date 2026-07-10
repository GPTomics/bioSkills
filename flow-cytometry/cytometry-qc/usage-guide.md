# Cytometry QC - Usage Guide

## Overview
Quality control detects and removes acquisition artifacts before any downstream analysis. The Time parameter is the master axis: flow-rate instability, signal drift, clogs, and CyTOF sensitivity decay all show up as signal-versus-time anomalies, which is what flowAI, flowCut, flowClean, and PeacoQC detect. This skill covers tool selection, the mandatory QC ordering (compensate -> transform -> margins -> time QC -> gating), CyTOF-specific checks, the instrument calibration/standardization machinery (MESF, CS&T, peak-2) that makes cross-study comparison possible, and batch-level outlier flagging.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowAI', 'PeacoQC', 'flowDensity', 'CATALYST'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run automated QC on my FCS files and tell me what was removed"
- "Remove margin events, then clean unstable time segments with PeacoQC"
- "Check my CyTOF run for signal drift and gate single nucleated cells"
- "Flag outlier samples across my batch"

## Example Prompts
### Event-level cleaning
> "Run flowAI on these files but it's removing too much - it's a slow acquisition, so tune the flow-rate sensitivity and show me the report before committing."
> "Remove margin/boundary events first, then apply PeacoQC across the marker channels - and explain why margins have to go before the density step."

### CyTOF and drift
> "This is a multi-day CyTOF study - check EQ-bead signal versus time for drift and gate intercalator-positive single cells by Event_length and Gaussian parameters."
> "Detect per-channel signal drift over acquisition time and tell me which markers exceed 10% change."

### Calibration and batch QC
> "Convert my PE MFI to antibodies-bound-per-cell using Quantibrite beads."
> "Build a per-sample QC table (events, flow-rate stability, marker medians) and flag outliers by MAD."

## What the Agent Will Do
1. Confirm data is compensated and transformed (QC on raw data misbehaves).
2. Remove margin/boundary events before any density-based step.
3. Run time-based cleaning with the appropriate tool (flowAI/flowCut/PeacoQC) and the correct, non-inverted parameters.
4. Apply instrument-specific checks (dead-cell dye for flow; DNA/Event_length/Gaussian + bead drift for CyTOF).
5. Summarize per-sample metrics and flag outliers; report what was removed.

## Tips
- Order matters: compensate -> transform -> margins -> time-based QC -> gating -> normalization.
- flowAI checks are FR (flow rate) / FS (signal) / FM (dynamic range) - FM is not "flow"; flowAI is known to be aggressive.
- PeacoQC `MAD` and `IT_limit` are LESS strict when higher (counterintuitive); run `RemoveMargins()` first.
- PeacoQC is the only cleaning tool validated across flow, mass, and spectral.
- flowClean needs ~30,000+ events and writes a "GoodVsBad" parameter you must gate on.
- For CyTOF, EQ-bead-median-vs-Time is the single most informative drift readout (see bead-normalization).
- Report dead-cell % and outlier samples; don't silently auto-exclude a whole sample without inspection.
- Cross-study MFI comparison requires calibration (MESF/ERF beads, CS&T, standardized voltages via peak-2) - record voltages per MIFlowCyt.

## Related Skills
- compensation-transformation - Compensate/transform before time-based QC
- doublet-detection - Singlet discrimination after QC
- bead-normalization - EQ-bead drift correction for CyTOF
- clustering-phenotyping - Cluster only QC-passed events
- experimental-design/batch-design - Anchor/reference-sample design for batch QC
