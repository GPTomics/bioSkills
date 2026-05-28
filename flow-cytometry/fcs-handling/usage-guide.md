# FCS Handling - Usage Guide

## Overview
FCS (Flow Cytometry Standard) is the file format for conventional, spectral, and mass cytometry (CyTOF) data. This skill covers reading and writing FCS files, mapping detector channels to antibodies, understanding the FCS 2.0/3.0/3.1/3.2 keyword internals that determine how the data should be read, and choosing between the R (flowCore/flowWorkspace/CytoML) and Python (FlowKit/readfcs) ecosystems. The load step's defaults (silent log-linearization, range truncation) materially change downstream numbers, so they are treated as deliberate decisions, not boilerplate.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowWorkspace', 'CytoML'))

# Python (optional - for FlowKit workflows or the scanpy/AnnData bridge)
pip install flowkit readfcs
```

## Quick Start
Tell your AI agent what you want to do:
- "Read my FCS files into R without any transformation"
- "Map the detector channels in this FCS file to my antibody panel"
- "Load this FCS file into Python as an AnnData object"
- "Import my FlowJo workspace into R"

## Example Prompts
### Reading and inspecting
> "Load all FCS files from this directory into a flowSet, read them truly raw (no linearization, no range truncation), and show me the channel-to-antibody mapping."
> "This is a CyTOF file with metal channels - read it and tell me which channels are markers vs Time/Event_length/beads."

### Compensation keyword and reader choice
> "Pull the recorded spillover matrix out of this FCS file's keywords - it may be stored as $SPILLOVER, SPILL, or $COMP."
> "I want to do the rest of my analysis in scanpy - load these FCS files into AnnData with readfcs."

### Workspaces and writing
> "Import my FlowJo .wsp workspace and its gated populations into a GatingSet."
> "Rename the channels to the antibody names, attach my sample metadata table, and write out the result."

## What the Agent Will Do
1. Choose a reader based on instrument and downstream ecosystem (flowCore for R/Bioconductor, FlowKit/readfcs for Python).
2. Read with deliberate settings (`transformation=FALSE`, `truncate_max_range=FALSE`) so values are genuinely raw.
3. Build the detector-to-antibody map from `pData(parameters())` (or `pnn`/`pns` labels in FlowKit).
4. Resolve the compensation keyword across the $SPILLOVER / SPILL / $COMP conventions.
5. Optionally rename channels, attach sample metadata, subset, or write FCS.

## Tips
- `read.FCS()` defaults to `transformation="linearize"` - this silently applies $PnE scaling. Use `transformation=FALSE` for raw values that match other tools.
- `truncate_max_range=FALSE` preserves out-of-range events (common on CyTOF and digital instruments).
- Use `alter.names=TRUE` if channel names with hyphens/spaces break formulas or gating.
- FlowJo/Cytobank/Diva workspaces are parsed by CytoML (`flowjo_to_gatingset`), not flowWorkspace; only `.wsp` (FlowJo 10+), not legacy `.jo`.
- `flowFrame` = single file; `flowSet` = collection with shared `pData`. CATALYST/diffcyt require a `pData`/metadata table keyed by sample name.
- -A (Area) is the analysis channel; -H (Height) and -W (Width) drive doublet discrimination.

## Related Skills
- compensation-transformation - Compensate and transform after loading
- cytometry-qc - Assess acquisition quality on the loaded data
- gating-analysis - Define populations from the loaded GatingSet
- clustering-phenotyping - Unsupervised analysis of the event matrix
- single-cell/data-io - readfcs bridges FCS to the AnnData/scanpy ecosystem
- imaging-mass-cytometry/data-preprocessing - Shared metal-channel and FCS conventions
