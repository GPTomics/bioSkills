# FCS Handling - Usage Guide

## Overview
FCS (Flow Cytometry Standard) is the standard file format for cytometry data. flowCore provides comprehensive tools for reading, writing, and manipulating FCS files.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install('flowCore')
```

## Quick Start
Tell your AI agent what you want to do:
- "Read my FCS files into R for analysis"
- "Combine multiple FCS files into a flowSet"
- "Extract the expression matrix from this FCS file"

## Example Prompts
### Reading Data
> "Load all FCS files from this directory into a flowSet"
> "Read this FCS file and show me the channel names"
> "Extract raw expression data without any transformation"

### Writing Data
> "Save this flowFrame as a new FCS file"
> "Export my gated population to a new FCS file"

### Metadata
> "Show me the acquisition keywords from this FCS file"
> "Rename the channels to match my antibody panel"
> "Add sample metadata to my flowSet"

## What the Agent Will Do
1. Load flowCore library
2. Read FCS files with appropriate settings (transformation=FALSE, truncate_max_range=FALSE)
3. Create flowFrame (single file) or flowSet (multiple files) objects
4. Extract expression matrices and parameter metadata as needed
5. Optionally rename channels or add sample annotations

## Tips
- Use `transformation=FALSE` to get raw data without default transformations
- Use `truncate_max_range=FALSE` to preserve values beyond detector range
- flowFrame = single FCS file, flowSet = collection of samples
- -A (Area) parameters are most common, -H (Height) and -W (Width) are for doublet detection
- Check the Time parameter for acquisition issues

## Key Objects

| Object | Description |
|--------|-------------|
| flowFrame | Single FCS file with expression matrix + metadata |
| flowSet | Collection of flowFrames with shared pData |

## Common Parameters

| Type | Examples |
|------|----------|
| Scatter | FSC-A, FSC-H, SSC-A, SSC-H |
| Time | Time |
| Fluorescence | FL1-A, FITC-A, PE-A |
| Mass (CyTOF) | Ir191Di, Yb176Di |

## References
- FCS specification: https://isac-net.org/
- flowCore: doi:10.1186/1471-2105-10-106
