# FCS Handling Usage Guide

## Overview

FCS (Flow Cytometry Standard) is the standard file format for cytometry data. flowCore provides comprehensive tools for reading, writing, and manipulating FCS files.

## FCS File Structure

- **Header**: File format version, text/data offsets
- **Text segment**: Parameter names, acquisition info
- **Data segment**: Event measurements

## Key Objects

### flowFrame
Single FCS file. Contains:
- Expression matrix (events x parameters)
- Parameter metadata
- Acquisition keywords

### flowSet
Collection of flowFrames from multiple samples.
- Enables batch operations
- Stores sample metadata in pData

## Common Parameters

| Type | Examples |
|------|----------|
| Scatter | FSC-A, FSC-H, SSC-A, SSC-H |
| Time | Time |
| Fluorescence | FL1-A, FITC-A, PE-A |
| Mass (CyTOF) | Ir191Di, Yb176Di |

## -A vs -H vs -W

- **-A (Area)**: Total signal area, most common
- **-H (Height)**: Peak signal height
- **-W (Width)**: Signal width, for doublet exclusion

## Best Practices

1. Always check `transformation = FALSE` for raw data
2. Use `truncate_max_range = FALSE` to preserve high values
3. Store original channel names before renaming
4. Check time parameter for acquisition issues

## References

- FCS specification: https://isac-net.org/
- flowCore: doi:10.1186/1471-2105-10-106
