# Spatial Domain Detection - Usage Guide

## Overview

This skill covers identifying spatial domains and tissue regions in spatial transcriptomics data using methods that combine expression and spatial information.

## Prerequisites

```bash
pip install squidpy scanpy scikit-learn
```

## Quick Start

Tell your AI agent what you want to do:
- "Identify spatial domains in my tissue"
- "Cluster my spatial data considering location"

## Example Prompts

### Basic Domain Detection
> "Find spatial domains in my Visium data"

> "Cluster spots using spatial information"

### Combined Methods
> "Cluster using both expression and spatial location"

> "Run spatial-aware clustering"

### Domain Analysis
> "Find marker genes for each domain"

> "Annotate my spatial domains"

## What the Agent Will Do

1. Build spatial and/or expression neighbor graphs
2. Run clustering algorithm
3. Assign spots to domains
4. Optionally smooth domain boundaries
5. Return domain assignments and markers

## Tips

- **Spatial weight** - Balance spatial vs expression (0.2-0.5 typical for spatial weight)
- **Resolution** - Lower resolution = fewer, larger domains
- **Smoothing** - Post-hoc smoothing can clean up domain boundaries
- **Validation** - Check if domains match known tissue architecture
