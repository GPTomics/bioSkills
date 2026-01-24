# Matrix Operations - Usage Guide

## Overview

This skill covers balancing, normalizing, and transforming Hi-C contact matrices using cooler and cooltools.

## Prerequisites

```bash
pip install cooler cooltools numpy scipy
```

## Quick Start

Tell your AI agent what you want to do:
- "Balance my Hi-C matrix"
- "Compute observed/expected normalization"

## Example Prompts

### Balancing
> "Apply ICE balancing to my cooler file"

> "Normalize my Hi-C matrix"

### Transformations
> "Compute O/E matrix"

> "Log transform the contact matrix"

### Expected Values
> "Calculate distance-dependent expected values"

## What the Agent Will Do

1. Load the cooler file
2. Apply requested normalization
3. Compute expected values if needed
4. Return transformed matrix

## Tips

- **ICE balancing** - Standard normalization for Hi-C
- **O/E** - Reveals enrichment/depletion relative to expected
- **Log2(O/E)** - Symmetric around 0 for visualization
- **cis_only** - Usually balance cis contacts only
