# Primer Validation - Usage Guide

## Overview

This skill covers validating PCR primers for secondary structures, dimers, and other quality metrics using primer3-py thermodynamic calculations. Use this to check existing primers before ordering or troubleshoot PCR problems.

## Prerequisites

```bash
pip install primer3-py pandas
```

## Quick Start

Ask your AI agent:

- "Check my primers for dimers"
- "Validate this primer pair"
- "Does my primer form hairpins?"

## Example Prompts

### Individual Primer
> "Check this primer for secondary structures"

> "What's the hairpin Tm for ATGCGATCGATCGATCGATC?"

> "Validate my forward primer"

### Primer Pair
> "Check if my forward and reverse primers form dimers"

> "Validate this primer pair for PCR"

> "Are these two primers compatible?"

### Batch Validation
> "Validate all my primers and show results in a table"

> "Check these 10 primers for issues"

## What the Agent Will Do

1. Calculate Tm for the primer(s)
2. Check for hairpin formation
3. Check for self-dimer (homodimer) formation
4. For pairs, check heterodimer formation
5. Evaluate 3' end stability and specificity
6. Report any warnings or issues

## Tips

- **Hairpin Tm** - Should be >10C below your annealing temperature
- **Dimer Tm** - Homodimer and heterodimer Tm should be <40-45C
- **3' end** - Should have 1-2 G/C bases for specificity without being too stable
- **Tm matching** - Primer pairs should have Tms within 2C of each other
- **dG values** - More negative dG means more stable (problematic) structure
