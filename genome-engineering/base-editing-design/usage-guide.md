# Base Editing Design - Usage Guide

## Overview

Design guide RNAs for cytosine (CBE) and adenine (ABE) base editors to make precise C-to-T or A-to-G conversions without double-strand breaks.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a base editing guide for the C282Y mutation"
- "Design an ABE guide to correct this A>G variant"
- "Find the best CBE guide to convert C to T at this position"

## Example Prompts

### Cytosine Base Editing

> "Design a CBE guide to change C to T at position 50 in this sequence"

> "Find guides that place my target C in the optimal editing window"

### Adenine Base Editing

> "Design an ABE guide to correct this pathogenic A>G variant"

> "Convert the adenine at position 100 to guanine using base editing"

### Bystander Avoidance

> "Find CBE guides with no bystander Cs in the editing window"

> "Which guide has the fewest off-target edits?"

## What the Agent Will Do

1. Identify the target base and desired conversion
2. Search for PAMs that place target in editing window
3. Evaluate position-dependent editing efficiency
4. Identify bystander bases that may be edited
5. Score sequence context effects
6. Return ranked guide options with predictions

## Tips

- **Editing window** - CBE: positions 4-8; ABE: positions 4-7 (from PAM-distal end)
- **Peak efficiency** - Position 5-6 typically has highest editing rates
- **Bystanders** - Any C (CBE) or A (ABE) in the window may be edited
- **Context** - TC context is preferred for CBE (higher efficiency)
- **No bystanders** - Prioritize guides with target as only editable base
