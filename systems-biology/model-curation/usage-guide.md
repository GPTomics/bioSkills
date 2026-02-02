# Model Curation - Usage Guide

## Overview

Validate, gap-fill, and curate genome-scale metabolic models using memote for quality assessment and COBRApy for corrections.

## Prerequisites

```bash
pip install memote cobra
```

## Quick Start

Tell your AI agent what you want to do:
- "Validate my metabolic model with memote"
- "Gap-fill my model to grow on M9 media"
- "Find dead-end metabolites in my model"

## Example Prompts

### Quality Assessment

> "Run memote on my model and show the quality score"

> "Check which SBML compliance tests my model fails"

### Gap-Filling

> "My model can't grow - find reactions to add"

> "Gap-fill my model using the BiGG universal model"

### Debugging

> "Find all unbalanced reactions in my model"

> "Which metabolites are dead-ends?"

### Standards Compliance

> "Add SBO annotations to my model"

> "Fix GPR rules to standard format"

## What the Agent Will Do

1. Load the metabolic model
2. Run memote quality tests
3. Identify issues (dead-ends, unbalanced reactions, missing annotations)
4. Suggest or implement fixes
5. Re-test to confirm improvements
6. Export curated model

## Tips

- **Memote score** - Aim for >70% for publication-ready models
- **Dead-ends** - Often indicate missing transport or exchange reactions
- **Mass balance** - Check proton accounting (common source of imbalance)
- **GPR format** - Use 'and' for complexes, 'or' for isozymes
- **Iterative curation** - Run tests after each major change
