# Flux Balance Analysis - Usage Guide

## Overview

Predict metabolic fluxes and growth rates using constraint-based modeling with COBRApy.

## Prerequisites

```bash
pip install cobra
# For visualization:
pip install escher
```

## Quick Start

Tell your AI agent what you want to do:
- "Run FBA on the E. coli core model"
- "Predict growth rate on glucose minimal media"
- "Compare growth on different carbon sources"

## Example Prompts

### Basic Analysis

> "Run FBA on my metabolic model and show the growth rate"

> "What is the maximum theoretical yield of ethanol from glucose?"

### Media Conditions

> "Simulate growth on acetate instead of glucose"

> "What happens to growth under anaerobic conditions?"

### Flux Distributions

> "Show me the flux through glycolysis"

> "Run FVA to find which reactions are flexible"

### Production Analysis

> "Calculate the production envelope for lactate"

> "What is the growth-coupled production rate of acetate?"

## What the Agent Will Do

1. Load metabolic model (from file or BiGG database)
2. Configure media conditions and constraints
3. Run FBA optimization
4. Calculate growth rate and flux distribution
5. Optionally run FVA for flux ranges
6. Report key results and visualize if requested

## Tips

- **Model sources** - BiGG database has curated models for many organisms
- **Solver choice** - GLPK is free; CPLEX/Gurobi faster for large models
- **Growth units** - Growth rate is typically in h^-1 (per hour)
- **Uptake direction** - Negative lower bound = allowed uptake
- **pFBA** - Use parsimonious FBA for more realistic flux distributions
