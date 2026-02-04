# Longitudinal Monitoring - Usage Guide

## Overview
Track ctDNA dynamics over time for treatment response monitoring. Analyze tumor fraction trends, mutation clearance, and detect molecular relapse before clinical progression.

## Prerequisites
```bash
pip install pandas numpy scipy matplotlib
```

## Quick Start
Tell your AI agent what you want to do:
- "Track ctDNA levels across my serial samples"
- "Calculate mutation clearance kinetics"
- "Detect molecular relapse from tumor fraction trends"
- "Generate a treatment monitoring report"

## Example Prompts

### Tumor Fraction Tracking
> "Plot tumor fraction dynamics over treatment for this patient."

> "Calculate the log-fold change in tumor fraction from baseline."

### Mutation Tracking
> "Track these specific mutations across all timepoints."

> "Calculate the half-life of mutation clearance."

### Response Assessment
> "Has this patient achieved a 2-log reduction in ctDNA?"

> "Detect if there's evidence of molecular relapse."

## What the Agent Will Do
1. Compile serial ctDNA measurements
2. Calculate response metrics (nadir, fold-change)
3. Track individual mutation dynamics
4. Assess molecular response criteria
5. Detect rising ctDNA for relapse

## Tips
- Use consistent assay across timepoints for valid comparison
- Plot tumor fraction on log scale for exponential kinetics
- ctDNA often rises months before imaging progression
- Define response criteria based on assay LOD (e.g., 2-log reduction)
- Consider both tumor fraction and mutation clearance

## Related Skills
- ctdna-mutation-detection - Detect mutations to track
- tumor-fraction-estimation - Estimate tumor fraction per timepoint
- fragment-analysis - Complement with fragmentomics trends
