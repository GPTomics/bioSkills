# Multiplicity Control via Graphical Procedures - Usage Guide

## Overview

Implements multiplicity control for confirmatory clinical trials using closed testing, graphical procedures (Bretz-Maurer-Hommel), gatekeeping (parallel, serial, mixed), and step-down/step-up procedures (Holm, Hochberg, Hommel). Covers FDA Multiple Endpoints Final Guidance (October 2022) and the Goeman 2021 admissibility result that all valid procedures are closed tests.

## Prerequisites

```r
install.packages(c('gMCP', 'graphicalMCP', 'multcomp', 'gatekeeper'))
```

For Python (basic FDR/FWER methods):

```bash
pip install statsmodels
```

## Quick Start

Tell your AI agent what you want to do:
- "Design a graphical procedure for primary endpoint + 2 key secondary endpoints in gMCP"
- "Apply Holm step-down to my 5 endpoint p-values"
- "Choose between Hochberg and Holm given my endpoint correlation structure"
- "Set up serial gatekeeping (hierarchical) for co-primary then 3 secondary tiers"
- "Allocate 20% of alpha to a discovery subgroup analysis per Dane 2019"

## Example Prompts

### Graphical procedures

> "Pivotal trial with primary endpoint A and two key secondary endpoints B and C, both supporting labeling claims. Design a graph in gMCP: primary at full alpha; if rejected, alpha propagates 50/50 to B and C; if B rejects, weight returns to C, vice versa. Generate the gMCP code and matrix2graph specification."

> "My SAP has hierarchical primary (A) then key secondary (B) then exploratory secondaries (C, D). Implement as Bretz-Maurer 2009 graph with full sequential transitions."

### Gatekeeping

> "Co-primary endpoints A1 and A2 (both must reach significance). If both reject, test family of 3 secondary endpoints with Holm. Implement as mixed gatekeeping per Dmitrienko-Tamhane 2008."

> "Parallel gatekeeping: secondary family tested only if at least one of primary family rejects. Implement in gMCP."

### Choosing step-down vs step-up

> "Five correlated endpoint p-values: 0.018, 0.042, 0.038, 0.015, 0.029. Endpoints are positively correlated. Should I use Hochberg or Holm? Provide reasoning and both adjusted p-values."

> "My two endpoints are LDL-C and HDL-C, which move in opposite directions. Is Hochberg valid? Provide PRDS check and recommendation."

### FDA Multiple Endpoints Guidance

> "Justify my multiplicity strategy per FDA 2022 Multiple Endpoints Final Guidance. Distinguish co-primary (all-win), multiple primary (any-wins), and key secondary."

### Subgroup alpha allocation

> "Allocate alpha across primary endpoint (60%), key secondary (30%), and pre-specified discovery subgroup (10%) per Dane 2019 EFSPI white paper. Design the graph."

## What the Agent Will Do

1. Identify the hypothesis structure (primary, co-primary, multiple primary, key secondary, subgroup)
2. Choose the multiplicity procedure (graphical, gatekeeping, hierarchical, Holm, Hochberg, Hommel)
3. Verify PRDS assumption for Hochberg/Hommel via correlation reasoning
4. Implement the procedure in gMCP (graphical) or specify Holm/Hochberg/Hommel in statsmodels
5. Report adjusted p-values and which hypotheses are rejected at α
6. Generate graph visualisation for SAP appendix

## Tips

- **Closed testing is necessary, not just sufficient** (Goeman 2021). Every valid multiplicity procedure is a closed test in disguise.
- **FWER for confirmatory; FDR for exploratory** is the universal regulatory standard.
- **Hochberg requires PRDS** (Sarkar 1998); fails (Type-I inflated) under negative dependence. Fall back to Holm when PRDS unprovable.
- **Hommel uniformly dominates Hochberg** by 1-3% with no Type-I cost when both are valid; use Hommel by default.
- **Bonferroni loses 30-50% power vs Hommel** under positive dependence at m≈10. Always design a graph.
- **Co-primary vs multiple primary distinction matters**: co-primary doesn't split alpha (inflate n for joint power); multiple primary must split alpha.
- **Sensitivity analyses do NOT require multiplicity alpha** -- they're "what if" not "another claim."
- **Pre-specify the graph in SAP**, not at analysis time. Post-hoc weight tuning inflates Type-I.
- **Subgroup alpha budget ≤20%** per Dane 2019 EFSPI white paper; never allocate primary alpha to data-driven subgroups.
- **Win-ratio composite** (Pocock 2012) preserves a single test for multi-component endpoints; requires pre-specified component hierarchy.
- **`statsmodels.stats.multitest.multipletests` default is `method='hs'` (Holm-Sidak)** -- always specify `method='holm'`, `'hommel'`, or `'bonferroni'` explicitly.

## Related Skills

- clinical-biostatistics/trial-reporting - CONSORT 2025 multiplicity reporting
- clinical-biostatistics/subgroup-analysis - Subgroup alpha allocation
- clinical-biostatistics/power-and-sample-size - Power for co-primary
- clinical-biostatistics/adaptive-designs - Combination tests
- clinical-biostatistics/effect-measures - Effect reporting post-multiplicity
- experimental-design/multiple-testing - General FDR/FWER
