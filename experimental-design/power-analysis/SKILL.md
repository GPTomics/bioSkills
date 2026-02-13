---
name: bio-experimental-design-power-analysis
description: Calculates statistical power and minimum sample sizes for RNA-seq, ATAC-seq, and other sequencing experiments. Use when planning experiments, determining how many replicates are needed, or assessing whether a study is adequately powered to detect expected effect sizes.
tool_type: r
primary_tool: RNASeqPower
---

## Version Compatibility

Reference examples tested with: RNASeqPower 1.42+, pwr 1.3+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion("<pkg>")` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Power Analysis for Sequencing Experiments

**"How many replicates do I need for RNA-seq?"** → Calculate statistical power or minimum sample size given sequencing depth, biological variability, and expected effect size.
- R: `RNASeqPower::rnapower()`, `pwr::pwr.t.test()`

## Core Concept

Power = probability of detecting a true effect. Underpowered studies waste resources; overpowered studies are inefficient.

## RNA-seq Power Analysis

**Goal:** Determine whether a planned RNA-seq experiment has sufficient statistical power to detect biologically meaningful fold changes, or calculate the minimum sample size needed for a target power.

**Approach:** Provide sequencing depth, biological coefficient of variation, expected fold change, and significance level to rnapower, which uses a negative binomial model to compute power or required sample size.

```r
library(RNASeqPower)

# Typical parameters
# - depth: sequencing depth per sample (reads/gene)
# - cv: biological coefficient of variation (0.1-0.4 typical)
# - effect: fold change to detect (1.5 = 50% change)
# - alpha: significance level (0.05 standard)

# Calculate power for given sample size
rnapower(depth = 20, n = 3, cv = 0.4, effect = 2, alpha = 0.05)

# Calculate required samples for target power
rnapower(depth = 20, cv = 0.4, effect = 2, alpha = 0.05, power = 0.8)
```

## CV Guidelines

| Experiment Type | Typical CV | Notes |
|-----------------|------------|-------|
| Cell lines | 0.1-0.2 | Low variability |
| Inbred mice | 0.2-0.3 | Moderate |
| Human samples | 0.3-0.5 | High variability |
| Primary cells | 0.3-0.4 | Donor-dependent |

## ATAC-seq Power (ssizeRNA)

```r
library(ssizeRNA)

# For differential accessibility
size.zhao(m = 10000, m1 = 500, fc = 2, fdr = 0.05, power = 0.8,
          mu = 10, disp = 0.1)
```

## Quick Reference

| Effect Size | Recommended n (CV=0.4) |
|-------------|------------------------|
| 4-fold | 3 per group |
| 2-fold | 5-6 per group |
| 1.5-fold | 10-12 per group |
| 1.25-fold | 20+ per group |

## Related Skills

- experimental-design/sample-size - Detailed sample size calculations
- experimental-design/batch-design - Accounting for batch effects in design
- differential-expression/deseq2-basics - Running the actual DE analysis
