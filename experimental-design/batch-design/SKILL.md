---
name: bio-experimental-design-batch-design
description: Designs experiments to minimize and account for batch effects using balanced layouts and blocking strategies. Use when planning multi-batch experiments, assigning samples to sequencing lanes, or designing studies where technical variation could confound biological signals.
tool_type: r
primary_tool: sva
---

## Version Compatibility

Reference examples tested with: limma 3.58+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Batch Design and Mitigation

**"Design experiment to avoid batch effects"** → Plan sample-to-batch assignments that confound biology with technical variation, and apply correction methods post-hoc.
- R: `sva::ComBat()`, `limma::removeBatchEffect()`
- Python: `scanpy.pp.combat()` for single-cell data

## Core Principle

Batch effects are unavoidable. Good design makes them correctable.

## Design Rules

1. **Never confound batch with condition** - Each batch must contain all conditions
2. **Balance samples across batches** - Equal numbers per condition per batch
3. **Randomize within constraints** - Avoid systematic patterns
4. **Include controls** - Same samples across batches if possible

## Balanced Design Example

```r
# BAD: Confounded design
# Batch 1: All treated samples
# Batch 2: All control samples
# -> Cannot separate batch from treatment

# GOOD: Balanced design
# Batch 1: 3 treated, 3 control
# Batch 2: 3 treated, 3 control
# -> Batch effect can be estimated and removed
```

## Sample Assignment

```r
library(designit)

# Create balanced assignment
samples <- data.frame(
  sample_id = paste0('S', 1:24),
  condition = rep(c('ctrl', 'treat'), each = 12),
  sex = rep(c('M', 'F'), 12)
)

# Optimize batch assignment
batch_design <- osat(samples, batch_size = 8,
                     balance_cols = c('condition', 'sex'))
```

## Detecting Batch Effects

**Goal:** Identify hidden batch effects in expression data by estimating surrogate variables that capture unmodeled technical variation.

**Approach:** Fit a model matrix for the biological variable, estimate the number of surrogate variables using num.sv, then compute surrogate variables with sva for inclusion in downstream differential analysis.

```r
library(sva)

# From count matrix
mod <- model.matrix(~condition, colData)
mod0 <- model.matrix(~1, colData)

# Estimate number of surrogate variables (hidden batches)
n_sv <- num.sv(counts_normalized, mod)

# Estimate surrogate variables
svobj <- sva(counts_normalized, mod, mod0, n.sv = n_sv)
```

## Correction Methods

| Method | When to Use |
|--------|-------------|
| ComBat | Known batches, moderate effects |
| SVA | Unknown batches, exploratory |
| RUVseq | Using control genes |
| limma::removeBatchEffect | Visualization only |

## Documenting Design

Always record:
- Date of sample processing
- Reagent lot numbers
- Operator
- Equipment/lane assignments
- Any deviations from protocol

## Related Skills

- experimental-design/power-analysis - Account for batch in power calculations
- differential-expression/batch-correction - Correcting batch effects in analysis
- single-cell/batch-integration - scRNA-seq batch correction
