# Context-Specific Models - Usage Guide

## Overview

A context-specific model prunes a generic genome-scale model down to the metabolism active in a particular tissue, cell type, or condition, by integrating transcriptomic or proteomic data. The uncomfortable truth to internalize first: the extraction method, the on/off threshold, and (for some methods) the protected objective shape the resulting model more than the input data does. Systematic evaluations found no method reliably beats the others, and mRNA is a weak, sometimes sign-wrong proxy for flux. So a context-specific model is largely an artifact of analyst choices, and those choices - method, threshold strategy, objective - must be reported as first-class methods and stress-tested. COBRApy does not implement these algorithms; the real Python options are troppo and corda, and the reference multi-method implementations are the MATLAB COBRA Toolbox and RAVEN.

## Prerequisites

```bash
pip install cobra numpy pandas
pip install corda           # turnkey native-Python CORDA
# troppo (github.com/BioSystemsUM/troppo) for GIMME/iMAT/tINIT/FASTCORE in Python
# For iMAT/GIMME/INIT with the most validated implementation: MATLAB COBRA Toolbox / RAVEN
```

Inputs: a generic genome-scale model (e.g. Recon3D/Human1 for human, or a bacterial GEM) and expression data (bulk RNA-seq TPM, microarray, or proteomics) with gene IDs that match the model. For MADE, at least two conditions.

## Quick Start

Tell your AI agent:
- "Build a liver-specific model from my GTEx expression, and tell me which method fits non-growing tissue"
- "Use CORDA in Python to extract a context model from my confidence classes"
- "Map my TPM through the GPR rules to reaction activity scores"
- "Rebuild the model at three thresholds and show me which reactions are threshold-dependent"
- "Which extraction method needs an objective and which does not?"

## Example Prompts

### Method Choice
> "I have hepatocyte RNA-seq. Recommend an extraction method given hepatocytes do not proliferate, and explain why forcing a biomass objective would be wrong here."

### Threshold Sensitivity
> "Extract a context model at the 10th, 25th, and 50th expression percentile and report which pathways are stable versus which appear or vanish with the threshold."

### CORDA in Python
> "Translate my expression into CORDA confidence classes through the GPR and build a context-specific model with corda."

### GPR Mapping
> "Aggregate my per-gene expression to per-reaction activity using min-for-AND and max-for-OR, and flag reactions where the aggregation is ambiguous."

## What the Agent Will Do

1. Load the generic model and confirm gene IDs match the expression data namespace.
2. Choose an extraction method from the decision table (objective-required vs not, task-based, differential).
3. Map expression through the GPR (min for AND, max for OR) into reaction activity or confidence classes.
4. Pick and justify a threshold strategy (global/local/StanDep), then extract with troppo/corda (Python) or COBRA Toolbox/RAVEN (MATLAB).
5. Rebuild at 2-3 thresholds and report stable vs threshold-dependent content.
6. Hand the extracted model to FBA/essentiality, noting that pruning has no enzyme-capacity budget.

## Tips

- Expression is not flux: presence is a weak signal, absence a moderate one; never report a pruned/kept reaction as biological fact without a sensitivity check.
- The threshold is the highest-leverage decision. Sweep it and report which reactions flip; a pathway that appears only at one cutoff is a threshold artifact.
- Match method to biology: GIMME needs an objective (wrong for non-proliferating tissue); iMAT is objective-free; tINIT uses metabolic tasks; MADE needs >=2 conditions.
- Do not force biomass on a hepatocyte or neuron - define a maintenance/functional task or use an objective-free method.
- Use min-for-AND/max-for-OR to aggregate GPRs, but remember it discards all but the limiting/dominant gene's quantitative value.
- COBRApy has no `gimme()`/`imat()`; use troppo or corda in Python, or the COBRA Toolbox/RAVEN in MATLAB - the MATLAB implementations are the most validated.
- Proteomics is closer to flux capacity than mRNA but still not flux; for capacity effects (overflow/Warburg) use enzyme-constrained models (GECKO/sMOMENT).
- For single-cell input, scRNA-seq zeros are mostly technical dropout, which breaks the "absence is a strong constraint" assumption; aggregate to pseudobulk or metacells per cell type before extraction rather than thresholding individual cells.
- Report method, threshold strategy, and objective together; two analysts with different choices will get different tissue models from identical data.

## Related Skills

- systems-biology/flux-balance-analysis - Run FBA/FVA on the extracted context model
- systems-biology/gene-essentiality - Context-specific essentiality on the tissue model
- systems-biology/metabolic-reconstruction - The generic model these methods prune
- differential-expression/de-results - Bulk expression input for extraction
- single-cell/cell-annotation - Cell-type expression for cell-type-specific models
