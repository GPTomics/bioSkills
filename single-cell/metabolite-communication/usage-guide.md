# Metabolite Communication - Usage Guide

## Overview

Metabolite-mediated cell-cell communication infers which cell types exchange metabolites by scoring enzyme-to-sensor expression from scRNA-seq (MEBOCOST), with flux (scFEA), FBA-state (Compass), and neurotransmitter (NeuronChat) alternatives. Metabolite levels are never measured - they are inferred through a chain of enzyme expression -> flux -> level -> sensing - so this is the most speculative communication layer and is hypothesis-generation only, requiring metabolomics, mass-spectrometry imaging, or tracing to validate.

## Prerequisites

```bash
pip install mebocost scanpy anndata
# If the PyPI build is unavailable or stale, install from source:
# pip install git+https://github.com/kaifuchenlab/MEBOCOST
# MEBOCOST also needs its metabolite-enzyme-sensor database and a mebocost.conf
# pointing at those files (download from the MEBOCOST repository).
```

## Quick Start

Tell your AI agent what you want to do:
- "Score metabolite communication between my cell types with MEBOCOST"
- "Which cell types express the machinery to secrete lactate, and which sense it?"
- "Compare metabolic signaling between tumor and normal tissue"
- "Separate the high-confidence receptor calls from the bidirectional-transporter calls"

## Example Prompts

### Basic Analysis

> "Run MEBOCOST on my annotated scRNA-seq data and filter on permutation FDR"

> "Find significant metabolite sender-receiver pairs between cell types"

### Specific Metabolites

> "Which cells express the enzymes for prostaglandin E2 and which express its sensors?"

> "Show the inferred amino-acid signaling and flag the transporter-based calls"

### Comparative Analysis

> "Compare metabolite communication between treatment and control"

> "Find metabolic signaling differences in the tumor microenvironment"

### Validation Framing

> "List orthogonal experiments to confirm the top metabolite communications"

## What the Agent Will Do

1. Verify the data is log-normalized and uses gene symbols, and that ambient RNA was handled in preprocessing
2. Build a MEBOCOST object with cell-type labels and a config file pointing at the metabolite-sensor database
3. Run permutation inference and filter on the FDR, not the raw p-value
4. Summarize by metabolite and sender->receiver pair, separating lower-confidence transporter-based calls
5. State the enzyme->flux->level->sensing chain explicitly and frame every result as a hypothesis
6. Recommend orthogonal validation (metabolomics, MSI, isotope tracing, or enzyme/sensor perturbation)

## Tips

- **Double inference** - Metabolite levels are inferred from enzyme expression, never measured; report "machinery consistent with producing X", never "cell A produces X".
- **Filter on FDR** - Use `permutation_test_fdr`, not the raw permutation p-value; significance is statistical, not a measured concentration.
- **Capitalized columns** - The result table uses Sender, Receiver, Metabolite_Name, Sensor, Annotation, Commu_Score, Norm_Commu_Score, permutation_test_fdr.
- **Transporters are bidirectional** - A transporter "sensor" may export rather than import and may move several metabolites, so sender/receiver direction can be wrong; treat transporter calls as lower confidence.
- **Decontaminate first** - Ambient RNA inflates enzyme and sensor expression; run SoupX/DecontX/CellBender before inference.
- **Gene symbols and normalization** - MEBOCOST needs gene symbols (not Ensembl IDs) and log-normalized data, plus at least ~10 cells per group for stable statistics.
- **No spatial geometry** - Metabolites diffuse and degrade, but dissociated data has no coordinates; do not claim neighbor exchange without spatial metabolomics.
- **Pick the right tool** - MEBOCOST for enzyme-sensor crosstalk, scFEA for per-cell flux, Compass for metabolic-state comparison, NeuronChat for neural systems; verify the database version and required files against the installed docs.
- **Validation is non-optional** - Confirm with targeted metabolomics, MALDI/DESI imaging, isotope tracing, or enzyme/sensor knockout; expression-only metabolite CCC alone is a lead, not a finding.

## Related Skills

single-cell/cell-communication - Ligand-receptor CCC; the single-inference counterpart this skill mirrors at one extra remove
single-cell/cell-annotation - Cell-type labels define metabolite senders and receivers
single-cell/preprocessing - Log-normalization, gene-symbol mapping, and ambient-RNA decontamination happen here, before inference
metabolomics/pathway-mapping - Places inferred metabolites in pathway context and informs which to prioritize
metabolomics/isotope-tracing - Orthogonal flux validation that a producing cell actually makes the metabolite
systems-biology/flux-balance-analysis - Genome-scale FBA underlying Compass-style per-cell metabolic state
