# Spatial Multi-omics Integration - Usage Guide

## Overview

Integrate spatial RNA with a second modality -- protein, chromatin accessibility, or a histone mark -- on spatial CITE-seq, DBiT-seq, spatial-ATAC, spatial epigenome-transcriptome, or Visium CytAssist data. The governing decision is whether the modalities were co-measured on the SAME pixels (vertical integration, where joint same-cell methods like MOFA and WNN apply) or come from different serial sections (diagonal integration, where the right operation is spatial registration of DISTINCT cell populations, not joint modeling).

## Prerequisites

```bash
pip install muon mudata scanpy anndata
pip install mofapy2          # MOFA backend for muon.tl.mofa
pip install paste-bio        # serial-section registration (diagonal)
```

## Quick Start

- "Integrate my spatial RNA and protein from spatial CITE-seq"
- "Do MOFA on my DBiT-seq RNA + protein pixels"
- "My ATAC and RNA are on adjacent sections -- how do I integrate them?"
- "Register two serial sections into a common coordinate frame"
- "Is my Visium CytAssist protein panel safe to call markers absent?"

## Example Prompts

### Vertical (same-pixel co-profiling)

> "I have spatial CITE-seq with whole-transcriptome RNA and an ADT protein panel on the same pixels -- build a joint embedding and find shared factors."

> "Run MOFA on my DBiT-seq RNA + protein and tell me which factors are RNA-driven versus protein-driven."

> "Cluster my Visium CytAssist Gene+Protein spots on a joint RNA+protein representation, accounting for the spot being multiple cells."

### Diagonal (serial sections, different cells)

> "My chromatin accessibility was measured on one section and RNA on the adjacent section -- align them into a common coordinate frame."

> "Register my two serial spatial sections with PASTE and warn me about what a registered map does and does not mean."

### Panel and regime checks

> "Before I report that a cell type is absent in my spatial protein data, check whether its markers were even in the antibody panel."

> "Tell me whether my platform co-captures both modalities on one pixel or runs adjacent sections, and which integration method that implies."

## What the Agent Will Do

1. Classify the integration regime: same-pixel (vertical) versus serial-section (diagonal).
2. For vertical data, assemble a MuData with both modalities indexed on identical pixels and intersect observations.
3. Normalize each modality in its own statistics (log1p for RNA counts, CLR/arcsinh for the targeted protein panel).
4. Run a joint method (MOFA factors or WNN weighting) and inspect per-modality variance or weights.
5. For diagonal data, compute an optimal-transport or diffeomorphic registration and stack sections on a shared frame, flagging that the result is a coordinate map of different cells.
6. Disclose the pixel-is-not-a-cell and bounded-panel caveats in any biological conclusion.

## Tips

- Serial-section modalities are different cells -- aligning them is registration, never same-cell coupling. Reserve MOFA/WNN for same-pixel data.
- When a vendor's "multiomics" product runs adjacent sections rather than one pixel, default to diagonal and register.
- Every spatial multi-omics platform measures pixels or spots (10-55 um) or sub-micron DNBs; "single-cell multi-omics" is a downstream binning claim, not a measurement.
- A multi-cell pixel mixing two cell types can manufacture apparent within-cell cross-modal coupling -- bin or segment explicitly and say so.
- The protein/ATAC/histone-mark side is a targeted panel; absence is uninformative, exactly like a targeted RNA panel.
- Inspect per-modality variance explained: a factor that is ~100% one modality is not multi-omic structure.
- Model the protein panel as intensity (CLR/arcsinh across pixels), not as RNA counts; see spatial-proteomics for the full intensity-not-counts handling.

## Related Skills

- spatial-transcriptomics/spatial-proteomics - protein-intensity (not count) handling, arcsinh, segmentation-dominated panels for the protein side
- spatial-transcriptomics/spatial-data-io - load each modality with the correct reader before assembling a MuData
- single-cell/multimodal-integration - WNN/totalVI/MOFA+ mechanics, ADT denoising, and the anchor-structure fork in non-spatial data
- multi-omics-integration/mofa-integration - MOFA factor interpretation and likelihood choice across modalities
