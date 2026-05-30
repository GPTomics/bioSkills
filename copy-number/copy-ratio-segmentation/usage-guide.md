# Copy-Ratio Segmentation Usage Guide

## Overview

Segmentation turns a noisy per-bin read-depth profile into clean copy-number segments, and it has two error-prone stages. First, normalization must remove GC, mappability, and replication-timing bias so the only remaining variation is copy number; GC correction alone leaves the replication-timing "wave artifact" behind. Second, the segmentation algorithm (CBS, HMM, HaarSeg, fused lasso) must partition the profile, and the choice has a predictable, depth- and event-size-dependent bias. This skill covers bias correction, algorithm selection, parameter tuning, and the diploid-baseline trap that can invert every call.

## Prerequisites

```bash
R -e "BiocManager::install('DNAcopy')"          # CBS reference implementation
pip install numpy pandas statsmodels            # GC loess normalization
R -e "BiocManager::install('QDNAseq')"          # optional: GC/mappability normalization
```

Inputs: a per-bin depth or copy-ratio table (chrom, position, depth or log2), optionally with GC content and mappability per bin; for bias removal beyond GC, a matched normal or a panel of normals.

## Quick Start

Tell the AI agent what to do:
- "GC-correct this depth profile and segment it with CBS"
- "Choose a segmentation algorithm for my 3x shallow WGS sample"
- "Diagnose why my segmentation produced hundreds of tiny segments"
- "Explain why a wave artifact remains after GC correction and how to remove it"
- "Tune CBS alpha and undo.SD to stop oversegmentation"

## Example Prompts

### Normalization and segmentation

> "Normalize this WGS depth profile for GC bias, then segment with CBS using DNAcopy and a noise-scaled undo threshold to avoid oversegmentation."

> "My exome panel still shows a wavy baseline after GC correction. Explain why a panel of normals is needed and what bias GC correction cannot remove."

### Algorithm choice

> "I have 3x shallow WGS. Decide between CBS, HMM, and HaarSeg and justify the choice with the depth-dependent precision/recall trade-off."

> "My focal amplifications are being blurred by HMM segmentation. Recommend a fix."

### Diagnosis

> "My segmentation shattered into hundreds of micro-segments. Diagnose oversegmentation and tune the parameters."

> "Every call in this whole-genome-doubled tumor has the wrong sign. Diagnose the diploid-baseline centering trap."

## What the Agent Will Do

1. Identify which biases are present (GC, mappability, replication timing, capture)
2. Normalize: GC loess for WGS, panel of normals for exome/panel capture bias
3. Select the segmentation algorithm from sequencing depth and target event size
4. Run CBS, HMM, or HaarSeg with parameters tuned to the noise level
5. Guard against oversegmentation (alpha, undo.SD) and gap-bridging artifacts
6. Verify the diploid baseline is correctly anchored before sign-dependent interpretation
7. Reconcile algorithm disagreements and report arm vs focal confidence separately

## Tips

- GC correction does not remove the replication-timing wave artifact; only a matched normal or panel of normals does.
- For hybrid capture, per-probe capture bias is the dominant effect; a panel of normals is essential, not optional.
- CBS is precise on focal events but its recall collapses below ~5x depth; HMM and HaarSeg are depth-robust but blur small segments.
- Oversegmentation propagates into every downstream analysis, including copy-number signatures; tighten alpha and undo.SD or denoise the input.
- Never depth-center an aneuploid or whole-genome-doubled genome; anchor the diploid baseline with an allele-specific ploidy estimate.
- Segment per chromosome arm so CBS does not bridge centromere/telomere gaps.
- Fix random seeds for HMM; Baum-Welch EM only finds local optima.

## Related Skills

- copy-number/cnvkit-analysis - Read-depth caller exposing segmentation choices
- copy-number/gatk-cnv - Tangent normalization and ModelSegments segmentation
- copy-number/allele-specific-copy-number - ASPCF joint logR+BAF segmentation
- copy-number/recurrent-cnv - Copy-number signatures sensitive to segmentation quality
- copy-number/cnv-visualization - Visual diagnosis of oversegmentation and baseline shift
- genome-intervals/coverage-analysis - Per-bin depth computation upstream of segmentation
