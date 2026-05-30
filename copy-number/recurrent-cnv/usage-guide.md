# Recurrent and Driver CNV Usage Guide

## Overview

A copy number change in one tumor is an observation; one recurring across many tumors beyond chance is evidence of selection. GISTIC2 separates recurrent driver events from passengers by modeling a background rate and scoring each locus by frequency and amplitude (the G-score), with separate background estimation for focal and arm-level events. Copy-number signatures decompose the genome-wide alteration pattern into the mutational processes that generated it. This skill covers GISTIC2, driver-gene localization from recurrence peaks, and the Steele 2022 and Drews 2022 copy-number signature frameworks, including the caveats that make their output easy to misread.

## Prerequisites

```bash
# GISTIC 2.0: download the compiled binary + MATLAB Compiler Runtime from the Broad
# Institute; there is no conda/pip package. Verify the -refgene .mat matches the genome build.
R -e "BiocManager::install('CINSignatureQuantification')"   # Drews 2022 signatures
pip install SigProfilerAssignment                            # Steele 2022 COSMIC CN
```

Inputs: a pooled cohort segment file (6 columns: sample, chrom, start, end, num_markers, seg.mean), diploid-centered; for signatures, absolute allele-specific copy number from ASCAT/Sequenza/FACETS.

## Quick Start

Tell the AI agent what to do:
- "Run GISTIC2 on my cohort segment file to find recurrent focal amplifications"
- "Localize the driver gene inside this wide GISTIC peak"
- "Separate focal drivers from arm-level passenger events in my cohort"
- "Quantify copy-number signatures for my tumor cohort"
- "Explain why I cannot compare GISTIC q-values between my two cohorts"

## Example Prompts

### Recurrence and drivers

> "Pool these CNVkit segments into a cohort seg file and run GISTIC2 with high peak-boundary confidence, then list focal amplification peaks and their candidate driver genes."

> "This GISTIC peak is 3 Mb wide and contains 25 genes. Localize the likely driver using the COSMIC Cancer Gene Census and explain why the peak gene may be a passenger."

### Focal vs broad

> "Distinguish the focal driver events from arm-level passenger events in my cohort and explain how Ziggurat deconstruction and arm-peeling make that separation."

### Signatures

> "Quantify copy-number signatures for my cohort with the Drews CINSignatures framework using absolute copy number from ASCAT."

> "Explain why three of the Steele 2022 copy-number signatures are oversegmentation artifacts and how that affects my interpretation."

## What the Agent Will Do

1. QC and diploid-center the cohort segmentation before pooling
2. Run GISTIC2 with focal and broad analysis and high peak-boundary confidence
3. Localize candidate driver genes from recurrence peaks plus known-driver evidence
4. Report recurrence frequency (portable) alongside q-values (cohort-size dependent)
5. Quantify copy-number signatures from absolute CN with one consistent framework
6. Flag oversegmentation-artifact signatures and caller-sensitivity caveats

## Tips

- GISTIC q-values fall as cohort size rises; compare recurrence frequency across cohorts, not q-values.
- The seg file must be diploid-centered before pooling; a mis-centered profile inverts every call before GISTIC runs.
- Oversegmented input creates spurious narrow peaks; QC the segmentation first.
- A wide GISTIC peak localizes a region, it does not nominate a gene; intersect with known drivers, expression, and dependency data.
- Copy-number signatures require absolute allele-specific copy number, not relative log2.
- Pick one signature framework (Steele/COSMIC or Drews/CINSignatures) and do not mix exposures; both are caller-sensitive.
- GISTIC 2.0 is a frozen MATLAB-compiled binary needing the MCR runtime; no package.

## Related Skills

- copy-number/allele-specific-copy-number - Absolute CN for GISTIC and signatures
- copy-number/copy-ratio-segmentation - Segmentation quality controlling GISTIC peaks
- copy-number/cnv-annotation - Annotating peaks with genes and driver roles
- copy-number/focal-amplification-ecdna - Resolving focal amplicon architecture
- copy-number/cnv-visualization - Cohort heatmaps of recurrent CNV
- pathway-analysis/go-enrichment - Pathway context for recurrently altered genes
