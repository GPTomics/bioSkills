# Loop Calling - Usage Guide

## Overview

This skill detects focal chromatin loops (point interactions / corner-dots) in balanced Hi-C and Micro-C contact maps and validates a loop set. It covers de-novo calling with cooltools `dots` (a pure-CPU HiCCUPS reimplementation using four local-background kernels and lambda-chunked FDR), chromosight (template-correlation), and Mustache (scale-space blob detection), plus aggregate peak analysis (APA) via cooltools `pileup`. The decision that gates everything is depth: de-novo loop calling realistically needs 5-10kb resolution (hundreds of millions to billions of valid pairs); below that, the right move is APA on a known anchor set rather than de-novo calling. The skill stresses that de-novo calling DISCOVERS while APA CONFIRMS, that callers disagree (so consensus plus convergent-CTCF support is how loops are trusted), and that HiChIP/PLAC-seq data needs FitHiChIP/MAPS, not dots.

## Prerequisites

```bash
pip install cooler cooltools bioframe chromosight matplotlib
# Mustache (scale-space loop caller):
pip install mustache-hic
# Optional GPU HiCCUPS path on .hic:
# Juicer Tools (requires CUDA)
```

The cooler must be balanced before loop calling (`cooler balance matrix.mcool::/resolutions/10000`). For an `.mcool`, always pass a single-resolution URI (`file.mcool::/resolutions/10000`), not the bare file.

## Quick Start

Tell your AI agent what you want to do:
- "Call chromatin loops from my balanced cooler at 10kb"
- "Is my map deep enough to call loops de-novo, or should I run APA on known anchors?"
- "Run aggregate peak analysis on my loop set and report the APA score"
- "Validate my loops against convergent CTCF motifs"
- "Compare loop strength between my two conditions"

## Example Prompts

### De-novo loop calling
> "I have a balanced .mcool from a deep in-situ Hi-C library. Call loops at 10kb with cooltools dots: build chromosome-arm regions, compute expected_cis on those arms, then run dots with the same view, max_loci_separation 10Mb, and the default lambda-chunked FDR. Report how many loops and their size distribution."

> "Run chromosight detect for loops on my 5kb cooler between 20kb and 2Mb separation with a Pearson cutoff of 0.4, then also detect borders and stripes so I can keep them separate from the loop set."

### Depth decision and APA confirmation
> "My library only has ~40 million valid pairs. Don't de-novo call. Instead build anchor pairs from my CTCF and RAD21 ChIP-seq peaks and run cooltools pileup as APA, then give me the APA score as center pixel over the lower-left corner block."

> "Pile up my called loop set with cooltools pileup using expected so the snippets are O/E, average across the stack, and tell me whether the aggregate center dot is clean enough that the set is mostly real."

### Validation and differential
> "Intersect the loops called by cooltools dots and Mustache, then check how many consensus loops have convergent CTCF motifs at both anchors."

> "Build a union anchor set across my WT and cohesin-degron conditions, quantify each loop's strength per condition with chromosight quantify, and report which loops weaken after cohesin loss."

## What the Agent Will Do

1. Ask how deep the map is and pick the path: deep -> de-novo calling; shallow -> APA on known anchors.
2. Confirm the cooler is balanced and select a single-resolution URI (5-10kb for de-novo).
3. Build a chromosome-arm `view_df` and compute `expected_cis` on it.
4. Run `cooltools.dots` (or chromosight / Mustache) with the same view, keeping FDR control via lambda-chunking.
5. Validate the call set: consensus across >=2 callers, convergent-CTCF support, and APA of the full set with a corner control.
6. For differential, quantify a union anchor set per condition and test the strength delta.

## Tips

- Decide on depth first - de-novo callers return noise below ~5-10kb resolution; shallow maps get APA on known anchors, not de-novo calls.
- Balance the cooler before anything - `dots` and `pileup` read the `weight` column and raw counts are unsupported.
- Reuse the SAME `view_df` for `expected_cis` and `dots`/`pileup`, or the shapes will not align.
- An enriched APA confirms a SET on average; it never validates an individual loop and is meaningless without a corner control.
- Callers disagree (Forcato 2017) - trust consensus across >=2 tools plus convergent-CTCF or ChIA-PET/HiChIP support, not a single list.
- Stripes are NOT loops - the horizontal/vertical kernels exist to suppress them; detect stripes separately with chromosight stripes_left/stripes_right.
- For HiChIP / PLAC-seq / PCHi-C use FitHiChIP / MAPS / HiC-DC+, not dots - the protein-anchored coverage bias breaks the Hi-C null.
- For sub-5kb Micro-C, shrink HiCCUPS kernels or prefer Mustache/chromosight, which adapt to resolution more gracefully.

## Related Skills

- hic-data-io - Load and access the cooler files this skill calls loops on
- matrix-operations - Balancing and expected/O/E that dots and pileup depend on
- hic-visualization - Render called loops and APA pileups on the heatmap
- hic-differential - Bin/compartment-level differential (the regime loops are NOT in)
- tad-detection - TAD corners vs point loops; the lower-left background separates them
- chip-seq/peak-calling - CTCF/cohesin peaks to anchor and validate loops; HiChIP peak context
- chip-seq/peak-annotation - Annotate loop anchors with TF/CTCF peaks
- atac-seq/enhancer-gene-linking - E-P contacts complementing loop calls
- genome-intervals/overlap-significance - Permutation test for anchor/feature enrichment
