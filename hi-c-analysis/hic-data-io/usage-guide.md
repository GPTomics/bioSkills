# Hi-C Data I/O - Usage Guide

## Overview

This skill covers loading, converting, and manipulating Hi-C contact matrices in the cooler ecosystem (.cool/.mcool/.scool) and the Juicer .hic format. The central idea is that cooler stores raw observed counts and computes everything else on the fly, whereas .hic bakes normalization and expected vectors into a sealed binary; conversion between them is therefore never lossless in both directions. The skill teaches the single-resolution mcool URI, the divisive-vs-multiplicative weight-naming rule that silently corrupts balanced values, what survives .hic conversion (FRAG matrices and some norms do not), how to coarsen correctly by summing raw counts and re-balancing, and the chrom-naming and bin-table provenance that decide whether two coolers can even be compared.

## Prerequisites

```bash
pip install cooler hic2cool hictk bioframe numpy pandas
# hictk also ships as a conda package and a static binary:
conda install -c bioconda hictk cooler
```

A `.cool` must be balanced (carry a `weight` column) before most downstream analysis. Balancing mechanics live in matrix-operations; this skill loads, converts, and inspects.

## Quick Start

Tell your AI agent what you want to do:
- "Load the 10kb resolution from my mcool and pull out chr1"
- "Convert this .hic file to a multi-resolution mcool, keeping the KR norms"
- "Build a cooler from my pairs file at 10kb"
- "Why is my balanced matrix all NaN?"
- "My fetch('1') returns nothing - what is wrong?"
- "Coarsen this 10kb cooler to 100kb correctly"

## Example Prompts

### Loading and selecting a resolution
> "I have an mcool. List its resolutions, open the 25kb level, tell me whether it is balanced, and return the balanced chr2 matrix as a numpy array."

> "Fetch the inter-chromosomal (trans) submatrix between chr1 and chr8 from my balanced cooler."

### Converting formats
> "Convert in.hic to an mcool with all resolutions using hic2cool, and make sure the Juicer norm vectors come through as divisive so balancing is correct."

> "This .hic is large - use hictk to convert just the 5kb and 10kb resolutions to a cooler, and warn me if any matrices are FRAG-binned and would be lost."

### Building from pairs or a matrix
> "Build a 10kb cooler from my flipped, deduped .pairs file using cooler cload, with hg38 chromsizes."

> "I have an in-memory numpy contact matrix for chr1 - write it to a cooler without a slow per-bin-pair loop."

### Debugging provenance
> "My two coolers give different balanced values from the same .hic - figure out if it is a hic2cool version mismatch."

> "Coarsen my 10kb cooler to 50kb and 100kb the correct way (raw sum then re-balance), not by summing balanced pixels."

## What the Agent Will Do

1. Resolve the file to a single-resolution cooler URI (`file.mcool::/resolutions/<bp>`), listing resolutions if needed.
2. Inspect the bin table for a `weight` column and report whether the chosen resolution is balanced.
3. Fetch raw or balanced pixels, handling the KR/VC/VC_SQRT divisive auto-rule and accepting masked-bin NaNs as correct.
4. For conversion, pick hic2cool (norm-faithful, BP only) or hictk (fast, no FRAG/asymmetric) and flag what will not survive.
5. For building, binnify the assembly and write the cooler vectorized (or via `cooler cload` for pairs).
6. For coarsening, use `zoomify`/`coarsen` so raw counts are summed and balancing is re-run per resolution.
7. Check chrom naming and bin-table identity before any cross-cooler comparison.

## Tips

- Always use the `::/resolutions/<bp>` URI; never hand a bare `.mcool` to a downstream tool.
- The `weight` column name is not cosmetic: cooler `weight` is multiplicative; KR/VC/VC_SQRT are auto-divisive. Do not rename Juicer norms to `weight`.
- Record the hic2cool version - the 0.5.0 boundary flips how norms are stored.
- All-NaN balanced rows on a balanced file are masked low-coverage bins, not a bug; an outright error means the file is unbalanced.
- Balanced pixels cannot be summed to coarsen - sum raw counts and re-ICE per resolution (`cooler zoomify --balance`).
- FRAG-binned `.hic` matrices do not survive any converter; re-bin in base pairs from the pairs.
- `chr1` vs `1` silently zeros every join; pin one chrom naming and one assembly across the cooler, FASTA, tracks, and blacklist.

## Related Skills

- contact-pairs - Produces the flipped/deduped .pairs that cooler cload bins into a matrix
- matrix-operations - Balancing (ICE/KR), expected, and O/E that operate on the loaded cooler
- hic-visualization - Render the matrices loaded here
- compartment-analysis - Consumes the balanced cooler at compartment resolution
- read-alignment/bwa-alignment - Produces the BAM that pairtools converts to .pairs
- genome-intervals/bed-file-basics - Chrom-naming and BED handling for blacklists/annotation joins
- single-cell/scatac-analysis - Single-cell chromatin context for .scool / scHi-C
