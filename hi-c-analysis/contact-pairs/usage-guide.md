# Contact Pairs Processing - Usage Guide

## Overview

This skill turns raw Hi-C, Micro-C, or Omni-C FASTQ into a deduplicated, filtered `.pairs` file with pairtools, and - just as importantly - decides whether the library worked. A clean `.pairs` file is a list of 5'-canonicalized, single-ligation, deduplicated, both-uniquely-mapped contacts; every step strips a specific proximity-ligation artifact. The skill covers the `bwa mem -SP5M` alignment idiom (mates aligned independently through ligation junctions), the pairtools parse/sort/dedup/select/stats/restrict/phase pipeline, the walks-policy and parse-vs-parse2 choices, restriction-fragment handling for single- and multi-enzyme protocols, allele-specific phasing, and the library-QC readout (cis/trans ratio, orientation balance, duplicate complexity) that confirms the experiment succeeded before any matrix is built.

## Prerequisites

```bash
conda install -c bioconda pairtools bwa bwa-mem2 chromap samtools cooler bcftools
# Optional complexity-curve estimation:
conda install -c bioconda preseq
```

A reference FASTA (indexed with `bwa index` or `bwa-mem2 index`) and a `chrom.sizes` file are required. For allele-specific work, a phased VCF and a diploid reference are also needed.

## Quick Start

Tell your AI agent what you want to do:
- "Align my Hi-C FASTQ and turn it into a deduplicated pairs file"
- "Run pairtools parse, sort, dedup, and select on this BAM"
- "Tell me if my Hi-C library is good - check cis/trans and duplicates"
- "Process Micro-C reads without applying a 1kb distance cut"
- "Handle Pore-C multi-way contacts instead of collapsing them to pairwise"
- "Set up allele-specific phasing into two coolers"

## Example Prompts

### Standard processing
> "I have Hi-C R1/R2 FASTQ and an hg38 reference. Align with bwa mem -SP5M, parse to pairs with pairtools, sort, dedup keeping both UU and UC, and write a valid pairs file ready for cooler cload."

### Library QC decision
> "Run pairtools stats on my deduplicated pairs and interpret it: what fraction is long-range cis, what is the trans noise floor, do the orientations converge to 25% above 1kb, and is the duplicate rate telling me to sequence deeper or stop?"

### Protocol-specific handling
> "This is a Micro-C library - process it the restriction-agnostic way and do NOT apply Hi-C's 1kb min-distance cut, since the nucleosome ladder below 1kb is the signal I want."

> "These are Arima dual-enzyme Hi-C reads and I need fragment-level assignment - make sure the digest file encodes all four junction motifs, not just DpnII."

### Multi-way and allele-specific
> "My data is Pore-C concatemers - use a walks policy that keeps the multi-way contacts instead of collapsing them to pairwise."

> "Set up allele-specific Hi-C: build a diploid reference from my phased VCF, align with suboptimal hits, run pairtools phase, and split into haplotype-1, haplotype-2, and trans pairs for two separate coolers."

### Choosing the fast path
> "I have hundreds of Hi-C libraries to scan quickly - use chromap --preset hic for integrated alignment, dedup, and pairs output instead of the bwa-mem2 + pairtools chain."

## What the Agent Will Do

1. Align R1 and R2 as independent single-end reads with `bwa mem -SP5M` (or bwa-mem2, or chromap --preset hic) so in-read ligation junctions are reconstructed and long-range/trans contacts are preserved.
2. Parse alignments into a 5'-canonical `.pairs` file, choosing the walks policy by protocol (5unique for standard pairwise Hi-C, all/parse2 for multi-way).
3. Sort (flip to upper-triangular, 5'-canonical), then deduplicate, separating optical from PCR duplicates with by-tile statistics.
4. Select valid pairs by type (keep UU and rescued UC), MAPQ, and a protocol-appropriate distance cut.
5. Run `pairtools stats` and interpret the cis/trans ratio, orientation-vs-distance balance, and duplicate-complexity curve to judge whether the library worked.
6. Optionally assign restriction fragments (for sub-kb/capture/bench QC) or phase pairs to haplotypes (for allele-specific folding).

## Tips

- Never drop `-SP5` from the aligner flags - `-SP` aligns mates independently (proper-pair logic destroys long-range/trans contacts) and `-5` anchors pairtools' 5' convention. `-M` is legacy compatibility only.
- Keep both UU and UC; selecting only UU silently discards every rescued chimeric ligation.
- Sort before dedup - dedup on a non-canonical file under-collapses and makes the library look better than it is.
- Read % long-range cis as the one-number quality metric and trans as the noise floor, but never apply a human trans threshold to a small (microbial/yeast) genome.
- Derive the min-distance cut from the orientation-vs-distance plot, not a hardcoded 1kb - Micro-C signal lives below 1kb.
- Use `--bytile-dups` to separate optical from PCR duplicates; only the PCR fraction reflects complexity.
- Modern pipelines skip fragment filtering on purpose; reach for `pairtools restrict` only for sub-kb, capture-Hi-C, or bench-QC work, and encode all enzyme motifs for Arima/Hi-C 3.0.

## Related Skills

- hic-data-io - Bins the deduplicated valid pairs into a .cool/.mcool matrix
- matrix-operations - Balancing and O/E that the binned pairs feed into
- hic-visualization - Render contact maps from the resulting cooler
- read-alignment/bwa-alignment - Aligner upstream; this skill adds the Hi-C `-SP5M` idiom
- alignment-files/duplicate-handling - General duplicate-marking context for the pairtools dedup step
- genome-intervals/bed-file-basics - Coordinate/digest BED handling for restriction fragments and anchors
- genome-assembly/scaffolding - Same Hi-C reads used to order contigs into chromosomes
