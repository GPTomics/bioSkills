# UMI Processing - Usage Guide

## Overview
Unique Molecular Identifiers (UMIs) are random sequences attached to molecules before PCR, enabling collapse by (mapping coordinate + UMI) to count ORIGINAL molecules and distinguish PCR duplicates from biological ones. umi_tools counts molecules (directional dedup); fgbio builds error-corrected single-strand or duplex consensus reads for rare-variant detection. The pipeline order is fixed: extract the UMI on the FASTQ, align, then dedup or call consensus (duplicate identity needs mapping coordinates).

## Prerequisites
```bash
conda install -c bioconda umi_tools fgbio samtools star subread
```

## Quick Start
Tell your AI agent what you want to do:
- "Extract UMIs from read 1 and deduplicate after alignment"
- "Count molecules in my UMI-tagged targeted panel"
- "Build duplex consensus reads for my ctDNA library"
- "Do I need to dedup my standard RNA-seq?" (answer: only if it has UMIs)

## Example Prompts

### Molecule counting
> "Extract 8 bp UMIs from read 1, align with STAR, and deduplicate with the directional method"

> "Deduplicate my aligned BAM using UMI information and report the rate"

### Single-cell
> "Process my 10x library: 16 bp cell barcode + 12 bp UMI on read 1"

> "I already have a CellRanger matrix; should I re-deduplicate?" (answer: no)

### Error correction
> "Build single-strand consensus reads for my targeted panel"

> "Build duplex consensus reads and filter them for ctDNA rare-variant detection"

## What the Agent Will Do
1. Extract the UMI from the read into the header (or RX tag) before alignment
2. Align with UMI preserved, sort, and index
3. Count molecules with umi_tools directional dedup, OR group + call consensus with fgbio
4. For ctDNA, group with the paired strategy and call duplex consensus, then FilterConsensusReads
5. Report dedup/consensus statistics and the edit-distance diagnostic

## UMI pattern syntax (umi_tools extract)
| Pattern | Description |
|---------|-------------|
| `N` | UMI base |
| `C` | Cell barcode base |
| `X` | Fixed/known base reattached to the read (true discard uses the regex `(?P<discard_N>...)` group) |
| `NNNNNNNN` | 8 bp UMI |
| `CCCCCCCCCCCCCCCCNNNNNNNNNNNN` | 10x 3' v3: 16 bp cell barcode + 12 bp UMI |

## Tips
- Extract UMIs before alignment; they must leave the aligned sequence (header or RX tag), or the aligner mismaps them.
- Directional is the default and the right choice; it models UMI errors via the count-gradient rule. `cluster` over-merges and `unique` over-counts, so neither is a "high-error" upgrade.
- Do NOT re-deduplicate CellRanger/STARsolo output; it is already UMI-collapsed.
- Do NOT dedup non-UMI bulk RNA-seq at all; duplicate coordinates there are biological signal.
- Single-strand consensus halves errors; duplex consensus is required for sub-0.1% VAF (ctDNA/MRD) and costs ~2x raw reads.
- Always run FilterConsensusReads after calling consensus.
- Deep amplicon needs longer UMIs because every molecule shares coordinates, so the UMI alone must separate them.

## Resources
- [umi_tools Documentation](https://umi-tools.readthedocs.io/)
- [fgbio Documentation](https://fulcrumgenomics.github.io/fgbio/)
- [umi_tools Publication](https://doi.org/10.1101/gr.209601.116)

## Related Skills
read-qc/fastp-workflow - UMI extraction folded into preprocessing
read-qc/rnaseq-qc - Why non-UMI bulk RNA-seq must NOT be deduplicated
alignment-files/duplicate-handling - Coordinate dedup for non-UMI DNA
single-cell/preprocessing - scRNA-seq UMI matrices and downstream
liquid-biopsy/ctdna-mutation-detection - Duplex consensus for rare-variant detection
