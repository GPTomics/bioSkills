# CRISPResso2 Usage Guide

## Overview

CRISPResso2 analyzes amplicon sequencing data from CRISPR editing experiments. It quantifies indels, HDR efficiency, and characterizes editing outcomes.

## Analysis Types

| Tool | Use Case |
|------|----------|
| CRISPResso | Single amplicon |
| CRISPRessoBatch | Multiple samples, same amplicon |
| CRISPRessoPooled | Multiple amplicons per sample |
| CRISPRessoWGS | Whole genome off-target |
| CRISPRessoCompare | Compare conditions |

## Key Parameters

### --quantification_window_size
Window around cut site for quantifying edits. Default: 1.

### --quantification_window_center
Offset from cut site. Default: -3 (3bp 5' of PAM).

### --exclude_bp_from_left/right
Ignore bases at amplicon ends (primer artifacts).

## Editing Types

### NHEJ (Indels)
- Insertions and deletions from non-homologous end joining
- Quantified in "Modified" percentage

### HDR
- Homology-directed repair using template
- Requires --expected_hdr_amplicon_seq

### Base Editing
- C>T or A>G conversions
- Use --base_editor_output

### Prime Editing
- Precise insertions/deletions
- Use prime editing parameters

## Quality Considerations

- Min read depth: 1000x recommended
- Check mapping percentage (>80%)
- Review indel size distribution
- Validate with negative controls

## References

- CRISPResso2: doi:10.1038/s41587-019-0032-3
- Documentation: https://crispresso.pinellolab.partners.org/
