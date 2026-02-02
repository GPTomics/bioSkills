# CRISPR Editing Pipeline - Usage Guide

## Overview

This workflow designs complete CRISPR editing experiments from target selection to delivery-ready constructs. It supports multiple editing strategies: knockouts (frameshift), base editing (CBE/ABE), and knockins (HDR). Each path includes guide design, off-target assessment, and construct generation.

## Prerequisites

```bash
pip install crisprscan biopython pandas numpy matplotlib primer3-py

conda install -c bioconda cas-offinder
```

**Optional for advanced scoring:**
- CRISPRscan web tool (https://www.crisprscan.org/)
- Benchling (commercial, integrates multiple algorithms)

## Quick Start

Tell your AI agent what you want to do:
- "Design CRISPR guides to knock out BRCA1"
- "Create a base editing experiment to introduce a stop codon"
- "Design HDR template to tag MYC with GFP"
- "Find the best guides for editing exon 3 of TP53"

## Example Prompts

### Knockout
> "Design guides to knock out gene X with frameshift"

> "Find CRISPR targets in the first coding exon of KRAS"

### Base editing
> "Design CBE guides to convert C to T at position 234"

> "Create an ABE experiment to correct the G>A mutation"

### Knockin
> "Design HDR template to insert FLAG tag at the N-terminus"

> "Create donor template for CRISPR knockin with 800bp homology arms"

### Off-target analysis
> "Check off-targets for this guide sequence"

> "Find guides with high specificity for EGFR"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Target sequence | FASTA/string | Gene or region to edit |
| Edit type | String | knockout, base_edit, or knockin |
| Insert sequence (knockin) | FASTA | Sequence to insert |
| Target position (base edit) | Integer | Position of base to edit |

## What the Agent Will Do

1. **Guide Design** - Identifies all PAM sites and scores for on-target activity
2. **Off-Target Assessment** - Searches genome for potential off-target sites
3. **Strategy Selection** - Branches based on edit type
4. **Construct Design** - Generates guide sequences, donors, or editing windows
5. **Validation Design** - Creates primers for confirming edits
6. **Visualization** - Maps guides onto gene structure

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Activity score threshold | 0.6 | Minimum guide activity |
| Specificity score threshold | 0.7 | Maximum off-target penalty |
| GC content range | 40-70% | Optimal guide composition |
| Homology arm length | 800bp | For plasmid HDR donors |
| Base editing window | positions 4-8 | Optimal for BE3/BE4 |

## Editing Strategy Decision Tree

```
What type of edit do you need?
    |
    +-- Complete gene disruption --> Knockout (Cas9 + NHEJ)
    |
    +-- Single base change
    |       |
    |       +-- C>T or G>A --> CBE (BE3/BE4)
    |       |
    |       +-- A>G or T>C --> ABE (ABE7.10/ABE8)
    |
    +-- Insert sequence (tag, reporter) --> HDR knockin
    |
    +-- Complex edit (insertion + deletion) --> Prime editing
```

## Tips

- **Target early exons**: For knockouts, target early in the coding sequence
- **Avoid repetitive regions**: High off-target risk in repeats
- **Verify PAM orientation**: Ensure NGG is on the correct strand
- **Test multiple guides**: Design 3-5 guides per target for backup
- **Consider Cas9 variants**: Use high-fidelity variants for sensitive applications
- **HDR timing matters**: Co-deliver donor with RNP for best efficiency
- **ssODN for small edits**: Single-strand oligos work well for <50bp changes
