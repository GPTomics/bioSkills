# genome-engineering

## Overview

Design CRISPR guides, predict off-targets, and create templates for genome editing experiments including Cas9/Cas12a knockouts, prime editing, base editing, and HDR knock-ins.

**Tool type:** mixed | **Primary tools:** crisprscan, Cas-OFFinder, PrimeDesign, primer3-py

## Skills

| Skill | Description |
|-------|-------------|
| grna-design | Design guide RNAs with activity scoring using CRISPRscan, CHOPCHOP |
| off-target-prediction | Predict off-target sites with Cas-OFFinder, CFD scores |
| prime-editing-design | Design pegRNAs for prime editing with PrimeDesign |
| base-editing-design | Design CBE/ABE guides with BE-Designer, BE-Hive |
| hdr-template-design | Design HDR donor templates with primer3-py |

## Example Prompts

- "Design guides to knock out BRCA1"
- "Find the best sgRNAs targeting exon 3 of TP53"
- "Check off-target sites for my guide sequence ATCGATCGATCGATCGATCG"
- "Design a pegRNA to correct the G551D mutation in CFTR"
- "Create a base editing guide for the C282Y mutation"
- "Design an HDR template to add a GFP tag to MYC"

## Requirements

```bash
pip install crisprscan primer3-py biopython
# Cas-OFFinder: download from http://www.rgenome.net/cas-offinder/
```

## Related Skills

- **crispr-screens** - Analyze pooled CRISPR screens after editing
- **primer-design** - General PCR primer design with primer3
- **variant-calling** - Detect edits in sequencing data
