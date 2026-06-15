# genome-engineering

## Overview

Design CRISPR edits and their safety checks: guide RNAs for knockout, off-target prediction, base editing (CBE/ABE), prime editing (pegRNAs), and HDR donor templates for knock-ins. The category covers the design side -- before the experiment -- and routes data analysis after the experiment to crispr-screens.

**Tool type:** mixed | **Primary tools:** CRISPOR, Cas-OFFinder, PrimeDesign, BE-Hive, primer3-py

## Skills

| Skill | Description |
|-------|-------------|
| grna-design | Design and outcome-rank sgRNAs for Cas9/Cas12a knockout (context-valid on-target scoring, frameshift/exon biology) |
| off-target-prediction | Nominate and assess off-targets (Cas-OFFinder/CFD, variant-aware, empirical assays, high-fidelity nucleases) |
| base-editing-design | Design CBE/ABE guides for transitions by window/bystander purity, editor variant, and the three off-target classes |
| prime-editing-design | Design pegRNA panels and choose the PE2/PE3b/PE5max/PE7 system (PBS/RTT, MMR evasion, epegRNA) |
| hdr-template-design | Design knock-in donors (ssODN/lssDNA/dsDNA/AAV) with a codon-checked blocking mutation |

## Example Prompts

- "Design guides to knock out BRCA1 and rank them by frameshift outcome"
- "Check off-target sites for my guide, including bulges, and tell me if it's therapeutic-grade"
- "I want a G-to-A correction -- base, prime, or HDR, and why?"
- "Design a pegRNA panel to correct the G551D mutation in CFTR"
- "Design a base-editing guide to install a stop codon without a double-strand break"
- "Design an HDR donor to add a GFP tag to MYC, with a blocking mutation"

## Requirements

```bash
pip install biopython primer3-py pandas
# CRISPOR (web/CLI) for context-valid on-target + off-target nomination
# Cas-OFFinder: standalone binary (https://github.com/snugel/cas-offinder) + a reference genome FASTA
# PrimeDesign: Docker (pinellolab/primedesign); BE-Hive / PRIDICT / DeepPrime: web or clone
```

## Related Skills

- **crispr-screens** - Analyze pooled screens and amplicon editing outcomes after the experiment
- **primer-design** - PCR/qPCR primer design, validation, and genome specificity (genotyping amplicons)
- **variant-calling** - Detect and annotate edits in sequencing data
- **genome-intervals** - Exon/feature coordinates to restrict guide placement
