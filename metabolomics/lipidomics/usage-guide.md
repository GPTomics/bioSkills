# Lipidomics Usage Guide

## Overview

Lipidomics is a combinatorial-structure problem wearing a quantification problem's clothes. This skill keeps the agent honest about two things software routinely overstates: the structural-resolution level a lipid name actually carries (the shorthand separator `space` -> `_` -> `/` -> `(9Z)` encodes what was measured, and sn-position is almost never measured), and quantification (one isotope-labeled internal standard per class is non-negotiable because ESI response is head-group-dependent). It guards against in-source-fragment phantom lyso-lipids, sn over-claims, ether/plasmalogen ambiguity, and invalid cross-class comparisons.

## Prerequisites

```bash
# R (primary)
# BiocManager::install("lipidr")

# Nomenclature parsing
pip install pygoslin
```

Conceptual prerequisites: a quantified lipid table (Skyline, MS-DIAL, or LipidSearch export) with sample group annotations, and knowledge of which internal standards were spiked and at what stage. To trust any sn or double-bond claim, the agent needs to know whether EAD/OzID/PB/UVPD data were acquired - under routine CID they were not.

## Quick Start

Tell your AI agent what you want to do:
- "Load my lipidomics table, normalize within class, and find lipids changing between groups"
- "Canonicalize these lipid names and downgrade any that over-claim sn-position"
- "Check whether my elevated LPC signal is real or an in-source fragment of PC"
- "Design an internal-standard strategy for class-based quantification"

## Example Prompts

### Annotation honesty
> "Parse these lipid names through Goslin and report the structural-resolution level each one actually claims."
> "Re-emit my LipidSearch output at molecular-species level since we only ran CID - drop the sn slashes."
> "Flag any plasmalogen (P-) call that lacks vinyl-ether diagnostic evidence."

### Quantification
> "Normalize each lipid class to its matched internal standard, not a single global standard."
> "Is comparing PE to PC abundance valid in my dataset given the standards I used?"
> "Set up class-based quantification using my EquiSPLASH internal standards."

### Differential and enrichment analysis
> "Run differential lipid analysis between treatment and control and make a class-faceted volcano plot."
> "Test whether any lipid class, chain length, or unsaturation pattern is enriched among the changed lipids."

### Artifact triage
> "My LPC pool is unexpectedly high - check retention-time co-elution against the parent PCs."
> "This apparent odd-chain PC 33:1 - is it real or an isotope/in-source artifact?"

## What the Agent Will Do

1. Canonicalize names through Goslin and assign each an honest resolution level, defaulting down where evidence is absent.
2. Set up class-based internal-standard normalization (one IS per class) or PQN, and state semi-quant vs absolute.
3. Run `de_analysis` with an explicit contrast and class/chain-aware visualization.
4. Run lipid set enrichment (class / chain length / unsaturation) with `lsea`.
5. Triage in-source fragments, odd-chain sums, and ether/plasmalogen calls before reporting.
6. Export results with class, chain, and resolution-level metadata.

## Tips

- The separator is the claim: `space` = sum composition, `_` = chains known (sn unresolved), `/` = sn-resolved, `(9Z)` = double-bond position+geometry. Default down when in doubt.
- One isotope-labeled internal standard per class, spiked before extraction. Never quantify a class with another class's standard.
- An "LPC" eluting at a PC's retention time is an in-source fragment, not biology. Shotgun has no retention-time axis to run this test.
- "Number of lipids identified" is a vanity metric unless paired with the resolution level and confidence grade.
- A single-software species-level annotation needs orthogonal validation (ECN/retention-time consistency, a second adduct/polarity, CCS, or manual MS/MS) before it is trustworthy.
- Treat isomer resolution (sn, C=C) as an emerging capability; verify what the installed instrument/method actually supports rather than assuming.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream peak detection and feature extraction
- metabolomics/msdial-preprocessing - MS-DIAL alignment and deconvolution upstream of lipid annotation
- metabolomics/metabolite-annotation - General (non-lipid) annotation and confidence levels
- metabolomics/normalization-qc - Sample normalization and QC framing
- metabolomics/statistical-analysis - Multivariate stats on the lipid abundance matrix
