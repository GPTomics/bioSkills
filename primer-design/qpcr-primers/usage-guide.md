# qPCR Primer and Probe Design - Usage Guide

## Overview

This skill co-designs qPCR/RT-qPCR primers and hydrolysis (TaqMan) or molecular-beacon probes with primer3-py. The central idea it enforces: a qPCR assay is a quantitative measurement device, so its validity rests on amplification efficiency (90-110%, standard-curve slope -3.6 to -3.1) and a single product, and every qPCR-specific constraint (short 70-150 bp amplicon, tight Tm, zero dimers, probe Tm offset) exists to protect the efficiency the 2^-ddCq / Pfaffl math assumes. It covers the coupled probe rules (probe Tm 8-10 C above primers, no 5' G, C-rich strand, enforced with PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME=HNNNN), gDNA exclusion by exon-junction spanning or intron flanking and why processed pseudogenes defeat it, SYBR melt-curve QC, and reference-gene validation. Use it for TaqMan/SYBR assays, exon-spanning expression primers, and multiplex panels; route genome specificity to primer-specificity and dimer checks to primer-validation.

## Prerequisites

```bash
pip install primer3-py
```

- Input is the target sequence (cDNA for expression assays; mark the splice-junction position if cDNA-specific).
- Probe design needs PRIMER_PICK_INTERNAL_OLIGO=1 and raised PRIMER_INTERNAL_* Tm; primer3's internal Tm defaults equal the primer Tm, so they must be set 8-10 C higher.
- primer3 has no dedicated no-5'-G option; enforce it with PRIMER_INTERNAL_MUST_MATCH_FIVE_PRIME='HNNNN' (IUPAC H = not G) or a post-hoc filter.
- Design-level gDNA exclusion is leaky: processed pseudogenes defeat junction-spanning, so a genome specificity check, DNase, and a no-RT control remain mandatory (primer-specificity).
- A passing design is not a validated assay; a standard curve (efficiency, R^2) and controls (NTC, no-RT) are still required (MIQE).

## Quick Start

Tell your AI agent what you want to do:
- "Design a TaqMan assay for this cDNA: primers around 60 C, probe 8-10 C higher, amplicon under 150 bp, probe not starting with G"
- "Design SYBR Green primers for a 100 bp amplicon with minimal dimer potential"
- "Make exon-junction-spanning primers so I do not amplify genomic DNA"
- "Pick stable reference genes for my qPCR normalization"

## Example Prompts

### TaqMan assay
> "Design a TaqMan assay for my cDNA target: a Tm-matched primer pair near 60 C plus a hydrolysis probe with Tm about 68-70 C, amplicon under 150 bp, and make sure the probe does not start with G and sits on the C-rich strand."

### SYBR Green assay
> "Design SYBR Green primers for a 90 bp amplicon, tight on 3'-end cross-complementarity to avoid primer-dimers, and remind me to inspect the melt curve for a single peak."

### cDNA-specific (avoid gDNA)
> "Design exon-junction-spanning primers for this transcript so they will not amplify unspliced genomic DNA, and check whether this gene has a processed pseudogene that would defeat that."

### Multiplex and normalization
> "Design qPCR assays for my target and two reference genes with matched Tm and efficiency, screen all primers and probes for cross-dimers, and tell me how to validate the reference genes with geNorm."

## What the Agent Will Do

1. Load the target sequence and confirm the assay type (TaqMan vs SYBR, expression vs genomic).
2. Co-design primers and (for TaqMan) an internal probe: short amplicon, Tm-matched primers, probe Tm raised 8-10 C, no 5'-G via HNNNN.
3. For expression assays, add exon-junction spanning and flag the processed-pseudogene caveat.
4. Report the primer/probe set with Tm offsets and amplicon size, and note that SYBR needs a melt-curve check.
5. Route the chosen pair to primer-specificity (genome/pseudogene) and primer-validation (dimers), and state the standard-curve/efficiency and control requirements (MIQE).
6. For multiplex, check all primer/probe pairs for cross-dimers and recommend matched efficiencies and primer-limiting.

## Tips

- Keep the amplicon 70-150 bp; long amplicons lose efficiency and bend the standard curve.
- Raise the probe Tm 8-10 C above the primers so it is bound when the polymerase cleaves it; the internal defaults equal the primer Tm and must be changed.
- Forbid a 5'-G on the probe and prefer the C-rich strand; a 5'-G quenches the reporter.
- Do not trust exon-junction spanning alone for gDNA exclusion; check the genome for pseudogenes and keep DNase plus a no-RT control.
- For SYBR, the melt curve is the specificity readout; a low-Tm shoulder means primer-dimer.
- Validate reference genes with geNorm/NormFinder in the actual experimental conditions; never normalize to a single unvalidated gene.
- Efficiency is part of the quantification equation; run a standard curve and use Pfaffl if target and reference efficiencies are not matched.

## Related Skills

- primer-basics - Design fundamentals, Tm matching, and the constraint model
- primer-validation - Dimers/hairpins of primers and probe at reaction conditions
- primer-specificity - Genome/pseudogene specificity and gDNA exclusion checking
- sequence-manipulation/transcription-translation - Work with cDNA and reading frames
- differential-expression/deseq2-basics - Downstream analysis qPCR validates against
