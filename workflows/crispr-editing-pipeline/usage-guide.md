# CRISPR Editing Pipeline - Usage Guide

## Overview

Orchestrates a complete CRISPR editing experiment from target gene to delivery-ready, validatable constructs. It sequences guide design, off-target assessment, edit-modality selection (knockout, base editing, prime editing, HDR knock-in), and template/donor design, applying a QC checkpoint at every handoff. The workflow's value is the order of operations, the pivotal edit-modality decision, and the cross-cutting traps; it routes each step's mechanics to the five genome-engineering skills and does not re-implement their scoring.

## Prerequisites

```bash
pip install biopython pandas matplotlib
# CRISPOR (web/CLI) for context-valid on-target + off-target; Cas-OFFinder / CRISPRme for
#   off-target nomination; BE-Hive for base-editor outcomes; PRIDICT/DeepPrime for pegRNAs;
#   CRISPResso2 for amplicon validation. Each step's detail lives in its genome-engineering skill.
```

## Quick Start

Tell your AI agent what you want to do:
- "Design a complete CRISPR experiment to knock out BRCA1"
- "Plan an end-to-end base-editing experiment to install a stop codon"
- "Design an HDR knock-in to tag MYC with GFP, with validation"
- "Walk me through guide -> off-target -> donor for correcting a point mutation"

## Example Prompts

### Knockout
> "Design guides to knock out KRAS, rank them by frameshift outcome, check off-targets, and give me validation primers."

### Base editing
> "Plan a CBE experiment to install a premature stop in this gene without a double-strand break, and report the bystander/purity profile."

### Knock-in
> "Design the full workflow to insert a FLAG tag at the C-terminus -- guide, off-target check, donor with a blocking mutation, and genotyping primers."

### Modality choice
> "I want to correct a G>A point mutation -- which modality (base, prime, or HDR), and why?"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Target sequence | FASTA/string | Gene or region to edit (with reading frame for coding edits) |
| Edit goal | String | knockout, transition, transversion/indel, or insertion |
| Insert sequence (knock-in) | FASTA | Sequence to insert |
| Target position (point edit) | Integer | Position of the base to change |
| Cell type / delivery | String | Sets valid on-target model, hard filters, and donor format |

## What the Agent Will Do

1. Stage 1 -- design guides and shortlist 3-6 by predicted frameshift outcome in an early constitutive exon (-> grna-design).
2. Stage 2 -- assess off-targets, escalating predicted -> detected -> validated; variant-aware for therapeutics (-> off-target-prediction).
3. Decide the edit modality (knockout / base / prime / HDR) from the decision tree.
4. Stage 3 -- design the modality-specific construct (frameshift guide; base-editor window+purity; pegRNA panel; HDR donor with a codon-checked blocking mutation).
5. Stage 4 -- design validation/genotyping primers and quantify outcomes by amplicon deep sequencing, reporting purity and the limit of detection.

## Editing Strategy Decision Tree

```
What edit do you need?
  |
  +-- Gene knockout (any frameshift) --> nuclease + NHEJ (grna-design)
  +-- Knockout without a DSB / non-dividing / multiplex --> base-editor stop or splice disruption
  +-- C*G->T*A or A*T->G*C transition --> base editing (CBE/ABE)
  +-- C->G transversion --> CGBE
  +-- Other transversion / small indel / combined --> prime editing
  +-- Tag / reporter / allele replacement (cycling cells) --> HDR knock-in
  +-- Large insertion / post-mitotic cells --> HDR (AAV/HITI) or PE+integrase (PASTE/twinPE)
```

## Tips

- The pivotal decision is the edit modality -- choose it first; it reframes every downstream step.
- On-target activity and off-target specificity are separate axes -- never pick a guide on activity alone.
- Efficient editing is not knockout -- rank by frameshift fraction, target an early constitutive NMD-competent exon, and verify protein.
- For base editing, report product purity (target edit and no bystander), not a lone efficiency number.
- An HDR donor without a codon-checked blocking mutation self-destructs -- the edit is re-cut and reads out as a plain indel.
- Off-target claims are predicted, detected, or validated -- keep them distinct and state the limit of detection.
- Each step has a dedicated skill with the full nuance; this workflow is the map and the checkpoints.

## Related Skills

- genome-engineering/grna-design - Guide design and outcome-aware knockout ranking
- genome-engineering/off-target-prediction - Specificity assessment and the evidence ladder
- genome-engineering/base-editing-design - CBE/ABE window, bystander purity, off-target classes
- genome-engineering/prime-editing-design - pegRNA panel design and PE system selection
- genome-engineering/hdr-template-design - Donor format and codon-checked blocking mutation
- crispr-screens/crispresso-editing - Quantify and validate editing outcomes from amplicon reads
- crispr-screens/library-design - Scale single-gene design to a pooled screen
- primer-design/primer-basics - Design the genotyping/amplicon primers around the edit
- primer-design/primer-specificity - Confirm the genotyping amplicon is unique near paralogs/pseudogenes
