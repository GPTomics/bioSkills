# Guide RNA Design - Usage Guide

## Overview

Designs and ranks CRISPR guide RNAs for Cas9/Cas12a gene knockout. The skill scans a target for the chosen nuclease's PAM, applies hard filters (Pol-III TTTT terminator, 5' G, GC), ranks on-target activity with the model that matches the delivery context, picks the cut site by exon/transcript biology, and predicts the indel/frameshift outcome -- because the indel *distribution*, not the cutting rate, decides whether a gene is knocked out. On-target activity and off-target specificity are two separate axes; this skill owns the activity/site-selection axis and hands off specificity to off-target-prediction.

## Prerequisites

```bash
pip install biopython
# CRISPOR (the de facto standard aggregator) is web (crispor.tefor.net) or a CLI clone;
#   it selects the context-appropriate on-target score and nominates off-targets per genome.
# Outcome prediction: inDelphi, FORECasT, Lindel are web servers / authors' code.
```

There is no maintained `pip install crisprscan` package with a generic `.score()` API -- CRISPRscan is a web tool/model. Route real on-target scoring to CRISPOR rather than hand-rolling a scoring matrix.

## Quick Start

Tell your AI agent what you want to do:
- "Design guides to knock out BRCA1"
- "Find SpCas9 guides in an early constitutive exon of TP53 and rank them"
- "My guides are for an in-vitro-transcribed RNP injection -- which on-target score applies?"
- "Design Cas12a guides for this AT-rich target"
- "Which of these guides has the most frameshift-prone repair outcome?"

## Example Prompts

### Knockout guide design

> "I have the coding sequence of KRAS. Find SpCas9 NGG guides in the first two constitutive coding exons, drop any with a TTTT or GC outside 40-70%, and give me the top 5 ranked so I can run them through CRISPOR for on-target scoring."

> "Design guides to knock out MYC, but explain which exon to target and why -- I keep getting edited cells with no phenotype."

### Context and nuclease choice

> "My guides will be expressed from a U6 lentiviral vector in a pooled screen -- which on-target model is valid, and should I tile the functional domain instead of the 5' exon?"

> "This target is AT-rich and has no good NGG. What are my options?"

> "I need AAV delivery in vivo -- design SaCas9-compatible guides and note the PAM change."

### Outcome and validation

> "Rank these candidate guides by frameshift probability, not just cutting efficiency."

> "I confirmed a frameshift but the protein is still there -- what could be happening?"

## What the Agent Will Do

1. Establish the delivery context (nuclease, U6/lentiviral vs in-vitro T7/RNP, cell type, genetic background) -- this sets which PAMs exist, which on-target score is valid, and which hard filters apply.
2. Enumerate candidate spacers adjacent to valid PAMs on both strands.
3. Apply hard filters (reject TTTT for U6/H1; prepend a 5' G if absent; flag GC extremes).
4. Route on-target ranking to CRISPOR (context-appropriate model), treating the score as a shortlist signal, not an oracle.
5. Choose the cut site by exon biology: early, constitutive, NMD-competent, away from splice sites; tile the conserved domain for screens.
6. Predict the editing outcome (Bae out-of-frame score; inDelphi/FORECasT/Lindel) and prefer frameshift-rich, dominant-outcome guides.
7. Return 3-6 candidates and hand off to off-target-prediction before ordering.

## Tips

- On-target scores are weak (Spearman ~0.4 across context) and context-locked -- use the model that matches the delivery context, rank to a shortlist, and design 3-6 guides.
- Efficient editing is not knockout: ~1/3 of indels are in-frame and ~1/3 of verified frameshift KOs retain protein -- rank by frameshift fraction and verify at the protein level.
- "Early exon" is a half-truth: early but not the start ATG, constitutive in the cell type of interest, not the last exon (escapes NMD), and away from splice junctions.
- Reaching for SpRY/xCas9 to hit a no-NGG site costs activity and specificity -- default to WT-SpCas9-NGG and escalate only when needed.
- In a non-reference background (patient/hybrid/cancer line), screen for SNPs under the guide/PAM or one allele drops out silently.
- Cas12a's self-processing crRNA array makes multi-gene knockout easy from one transcript.

## Related Skills

- off-target-prediction - Check genome-wide specificity after on-target design
- base-editing-design - DSB-free knockout via premature stop or splice disruption
- prime-editing-design - Scarless small edits without a double-strand break
- hdr-template-design - Design the donor for a precise knock-in
- crispr-screens/library-design - Pool guides into a screening library
- crispr-screens/crispresso-editing - Quantify editing outcomes from amplicon reads
- primer-design/primer-basics - Validation/genotyping primers around the cut
- primer-design/primer-specificity - Confirm genotyping primers are unique near paralogs/off-targets
- genome-intervals/gtf-gff-handling - Exon coordinates to restrict guide placement
