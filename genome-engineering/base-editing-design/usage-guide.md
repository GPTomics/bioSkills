# Base Editing Design - Usage Guide

## Overview

Designs cytosine (CBE, C-to-T) and adenine (ABE, A-to-G) base-editor guides for transition mutations without a double-strand break. The skill writes the edit as a base-pair change to pick the editor family, positions the target base at the activity peak of the window (~positions 5-7), minimizes in-window bystander edits for product purity, reads dinucleotide context, selects the editor variant, and frames the three off-target classes and knockout-by-base-editing. Its discipline is that the job is product *purity* (the fraction of alleles with the target edit and no bystander), not editing efficiency.

## Prerequisites

```bash
pip install biopython
# BE-Hive: web tool (crisprbehive.design) + clone-only repos (maxwshen/be_predict_efficiency,
#   be_predict_bystander) with old pinned deps (scikit-learn 0.20.3, BioPython 1.73) -- NOT pip.
# BE-Designer / BE-Analyzer: web (rgenome.net). CRISPResso2: amplicon editing-spectrum analysis.
```

## Quick Start

Tell your AI agent what you want to do:
- "Design a CBE guide to install the C282Y correction, with no bystanders"
- "I want a G-to-A change here -- which base editor do I use?"
- "Knock out this gene without a double-strand break"
- "Should I use CBE or ABE for this edit, and which is cleaner for off-targets?"
- "My guide is 80% efficient but the genotypes are a mess -- why?"

## Example Prompts

### Choosing the editor and window

> "I want to install the A*T->G*C change for a pathogenic variant. Confirm this is an ABE job, find a PAM that puts the A at position 5-7, and tell me if there are bystander A's I should worry about."

> "I need a G-to-A edit -- explain why this is a CBE job and design the guide on the right strand."

### Purity and variant choice

> "Find CBE guides where the target C is the only editable C in the window, and exploit GC-context disfavor for any bystanders."

> "We're switching from ABE7.10 to ABE8e for activity -- what happens to my bystander profile?"

### Knockout and off-targets

> "Knock out this gene in non-dividing cells without a DSB -- design a CRISPR-STOP or splice-disruption edit."

> "How do I assess off-targets for a base editor, and which class matters for CBE vs ABE safety?"

## What the Agent Will Do

1. Write the edit as a base-pair change to pick the family (C*G->T*A = CBE; A*T->G*C = ABE; C->G = CGBE; else prime editing).
2. Find a PAM that lands the target base at the window peak (~5-7), checking both strands.
3. List in-window bystander C's/A's and read the 5' context of each (TC favored / GC disfavored for APOBEC1).
4. Select the editor variant (BE4max/ABEmax default; ABE8e for activity at a purity cost; SECURE/TadCBE for off-target; SpG/SpRY when no NGG positions the base).
5. Report the predicted genotype spectrum (BE-Hive/DeepBE as a relative ranker), not a lone efficiency number.
6. For knockout, design an early premature stop (CRISPR-STOP/iSTOP) or splice-site disruption.
7. Frame the three off-target classes and route class-1 to off-target-prediction; recommend NGS validation.

## Tips

- Rank by purity (target edit and no bystander), not by editing efficiency -- "80% edited" can be a genotype soup.
- Write the edit as a base-pair change first; a "G->A" edit is a CBE job (C->T on the complement), and getting the strand wrong is the most common editor-choice error.
- The window is an activity gradient, not a box, and its width is editor-specific -- ABE8e is wider (~3-11), so re-check bystanders after any editor switch.
- Three off-target classes: class 1 (Cas-dependent, guide-directed) is what predictors/GUIDE-seq cover; classes 2 (Cas-independent DNA) and 3 (RNA) are guide-invisible and solved only by choosing a cleaner editor (SECURE/TadCBE); CBE was historically dirtier than early ABE.
- Knockout without a DSB: prefer an early premature stop or splice-site disruption; the KO is only as complete as the editing, so verify protein.
- BE-Hive is a web tool / clone, not a pip install -- predictions are relative rankers in their training regime and say nothing about classes 2-3; always validate by NGS.

## Related Skills

- grna-design - Cas-domain on-target/class-1 guide scoring this skill builds on
- prime-editing-design - For transversions (except C->G), indels, and multi-base replacements
- off-target-prediction - Owns the class-1 Cas-dependent off-target assays/prediction
- crispr-screens/base-editing-analysis - Quantify editing outcomes from NGS after the experiment
- crispr-screens/crispresso-editing - Amplicon editing-spectrum quantification
- variant-calling/variant-annotation - Annotate the consequence of the installed edit
