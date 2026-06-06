# Metabolite Annotation Usage Guide

## Overview

Metabolite annotation turns untargeted LC-MS/MS features (m/z, RT, MS/MS) into named compounds, but the central job is honesty: every name must carry a confidence level. This skill guards against the field's recurring failures -- reporting a database hit as an identification, treating a high cosine score as proof, claiming a specific isomer MS/MS cannot resolve, and letting ambiguous annotations poison downstream pathway analysis.

## Prerequisites

```bash
pip install matchms
# SIRIUS 6 (in-silico formula/structure/class): https://v6.docs.sirius-ms.io/ (free academic account/license)
# MetFrag (transparent in-silico fragmenter): MetFragCommandLine.jar
```

Conceptual prerequisites: a feature table from upstream preprocessing (metabolomics/xcms-preprocessing or metabolomics/msdial-preprocessing) with ion families already collapsed (adducts, isotopes, in-source fragments); the ion mode and expected adducts; and the distinction between annotation (a hypothesis with a level) and identification (Level 1, an in-house standard).

## Quick Start

Tell your AI agent what you want to do:
- "Match my MS/MS spectra against a reference library and report the matched-peak count, not just the score"
- "Run SIRIUS for molecular formula and compound class on features with no library spectrum"
- "Assign a defensible MSI/Schymanski confidence level to each annotation given the evidence I have"
- "Tell me which annotations are safe to carry into pathway analysis and which are Level 3 hypotheses"

## Example Prompts

### Spectral Library Matching
> "Score my query MS/MS against MassBank using modified cosine with a 0.7 score and 6-peak floor."
> "Use spectral entropy similarity for identity matching and flag anything below the 0.75 natural-products threshold."

### In-silico Annotation
> "Run the SIRIUS 6 subcommand chain for formulas, fingerprints, structures against the bio database, and CANOPUS."
> "I trust formula more than structure -- report the ZODIAC formula and only call a structure confident if COSMIC FDR is set."

### Confidence Assignment
> "Given a library match but no in-house standard, what Schymanski level is this and why?"
> "Collapse my evidence set into a single confidence level and explain what would promote it."

### Avoiding Over-claiming
> "Check whether any of these annotations claim a specific isomer that MS/MS cannot resolve."
> "Before pathway analysis, flag features whose ambiguous candidate sets would inflate enrichment."

## What the Agent Will Do

1. Confirm ion families are collapsed and the adduct is assigned before inferring any neutral mass.
2. Match MS/MS against a library (matchms), enforcing both a score and a matched-peak floor.
3. Run in-silico tools (SIRIUS formula/structure/class, MetFrag) for features without library spectra.
4. Treat formula and class as more trustworthy than top-1 structure; require COSMIC FDR for confident structures.
5. Assign each annotation an MSI/Schymanski level, capping at Level 2 without an in-house standard.
6. Carry annotation uncertainty (candidate sets, levels) into downstream analysis rather than laundering it into facts.

## Tips

- A name with no confidence level is incomplete -- reject annotations that do not state one.
- High cosine on few peaks is noise; always pair the score with a matched-peak count.
- The same feature maps to many isomers SIRIUS/CSI:FingerID cannot separate -- a top-1 structure is Level 3 by default.
- Only an in-house authentic standard, same method, earns Level 1; literature/external RT and great in-silico scores do not.
- Database choice biases the answer: the truly novel metabolite is absent from your database by definition.
- Network/analogue edges mean "related to," not "is" -- propagated names are Level 3 scaffold hypotheses.

## Related Skills

- metabolomics/xcms-preprocessing - Upstream feature extraction
- metabolomics/msdial-preprocessing - Alternative feature extraction and deconvolution
- metabolomics/pathway-mapping - Downstream enrichment that must respect these confidence levels
- metabolomics/lipidomics - Lipid-specific annotation and structural resolution
- proteomics/spectral-libraries - Related spectral-matching concepts
