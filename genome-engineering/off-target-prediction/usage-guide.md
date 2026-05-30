# Off-Target Prediction - Usage Guide

## Overview

Nominates and assesses CRISPR off-target sites genome-wide. The skill enumerates candidate sites by mismatch and bulge tolerance (Cas-OFFinder), ranks them with the published CFD score, runs variant-aware screening against population/individual genomes (CRISPRme), frames the empirical genome-wide discovery assays (GUIDE-seq, CIRCLE-seq, CHANGE-seq, DISCOVER-seq), and helps choose a high-fidelity nuclease. Its central discipline is the evidence ladder: a site is *predicted* (in-silico), *detected* (an unbiased assay shows the nuclease acts there), or *validated* (confirmed editing at a measured allele frequency) -- and an in-silico list is a hypothesis, never a safety verdict.

## Prerequisites

```bash
pip install pandas
# Cas-OFFinder: standalone binary from https://github.com/snugel/cas-offinder (native
#   DNA/RNA bulge support from v3.0.0); needs a directory of reference-genome FASTA.
# CFD tables: mismatch_score.pkl + pam_scores.pkl ship with CRISPOR / the Doench 2016 code
#   -- load them, do not hand-type CFD values.
# Variant-aware nomination: CRISPRme (web or CLI) against gnomAD/1000G + an individual genome.
# Empirical assays (GUIDE-seq, CIRCLE-seq, etc.) are wet-lab methods; this skill plans/interprets them.
```

## Quick Start

Tell your AI agent what you want to do:
- "Check this guide for off-targets genome-wide, including bulges"
- "Rank my three candidate guides by specificity"
- "I'm designing a therapeutic guide -- how should I screen off-targets across ancestries?"
- "An assay found an off-target my in-silico search missed -- why?"
- "Which high-fidelity Cas9 should I use for RNP delivery?"

## Example Prompts

### In-silico nomination and ranking

> "I have guide GGCCGACCTGTCGCTGACGC. Build a Cas-OFFinder search against hg38 with up to 4 mismatches and 2 nt DNA/RNA bulges, then rank the hits by CFD and flag any with a low-mismatch hit in a coding gene."

> "Compare the aggregate specificity of these three guides and tell me which to take forward."

### Therapeutic-grade assessment

> "This is a clinical-stage guide. Lay out a variant-aware off-target screen against gnomAD plus the patient genome, the empirical assays I should run, and how to validate."

> "We'll deliver as RNP -- which high-fidelity nuclease keeps its specificity, and what's the catch?"

### Interpretation

> "GUIDE-seq in K562 and CIRCLE-seq give different site lists -- which do I believe?"

> "Our amplicon panel was clean. Is the guide safe?"

## What the Agent Will Do

1. Ask what's at stake (research knockout vs therapeutic) -- this sets the depth.
2. Enumerate candidate sites with Cas-OFFinder, including bulges and a relaxed PAM (NRG).
3. Rank sites with the published CFD tables and compute an aggregate specificity to compare candidate guides (not as a safety verdict).
4. For therapeutic guides, run variant-aware nomination (CRISPRme vs gnomAD + individual).
5. Recommend a high-fidelity nuclease appropriate to the delivery format, with the on-target caveat.
6. Plan empirical discovery (>=2 orthogonal assays) and amplicon validation, plus a structural readout for large deletions/translocations.
7. Report each claim at the right rung (predicted/detected/validated) with a stated LoD.

## Tips

- An in-silico mismatch list is a hypothesis generator, not a verdict -- the failure mode (bulges, chromatin, sub-LoD editing) is silent, so a clean report is not proof of safety.
- CFD is a SpCas9/NGG relative ranker, not a probability -- compare guides, don't read an absolute CFD as safe/unsafe; load the published tables rather than hand-typing them.
- The empirical assays disagree by design: in-vitro (CIRCLE/CHANGE-seq) over-call, cell-based (GUIDE/DISCOVER-seq) under-call rare sites and are cell-type-specific -- triangulate.
- "Switch to a high-fidelity nuclease" is often the biggest lever, but plasmid-tuned variants can lose their edge as RNP (use HiFi Cas9/Sniper-Cas9), and the on-target tax is guide-dependent.
- For human therapeutics, reference-only screening is a liability -- screen variant-aware across ancestries (a common SNP can create a real off-target).
- "Not detected" means "below the limit of detection," never "zero" -- always state the LoD.
- Amplicon panels miss large deletions (PCR dropout) and translocations -- add a junction-capture readout.

## Related Skills

- grna-design - Design and on-target-score guides before the specificity check
- base-editing-design - Owns the deaminase (Cas-independent) DNA/RNA off-target classes
- prime-editing-design - pegRNA off-target and PE3 nicking-guide specificity
- crispr-screens/crispresso-editing - Quantify and validate editing at candidate sites
- variant-calling/variant-annotation - Annotate whether off-targets hit genes/pathogenic loci
- genome-intervals/bed-file-basics - Intersect off-target sites with genomic features
- database-access/ncbi-datasets-cli - Download the reference genome for the search
