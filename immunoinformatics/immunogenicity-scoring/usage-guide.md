# Immunogenicity Scoring - Usage Guide

## Overview

Rank and prioritize neoantigen/epitope candidates by likely T-cell response, with the discipline this least-solved layer demands: annotate features transparently (NeoFox), rank within a patient rather than thresholding absolutely, and treat the output as an ordered, uncertainty-stated shortlist for validation, never as a verdict.

## Prerequisites

```bash
pip install neofox          # feature-annotation engine (~16 published features)
pip install pvactools       # rule-based tiering + within-tier ranking
# PRIME 2.x and MixMHCpred v3.0+ are standalone (Gfeller lab); BigMHC is a PyTorch repo
# (separate -m el and -m im heads). NeoFox also wraps NetMHCpan/NetMHCIIpan/MixMHC(2)pred.
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate these neoantigen candidates with NeoFox features and rank them for one patient"
- "Compute agretopicity defensively and flag anchor or unstable-WT artifacts"
- "Order my candidate list by presentation, abundance, and foreignness"
- "Is this a within-patient ranking or an absolute immunogenicity claim?"

## Example Prompts

### Ranking

> "Filter on expression and clonality, then rank these candidates within the patient keeping features visible"

> "Run PRIME2.0 and BigMHC-IM and report them alongside the NeoFox features, not collapsed into one score"

### Quality Features

> "Compute the fitness-model quality (amplitude x foreignness) for these neoantigens"

> "Which top candidates are high-DAI artifacts from anchor mutations or barely-presented wild-type?"

### Caveats

> "Can I say this neoantigen is more immunogenic than another patient's top candidate?"

> "These are CD4/class II candidates - how confident should the immunogenicity ranking be?"

## What the Agent Will Do

1. Apply the non-negotiable expression and clonality filters before any ranking
2. Annotate the published feature panel with NeoFox (presentation, agretopicity, foreignness, dissimilarity, PRIME, PHBR)
3. Compute agretopicity defensively, flagging anchor inflation and WT-denominator instability
4. Rank within the patient by presentation + abundance + quality, keeping features side by side
5. Use pVACtools tiers for an auditable, trap-quarantining default and hand off to human curation
6. State that scores are within-patient only and that the shortlist needs functional validation

## Tips

- **Rank, don't threshold** - scores are calibrated within tool/allele/patient; no absolute go/no-go
- **Don't stack into one number** - a confident 3-decimal composite hides fragile, correlated, biased features
- **Presentation carries the signal** - TESLA found dedicated immunogenicity scores weak; binding/stability/abundance did the work
- **DAI traps** - anchor mutations inflate it; a barely-presented WT makes the ratio noise; inspect both
- **Negative-set skepticism** - ask how a tool defined non-immunogenic before trusting its AUROC
- **CD4 is a frontier** - class II immunogenicity is even less solved; flag it as unproven
- **ImmunoBERT is presentation** - do not use it as an immunogenicity predictor

## Related Skills

- immunoinformatics/neoantigen-prediction - produces the candidate list this skill ranks
- immunoinformatics/mhc-binding-prediction - the presentation features that carry most of the signal
- immunoinformatics/mhc-class-ii-prediction - CD4 immunogenicity, the under-served frontier
- immunoinformatics/epitope-prediction - epitope candidates feeding the ranking
