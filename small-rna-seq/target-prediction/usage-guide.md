# miRNA Target Prediction - Usage Guide

## Overview

Predict and prioritize miRNA target genes using seed-based tools (miRanda, TargetScan, miRDB) and experimentally validated databases (miRTarBase, multiMiR). The governing reality is that a predicted target is a hypothesis, not a finding: seed prediction is roughly 50% false-positive even for conserved sites, a 6mer seed occurs by chance about once per 4 kb, and a single miRNA represses most real targets only modestly. Confidence comes from climbing an evidence ladder - conservation, validated CLIP/reporter data, and above all intersecting predictions with inversely-correlated mRNA differential expression from the same samples - not from stacking more seed-based predictors, which is pseudo-replication.

## Prerequisites

```bash
conda install -c bioconda miranda
pip install pandas biopython gseapy
# multiMiR (unifies predicted + validated) is an R/Bioconductor package
```

## Quick Start

Tell your AI agent:
- "Predict targets for my DE miRNAs, then keep only the functional ones"
- "Rank TargetScan targets by weighted context++ score"
- "Get experimentally validated targets from miRTarBase, weighted by evidence type"
- "Intersect predicted targets with my anti-correlated mRNA-seq"
- "Why shouldn't I just run GO enrichment on the predicted target list?"

## Example Prompts

### Prediction and Ranking

> "Run miRanda with strict seed pairing for hsa-miR-21-5p against my 3' UTRs"

> "Rank this miRNA's TargetScan targets by weighted context++ score"

> "Get miRDB targets with score >= 80"

### Evidence and Validation

> "Pull validated targets from miRTarBase and separate strong from less-strong evidence"

> "Use multiMiR to combine predicted and validated targets in one query"

> "Which of my predicted targets are supported by AGO-CLIP?"

### The Confidence Move

> "Intersect predicted targets of my UP miRNAs with DOWN mRNAs from matched RNA-seq"

> "My target enrichment says 'cancer pathways' for every miRNA - how do I fix the circularity?"

## What the Agent Will Do

1. Generate candidate targets by seed complementarity and thermodynamics (miRanda) or database lookup (TargetScan/miRDB)
2. Rank by the appropriate score (weighted context++, mirSVR, miRDB score)
3. Anchor in experimentally validated interactions (miRTarBase strong, CLIP/CLASH)
4. Intersect with inversely-correlated mRNA/protein DE from matched samples - the key confidence filter
5. Only then run functional enrichment, on an evidence-filtered list, to avoid circular inflation

## Tips

- A predicted target is a hypothesis for the wet lab, not a result - seed prediction is ~50% false-positive even for conserved sites
- Five seed-based tools agreeing is pseudo-replication; require an orthogonal evidence tier (CLIP, validated, expression) instead
- The single most useful move is intersecting predictions with anti-correlated mRNA DE from the same samples
- Rank TargetScan by `weighted context++ score` (more negative = stronger); it is relative within a miRNA, not an absolute cutoff
- The miRNA column in the TargetScan context-scores file is `Mirbase ID`, not `miRNA family`
- miRTarBase separates strong (reporter/western/qPCR) from less-strong (CLIP/NGS) evidence - weight accordingly
- Enrichment on an unfiltered predicted target list is circular and inflates big pathways; filter by evidence first
- miRNAs also repress translationally, so some real targets do not drop at the mRNA level (ribosome profiling/proteomics catches those)
- Seed-only prediction also misses the bulged/seedless/3'-compensatory sites (a large fraction of real CLASH interactions are noncanonical, Helwak 2013), so a clean seed list is incomplete, not just imprecise
- AGO-CLIP evidence is cell-type/state-specific; a peak from another tissue is weak support for yours
- Treat ceRNA/sponge mechanisms skeptically - stoichiometry is usually too low to derepress targets

## Related Skills

- differential-mirna - Source of DE miRNAs
- pathway-analysis/go-enrichment - Enrich an evidence-filtered target list
- database-access/entrez-fetch - Fetch UTR/gene sequences
- clip-seq/ago-clip-mirna-targets - Direct AGO-CLIP/CLASH target evidence
