# Ortholog Inference (Database Access) Usage Guide

## Overview

Pull pre-computed ortholog calls from public databases (OrthoDB, Ensembl Compara, OMA, eggNOG, PANTHER, KEGG Orthology, HomoloGene). This skill is the database-access view: how to query each resource programmatically, what their confidence semantics mean (and why they're not comparable across resources), the orthology-conjecture caveat, when to defect to de novo computation in `comparative-genomics/ortholog-inference`, and how to handle symbol-renaming, rate limits, and stale snapshots.

## Prerequisites

```bash
pip install requests pandas
```

No API keys required for the listed resources (as of 2026). Rate limits: Ensembl REST is 15 req/sec / 55K/hour; others vary.

## Quick Start

- "Get the mouse ortholog of human BRCA1 from Ensembl Compara, with confidence score"
- "Pull orthologs of TP53 across all vertebrate species in OrthoDB"
- "Compare Compara, OMA, and OrthoDB calls for the same gene and flag disagreements"
- "Batch-fetch mouse orthologs for 500 human genes, respecting Compara's 15-req/sec limit"
- "Resolve the human/yeast ortholog of MARCHF1 -- remember the symbol was MARCH1 before 2020"

## Example Prompts

### Single-gene cross-species lookup

> "Get the high-confidence 1:1 ortholog of human BRCA1 in mouse via Ensembl Compara REST. Filter type='ortholog_one2one' and confidence=1. Return the Ensembl Mouse Gene ID plus identity percentages."

### Cross-resource consensus

> "For TP53, pull orthologs in mouse from Ensembl Compara, OMA, and OrthoDB. Convert all to a common namespace (UniProt accession) and intersect. Report which resources agree on the 1:1 call vs which include co-orthologs."

### Batch with rate-limit awareness

> "I have 500 human gene symbols. Get their zebrafish orthologs from Ensembl Compara. Use sleep=0.07 between calls to stay under 15 req/sec; handle 429 with Retry-After. Cache the results in a DataFrame."

### Symbol renaming

> "Get the mouse ortholog of MARCH1. Don't use the symbol directly -- it was renamed to MARCHF1 in 2020. First resolve via Ensembl /lookup/symbol/human/MARCHF1 to get the Ensembl Gene ID, then use the ID-based homology endpoint."

### Functional annotation transfer

> "I have an unannotated bacterial proteome. Annotate via eggNOG-mapper (not the REST API -- mapper is the right tool for batches). Cross-validate single-copy orthologs from eggNOG with KEGG Orthology assignments."

### When to give up on pre-computed

> "My species isn't in any ortholog database (just-sequenced non-model). Pre-computed resources won't have it. Switch to comparative-genomics/ortholog-inference for de novo OrthoFinder on my proteomes vs reference proteomes."

## What the Agent Will Do

1. Pick the right resource based on the (source species, target species, scale) tuple:
   - Vertebrates, 1-100 genes: Ensembl Compara
   - All species, broad coverage: OrthoDB
   - High precision, conservative: OMA
   - Functional groups: eggNOG
   - Pathway-centric: KEGG Orthology
2. Resolve gene symbols to canonical IDs before symbol-based lookups (HGNC symbols are unstable).
3. Use the right confidence field per resource; don't compare confidence across resources.
4. For batches >100, throttle to per-resource rate limits; respect Retry-After.
5. For batches >5000, use BioMart bulk export instead of REST loops.
6. When resources disagree, intersect; surface disagreement as a signal not an error.
7. Document the resource + release version for reproducibility.
8. Recommend de novo computation when species isn't in any database.

## Tips

- Ensembl REST has a 15 req/sec rate limit (55K/hour). Sleep 0.07s between calls; check `Retry-After` on 429.
- HomoloGene is frozen at 2014 -- use only for legacy comparisons; switch live workflows to Compara or OrthoDB.
- For thousands of genes, BioMart bulk export beats REST loops by 100x -- see `biomart-queries`.
- HGNC gene-symbol renames break symbol-based queries. The big ones: MARCH1->MARCHF1, SEPT1->SEPTIN1 (Excel autocorrect drove the renaming in 2020). Always resolve to Ensembl Gene IDs first.
- KEGG license is academic-free, commercial-paid. eggNOG and OrthoDB licenses are more permissive for commercial use.
- Resource disagreement is informative. If Ensembl says 1:1 but OMA gives no call, the gene tree is probably ambiguous -- flag this for review.
- The orthology conjecture (orthologs share function) is supported but weakly (Studer & Robinson-Rechavi 2009). For drug-target or clinical questions, cross-validate orthology with shared catalytic residues or expression conservation.
- eggNOG-mapper (Cantalapiedra 2021) is the right tool for batch annotation; the eggNOG REST API is for ad hoc lookup only.
- KEGG returns text TSV by default (not JSON) -- parse with `line.split('\t')`.

## Related Skills

- comparative-genomics/ortholog-inference - De novo orthology (OrthoFinder, OMA standalone, SonicParanoid)
- ensembl-rest - Broader Ensembl REST workflows
- biomart-queries - Bulk export of ortholog tables via BioMart
- uniprot-access - Namespace conversion for OMA queries
- pathway-analysis/kegg-pathways - KO to pathway mapping
