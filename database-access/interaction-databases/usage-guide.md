# Interaction Databases Usage Guide

## Overview

Query protein-protein and gene interaction databases (STRING, BioGRID, IntAct, SIGNOR, Reactome, HuRI, HuMAP, OmniPath). Encodes the decision matrix (which resource for which question), STRING v12 channel semantics and confidence tiers, SIGNOR as the only major signed/directed signaling resource, OmniPath as the modern meta-database, BioGRID's HT-vs-LT distinction, per-resource license constraints (commercial restrictions on ConsensusPathDB and PhosphoSitePlus), and the STRING version-pinning trap (v11.5 deprecated 2023).

## Prerequisites

```bash
pip install requests pandas networkx
# Optional:
pip install omnipath              # OmniPath Python client (cleaner than raw REST)
# R alternatives: STRINGdb, OmnipathR (Bioconductor)
```

BioGRID requires a free API key from `https://webservice.thebiogrid.org/`. All other listed resources are key-free.

## Quick Start

- "Get STRING v12 interactions for [TP53, MDM2, BRCA1] at high confidence (700+); show per-channel scores"
- "Filter BioGRID interactions to physical, low-throughput only (HT screens have higher FP rates)"
- "Get signed signaling interactions from SIGNOR for MAPK pathway -- with direction and mechanism"
- "Aggregate STRING + OmniPath + BioGRID into one network; track per-edge provenance"
- "Use OmniPath with license=commercial to restrict to permissive-license sources"

## Example Prompts

### Decision: which resource

> "I want only physically interacting proteins for a 10-gene set, publication-grade. Don't use STRING's combined score (which mixes textmining). Use IntAct or BioGRID filtered to LT physical systems (Affinity Capture-MS, Two-hybrid, Reconstituted Complex)."

### STRING confidence + channels

> "Get STRING v12 interactions for my gene list at score 700. Return the full per-channel breakdown (escore, dscore, tscore, ascore, etc.) so I can see whether each edge is supported by experiments vs textmining vs coexpression."

### Signed signaling

> "I need signed, directed signaling for a TP53 pathway analysis. STRING and BioGRID won't do -- they're undirected. Use SIGNOR, which has direction, effect (up/down-regulates), and mechanism (phosphorylation, etc.)."

### Multi-resource aggregation

> "Build a union network across STRING (high confidence), OmniPath, and BioGRID LT physical. Track sources per edge; flag edges supported by multiple resources as highest confidence."

### License-aware OmniPath

> "I'm building a commercial product. Query OmniPath with license=commercial to filter to commercially-permissive sources only. Audit each source's license before adding to the pipeline."

### Symbol disambiguation

> "Resolve gene symbols to UniProt accessions before querying. MARCH1 was renamed to MARCHF1 in 2020 by HGNC. Some interaction resources updated; some haven't."

## What the Agent Will Do

1. Pick the right resource based on the question (physical/functional, signed/unsigned, HT/LT, species).
2. For STRING, set `caller_identity` and pin to v12 URL; pick confidence threshold from the use case (700+ for publication).
3. For BioGRID, filter `THROUGHPUT='Low Throughput'` and physical experimental systems for high-quality calls.
4. For SIGNOR, preserve direction and mechanism; use DiGraph not Graph.
5. For OmniPath, preserve `sources` and `references` per edge for provenance.
6. Aggregate across resources for consensus; flag multi-source edges as higher confidence.
7. Document license constraints for any commercial pipeline (avoid ConsensusPathDB / PhosphoSitePlus without legal review).
8. Resolve gene symbols via UniProt or HGNC IDs before querying (HGNC rename instability).

## Tips

- STRING URL is version-pinned: `version-12-0.string-db.org/api` as of 2024. Older `version-11-5` URLs were deprecated 2023.
- STRING `caller_identity` is requested for usage attribution and rate-limit allocation; non-compliant clients get throttled first.
- STRING combines 7 channels into a single score by default. Treating combined score as "physical interaction" is the most common misuse. Filter to `escore` (experiments channel) for physical-only.
- BioGRID's `THROUGHPUT` flag is the most important quality filter -- LT (low-throughput) is curated and high-confidence; HT (high-throughput screens) has higher FP rates.
- SIGNOR is **the only major curated database with signed and mechanism-typed interactions**. For any signaling pathway model, SIGNOR is essential.
- OmniPath (Türei 2021) aggregates 100+ sources with provenance -- the modern one-stop for "give me everything available for X".
- Reactome is gold-standard for human pathway curation; species coverage outside human is limited.
- HuRI (binary Y2H interactome) and HuMAP (AP-MS complexes) are reference human-specific resources.
- For commercial use: stick to STRING + BioGRID + IntAct + SIGNOR + Reactome + HuRI/HuMAP + OmniPath; ConsensusPathDB is academic-only; PhosphoSitePlus requires paid commercial license.
- For directional resources (SIGNOR, OmniPath), use `nx.DiGraph` not `nx.Graph` to preserve direction.

## Related Skills

- uniprot-access - Resolve symbols to UniProt accessions
- ensembl-rest - Cross-reference Ensembl IDs in network nodes
- gene-regulatory-networks/coexpression-networks - Co-expression as complement to PPI
- pathway-analysis/go-enrichment - Functional enrichment of network genes
- pathway-analysis/reactome-pathways - Pathways alongside Reactome interactions
- data-visualization/network-visualization - Network rendering
