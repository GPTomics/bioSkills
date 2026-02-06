# Interaction Databases - Usage Guide

## Overview

This skill enables AI agents to query protein-protein interaction (PPI) databases including STRING, BioGRID, IntAct, and OmniPath. Retrieves interaction partners, confidence scores, and functional enrichment for gene sets, then converts results into NetworkX graphs for downstream analysis.

## Prerequisites

```bash
pip install requests pandas networkx
```

Optional for specific databases:
```bash
# OmniPath Python client (alternative to REST)
pip install omnipath

# STRINGdb R package (alternative to REST)
# install.packages('BiocManager')
# BiocManager::install('STRINGdb')
```

- BioGRID requires a free API key from [thebiogrid.org](https://wiki.thebiogrid.org/doku.php/biogridrest)
- STRING, IntAct, and OmniPath are key-free

## Quick Start

Tell your AI agent what you want to do:

- "Get all protein interactions for TP53 from STRING"
- "Build a PPI network for my list of DE genes"
- "Find high-confidence interactions between these kinases"
- "Query BioGRID for physical interactions with BRCA1"
- "Combine STRING and OmniPath interactions into one network"

## Example Prompts

### STRING Queries

> "Get STRING interactions for TP53, MDM2, BRCA1, ATM, and CHEK2 with high confidence"

> "Download the STRING network image for my DNA damage response genes"

> "Run STRING enrichment analysis on my upregulated gene list"

### BioGRID Queries

> "Find all low-throughput physical interactions for MYC in BioGRID"

> "Get BioGRID interactions for EGFR filtered to co-immunoprecipitation experiments"

### Multi-Database

> "Query both STRING and OmniPath for my gene list and merge the results"

> "Build a combined PPI network from STRING, BioGRID, and IntAct for these 50 genes"

### Network Construction

> "Convert my STRING interactions to a NetworkX graph and find the hub genes"

> "Build a PPI network and compute centrality measures for my gene list"

## What the Agent Will Do

1. Resolve gene identifiers to database-specific IDs
2. Query one or more interaction databases via REST APIs
3. Filter interactions by confidence score and evidence type
4. Convert results to a NetworkX graph
5. Compute network statistics (degree, clustering, components)
6. Export the network for visualization or downstream analysis

## Tips

- **Score threshold 700** - Use high-confidence STRING scores (700+) for publication networks; use 400 for exploratory analysis where recall matters
- **BioGRID API key** - Free but required; register at thebiogrid.org for access
- **Low-throughput evidence** - Filter BioGRID to low-throughput experiments for higher-quality physical interactions
- **OmniPath for signaling** - OmniPath curates directed signaling interactions, ideal for pathway reconstruction
- **Multi-database consensus** - Interactions found in multiple databases are more reliable; aggregate and filter by source count
- **Species codes** - STRING uses NCBI taxonomy IDs: 9606 (human), 10090 (mouse), 7955 (zebrafish), 7227 (fly)
- **Rate limiting** - Add delays between batch queries to STRING and BioGRID; process gene lists in chunks of 50-100

## Related Skills

- database-access/uniprot-access - Protein annotations and cross-references
- pathway-analysis/go-enrichment - Functional enrichment of network genes
- gene-regulatory-networks/coexpression-networks - Co-expression network construction
- data-visualization/network-visualization - Visualize interaction networks
