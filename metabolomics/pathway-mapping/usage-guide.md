# Pathway Mapping - Usage Guide

## Overview

Pathway mapping places differential metabolites in biological context, identifying affected metabolic processes using enrichment analysis and pathway topology.

## Prerequisites

```bash
# R packages
install.packages("MetaboAnalystR")
BiocManager::install("clusterProfiler")

# Python
pip install requests pandas  # For API access

# Web tool: https://www.metaboanalyst.ca/
```

## Quick Start

Tell your AI agent what you want to do:
- "Map my significant metabolites to KEGG pathways"
- "Run pathway enrichment analysis on differential metabolites"

## Example Prompts

### Over-Representation Analysis
> "Test which KEGG pathways are enriched in my list of significant metabolites"
> "Run hypergeometric test for pathway enrichment using HMDB IDs"

### Quantitative Enrichment
> "Perform GSEA-style enrichment using metabolite fold changes"
> "Run topology-based pathway analysis weighting by network centrality"

### ID Conversion
> "Convert my metabolite names to KEGG compound IDs for pathway analysis"
> "Map HMDB IDs to KEGG pathway terms"

### Visualization
> "Create a pathway bubble plot showing enrichment and impact scores"
> "Highlight my differential metabolites on the KEGG glycolysis pathway"

## What the Agent Will Do

1. Convert metabolite IDs to database format (KEGG, HMDB)
2. Map metabolites to pathway databases
3. Calculate enrichment statistics
4. Compute pathway impact scores (topology)
5. Generate pathway visualizations
6. Export enriched pathways with statistics

## Tips

- Convert metabolite IDs before analysis (use MetaboAnalyst or UniChem)
- Use multiple databases (KEGG, Reactome, SMPDB) for comprehensive coverage
- Report both enrichment p-value and pathway impact score
- Consider pathway size when interpreting results
- Validate computational findings with biological knowledge

## Analysis Types

| Method | Input | Question |
|--------|-------|----------|
| ORA | Metabolite list | Are pathways over-represented? |
| QEA | Metabolites + values | Are pathways affected overall? |
| Topology | Metabolites + network | Which central metabolites affected? |

## Key Databases

| Database | Content | ID Types |
|----------|---------|----------|
| KEGG | Metabolic pathways | C-numbers |
| Reactome | All pathways | ChEBI |
| SMPDB | Small molecule | HMDB |
| BioCyc | Multi-organism | BioCyc IDs |

## Interpretation

| Metric | Threshold | Meaning |
|--------|-----------|---------|
| FDR | < 0.05 | Statistically significant |
| Impact | > 0.1 | Biologically relevant |
| Hits/Total | Higher = better | Pathway coverage |

## References

- MetaboAnalyst: doi:10.1093/nar/gkz240
- KEGG: doi:10.1093/nar/gkaa970
- SMPDB: doi:10.1093/nar/gkab1086
