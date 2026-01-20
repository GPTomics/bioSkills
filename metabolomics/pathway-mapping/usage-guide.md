# Pathway Mapping Usage Guide

## Overview

Pathway mapping places differential metabolites in biological context, identifying affected metabolic processes.

## Analysis Types

### Over-Representation Analysis (ORA)
- Input: List of significant metabolites
- Tests: Hypergeometric/Fisher's exact
- Question: Are any pathways enriched?

### Quantitative Enrichment (QEA)
- Input: Metabolites with values (FC, concentration)
- Tests: Globaltest, GSEA-like
- Question: Are pathways overall affected?

### Topology-Based
- Considers pathway structure
- Central metabolites weighted more
- More specific results

## Databases

| Database | Content | ID Types |
|----------|---------|----------|
| KEGG | Metabolic pathways | C-numbers |
| Reactome | All pathways | ChEBI |
| SMPDB | Small molecule | HMDB |
| BioCyc | Multi-organism | BioCyc IDs |

## ID Conversion

Metabolites need standard IDs:
1. HMDB ID (Human Metabolome Database)
2. KEGG compound ID (C-number)
3. ChEBI ID
4. PubChem CID

Use MetaboAnalyst or UniChem for conversion.

## Interpretation

- FDR < 0.05: Statistically significant
- Impact > 0.1: Biologically relevant (topology)
- Hits/Total: Coverage of pathway

## Best Practices

1. Use multiple databases
2. Consider pathway size
3. Report both ORA and topology
4. Validate with biology

## References

- MetaboAnalyst: doi:10.1093/nar/gkz240
- KEGG: doi:10.1093/nar/gkaa970
- SMPDB: doi:10.1093/nar/gkab1086
