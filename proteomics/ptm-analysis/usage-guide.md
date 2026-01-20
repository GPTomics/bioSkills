# Post-Translational Modification Analysis

## Overview

PTM analysis identifies and quantifies chemical modifications on proteins that regulate function, localization, and interactions.

## Common PTMs

| Modification | Sites | Mass Shift | Enrichment |
|--------------|-------|------------|------------|
| Phosphorylation | S, T, Y | +79.97 | TiO2, IMAC |
| Acetylation | K, N-term | +42.01 | Anti-acetyl antibody |
| Methylation | K, R | +14.02 | Anti-methyl antibody |
| Ubiquitination | K | +114.04 (GG) | Anti-K-GG antibody |
| Glycosylation | N, S, T | Variable | Lectin enrichment |

## Workflow

1. **Enrichment**: Concentrate modified peptides
2. **LC-MS/MS**: Acquire fragmentation spectra
3. **Database search**: Include variable modifications
4. **Site localization**: Score confidence in site assignment
5. **Quantification**: Site-level abundance
6. **Normalization**: Often adjust for protein-level changes

## Site Localization

Localization probability indicates confidence in exact modification site:
- **>0.75**: Class I (confident)
- **0.50-0.75**: Class II (probable)
- **<0.50**: Class III (ambiguous)

Tools: MaxQuant (localization prob), Ascore, PhosphoRS

## Motif Analysis

Look for enriched sequence patterns around modification sites:
- **Kinase motifs**: [RK]xx[ST] for basophilic kinases
- **Proline-directed**: [ST]P for CDKs, MAPKs
- **Tools**: motif-x, pLogo, icelogo

## Normalization Strategies

1. **Simple**: Median centering of PTM sites
2. **Protein-adjusted**: Normalize PTM to protein abundance
3. **MSstatsPTM**: Statistical framework for PTM vs protein

## Key Databases

- **PhosphoSitePlus**: Curated PTM database
- **UniProt**: Annotated modifications
- **Phospho.ELM**: Kinase-substrate relationships
