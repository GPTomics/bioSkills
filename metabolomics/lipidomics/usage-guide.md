# Lipidomics Usage Guide

## Overview

Lipidomics is a specialized branch of metabolomics focused on the comprehensive analysis of lipids. Lipids have unique characteristics requiring specialized annotation and analysis approaches.

## When to Use This Skill

- Lipid-focused studies (membranes, signaling, energy storage)
- Shotgun lipidomics (direct infusion MS)
- LC-MS lipidomics (chromatographic separation)
- Lipid biomarker discovery
- Membrane composition analysis

## Lipid Classes

| Class | Abbreviation | Examples |
|-------|--------------|----------|
| Fatty Acyls | FA | Palmitic acid (FA 16:0) |
| Glycerolipids | GL | Triacylglycerol (TG) |
| Glycerophospholipids | GP | PC, PE, PS, PI, PG |
| Sphingolipids | SP | Ceramide (Cer), Sphingomyelin (SM) |
| Sterol lipids | ST | Cholesterol, CE |
| Prenol lipids | PR | Coenzyme Q |

## Lipid Nomenclature

Standard shorthand notation: `Class(carbon:double_bonds)`

Examples:
- `PC(34:1)` - Phosphatidylcholine with 34 carbons and 1 double bond
- `TG(52:2)` - Triacylglycerol with 52 carbons and 2 double bonds
- `Cer(d18:1/16:0)` - Ceramide with specific chains

## Analysis Workflow

1. **Data acquisition** - LC-MS or shotgun MS
2. **Peak detection** - XCMS, MS-DIAL, LipidSearch
3. **Lipid annotation** - LipidMaps, in-house database
4. **Normalization** - Internal standards, PQN
5. **Statistical analysis** - Univariate/multivariate
6. **Pathway interpretation** - KEGG lipid pathways

## Normalization Strategies

### Internal Standards
- Use class-matched internal standards (e.g., PC-d7 for PC)
- One IS per lipid class minimum
- Calculate ratios to IS

### Sample Normalization
- Total ion current (TIC)
- Probabilistic quotient normalization (PQN)
- Median centering

## Common Issues

### Isomer separation
- LC may not separate all isomers
- Consider reporting sum composition

### Ion suppression
- Matrix effects vary by lipid class
- Use multiple internal standards

### Annotation confidence
- Level 1: MS/MS match + RT
- Level 2: MS/MS match only
- Level 3: Mass match only

## Databases

- **LipidMaps** - Comprehensive lipid database
- **SwissLipids** - Curated lipid structures
- **LipidBlast** - In-silico MS/MS library

## References

- lipidr: doi:10.1093/bioinformatics/btaa706
- LipidMaps: doi:10.1093/nar/gkl838
- MS-DIAL lipidomics: doi:10.1038/nmeth.4512
