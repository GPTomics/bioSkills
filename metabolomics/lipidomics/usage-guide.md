# Lipidomics - Usage Guide

## Overview

Lipidomics is a specialized branch of metabolomics focused on comprehensive lipid analysis. Lipids have unique characteristics requiring specialized annotation, normalization, and class-based analysis approaches.

## Prerequisites

```bash
# R packages
BiocManager::install("lipidr")
install.packages("ggplot2")

# Python
pip install pandas numpy scipy

# Software: MS-DIAL, LipidSearch
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze my lipidomics data by lipid class"
- "Normalize to internal standards and compare lipid profiles"

## Example Prompts

### Data Processing
> "Import my MS-DIAL lipidomics output and parse lipid annotations"
> "Extract lipid class, carbon chain length, and double bond information from annotations"

### Normalization
> "Normalize each lipid class to its matched internal standard"
> "Apply PQN normalization and calculate class-level summaries"

### Class Analysis
> "Compare lipid class abundances between treatment groups"
> "Create a lipid class composition bar plot for each sample group"
> "Test for changes in sphingolipids between disease and control"

### Chain Analysis
> "Analyze the distribution of chain lengths and saturation by class"
> "Compare fatty acid profiles between groups"

### Visualization
> "Create a lipid bubble plot showing class, saturation, and fold change"
> "Generate a heatmap of individual lipid species grouped by class"

## What the Agent Will Do

1. Parse lipid annotations (class, chains, saturation)
2. Normalize to internal standards by class
3. Aggregate by lipid class or species
4. Run statistical comparisons
5. Create lipid-specific visualizations
6. Export results with lipid metadata

## Tips

- Use class-matched internal standards (e.g., PC-d7 for PC class)
- Report sum composition when isomers cannot be separated
- Include at least one internal standard per lipid class
- Consider ion suppression varies by lipid class
- LipidMaps shorthand notation is standard (e.g., PC 34:1)

## Lipid Classes

| Class | Abbreviation | Examples |
|-------|--------------|----------|
| Glycerophospholipids | GP | PC, PE, PS, PI |
| Sphingolipids | SP | Cer, SM |
| Glycerolipids | GL | TG, DG |
| Sterol lipids | ST | CE, Cholesterol |
| Fatty acyls | FA | Free fatty acids |

## Nomenclature

Standard shorthand: `Class(carbons:double_bonds)`
- `PC(34:1)` - Phosphatidylcholine, 34 carbons, 1 double bond
- `TG(52:2)` - Triacylglycerol, 52 carbons, 2 double bonds
- `Cer(d18:1/16:0)` - Ceramide with specific chains

## Databases

- **LipidMaps** - Comprehensive lipid database
- **SwissLipids** - Curated structures
- **LipidBlast** - In-silico MS/MS library

## References

- lipidr: doi:10.1093/bioinformatics/btaa706
- LipidMaps: doi:10.1093/nar/gkl838
- MS-DIAL lipidomics: doi:10.1038/nmeth.4512
