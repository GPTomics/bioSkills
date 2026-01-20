# Diversity Analysis

## Overview

Diversity analysis characterizes microbial community structure through alpha (within-sample) and beta (between-sample) diversity metrics.

## Alpha Diversity Metrics

| Metric | Measures | Interpretation |
|--------|----------|----------------|
| Observed | Richness | Number of taxa |
| Chao1 | Richness | Estimated total richness |
| Shannon | Richness + Evenness | Higher = more diverse |
| Simpson | Dominance | Higher = more even |
| Faith's PD | Phylogenetic diversity | Requires tree |

## Beta Diversity Metrics

| Metric | Type | When to Use |
|--------|------|-------------|
| Bray-Curtis | Abundance-based | General purpose |
| Jaccard | Presence/absence | Focus on composition |
| UniFrac (weighted) | Phylogenetic | Abundant taxa matter |
| UniFrac (unweighted) | Phylogenetic | Rare taxa matter |
| Aitchison | Compositional | Compositionally-aware |

## Workflow

1. **Rarefy** (optional) - Normalize sequencing depth
2. **Alpha diversity** - Calculate per-sample metrics
3. **Test differences** - Kruskal-Wallis, Wilcoxon
4. **Beta diversity** - Calculate distance matrix
5. **Ordination** - PCoA, NMDS visualization
6. **PERMANOVA** - Test group differences

## Rarefaction Debate

**For rarefaction:**
- Controls for sequencing depth
- Traditional approach
- Required for some metrics

**Against rarefaction:**
- Discards data
- Modern methods handle variable depth
- Use variance-stabilizing transforms instead

## Statistical Tests

### Alpha Diversity
- **Kruskal-Wallis**: Non-parametric, >2 groups
- **Wilcoxon**: Non-parametric, 2 groups
- **ANOVA**: Parametric (check normality)

### Beta Diversity
- **PERMANOVA** (adonis2): Group differences
- **ANOSIM**: Alternative to PERMANOVA
- **betadisper**: Test dispersion homogeneity

## Interpretation Guidelines

### PERMANOVA R2
- 0.1-0.2: Small effect
- 0.2-0.4: Moderate effect
- >0.4: Large effect

### NMDS Stress
- <0.1: Excellent
- 0.1-0.2: Good
- >0.2: Poor (try more dimensions)
