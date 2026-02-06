---
name: bio-gene-regulatory-networks-differential-networks
description: Compare gene regulatory and co-expression networks between biological conditions to identify rewired regulatory relationships using DiffCorr. Detects gained, lost, and reversed gene-gene correlations between conditions. Use when comparing co-expression networks between disease vs control, treatment conditions, or developmental stages.
tool_type: r
primary_tool: DiffCorr
---

# Differential Networks

Compare co-expression and regulatory networks between biological conditions to identify rewired gene-gene relationships.

## DiffCorr Workflow

DiffCorr uses Fisher's z-transform to test whether the correlation between two genes differs significantly between two conditions.

### Required Libraries

```r
library(DiffCorr)
library(igraph)
```

### Input Preparation

```r
# Separate expression matrices by condition (genes as columns, samples as rows)
expr_all <- read.csv('normalized_counts.csv', row.names = 1)
sample_info <- read.csv('sample_info.csv', row.names = 1)

expr_condition1 <- t(expr_all[, sample_info$condition == 'control'])
expr_condition2 <- t(expr_all[, sample_info$condition == 'disease'])

# Filter to variable genes (top 2000-5000 most variable)
gene_vars <- apply(expr_all, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:3000]
expr_condition1 <- expr_condition1[, top_genes]
expr_condition2 <- expr_condition2[, top_genes]
```

### Run DiffCorr Analysis

```r
# Compute correlation matrices
cor1 <- cor(expr_condition1, method = 'pearson')
cor2 <- cor(expr_condition2, method = 'pearson')

n1 <- nrow(expr_condition1)
n2 <- nrow(expr_condition2)

# Fisher z-transform test for differential correlation
# Compares correlation coefficients between two groups
result <- comp.2.cc.fdr(
    output.file = 'diffcorr_results.txt',
    data1 = expr_condition1,
    data2 = expr_condition2,
    threshold = 0.05  # FDR threshold
)
```

### Parse and Classify Results

```r
diffcorr <- read.delim('diffcorr_results.txt')

# Classify differential correlations
classify_edge <- function(cor1, cor2, threshold = 0.3) {
    if (abs(cor1) < threshold & abs(cor2) >= threshold) return('gained')
    if (abs(cor1) >= threshold & abs(cor2) < threshold) return('lost')
    if (cor1 > threshold & cor2 < -threshold) return('reversed')
    if (cor1 < -threshold & cor2 > threshold) return('reversed')
    return('unchanged')
}

diffcorr$edge_type <- mapply(classify_edge, diffcorr$cor1, diffcorr$cor2)
table(diffcorr$edge_type)

# Focus on significant rewired edges
significant <- diffcorr[diffcorr$p.adj < 0.05, ]
rewired <- significant[significant$edge_type != 'unchanged', ]
print(paste('Significant rewired edges:', nrow(rewired)))
```

### Identify Rewired Hub Genes

```r
# Genes involved in most differential correlations
gene_rewiring <- c(rewired$gene1, rewired$gene2)
rewiring_counts <- sort(table(gene_rewiring), decreasing = TRUE)
top_rewired_genes <- head(rewiring_counts, 20)
print(top_rewired_genes)

# These genes have the most condition-specific regulatory changes
```

## DGCA (Alternative)

DGCA (Differential Gene Correlation Analysis) was archived from CRAN in May 2024. It is still available via GitHub.

```r
# Install from GitHub (archived from CRAN 2024-05)
devtools::install_github('andymckenzie/DGCA')
library(DGCA)

# Run differential correlation analysis
dgca_results <- ddcorAll(
    inputMat = expr_all,
    design = sample_info$condition,
    compare = c('control', 'disease'),
    adjust = 'BH',         # multiple testing correction
    nPerm = 0,             # 0 for analytical p-values (faster)
    corrType = 'pearson'
)

# Filter significant results
sig_dgca <- dgca_results[dgca_results$pValDiff_adj < 0.05, ]

# Classify by correlation change category
# DGCA provides: +/+, +/-, -/+, -/-, +/0, 0/+, -/0, 0/-, 0/0
table(sig_dgca$Classes)
```

## Python NetworkX Approach

```python
import pandas as pd
import numpy as np
from scipy import stats
import networkx as nx

def build_correlation_network(expr_matrix, threshold=0.6):
    '''Build network from correlation matrix, keeping edges above threshold.'''
    cor_matrix = expr_matrix.corr(method='pearson')
    genes = cor_matrix.columns.tolist()
    G = nx.Graph()
    G.add_nodes_from(genes)

    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            r = cor_matrix.iloc[i, j]
            if abs(r) >= threshold:
                G.add_edge(genes[i], genes[j], weight=r)
    return G


def fisher_z_test(r1, n1, r2, n2):
    '''Fisher z-transform test for difference between two correlations.'''
    z1 = np.arctanh(np.clip(r1, -0.9999, 0.9999))
    z2 = np.arctanh(np.clip(r2, -0.9999, 0.9999))
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z_stat = (z1 - z2) / se
    p_value = 2 * stats.norm.sf(abs(z_stat))
    return z_stat, p_value


def differential_network_analysis(expr_cond1, expr_cond2, fdr_threshold=0.05):
    '''Compare correlation networks between two conditions.'''
    genes = expr_cond1.columns.tolist()
    n1, n2 = len(expr_cond1), len(expr_cond2)
    cor1 = expr_cond1.corr()
    cor2 = expr_cond2.corr()

    results = []
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            r1, r2 = cor1.iloc[i, j], cor2.iloc[i, j]
            z_stat, p_val = fisher_z_test(r1, n1, r2, n2)
            results.append({
                'gene1': genes[i], 'gene2': genes[j],
                'cor_cond1': r1, 'cor_cond2': r2,
                'z_stat': z_stat, 'p_value': p_val
            })

    df = pd.DataFrame(results)

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    _, df['p_adj'], _, _ = multipletests(df['p_value'], method='fdr_bh')

    # Classify edges
    def classify(row, threshold=0.3):
        r1, r2 = row['cor_cond1'], row['cor_cond2']
        if abs(r1) < threshold and abs(r2) >= threshold:
            return 'gained'
        if abs(r1) >= threshold and abs(r2) < threshold:
            return 'lost'
        if r1 > threshold and r2 < -threshold:
            return 'reversed'
        if r1 < -threshold and r2 > threshold:
            return 'reversed'
        return 'unchanged'

    df['edge_type'] = df.apply(classify, axis=1)
    return df


# Usage
expr_all = pd.read_csv('normalized_counts.csv', index_col=0).T
sample_info = pd.read_csv('sample_info.csv', index_col=0)

expr_ctrl = expr_all.loc[sample_info[sample_info['condition'] == 'control'].index]
expr_disease = expr_all.loc[sample_info[sample_info['condition'] == 'disease'].index]

# Filter to top variable genes
gene_vars = expr_all.var()
top_genes = gene_vars.nlargest(3000).index.tolist()
expr_ctrl, expr_disease = expr_ctrl[top_genes], expr_disease[top_genes]

diff_results = differential_network_analysis(expr_ctrl, expr_disease)
significant = diff_results[diff_results['p_adj'] < 0.05]
print(significant['edge_type'].value_counts())
```

### Visualize Differential Network

```python
import matplotlib.pyplot as plt

rewired = significant[significant['edge_type'] != 'unchanged']

G_diff = nx.Graph()
color_map = {'gained': '#2ca02c', 'lost': '#d62728', 'reversed': '#9467bd'}

for _, row in rewired.iterrows():
    G_diff.add_edge(row['gene1'], row['gene2'],
                    color=color_map[row['edge_type']], weight=abs(row['z_stat']))

edge_colors = [G_diff[u][v]['color'] for u, v in G_diff.edges()]

fig, ax = plt.subplots(figsize=(12, 12))
pos = nx.spring_layout(G_diff, k=2, seed=42)
nx.draw_networkx(G_diff, pos, edge_color=edge_colors, node_size=50,
                 font_size=6, width=0.5, alpha=0.7, ax=ax)

# Legend
import matplotlib.patches as mpatches
legend_handles = [mpatches.Patch(color=c, label=l) for l, c in color_map.items()]
ax.legend(handles=legend_handles, loc='upper left')
ax.set_title('Differential co-expression network')
plt.savefig('differential_network.pdf', bbox_inches='tight')
```

## Statistical Considerations

| Consideration | Guideline | Rationale |
|---------------|-----------|-----------|
| Minimum samples per group | 20 recommended, 15 floor | Correlation estimates unstable below 15 |
| Multiple testing | BH FDR correction | Thousands of gene pairs tested |
| Correlation method | Pearson (default) | Spearman for non-linear; Pearson is standard |
| Gene filtering | Top 2000-5000 variable | Reduces test burden, focuses on informative genes |
| Effect size filter | abs(delta_r) > 0.3 | Avoid reporting trivially different correlations |

## Related Skills

- coexpression-networks - Build WGCNA networks for individual conditions
- scenic-regulons - TF regulon inference from scRNA-seq
- differential-expression/deseq2-basics - DE analysis to complement network rewiring
