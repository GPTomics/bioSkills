# Pathway Analysis

Functional enrichment and pathway analysis using R/Bioconductor (clusterProfiler, ReactomePA, rWikiPathways) for Gene Ontology, KEGG, Reactome, WikiPathways, and Gene Set Enrichment Analysis.

## Overview

This category covers functional annotation of gene lists: over-representation analysis (ORA) tests whether specific biological terms are enriched in a gene set, while Gene Set Enrichment Analysis (GSEA) uses ranked gene lists to find coordinated changes. Both approaches help interpret experimental results in biological context. Multiple pathway databases are supported for comprehensive analysis.

**Tool type:** `r`
**Primary tools:** clusterProfiler, ReactomePA, rWikiPathways, enrichplot, org.Hs.eg.db (and other OrgDb packages)

## Skills

| Skill | Description |
|-------|-------------|
| [go-enrichment](go-enrichment/) | Gene Ontology over-representation analysis with enrichGO |
| [kegg-pathways](kegg-pathways/) | KEGG pathway enrichment with enrichKEGG and enrichMKEGG |
| [reactome-pathways](reactome-pathways/) | Reactome pathway enrichment with ReactomePA |
| [wikipathways](wikipathways/) | WikiPathways enrichment with enrichWP and rWikiPathways |
| [gsea](gsea/) | Gene Set Enrichment Analysis with gseGO, gseKEGG |
| [enrichment-visualization](enrichment-visualization/) | Dot plots, bar plots, enrichment maps, cnetplots, GSEA plots |

## Workflow

```
Differential Expression Results
    |
    v
[Gene List] ---------> Significant genes (e.g., padj < 0.05, |log2FC| > 1)
    |
    v
[go-enrichment] -----> GO terms (BP, MF, CC)
    |
    v
[kegg-pathways] -----> KEGG pathways and modules
    |
    v
[reactome-pathways] -> Reactome signaling pathways
    |
    v
[wikipathways] ------> Community-curated pathways
    |
    v
[enrichment-visualization] -> Interpret and present results
```

**Alternative: GSEA workflow (uses all genes ranked by statistic)**

```
Differential Expression Results
    |
    v
[Ranked Gene List] --> All genes ranked by log2FC or -log10(p)*sign(FC)
    |
    v
[gsea] --------------> gseGO, gseKEGG
    |
    v
[enrichment-visualization] -> GSEA plots, ridge plots
```

## Analysis Types

| Type | Input | Function | Use Case |
|------|-------|----------|----------|
| Over-Representation | Gene list | enrichGO, enrichKEGG, enrichPathway, enrichWP | Significant gene set |
| GSEA | Ranked genes | gseGO, gseKEGG, gsePathway, gseWP | All genes with statistics |
| Module Enrichment | Gene list | enrichMKEGG | KEGG modules |

## Pathway Databases

| Database | Package | Organisms | License | Focus |
|----------|---------|-----------|---------|-------|
| GO | clusterProfiler | All with OrgDb | Open | Function/process |
| KEGG | clusterProfiler | 4000+ | Commercial | Metabolic/signaling |
| Reactome | ReactomePA | 7 | Open | Signaling detail |
| WikiPathways | rWikiPathways | 30+ | CC0 | Community/disease |

## Organism Support

clusterProfiler supports 4000+ species via KEGG and requires OrgDb packages for GO analysis.

**Common OrgDb packages:**

| Organism | Package |
|----------|---------|
| Human | org.Hs.eg.db |
| Mouse | org.Mm.eg.db |
| Rat | org.Rn.eg.db |
| Zebrafish | org.Dr.eg.db |
| Fruit fly | org.Dm.eg.db |
| C. elegans | org.Ce.eg.db |
| Yeast | org.Sc.sgd.db |
| Arabidopsis | org.At.tair.db |

## Example Prompts

### GO Enrichment
- "Run GO enrichment on my differentially expressed genes"
- "Find enriched biological processes for these genes"
- "What molecular functions are over-represented in my gene list?"

### KEGG Pathways
- "Find enriched KEGG pathways for my gene set"
- "What pathways are active in my differentially expressed genes?"
- "Run KEGG module enrichment analysis"

### Reactome Pathways
- "Run Reactome pathway enrichment on my genes"
- "Find enriched Reactome pathways for my DEGs"
- "What signaling pathways are enriched in Reactome?"

### WikiPathways
- "Run WikiPathways enrichment analysis"
- "Find community-curated pathways for my gene list"
- "Search WikiPathways for cancer-related pathways"

### GSEA
- "Run GSEA on my ranked gene list"
- "Perform gene set enrichment analysis using GO terms"
- "Run GSEA with KEGG pathways"

### Visualization
- "Create a dot plot of my enrichment results"
- "Make an enrichment map showing term relationships"
- "Show a gene-concept network for top pathways"
- "Create a GSEA running score plot"

## Requirements

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'enrichplot', 'org.Hs.eg.db'))

# For Reactome and WikiPathways:
BiocManager::install(c('ReactomePA', 'rWikiPathways'))

# For additional organisms:
BiocManager::install('org.Mm.eg.db')  # Mouse
```

## Key Functions

| Function | Package | Purpose |
|----------|---------|---------|
| enrichGO | clusterProfiler | GO over-representation |
| enrichKEGG | clusterProfiler | KEGG pathway enrichment |
| enrichMKEGG | clusterProfiler | KEGG module enrichment |
| gseGO | clusterProfiler | GO GSEA |
| gseKEGG | clusterProfiler | KEGG GSEA |
| setReadable | clusterProfiler | Convert IDs to symbols |
| bitr | clusterProfiler | ID conversion |
| dotplot | enrichplot | Dot plot visualization |
| barplot | enrichplot | Bar plot visualization |
| cnetplot | enrichplot | Gene-concept network |
| emapplot | enrichplot | Enrichment map |
| gseaplot2 | enrichplot | GSEA running score plot |
| ridgeplot | enrichplot | Ridge plot for GSEA |
| enrichPathway | ReactomePA | Reactome ORA |
| gsePathway | ReactomePA | Reactome GSEA |
| viewPathway | ReactomePA | Open pathway in browser |
| enrichWP | clusterProfiler | WikiPathways ORA |
| gseWP | clusterProfiler | WikiPathways GSEA |
| listPathways | rWikiPathways | List available pathways |

## Notes

- **enrichKEGG has no readable parameter** - use `setReadable()` separately
- **GSEA requires ranked gene list** - named vector sorted by statistic
- **OrgDb needed for GO** - KEGG uses online database
- **p-value adjustment** - BH (Benjamini-Hochberg) is default
- **keyType flexibility** - enrichGO accepts ENSEMBL, SYMBOL, ENTREZID, etc.

## Related Skills

- **differential-expression** - Generate gene lists and statistics for enrichment
- **single-cell** - Marker genes can be analyzed with pathway enrichment
- **database-access** - Fetch gene annotations from NCBI

## References

- [clusterProfiler Book](https://yulab-smu.top/biomedical-knowledge-mining-book/)
- [clusterProfiler Bioconductor](https://bioconductor.org/packages/clusterProfiler/)
- [enrichplot Bioconductor](https://bioconductor.org/packages/enrichplot/)
- [ReactomePA Bioconductor](https://bioconductor.org/packages/ReactomePA/)
- [rWikiPathways Bioconductor](https://bioconductor.org/packages/rWikiPathways/)
- [Reactome Website](https://reactome.org/)
- [WikiPathways Website](https://www.wikipathways.org/)
