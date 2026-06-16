# Expression to Pathways - Usage Guide

## Overview

This workflow orchestrates the path from differential expression results to redundancy-collapsed functional enrichment. It does not re-teach each enrichment method; it routes the meta-decisions - ORA vs GSEA, the background universe, per-method gene-ID conversion, and the live-vs-local database caveat - and hands off to the pathway-analysis skills that own each method. The first decision is the generation (ORA or GSEA), set by whether a per-gene ranking exists for all genes or only a pre-selected list.

## Prerequisites

```r
BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'enrichplot'))
BiocManager::install(c('ReactomePA', 'rWikiPathways'))
```

Conceptual prerequisites:
- The input is either a gene LIST (for ORA) or a NAMED numeric vector sorted in DECREASING order (for GSEA). The ranking metric (signed test statistic, DESeq2 Wald `stat`, or `-sign(log2FC)*log10(p)`) determines the GSEA result.
- The background universe for ORA is the genes that were testable (entered the DE test), not the whole genome. This is the single most common published enrichment error.
- Gene IDs must match the method: OrgDb keyType for enrichGO, 'kegg'/'ncbi-geneid' for enrichKEGG, ENTREZ for ReactomePA and enrichWP. Convert with bitr/bitr_kegg.
- KEGG and WikiPathways query LIVE databases - they need internet at runtime and are not reproducible across data releases (pin the access date). GO and Reactome read local Bioconductor annotation and are reproducible given the release.
- The DE list and the ranking statistic come from differential-expression; this workflow only shapes them into enrichment inputs.

## Quick Start

Tell your AI agent what you want to do:
- "Take my DESeq2 results and find enriched pathways"
- "Should I run ORA or GSEA on these DE results, then run the right one"
- "Run GO, KEGG, and Reactome enrichment on my significant genes with the correct background"
- "Make a redundancy-collapsed dot plot and enrichment map of the results"

## Example Prompts

### Deciding ORA vs GSEA
> "I have a full DESeq2 result with the Wald statistic for every tested gene and no obvious significance cutoff. Decide whether ORA or GSEA is appropriate and run whichever fits on GO biological processes, using a named decreasing vector if it is GSEA."

### ORA with a defensible background
> "I have 240 significant genes (padj < 0.05, |log2FC| > 1) as Ensembl IDs and the full set of ~13,000 expressed genes. Run GO BP and KEGG over-representation using the expressed genes as the background, convert IDs to what each method needs, and simplify the redundant GO terms."

### End-to-end DE to pathways
> "Take my edgeR results to enriched pathways: convert IDs, run GO and Reactome ORA plus GSEA, collapse redundant terms, and give me a dot plot and an enrichment map with the top 20 terms."

### Multi-condition comparison
> "I have significant gene lists for three treatments. Compare KEGG enrichment across them in one model and make a faceted dot plot - do not compare p-values from separate runs."

## Input Requirements

| Input | Format | Required for |
|-------|--------|--------------|
| Gene list | Character vector of IDs | ORA |
| Ranked genes | Named numeric vector, sorted decreasing | GSEA |
| Background universe | Character vector of testable genes | ORA |
| Organism | OrgDb package / KEGG organism code | ID mapping and enrichment |

## What the Agent Will Do

1. Decide the generation - check whether a ranking exists for all genes (GSEA) or only a pre-selected list (ORA), and state the choice.
2. Prepare the input - build the significant list and/or the named decreasing ranked vector, and set the universe to the testable genes.
3. Convert gene IDs to the form each method needs and report the conversion rate.
4. Run the chosen enrichment - ORA (enrichGO/enrichKEGG/enrichPathway/enrichWP) and/or GSEA (gseGO/gseKEGG), flagging the live-DB steps.
5. Collapse redundancy with simplify/pairwise_termsim, then visualize deliberately (dotplot/emapplot/gseaplot2).
6. Record provenance - tool and database version/date, ranking metric, p-adjust method, and the universe - so the result is reproducible.

## ORA vs GSEA

| Feature | ORA | GSEA |
|---------|-----|------|
| Input | pre-selected gene list | all genes, ranked (named decreasing vector) |
| Cutoff | uses a DE threshold | no threshold needed |
| Information | binary (in / out) | uses magnitude and direction |
| Key parameter | the background universe | the ranking metric |
| Best for | unranked lists (modules, GWAS, screens) | full rankings; subtle coordinated shifts |

## Tips

- Generation first: if a ranking exists for all genes, prefer GSEA; reserve ORA for genuinely unranked lists. See the pathway-analysis README for why.
- Background universe: always set the universe to the genes that entered the DE test, not the genome. If the background is the whole genome, the analysis measures expression bias, not enrichment.
- Gene counts: ORA needs roughly 50-500 genes to be reliable; far fewer suggests switching to GSEA.
- GSEA ranking: use the Wald statistic (DESeq2), the moderated t (limma), or a signed p-value (edgeR), not a bare log2FC, which over-weights noisy low-count genes.
- GSEA reproducibility: set a seed before any gseGO/gseKEGG run, or permutation p-values drift between runs.
- ID conversion: deduplicate after bitr; a conversion rate below 85% usually means the wrong ID type or organism.
- Live databases: KEGG and WikiPathways query the current data release. Record the access date, and prefer local GO/Reactome when reproducibility is paramount.
- Redundancy: simplify GO terms and inspect the leading-edge or geneID core. If the same 3-5 genes explain the top 20 terms, that is one finding, not twenty.
- Multi-condition: use compareCluster, and compare NES (GSEA) rather than p-values across conditions; never compare raw p-values from separate enrichment runs.

## Related Skills

- pathway-analysis/go-enrichment - GO over-representation, background universe, redundancy reduction, length bias
- pathway-analysis/gsea - Ranked-list GSEA, named decreasing vector, ranking metric, leading edge, NES
- pathway-analysis/kegg-pathways - KEGG pathway/module enrichment, live DB, prokaryotes, multi-condition
- pathway-analysis/reactome-pathways - Reactome curated-pathway ORA and GSEA, ENTREZ IDs, reproducible local DB
- pathway-analysis/wikipathways - WikiPathways community-pathway enrichment, versioned GMT, broad species
- pathway-analysis/enrichment-visualization - Dot/bar/cnet/emap/GSEA plots and required pre-steps
- differential-expression/de-results - Source of the gene list and the ranking statistic
