# GO Enrichment - Usage Guide

## Overview
Gene Ontology over-representation analysis (ORA) tests whether a gene LIST is enriched for genes annotated to each GO term, using the one-sided hypergeometric (Fisher 2x2) test in clusterProfiler's enrichGO. The single decision that governs the result is the BACKGROUND universe: ORA is a competitive test comparing the foreground against the rest of the measured genome, so the universe IS the null distribution. Pick the genes that could have entered the list (the tested-gene set), not the whole genome, or every p-value is wrong in a direction that cannot be seen.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))
BiocManager::install('goseq')   # only for RNA-seq gene-length bias correction
```

Conceptual prerequisites:
- The input is a pre-selected gene LIST plus a matched BACKGROUND universe. A ranked vector of all genes with no cutoff is a GSEA input instead (see the gsea skill).
- The universe must be the genes that could have been selected (RNA-seq: genes that passed filtering / have a non-NA DESeq2 pvalue; proteomics: detected proteins; targeted panel: panel genes), mapped with the same ID call as the foreground. Omitting the universe defaults it to all annotated genes and inflates significance.
- enrichGO's source default is ont='MF', not 'BP' - set ont explicitly every time.
- pvalueCutoff filters the adjusted p (p.adjust), not the raw p; an empty result usually means the cutoff or the universe, not biology.
- GO annotation is local (org.*.eg.db + GO.db, pinned to the Bioconductor release), so a GO ORA is reproducible given the package versions - no internet needed at runtime.
- The DE list and the significance threshold that defines the foreground come from differential-expression; this skill consumes them.

## Quick Start
Tell your AI agent what you want to do:
- "Run GO enrichment on my differentially expressed genes with the tested genes as background"
- "Which biological processes are over-represented in this gene list?"
- "Simplify the redundant GO terms and show fold enrichment, not just p-values"
- "My data is RNA-seq - correct for gene length bias with GOseq"

## Example Prompts

### GO over-representation
> "I have 240 significant genes from a DESeq2 contrast (padj < 0.05, |log2FC| > 1) as gene symbols, and the ~13,000 expressed genes as the background. Run GO BP over-representation with the expressed genes as the universe, simplify the redundant terms, and give me the top 15 by adjusted p-value with fold enrichment."

> "Run GO enrichment for all three ontologies (BP, MF, CC) on this gene list and simplify each ontology separately."

### ORA vs GSEA choice
> "I have a full ranked DESeq2 result with the test statistic for every gene and no clear significance cutoff. Should I run ORA or GSEA, and run whichever is appropriate on GO terms."

### Direction and bias
> "Run GO enrichment separately for the upregulated and downregulated genes and report both."

> "This is RNA-seq and I am worried about long-gene bias - run length-corrected GO enrichment with GOseq."

### Custom gene sets
> "Run an over-representation test of my gene list against these MSigDB hallmark sets with my expressed genes as background."

## What the Agent Will Do
1. Load the DE results, extract the foreground (the hits) and the universe (the genes actually tested), and map both ID sets with the same bitr call, reporting the conversion rate.
2. Run enrichGO with ont set explicitly, the universe passed, and BH adjustment, returning an enrichResult with GeneRatio, BgRatio, p.adjust, and Count.
3. Reduce GO-DAG redundancy with simplify() per ontology, or use topGO if in-test decorrelation is wanted.
4. For RNA-seq, optionally correct gene-length bias with GOseq (named 0/1 vector, Wallenius PWF, then BH).
5. Report terms ranked with fold enrichment alongside p.adjust, treating the list as hypothesis generation.

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | GO term ID (GO:XXXXXXX) |
| Description | GO term name |
| GeneRatio | k/n = foreground genes in term / foreground genes annotated to any term |
| BgRatio | M/N = universe genes in term / universe genes annotated to any term |
| pvalue | Raw hypergeometric p-value |
| p.adjust | BH-adjusted p-value (this is what pvalueCutoff filters) |
| qvalue | Q-value |
| geneID | Genes in the term |
| Count | k = number of foreground genes in the term |

Fold enrichment = GeneRatio / BgRatio. Both denominators are restricted to ANNOTATED genes.

## Three Ontologies

| Ontology | Code | Description |
|----------|------|-------------|
| Biological Process | BP | what the genes participate in |
| Molecular Function | MF | biochemical activity (the enrichGO default) |
| Cellular Component | CC | where in the cell |

## Tips
- Always pass the universe (the tested genes, not the genome); this is the deepest ORA error. The genome is defensible only when every gene truly could have been detected.
- Read fold enrichment, not just p-values: a 2000-gene term beats a 12-gene term on p at a fraction of the effect size.
- simplify() works on one ontology at a time (semantic similarity is defined within a single DAG); run BP, MF, CC separately and simplify each. It does not de-redundify an ont='ALL' object.
- enrichGO's default ont is 'MF'; set ont explicitly to avoid silently testing the wrong ontology.
- pvalueCutoff filters the adjusted p; if no terms appear, set pvalueCutoff=1 and qvalueCutoff=1 to inspect everything before loosening real thresholds.
- After bitr(), deduplicate one-to-many maps and report the conversion rate; flag results when more than ~15% of genes are lost.
- Run separate enrichment on up- and down-regulated genes; a mixed-direction list can cancel and hide a real term.
- For RNA-seq, consider GOseq (Wallenius PWF); TMM/RPKM normalization does not remove the length/selection bias.
- The clusterProfiler `enrichment_force_universe` option keeps unannotated genes in the universe instead of intersecting with annotated genes; the intersection is usually the desired behavior for GO (an unannotated gene can never be a hit), so use the option only to match another tool's denominator.
- Treat enrichment as hypothesis generation, not validation; terms derived from a DE list cannot validate that same DE list.
- See the enrichment-visualization skill for dotplot, cnetplot, emapplot, and treeplot.

## Related Skills

- gsea - Ranked-list GSEA alternative when a full ranking exists
- kegg-pathways - KEGG pathway and module enrichment
- reactome-pathways - Reactome curated-pathway enrichment
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - Dot/bar/cnet/emap/tree plots of enrichment results
- differential-expression/de-results - Source of the gene list and the tested-gene universe
- database-access/entrez-fetch - Fetch gene annotations / ID maps from NCBI
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
