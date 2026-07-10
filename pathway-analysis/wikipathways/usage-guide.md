# WikiPathways Enrichment - Usage Guide

## Overview
WikiPathways is an open, community-curated pathway database (CC0 license, ~30+ species, no formal peer review) built on a wiki model. This skill runs over-representation analysis (`enrichWP`) and GSEA (`gseWP`) against it with clusterProfiler, and uses rWikiPathways to query the database and pin a dated GMT. The point that governs reproducibility: a WikiPathways result is a snapshot of a live, monthly-updated database - the same code returns different results months apart unless a dated release is pinned. WikiPathways is a complement to KEGG/Reactome, not a sole source.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('clusterProfiler', 'rWikiPathways', 'enrichplot', 'org.Hs.eg.db'))
```

Conceptual prerequisites:
- The input is either a thresholded gene LIST (ORA) or a NAMED Entrez vector sorted DECREASING by a ranking metric (GSEA). The DE list and the ranking statistic come from differential-expression.
- The WP GMT is Entrez-keyed: SYMBOL or ENSEMBL IDs silently overlap nothing. Convert with `bitr`/an OrgDb first.
- Set a defensible `universe` (the assayed/tested genes). The default `universe=NULL` makes the background all-WP-genes, which inflates significance.
- `enrichWP`/`gseWP`/`downloadPathwayArchive` need INTERNET at run time and are NOT reproducible across monthly releases unless a dated GMT is pinned.
- WikiPathways content is CC0 with no formal peer review; treat each hit as a community claim and corroborate against KEGG/Reactome.

## Quick Start
Tell your AI agent what you want to do:
- "Run WikiPathways enrichment on my significant genes"
- "Find disease-specific pathways in my gene list that KEGG and Reactome miss"
- "Run a reproducible WikiPathways analysis pinned to a dated release"

## Example Prompts

### Basic enrichment
> "I have ~200 significant human genes as symbols and the ~13,000 tested genes as the background. Convert them to Entrez, run WikiPathways over-representation with the tested set as the universe, and give me the top 15 pathways by adjusted p-value with the gene symbols readable."

### Reproducible analysis
> "Run WikiPathways enrichment but make it reproducible: pin a dated GMT release, split the term field into the WPID and name, run enrichment on the pinned sets, and tell me which release date to report in the methods."

### ORA vs GSEA
> "I have a full ranked DESeq2 result for every gene with no clear cutoff. Should I run WikiPathways ORA or GSEA, and run whichever is appropriate."

### Organism-specific
> "Run WikiPathways enrichment for zebrafish Entrez genes, and first confirm the exact organism string WikiPathways expects."

### Combining databases
> "Run enrichment against WikiPathways, KEGG, and Reactome and tell me which pathways are unique to WikiPathways versus shared."

### Exploring pathways
> "Search WikiPathways for cancer-related pathways and show me the last-edited date for the top hit before I trust it."

## What the Agent Will Do
1. Load DE results and extract the significant gene list (ORA) or build the named decreasing ranking vector (GSEA).
2. Convert gene IDs to Entrez with `bitr` and set the background universe to the tested genes.
3. For a reproducible run, pin a dated GMT with `downloadPathwayArchive(date=, format='gmt')`, split the `name%version%wpid%org` term field, and run `enricher`/`GSEA`; otherwise run `enrichWP`/`gseWP` and log that it used the `current/` release.
4. Make the result readable with `setReadable()` and report `p.adjust`/`qvalue` (not raw p).
5. Hand the result object to enrichment-visualization for plots and corroborate hits against KEGG/Reactome.

## Understanding Results

| Column | Description |
|--------|-------------|
| ID | WikiPathways stable ID (WP####) |
| Description | Pathway name |
| GeneRatio | Query genes in the pathway / query genes mapped to any pathway |
| BgRatio | Pathway genes in the universe / universe genes mapped |
| pvalue | Raw p-value |
| p.adjust | BH-adjusted p-value (report this, not raw p) |
| qvalue | q-value |
| geneID | Genes in the pathway (symbols after setReadable) |
| Count | Number of query genes in the pathway |

For GSEA results read `NES` (sign = direction along the ranking) and `core_enrichment` (the leading-edge genes).

## WikiPathways vs Other Databases

| Feature | WikiPathways | KEGG | Reactome |
|---------|--------------|------|----------|
| License | CC0 (fully open) | Restrictive (commercial bulk/API) | CC-BY / CC0 |
| Curation | Community wiki, no formal peer review | Largely automated KO reconstruction | Expert-curated and reviewed |
| Species | ~30+ | 4000+ (genome-derived) | ~15 (deep human) |
| Focus | Disease/drug + general | Metabolic/signaling | Reaction-level mechanism |
| Reproducibility | Pin a dated monthly GMT (live `current/` otherwise) | Live REST API (date-dependent) | Local reactome.db (version-pinned) |

## Tips
- Use the exact scientific name from `listOrganisms()` / `get_wp_organisms()` for the organism argument.
- Always convert to Entrez before `enrichWP`/`gseWP`; symbols and Ensembl give an empty result with no error.
- Always set a background universe (the tested genes); the default `universe=NULL` inflates significance.
- For anything reportable, pin a dated GMT with `downloadPathwayArchive(date='YYYYMMDD', organism=, format='gmt')` and report the date. `gson_WP()` only freezes the `current/` release within a session, not across time.
- Pass `format='gmt'` - the default is `gpml`, which `read.gmt` cannot read.
- Split the GMT term field on `%` into name/version/wpid/org, or the WPIDs and names stay buried in one column.
- WikiPathways has fewer total pathways than KEGG; best used as a complement, not a standalone source.
- Check a pathway's last-edited date and curation (`getPathwayInfo`) before relying on a single WP hit for a key conclusion.
- For maximum disease/process coverage that tolerates noise, consider PFOCR (a separate, figure-OCR'd resource), not `enrichWP`.
- Use `setReadable()` to convert Entrez IDs to gene symbols in the result.
- See the enrichment-visualization skill for `dotplot()`, `cnetplot()`, and other plots.

## Related Skills
- go-enrichment - Gene Ontology enrichment
- kegg-pathways - KEGG pathway enrichment
- reactome-pathways - Reactome pathway enrichment
- gsea - Gene Set Enrichment Analysis mechanics
- enrichment-visualization - Visualization of enrichment results
- differential-expression/de-results - Source of the gene list and ranking statistic
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
