---
name: bio-pathway-reactome
description: Tests a gene list or ranked gene vector for over-representation or coordinated shifts in Reactome's curated, peer-reviewed, reaction-level pathways using ReactomePA's enrichPathway (ORA) and gsePathway (GSEA), reading the local reactome.db so a run is reproducible given the Bioconductor release. Covers why Reactome's atomic unit is the REACTION and pathways are nested containers so a parent and child enrich on the same genes and double-count one signal, why only human is curated and every other species is orthology-inferred, why enrichPathway has NO keyType argument and returns nothing unless genes are ENTREZ (bitr first), and why viewPathway draws a LOCAL reaction network from a pathway NAME. Use when reaction-level granularity, peer-reviewed curation, or an offline-reproducible database is wanted; for comparative multi-sample or multi-omics analysis use ReactomeGSA. The DE list comes from differential-expression; plots from enrichment-visualization.
tool_type: r
primary_tool: ReactomePA
---

## Version Compatibility

Reference examples tested with: ReactomePA 1.54+, reactome.db 1.95+, clusterProfiler 4.18+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

reactome.db is a LOCAL Bioconductor annotation package pinned to the Bioconductor release, so enrichPathway/gsePathway are reproducible offline given the package version - unlike KEGG and WikiPathways, which query a live database. The reactome.org web AnalysisService tracks the current quarterly Reactome release and can therefore disagree with a local ReactomePA run (see Failure Modes).

# Reactome Pathway Enrichment

**"Which curated Reactome pathways does my gene list over-represent?"** -> Test each Reactome pathway for over-representation against a measured background, then deduplicate the hierarchy - because a Reactome result is one signal projected onto a tree of nested reactions, not a list of independent findings.
- R: `enrichPathway(gene_entrez, organism='human', universe=measured_entrez, readable=TRUE)`

Scope: ORA (enrichPathway) and GSEA (gsePathway) over Reactome reaction-rolled-to-pathway gene sets, the ENTREZ-only constraint, hierarchy deduplication, the human-curated-only species caveat, the local viewPathway reaction-network plot, and the ReactomeGSA comparative pointer. The hypergeometric test and background-universe theory -> go-enrichment. The GSEA running-sum engine and ranking metric -> gsea. The DE list / ranking statistic -> differential-expression/de-results. Dotplot/emapplot/cnetplot/gseaplot2 -> enrichment-visualization.

## The Single Most Important Modern Insight -- Reactome's Atomic Unit Is the Reaction and Pathways Are Nested Containers, So a Result Is One Signal Projected Onto a Tree, Not a List

Reactome is not a collection of pathway maps like KEGG. Its atomic unit is the ReactionlikeEvent - a single typed molecular transformation (binding, catalysis, transport, modification) with a PubMed citation - and pathways are containers that group reactions and sub-pathways into a deep event hierarchy (`TopLevelPathway -> Pathway -> sub-Pathway -> Reaction`). Two consequences define every decision in this skill:

1. **Granularity is the reason to choose Reactome AND the multiple-testing tax.** A Reactome hit is finer than a KEGG map - a specific, peer-reviewed, literature-grounded reaction - which is why it is worth using. But finer means MORE gene sets (many tiny leaf pathways plus a few huge top-level ones), so the multiple-testing burden is heavier and `minGSSize`/`maxGSSize` matter more than for GO.
2. **Nesting double-counts the signal.** A gene annotated to one leaf reaction is a member of that pathway AND every ancestor, so a parent and child enrich on the SAME genes. A live cell-cycle gene list returns "Cell Cycle Checkpoints" (parent), "G2/M Checkpoints", and "G1/S Transition" (children) stacked at the top - one signal, three "independent" small p-values. A Reactome table is therefore read by reasoning about the hierarchy (report the deepest significant node, ancestors as context), not by sorting on p.adjust.

Second load-bearing fact: **only human is curated; every non-human pathway is orthology-inferred from the human reactions, not independently curated.** Mouse projection is ~81% complete, exotic species far less, so a "Reactome mouse pathway" is a hypothesis from orthology that inherits human-curation gaps and misses mouse-specific biology. The honest output is "the deepest curated human pathway nodes my list over-represents, deduplicated against their ancestors, against a background of genes I actually measured."

## Tool Taxonomy

| Source / engine | Citation | Mechanism / role | When |
|-----------------|----------|------------------|------|
| Reactome database | Milacic 2024 *Nucleic Acids Res* 52:D672 | expert-authored, externally peer-reviewed, PubMed-cited reactions in a deep event hierarchy; CC0 | the gene-set source: reaction-level, reproducible, open |
| ReactomePA enrichPathway (ORA) | Yu & He 2016 *Mol BioSyst* 12:477 | one-sided hypergeometric test over reaction-rolled-to-pathway sets; local reactome.db; ENTREZ-only | a pre-selected gene LIST + a measured universe |
| ReactomePA gsePathway (GSEA) | Yu & He 2016 *Mol BioSyst* 12:477 | fgsea running-sum over a ranked vector; ENTREZ-only | all genes carry a statistic; distributed signal; no cutoff |
| ReactomePA viewPathway | Yu & He 2016 *Mol BioSyst* 12:477 | LOCAL ggraph reaction-network plot of ONE pathway by NAME | inspect the reactions/entities of a single hit, optionally colored by fold change |
| ReactomeGSA | Griss 2020 *Mol Cell Proteomics* 19:2115 | hosted AnalysisService client: comparative GSA / ssGSEA / PADOG, multi-omics, per-scRNA-cluster | BETWEEN-condition, multi-omics, or single-cell-cluster pathway comparison |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Pre-selected significant-gene list + a measured background | `enrichPathway(gene, universe=)` | ORA needs a list and the universe decides significance |
| All genes carry a DE statistic, cutoff would be arbitrary | `gsePathway(geneList)` -> gsea | running-sum over the full ranking; no cutoff |
| Genes are SYMBOL or ENSEMBL | `bitr(..., toType='ENTREZID')` FIRST | enrichPathway has no keyType; non-ENTREZ silently returns empty |
| Parent and child both enriched on the same genes | report the deepest significant node; ancestors = context | nesting double-counts; they are ONE finding |
| Compare pathways BETWEEN conditions / across omics / scRNA clusters | ReactomeGSA (`perform_reactome_analysis`/`analyse_sc_clusters`) | ReactomePA is single-list; the hosted service is comparative |
| Non-human within the 7 ReactomePA organisms | set `organism=`; flag results as orthology-inferred | the projection is a hypothesis, not curation |
| Species beyond the 7 (bacteria, plant, etc.) | web AnalysisService / ReactomeGSA, not ReactomePA | reactome.db maps only 7 organisms |
| Deeper metabolic coverage wanted | supplement with KEGG -> kegg-pathways | KEGG remains the deeper metabolic resource |
| The ORA-vs-GSEA decision itself, or null/benchmark theory | -> the category README | the cross-database method-selection fork lives there |
| The DE list / ranking statistic itself | -> differential-expression/de-results | upstream, not enrichment |

ReactomePA's `organism` accepts exactly seven values: human, rat, mouse, celegans, yeast, zebrafish, fly. This is a reactome.db mapping ceiling, NOT a Reactome ceiling - the database projects to ~14-20 species and the web AnalysisService covers them; do not conflate the two.

## Over-Representation Analysis (enrichPathway)

**Goal:** Find Reactome pathways over-represented in a significant-gene list, against the genes actually measured.

**Approach:** Convert the gene list to ENTREZ (mandatory - no keyType argument), pass the measured background as `universe`, run enrichPathway, then read the hierarchy rather than the raw row order.

```r
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)   # bitr

sig_entrez <- bitr(sig_symbols, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)$ENTREZID
universe   <- bitr(all_tested_symbols, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)$ENTREZID

ora <- enrichPathway(gene=sig_entrez, organism='human', universe=universe,
                     pvalueCutoff=0.05, qvalueCutoff=0.2,
                     minGSSize=10, maxGSSize=500, readable=TRUE)
# enrichResult columns: ID Description GeneRatio BgRatio RichFactor FoldEnrichment zScore
#                       pvalue p.adjust qvalue geneID Count
# FoldEnrichment and RichFactor are columns - read effect size there, do not compute GeneRatio/BgRatio by hand.
```

Without `universe`, the background is all ~11,200 Reactome-annotated ENTREZ genes (the live `BgRatio` denominator), not the genes measured, which over-states significance. `readable=TRUE` maps the `geneID` column back to symbols. Reading the result means deduplicating the hierarchy: identify the deepest significant pathway for each signal and note its ancestors as context, not as separate hits. ReactomePA has no `simplify()` equivalent (unlike GO's DAG), so this is a manual judgment call; the visual collapse (treeplot, emapplot) is owned by enrichment-visualization.

## GSEA (gsePathway)

**Goal:** Find Reactome pathways whose genes shift coordinately across the full ranking, without a significance cutoff.

**Approach:** Build a named numeric vector sorted decreasing by the ranking statistic with ENTREZ names, set a seed for permutation reproducibility, then run gsePathway and read the leading edge.

```r
gene_list <- de$stat                       # any per-gene statistic: t-stat, signed -log10 p, shrunken log2FC
names(gene_list) <- de$entrez              # names MUST be ENTREZ
gene_list <- sort(gene_list, decreasing=TRUE)

set.seed(123)                              # gsePathway permutes; fix the seed so p-values reproduce
gse <- gsePathway(geneList=gene_list, organism='human',
                  pvalueCutoff=0.05, pAdjustMethod='BH', verbose=FALSE)
# gseaResult columns: ID Description setSize enrichmentScore NES pvalue p.adjust qvalue rank leading_edge core_enrichment
```

GSEA uses the whole ranking, so the universe/background pitfall of ORA does not apply - but the hierarchy double-counting STILL does: a parent and child both score on the same leading-edge genes. The ranking metric IS the experiment (the same genes ranked differently give different leading edges); the metric choice is owned by gsea.

## viewPathway - the Reactome reaction-network plot

**Goal:** Draw the reactions and physical entities of ONE enriched pathway, optionally colored by fold change.

**Approach:** Pass the pathway NAME (the `Description`, NOT the R-HSA id) to viewPathway; it renders a LOCAL ggraph reaction network in the R graphics device. To open the actual web diagram, build the PathwayBrowser URL from the R-HSA id and browseURL it.

```r
top_name <- ora@result$Description[1]                       # the NAME, not $ID
viewPathway(top_name, organism='human', readable=TRUE, foldChange=gene_list)

# the interactive web diagram needs the R-HSA id, NOT viewPathway:
browseURL(paste0('https://reactome.org/PathwayBrowser/#/', ora@result$ID[1]))
```

`viewPathway`'s `keyType` controls the ID type of the `foldChange` names only (so a SYMBOL-named fold-change vector works with `keyType='SYMBOL'`); the ENTREZ-only constraint of enrichPathway does not extend to it. Route generic dotplot/emapplot/cnetplot/gseaplot2 to enrichment-visualization; viewPathway is owned here because it is Reactome-data-structure-specific.

## ReactomeGSA - comparative and multi-omics

**Goal:** Compare pathway activity BETWEEN conditions, across omics layers, or across single-cell clusters - which ReactomePA's single-list model cannot do.

**Approach:** ReactomeGSA is a separate Bioconductor client to Reactome's hosted AnalysisService; build a request, add datasets, and send it to the server (network required; results track the server's release, not a local reactome.db).

```r
library(ReactomeGSA)
req <- ReactomeAnalysisRequest(method='Camera')                  # or 'ssGSEA', 'PADOG'
req <- add_dataset(req, expression_values=expr_matrix, name='RNAseq',
                   type='rnaseq_counts', comparison_factor='condition',
                   comparison_group_1='A', comparison_group_2='B', sample_data=meta)
res <- perform_reactome_analysis(req)                            # sends to the server
pw  <- pathways(res)                                             # combined pathway table
sc  <- analyse_sc_clusters(seurat_obj, use_interactors=FALSE)    # per-cluster ssGSEA
```

Use ReactomePA for "is this one list over-represented / coordinately changed"; use ReactomeGSA for "which pathways DIFFER between conditions / omics / clusters".

## Per-Method Failure Modes

### SYMBOL or ENSEMBL passed to enrichPathway/gsePathway
**Trigger:** feeding a symbol or Ensembl vector because enrichGO accepted one. **Mechanism:** enrichPathway has NO keyType argument; the gene->pathway map is reactome.db's ENTREZ-keyed table, so non-ENTREZ ids match nothing. **Symptom:** zero rows on a clearly enriched list, no error. **Fix:** `bitr(..., toType='ENTREZID')` first; this is the #1 "why are my results empty" cause.

### No universe -> inflated significance
**Trigger:** calling enrichPathway without `universe=`. **Mechanism:** the background defaults to all ~11,200 Reactome-annotated genes, not the ~15,000 genes measured, shrinking every p-value. **Symptom:** implausibly significant pathways, BgRatio denominator ~11230. **Fix:** pass the measured ENTREZ set as `universe`.

### Reading the hierarchy as independent hits
**Trigger:** reporting the top-N rows by p.adjust as N findings. **Mechanism:** membership propagates up the event tree, so parent/child/sibling rows share genes. **Symptom:** "G1/S Transition", "Mitotic G1 phase and G1/S transition", and "S Phase" stacked at the top of one cell-cycle list. **Fix:** deduplicate to the deepest significant node per signal; note ancestors as context; treat BH over hierarchy rows as anti-conservative because the rows are not independent tests.

### viewPathway misuse
**Trigger:** `viewPathway('R-HSA-109582')` or expecting a browser to open. **Mechanism:** the first argument is `pathName` (the Description), and the function draws a LOCAL ggraph plot, not a browser window. **Symptom:** an error / empty plot from the id, or surprise that no browser opens. **Fix:** pass `ora@result$Description[i]`; for the web diagram `browseURL('https://reactome.org/PathwayBrowser/#/<R-HSA-id>')`.

### Trusting non-human Reactome as curated
**Trigger:** interpreting a mouse/rat/fly result as curated biology. **Mechanism:** only human is curated; all others are orthology-projected from the human reactions (mouse ~81% complete, less for distant species). **Symptom:** species-specific findings that have no human ortholog are simply absent, and projected hits inherit human-curation gaps. **Fix:** flag every non-human result as orthology-inferred and confirm species-specific hits independently.

### Assuming the R package and the web AnalysisService agree
**Trigger:** quoting "Reactome says pathway X, p=..." without naming the tool. **Mechanism:** ReactomePA pins to the installed reactome.db snapshot while the web tool tracks the current quarterly release, and their default backgrounds and identifier-projection differ. **Symptom:** a collaborator's reactome.org p-values differ from the local ones on the same list. **Fix:** state the tool, the release/reactome.db version, and the background; do not treat the two as interchangeable.

### Assuming arbitrary-organism support
**Trigger:** passing a bacterial or plant `organism`. **Mechanism:** reactome.db maps only 7 organisms for ReactomePA. **Symptom:** an unsupported-organism error. **Fix:** for species beyond the 7, use the web AnalysisService or ReactomeGSA.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `pvalueCutoff=0.05` | enrichPathway/gsePathway default | filters on p.adjust (BH) by default in the result; standard FDR gate |
| `qvalueCutoff=0.2` | enrichPathway default | secondary q-value gate on the ORA result |
| `pAdjustMethod='BH'` | enrichPathway default | Benjamini-Hochberg FDR; less conservative than Bonferroni, but anti-conservative across nested hierarchy rows |
| `minGSSize=10` | enrichPathway default | drop tiny leaf pathways (2-3 genes) that inflate false positives; matters MORE for Reactome's deep tree |
| `maxGSSize=500` | enrichPathway default | drop huge top-level pathways (e.g. "Signal Transduction") that always enrich and are uninformative |
| Reactome background ~11,200 | live BgRatio denominator | the ENTREZ genes with any Reactome annotation; the implicit universe if `universe=` is omitted |
| Mouse projection ~81% complete | Reactome inference docs | the fraction of human reactions projected to mouse by orthology; far lower for distant species |
| `set.seed(123)` before gsePathway | reproducibility | gsePathway permutes; without a fixed seed the permutation p-values drift between runs |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| enrichPathway returns 0 rows on a clear list | genes are SYMBOL/ENSEMBL, not ENTREZ | `bitr(..., toType='ENTREZID')` first (no keyType arg) |
| Implausibly significant pathways | no `universe=`, background is all ~11k Reactome genes | pass the measured ENTREZ set as `universe` |
| Top hits are parent/child of one pathway | hierarchy nesting double-counts the signal | report the deepest significant node; ancestors as context |
| `viewPathway('R-HSA-...')` errors or is empty | first arg is the NAME (Description), not the id | `viewPathway(ora@result$Description[i], ...)` |
| viewPathway did not open a browser | it draws a LOCAL ggraph plot | `browseURL('https://reactome.org/PathwayBrowser/#/<id>')` for the web diagram |
| Different p-values than reactome.org | release skew + different universe between local db and web service | name the tool, reactome.db version, and background |
| gsePathway results change each run | no `set.seed` before the permutation | set a fixed seed |
| Unsupported-organism error | organism outside the 7 reactome.db maps | use the web AnalysisService / ReactomeGSA |

## References

- Milacic M, Beavers D, Conley P, et al. 2024. The Reactome Pathway Knowledgebase 2024. *Nucleic Acids Res* 52:D672-D678.
- Yu G, He QY. 2016. ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. *Mol BioSyst* 12:477-479.
- Griss J, Viteri G, Sidiropoulos K, Nguyen V, Fabregat A, Hermjakob H. 2020. ReactomeGSA - Efficient Multi-Omics Comparative Pathway Analysis. *Mol Cell Proteomics* 19:2115-2125.

## Related Skills

- go-enrichment - The hypergeometric test and the background-universe problem
- gsea - The GSEA running-sum engine and ranking-metric choice
- kegg-pathways - KEGG pathway/module enrichment; deeper metabolic coverage
- wikipathways - WikiPathways community-pathway enrichment (also CC0)
- enrichment-visualization - Dot/bar/cnet/emap/tree/GSEA plots of enrichment results
- differential-expression/de-results - Source of the gene list and the ranking statistic
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
