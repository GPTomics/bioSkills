---
name: bio-pathway-enrichment-visualization
description: Turns an enrichResult or gseaResult from clusterProfiler/enrichplot into a figure that collapses or shows gene-set redundancy, using dotplot, barplot, cnetplot, emapplot, treeplot, ridgeplot, gseaplot2, and upsetplot. Covers why a default top-20 GO dotplot is one biological theme drawn twenty times (the DAG/nesting guarantees redundant overlapping terms), so the figure is a modeling choice between SHOWING redundancy (pairwise_termsim -> emapplot/treeplot) and DELETING it (simplify/REVIGO); why cnetplot/emapplot/treeplot need pairwise_termsim first; why enrichplot ships no barplot for gseaResult (a bar cannot carry a signed NES); why GeneRatio is not fold enrichment; and why showCategory silently truncates. Use when plotting ORA or GSEA results, collapsing redundant GO terms visually, encoding a dotplot, or building a publication enrichment figure. Statistics come from go-enrichment and gsea; generic ggplot -> data-visualization/ggplot2-fundamentals.
tool_type: r
primary_tool: enrichplot
---

## Version Compatibility

Reference examples tested with: enrichplot 1.30+, clusterProfiler 4.18+, ggplot2 3.5+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The single biggest hazard here is the cnetplot/emapplot/goplot API churn. enrichplot 1.25.5 (2024-10) moved these to the ggtangle backend and DROPPED several old arguments (`cex_label_category`, `cex_label_gene`, `circular`, `colorEdge`, `group`/`group_category`/`group_legend`). The examples target the post-churn API (1.30+), but installed bases span three argument generations - run `?cnetplot` / `?emapplot` and adapt rather than pinning one generation.

# Enrichment Visualization

**"Make a figure from my enrichment results"** -> Render an enrichResult or gseaResult with enrichplot, choosing how the gene-set REDUNDANCY is handled - because a raw top-N plot is one biological theme drawn N times, not N findings.
- R: `dotplot(ego, showCategory=20)`; redundancy as structure via `emapplot(pairwise_termsim(ego))`

Scope: turn an `enrichResult`/`gseaResult`/`compareClusterResult` into a figure, and decide whether to SHOW or DELETE redundancy. The ORA/GSEA statistics that produce the objects -> go-enrichment, gsea. The method-selection fork lives in the category README; this skill already explains the redundancy (DAG/nesting) it renders. `simplify()` existence (the GO-DAG dedup) -> go-enrichment. Generic ggplot2 grammar (scales, themes, faceting) -> data-visualization/ggplot2-fundamentals. Cytoscape UI mechanics -> data-visualization/network-visualization.

## The Single Most Important Modern Insight -- An Enrichment Figure Is a Modeling Choice, Not a Rendering of a Table

The default `dotplot(ego, showCategory=20)` is almost never twenty findings. The GO DAG and nested pathway databases guarantee that a real signal surfaces as a CLUSTER of near-identical overlapping terms driven by the same handful of genes: if "mitotic cell cycle" is enriched, then "cell cycle process," "cell cycle," "cell division," and a dozen ancestors and siblings enrich too. Sorting by p-value floats that redundant cluster to the top, so the figure shows ONE theme twenty times and crowds out the second and third themes entirely. The reader infers twenty independent findings; the figure lies by omission.

So the load-bearing question is never "which plotting function" but three decisions:

1. **How is the redundancy collapsed?** SHOW it as structure (`pairwise_termsim` -> `emapplot`/`treeplot`, or EnrichmentMap) so the cluster size conveys support, or DELETE it (`simplify()` for a shorter GO list, REVIGO for a flat-list treemap). Plotting raw top-20 with no collapse step is the error the whole skill exists to prevent.
2. **Is the direction kept?** GSEA results are SIGNED (NES > 0 = activated, NES < 0 = suppressed). Any GSEA figure that maps magnitude to a bar height or |NES|, or colors by a one-sided p-value ramp, silently merges activation and suppression. enrichplot deliberately ships NO barplot method for `gseaResult` for exactly this reason - a bar from zero cannot carry a sign.
3. **Does the caption admit truncation?** `showCategory=20` is a window, not a census. If 200 terms passed FDR it is a 10% sample chosen by whatever `orderBy` used. Report the total significant count and the similarity `method=`/`min_edge=` settings - two honest analysts get different emapplot modules from the same object.

## The Object Model -- What Gets Plotted

Every enrichplot function dispatches on the S4 class of its input, and the SAME function name encodes different things by class:

- `enrichResult` (ORA: enrichGO/enrichKEGG/enricher) - columns `ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count`.
- `gseaResult` (GSEA: gseGO/gseKEGG/GSEA) - columns `ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge, core_enrichment`, plus the `@geneList` slot (the ranked named vector that drove the analysis).
- `compareClusterResult` (compareCluster) - stacked results across gene lists; the substrate for faceted dotplots.

Encoding definitions (verified): **GeneRatio = k/n** (k = query genes annotated to the term, n = query genes mapped to any term; stored as the string `"k/n"`). **Count = k** (the numerator alone). **BgRatio = M/N** (M = universe genes annotated to the term, N = universe genes mapped). **Fold enrichment = (k/n)/(M/N)** = `GeneRatio/BgRatio`. GeneRatio is NOT effect size: a giant term (M=800) can post a large GeneRatio with trivial enrichment, while a small term (M=5, k=3) shows a modest GeneRatio but huge fold enrichment. The p-value, not GeneRatio, is the test statistic. dotplot can put `GeneRatio` OR `Count` on x; size = Count, color = p.adjust by default.

## Tool Taxonomy

| Plot / method | Encodes | Class | Redundancy handling | Direction-aware |
|---------------|---------|-------|---------------------|-----------------|
| dotplot | GeneRatio (x), Count (size), p.adjust (color) | ORA + GSEA | none (raw top-N) | only if x/color = NES |
| barplot | Count or GeneRatio (height), p.adjust (color) | ORA ONLY | none | no - misuse for GSEA |
| cnetplot | gene<->term bipartite net; item color = fold change | ORA + GSEA | shows shared genes (gene side) | yes (item color) |
| emapplot | term net; edge = gene overlap; clusters = redundant groups | ORA + GSEA | SHOWS redundancy (term side) | node color = p.adjust |
| treeplot | hierarchical Ward clusters of terms | ORA + GSEA | COLLAPSES into nCluster groups | node color = p.adjust |
| ridgeplot | leading-edge metric density per set | GSEA ONLY | per-set | YES (left/right shift) |
| gseaplot2 | running ES + hit ticks + ranked metric | GSEA ONLY | single / few sets | YES (peak sign) |
| upsetplot | gene-overlap combinations (ORA); per-set metric boxplots (GSEA) | ORA + GSEA | quantifies overlap | metric boxplots for GSEA |
| goplot | induced GO DAG subgraph | GO ONLY | exposes DAG nesting | no |
| heatplot | gene x term matrix, color by fold change | ORA + GSEA | flattened cnetplot | yes (fold change) |
| simplify() | semantic dedup of GO terms | GO ORA/GSEA | DELETES redundant terms (lives in go-enrichment) | n/a |
| REVIGO / EnrichmentMap | non-redundant subset / node-edge map | any list | DELETE / SHOW + annotate | EnrichmentMap by sign |

Citations: dotplot/barplot/cnet/tree/upset/goplot/heatplot are enrichplot, paper-of-record clusterProfiler 4.0 (Wu 2021 *The Innovation* 2:100141). emapplot reimplements EnrichmentMap (Merico 2010 *PLoS One* 5:e13984; protocol Reimand 2019 *Nat Protoc* 14:482). simplify/Wang use GOSemSim (Yu 2010 *Bioinformatics* 26:976). REVIGO (Supek 2011 *PLoS One* 6:e21800). ridgeplot/gseaplot2 display the ES Subramanian 2005 *PNAS* 102:15545 defined.

## Decision Tree by Intent

| Intent | Do this | Why / avoid |
|--------|---------|-------------|
| First look at ORA results | `dotplot(simplify(ego))` - collapse GO redundancy THEN dotplot | avoid raw `dotplot(ego, showCategory=20)` (redundant cluster floats up) |
| Many significant terms, show the structure | `pairwise_termsim()` -> `emapplot` (topology) or `treeplot` (named clusters) | the redundancy becomes the message, not hidden |
| Hundreds of sets, manuscript figure | EnrichmentMap (Cytoscape; Reimand 2019 protocol) -> data-visualization/network-visualization | a top-20 list is indefensible at that scale |
| Flat GO-ID + p-value list from a non-clusterProfiler tool | REVIGO (treemap / MDS) | external semantic collapse |
| GSEA overview, all sets | `ridgeplot(gse)` | direction + shape preserved; never a barplot of NES |
| GSEA, one pathway in detail | `gseaplot2(gse, geneSetID=1)` | the running ES; a single number hides the shape |
| Compare several pathways' running scores | `gseaplot2(gse, geneSetID=1:3)` | overlay in one panel |
| Which genes bridge multiple terms | `cnetplot` (<=5-8 terms) or `heatplot` | a 20-term cnetplot is a hairball |
| Need effect size, not GeneRatio | `dotplot(ego, x='FoldEnrichment')` or compute `GeneRatio/BgRatio` | a dot far right on GeneRatio is not strong over-representation |
| Compare conditions / gene lists | `dotplot(ck) + facet_grid(~Cluster)` on compareCluster | one model, faceted panels |
| term similarity for KEGG/Reactome/custom | `pairwise_termsim(x, method='JC')` (default) | Wang/Resnik need the GO DAG |
| term similarity for GO, want DAG-awareness | `pairwise_termsim(x, method='Wang', semData=godata(...))` | JC sees only gene overlap |
| The ORA/GSEA statistics themselves | -> go-enrichment, gsea | upstream, not visualization |

## Dotplot -- the Three-Channel Summary

`dotplot(object, x='geneRatio', color='p.adjust', showCategory=10, orderBy='x', label_format=30)`. The terms are ordered by `orderBy='x'` (the x variable), NOT by p-value, so by default the TOP dot is the highest GeneRatio, not the most significant. State the ordering or set it.

```r
dotplot(ego, showCategory = 20)                                       # x = GeneRatio, size = Count, color = p.adjust
dotplot(ego, x = 'FoldEnrichment', showCategory = 20)                 # effect size = (k/n)/(M/N), not GeneRatio
dotplot(gse, x = 'NES', showCategory = 20, color = 'p.adjust')        # signed GSEA summary (dotplot dispatches on gseaResult)
dotplot(gse, showCategory = 20, split = '.sign') + facet_grid(~.sign) # split GSEA up vs down
```

For a compareClusterResult, `dotplot.compareClusterResult` defaults `showCategory=5` per cluster and `includeAll=TRUE` (a term top-N in any cluster appears in every column).

## Barplot -- ORA Only

`barplot(height, x='Count', color='p.adjust', showCategory=8)`. There is NO `barplot` method for `gseaResult` (verified) - forcing a bar onto GSEA drops the NES sign. For signed GSEA use a NES dotplot, ridgeplot, or gseaplot2.

```r
barplot(ego, showCategory = 15)                          # height = Count, color = p.adjust
barplot(ego, x = 'GeneRatio', showCategory = 15)
```

## Show the Redundancy -- pairwise_termsim then emapplot/treeplot

**Goal:** Reveal that a block of near-identical enriched terms is one biological theme, by clustering terms on gene-set overlap and drawing the clusters.

**Approach:** Populate the term-similarity matrix first with `pairwise_termsim` (emapplot/treeplot READ `x@termsim` and do NOT compute it), then draw it as a force-directed map (emapplot, shows topology) or a deterministic Ward tree (treeplot, named clusters). The similarity `method=` and `min_edge=` are modeling choices that change the picture.

```r
ego_ts <- pairwise_termsim(ego)                          # JC (Jaccard on gene overlap), the default; any gene-set type
emapplot(ego_ts, showCategory = 30)                      # nodes = terms, edges = overlap >= min_edge (0.2), clusters = redundant groups
treeplot(ego_ts, showCategory = 30, nCluster = 5)        # deterministic Ward clustering into 5 labeled groups

# GO terms, DAG-aware similarity (Wang sees parent/child closeness even with modest gene overlap)
ego_ts <- pairwise_termsim(ego, method = 'Wang', semData = GOSemSim::godata('org.Hs.eg.db', ont = 'BP'))
```

`pairwise_termsim` `method` is exactly one of `{Resnik, Lin, Rel, Jiang, Wang, JC}`, default `JC`. Resnik/Lin/Rel/Jiang/Wang are GO-ONLY and need a `GOSemSimDATA` object; JC works for any gene-set type. Lower `min_edge` and everything connects to everything (the "if every node touches every node, the result IS redundant" diagnostic); raise it and only the strongest overlaps survive.

## Gene-Concept Network (cnetplot)

**Goal:** Show which genes are shared across enriched terms - the redundancy seen from the gene side - and their direction.

**Approach:** Draw a bipartite term-to-gene network, mapping the gene-node color to fold change. Keep to 5-8 terms or it collapses into a hairball. The ggtangle-era arguments differ from older tutorials - introspect before pinning args.

```r
cnetplot(ego, showCategory = 5)                          # ggtangle backend (enrichplot >= 1.25.5)
cnetplot(ego, showCategory = 5, foldChange = gene_list)  # gene color by fold change; node_label = 'all'|'category'|'item'|'none'
# OLDER installed versions used: cnetplot(ego, foldChange=fc, circular=TRUE, colorEdge=TRUE) -- those args were REMOVED; run ?cnetplot
```

## GSEA Plots -- ridgeplot, gseaplot2

`ridgeplot(gse, showCategory=30, fill='p.adjust', core_enrichment=TRUE, orderBy='NES')` draws, per set, a density of the `@geneList` metric values of its LEADING-EDGE genes. Shifted right = up-ranked, left = down-ranked, bimodal = the set straddles both extremes (often too broad). `ridgeplot` needs the `ggridges` package (an enrichplot Suggests-only dependency) or it errors. `gseaplot2(gse, geneSetID, subplots=1:3)` stacks the running ES curve, the hit ticks, and the ranked-metric profile; `geneSetID` is required and accepts an index, a vector (`1:3` to overlay), or an ID string.

```r
ridgeplot(gse, showCategory = 20)                        # direction + shape; the honest GSEA overview
gseaplot2(gse, geneSetID = 1:3, pvalue_table = TRUE)     # overlay three sets' running scores
```

## Specialized Views -- upsetplot, goplot, heatplot

```r
upsetplot(ego, n = 10)                                   # gene-overlap combinations across terms (gseaResult gives per-set metric boxplots)
goplot(ego)                                              # GO-ONLY: the induced DAG subgraph; needs the ggarchery package (enrichplot Suggests)
heatplot(ego, foldChange = gene_list, showCategory = 15) # gene x term matrix, color by direction; a flattened cnetplot
```

## All Outputs Are ggplot Objects

Every enrichplot function returns a ggplot object, so chain ggplot2 modifiers and save with `ggsave`. Generic grammar (themes, scales, faceting) lives in data-visualization/ggplot2-fundamentals.

```r
p <- dotplot(ego, showCategory = 20) + scale_color_viridis_c() + ggtitle('GO BP enrichment')
ggsave('fig.pdf', p, width = 10, height = 8)
```

## Per-Method Failure Modes

### Raw top-20 redundancy artifact
**Trigger:** `dotplot(ego, showCategory=20)` straight from enrichGO on GO results. **Mechanism:** the GO DAG guarantees a real signal surfaces as a nested cluster of overlapping terms driven by the same genes. **Symptom:** twenty bars/dots that are "cell cycle," "cell cycle process," "mitotic cell cycle," "cell division" - one theme repeated. **Fix:** `simplify()` for a shorter list, or `pairwise_termsim` -> `emapplot`/`treeplot` to show the structure, or REVIGO/EnrichmentMap.

### Missing pairwise_termsim
**Trigger:** `emapplot(ego)` or `treeplot(ego)` without the precursor. **Mechanism:** these read `x@termsim`, an empty slot until populated. **Symptom:** an error about a missing termsim slot, or an empty map. **Fix:** `ego_ts <- pairwise_termsim(ego)` first, every time.

### Barplot on gseaResult / dropped NES sign
**Trigger:** coercing a gseaResult to a data frame and bar/dot-plotting |NES| or a p-value ramp. **Mechanism:** a bar from zero is unsigned; |NES| merges activated and suppressed pathways. **Symptom:** a figure that hides that half the pathways are suppressed. **Fix:** there is deliberately no barplot for gseaResult; use a diverging color-by-NES dotplot, ridgeplot, or gseaplot2.

### GeneRatio read as effect size
**Trigger:** "term A has GeneRatio 0.6 so it is strongly over-represented." **Mechanism:** GeneRatio is k/n, not the fold enrichment (k/n)/(M/N). **Symptom:** a giant uninformative term ranked above a small specifically-enriched one. **Fix:** use `x='FoldEnrichment'` (or `GeneRatio/BgRatio`) when specificity is the point; report the p-value as the test statistic.

### Default-ordering misread
**Trigger:** reading the top dot of a default dotplot as "most significant." **Mechanism:** `orderBy='x'` orders by the x variable (GeneRatio), not p.adjust. **Symptom:** a low-significance high-GeneRatio term presented as the headline. **Fix:** order/color by p.adjust explicitly, or state the ordering in the caption.

### Over-trimmed showCategory
**Trigger:** `showCategory=20` when 200 terms passed FDR. **Mechanism:** showCategory truncates to a top-N window by whatever orderBy used. **Symptom:** a 10% sample read as the complete result. **Fix:** report the total significant count and selection criterion in the caption; the figure is a window, not a census.

### Pinned deprecated enrichplot args
**Trigger:** copying `circular=TRUE`, `colorEdge=TRUE`, `cex_label_gene=`, `cex_label_category=`, or `group_category=` from a pre-2024 tutorial. **Mechanism:** enrichplot 1.25.5+ moved cnet/emap/goplot to ggtangle and removed those arguments. **Symptom:** an unused-argument error or a silently ignored arg. **Fix:** `?cnetplot` / `?emapplot` and use the current arguments (`color_item`, `size_category`, `node_label`, `node_label_size`, `min_edge`).

### Wang/IC similarity on non-GO results
**Trigger:** `pairwise_termsim(kegg_result, method='Wang')`. **Mechanism:** Resnik/Lin/Rel/Jiang/Wang require the GO DAG and a GOSemSimDATA object. **Symptom:** an error or a meaningless similarity for KEGG/Reactome/custom sets. **Fix:** use `method='JC'` (gene overlap) for any non-GO gene set.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `showCategory = 10-30` | enrichplot defaults (dotplot 10, emapplot/treeplot 30) | more terms become unreadable; always report the total significant count alongside |
| `pairwise_termsim(method='JC')` default | enrichplot | Jaccard on gene overlap; works for any gene-set type; non-JC are GO-only |
| `simplify(cutoff=0.7)` | clusterProfiler / GOSemSim (Yu 2010 *Bioinformatics* 26:976) | semantic-similarity redundancy cutoff; lower keeps more terms (lives in go-enrichment) |
| `emapplot(min_edge=0.2)` | enrichplot | draw a term-term edge only above this overlap; if everything still connects, the result is redundant |
| `treeplot(nCluster=5, cluster_method='ward.D')` | enrichplot | deterministic Ward cut into 5 named groups; an explicit, reproducible alternative to emapplot's stochastic layout |
| cnetplot <=5-8 terms | enrichplot (showCategory default 5) | the bipartite layout hairballs past ~8 terms |
| diverging color centered at 0 for NES | Subramanian 2005 *PNAS* 102:15545 | NES is signed; a sequential p-value ramp hides activation vs suppression |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| emapplot/treeplot error about a missing termsim slot | skipped `pairwise_termsim()` | run `x <- pairwise_termsim(x)` first |
| `unused argument (circular = TRUE)` in cnetplot | pre-1.25.5 args under ggtangle backend | `?cnetplot`; use `color_item`/`node_label`/`size_category` |
| no applicable method for 'barplot' on gseaResult | GSEA has no barplot method by design | use a NES dotplot, ridgeplot, or gseaplot2 |
| top dot is not the most significant | default `orderBy='x'` orders by GeneRatio | order/color by p.adjust explicitly |
| dotplot terms all look modestly enriched | GeneRatio is not fold enrichment | `dotplot(ego, x='FoldEnrichment')` |
| Wang similarity errors on KEGG terms | IC/graph methods need the GO DAG | `pairwise_termsim(x, method='JC')` |
| gene labels are Entrez IDs not symbols | object not made readable | `setReadable(x, OrgDb, 'ENTREZID')` before plotting |
| two analysts get different emapplot modules | different `method=` / `min_edge=` | record both in the caption; the clustering is a choice |

## References

- Wu T, Hu E, Xu S, et al. 2021. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2:100141.
- Yu G, Li F, Qin Y, Bo X, Wu Y, Wang S. 2010. GOSemSim: an R package for measuring semantic similarity among GO terms and gene products. *Bioinformatics* 26:976-978.
- Supek F, Bosnjak M, Skunca N, Smuc T. 2011. REVIGO summarizes and visualizes long lists of gene ontology terms. *PLoS One* 6:e21800.
- Merico D, Isserlin R, Stueker O, Emili A, Bader GD. 2010. Enrichment Map: a network-based method for gene-set enrichment visualization and interpretation. *PLoS One* 5:e13984.
- Reimand J, Isserlin R, Voisin V, et al. 2019. Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap. *Nat Protoc* 14:482-517.
- Subramanian A, Tamayo P, Mootha VK, et al. 2005. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545-15550.

## Related Skills

- go-enrichment - Produces the enrichResult; owns simplify() the GO-DAG dedup
- gsea - Produces the gseaResult; owns the enrichment score and leading-edge concept
- kegg-pathways - KEGG enrichResult/gseaResult to plot (pathview pathway-diagram overlay lives there)
- reactome-pathways - Reactome enrichResult/gseaResult to plot
- wikipathways - WikiPathways enrichResult/gseaResult to plot
- data-visualization/ggplot2-fundamentals - Generic ggplot2 grammar for the returned objects
- workflows/expression-to-pathways - End-to-end DE-to-enrichment-to-figure pipeline
