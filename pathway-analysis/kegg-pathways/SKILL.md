---
name: bio-pathway-kegg-pathways
description: Tests gene lists, ranked vectors, and fold-change vectors against KEGG pathways and modules with clusterProfiler enrichKEGG/enrichMKEGG (ORA), gseKEGG (GSEA), and SPIA/graphite (signed-topology perturbation) in R. Owns the third pathway-analysis generation because KEGG ships signed directed signaling topology (KGML). Covers why a KEGG result is a timestamped join against a live REST API (irreproducible unless pinned with a gson snapshot, not the stale 2012 KEGG.db), why enrichKEGG keyType is kegg/ncbi-geneid not OrgDb ENSEMBL/SYMBOL (zero hits), why organism is a KEGG code (hsa, pae) with prokaryotic locus tags, and why SPIA works only on signaling maps. Use when finding enriched KEGG pathways or modules, scoring signed pathway perturbation, analyzing prokaryotes or non-model organisms via locus tags or KO, comparing conditions with compareCluster, or overlaying data with pathview. The hypergeometric universe lives in go-enrichment; the GSEA engine in gsea.
tool_type: r
primary_tool: clusterProfiler
---

## Version Compatibility

Reference examples tested with: clusterProfiler 4.18+, org.Hs.eg.db 3.18+, gson 0.1+ (snapshot pinning), SPIA 2.50+ and graphite 1.56+ (topology section).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

KEGG is a LIVE DATABASE, not a package. enrichKEGG/enrichMKEGG/gseKEGG query the KEGG REST API (https://rest.kegg.jp/) at call time, so the same code on the same genes returns DIFFERENT pathways months apart as KEGG updates. For any reported result, pin the release with a gson snapshot (below) and record the access date; `use_internal_data=TRUE` does NOT pin the current KEGG (it loads the deprecated 2012 KEGG.db).

# KEGG Pathway and Topology Enrichment

**"Which KEGG pathways are perturbed in my data?"** -> Join genes to KEGG's curated pathway/module gene sets (ORA or GSEA), or propagate fold-changes through KEGG's signed wiring (SPIA) - and pin the KEGG release, because the result is a timestamped query against a moving curation, not a fact about the biology.
- R: `enrichKEGG(gene, organism, keyType)` | `gseKEGG(geneList, organism)` | `spia(de, all, organism)`

Scope: KEGG-specific enrichment across all three generations - membership ORA (enrichKEGG/enrichMKEGG), ranked GSEA (gseKEGG), and signed-topology perturbation (SPIA/graphite). KEGG ID mapping (organism codes, keyType, bitr_kegg, prokaryotic locus tags, KO routing), reproducibility/pinning, and pathview map overlay live here. The hypergeometric test and the universe problem -> go-enrichment. The GSEA running-sum engine and ranking-metric choice -> gsea. Reactome/WikiPathways gene sets -> reactome-pathways, wikipathways. Generic dot/cnet/emap plots -> enrichment-visualization. The DE list and fold-changes -> differential-expression/de-results.

## The Single Most Important Modern Insight -- A KEGG Result Is a Timestamped Join Against a Moving, Partially-Paywalled Curation, Not a Fact About Biology

Two consequences follow, and both are invisible until someone reruns the analysis.

1. **The query is live, so the result is irreproducible unless the release is pinned.** enrichKEGG/gseKEGG/SPIA hit the KEGG REST API at call time; KEGG adds maps, re-annotates genes, and revises edges continuously, so identical code returns a different pathway list next quarter. The fix is a gson snapshot: `gson_KEGG('hsa')` downloads the current KEGG pathway/module sets into a GSON object, `write.gson()`/`read.gson()` persist it, and the generic `enricher(gene, gson=k)` / `GSEA(geneList, gson=k)` run frozen and offline against it. Record the access date. `use_internal_data=TRUE` is NOT this fix - it silently reaches for the deprecated 2012 `KEGG.db`, which is the wrong, stale snapshot.

2. **KEGG is the only mainstream database shipping signed, directed signaling topology (KGML), which is why this skill owns the third generation of pathway analysis.** ORA and GSEA treat a pathway as an unordered bag of exchangeable genes; SPIA asks a question they structurally cannot pose - given where each gene sits in the wiring and the sign of every edge, how perturbed is this pathway? That requires the topology only KEGG (and a few others via graphite) provides. The discipline: choose the generation by the question (membership? rank? signed perturbation?), match keyType/organism to the actual IDs (locus tags for bacteria, KO for non-model), set the universe to the genes that could have been called DE, and pin the release before publishing.

## Tool Taxonomy (KEGG Across the Three Generations)

| Method | Generation | Engine | Uses log2FC? | Uses topology/direction? | Suitable KEGG maps | Citation |
|--------|-----------|--------|--------------|--------------------------|--------------------|----------|
| enrichKEGG (ORA) | 1st (over-representation) | hypergeometric | no (gene list) | no | all | Wu 2021 *The Innovation* 2:100141; Kanehisa & Goto 2000 *Nucleic Acids Res* 28:27 |
| enrichMKEGG (ORA on modules) | 1st | hypergeometric | no | no | modules (M-numbers) | Wu 2021 *The Innovation* 2:100141 |
| gseKEGG (GSEA) | 2nd (functional class scoring) | fgsea running sum | yes (ranking) | no | all (as sets) | Wu 2021 *The Innovation* 2:100141; engine -> gsea |
| SPIA | 3rd (pathway topology) | pNDE (ORA) x pPERT (perturbation) -> pG | yes (named log2FC) | YES (signed KGML) | SIGNALING only | Tarca 2009 *Bioinformatics* 25:75; Draghici 2007 *Genome Res* 17:1537 |
| graphite + runSPIA | 3rd | SPIA over harmonized graphs | yes | YES | signaling (KEGG/Reactome) | Sales 2012 *BMC Bioinformatics* 13:20 |

The three-generations framing (ORA -> FCS -> pathway topology) is Khatri 2012 *PLoS Comput Biol* 8:e1002375; this skill is the KEGG instantiation of all three (the category README compares the generations across databases).

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Pre-selected gene list, "which KEGG pathways" | enrichKEGG (ORA), set the universe | no ranking available; membership test |
| All genes carry a DE statistic, no clear cutoff | gseKEGG -> gsea | uses the full ranking; no arbitrary cutoff |
| Want WHERE in a broad pathway the signal sits | enrichMKEGG (modules) | M-numbers are tighter functional units |
| Have named log2FC + want signed perturbation on a SIGNALING map | SPIA (or graphite + runSPIA) | propagates fold-changes through the wiring; uses direction |
| Metabolic-pathway question (glycolysis, TCA) | enrichKEGG / gseKEGG | metabolic maps are compound-mediated; SPIA is undefined there |
| Human / mouse / model eukaryote | bitr -> Entrez, keyType='ncbi-geneid' | KEGG gene ID == Entrez for these organisms |
| Bacterial / prokaryotic data | locus tags, keyType='kegg', NO OrgDb/bitr | bacterial KEGG IDs ARE locus tags; no org.*.eg.db exists |
| Non-model organism with no KEGG genome | map to KO, organism='ko' | the universal escape hatch into KEGG pathway space |
| Result must be reproducible / published | gson_KEGG snapshot + enricher/GSEA, record date | live unpinned queries drift; use_internal_data pins the WRONG 2012 db |
| Multiple conditions to compare side by side | compareCluster(fun='enrichKEGG') | one model, faceted dotplot; never compare raw p-values |
| Overlay per-gene data on the KEGG map image | pathview -> render | a KEGG-specific operation; generic plots -> enrichment-visualization |
| The DE list / fold-changes themselves | -> differential-expression/de-results | upstream, not enrichment |

## Prepare the Gene IDs (the Join That Decides Everything)

**Goal:** Get the query genes and the universe into the exact ID type KEGG expects for the organism, because every KEGG failure is a join failure.

**Approach:** For model eukaryotes convert SYMBOL/ENSEMBL to Entrez (KEGG's gene ID for hsa/mmu/rno) and pass keyType='ncbi-geneid'. For prokaryotes pass locus tags directly with keyType='kegg' and no OrgDb. Convert the universe the same way. Passing ENSEMBL/SYMBOL to enrichKEGG returns zero hits silently.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

de <- read.csv('de_results.csv')   # DE list source -> differential-expression/de-results
sig_symbols <- de$gene[de$padj < 0.05 & abs(de$log2FoldChange) > 1]   # padj is the DESeq2 adjusted-p column
sig_entrez  <- bitr(sig_symbols, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)$ENTREZID

# universe = genes that COULD have been called DE (non-NA test statistic), same ID type
universe <- bitr(de$gene[!is.na(de$pvalue)], fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)$ENTREZID
```

`bitr_kegg(geneID, fromType, toType, organism)` converts among KEGG's own ID flavors ('kegg', 'ncbi-geneid', 'ncbi-proteinid', 'uniprot') via the REST conv endpoint - use it when starting from UniProt or NCBI protein IDs. Check KEGG coverage of an organism with `search_kegg_organism('Pseudomonas aeruginosa', by='scientific_name')`.

## Run KEGG ORA (enrichKEGG / enrichMKEGG)

**Goal:** Find KEGG pathways (or modules) over-represented among the query genes relative to the measured universe.

**Approach:** Run enrichKEGG with the correct organism code, keyType, and an explicit universe; enrichKEGG has no `readable` argument, so translate the geneID column to symbols afterward with setReadable (eukaryotes only).

```r
kk <- enrichKEGG(gene=sig_entrez, organism='hsa', keyType='ncbi-geneid',
                 universe=universe, pvalueCutoff=0.05, pAdjustMethod='BH',
                 minGSSize=10, maxGSSize=500, qvalueCutoff=0.2)
kk <- setReadable(kk, OrgDb=org.Hs.eg.db, keyType='ENTREZID')   # eukaryotes only; no OrgDb -> keep raw IDs
head(as.data.frame(kk))   # ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count

mkk <- enrichMKEGG(gene=sig_entrez, organism='hsa', keyType='ncbi-geneid', universe=universe)   # KEGG MODULES (M-numbers)
```

Report `p.adjust`/`qvalue`, not raw `pvalue`. Fold enrichment = GeneRatio / BgRatio. enrichMKEGG tests smaller, sparser sets: higher resolution (which sub-process is hit) but lower power and many genes belong to no module.

## Run KEGG GSEA (gseKEGG)

**Goal:** Find KEGG sets whose genes shift coordinately across the full ranking, with no cutoff.

**Approach:** Build a named numeric vector sorted DECREASING by the ranking metric, fix the seed (gseKEGG defaults `seed=FALSE`), then run gseKEGG. The running-sum engine and the ranking-metric choice are owned by gsea; only the KEGG arguments (organism, keyType) are KEGG-specific.

```r
geneList <- de$log2FoldChange; names(geneList) <- de$entrez   # names = Entrez IDs
geneList <- sort(geneList[!is.na(geneList)], decreasing=TRUE)
set.seed(123)   # gseKEGG seed=FALSE by default; fix it so permutation p-values are reproducible
kk2 <- gseKEGG(geneList=geneList, organism='hsa', keyType='ncbi-geneid', minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)
```

## Run Signed-Topology Perturbation (SPIA) -- the Third Generation

**Goal:** Score how perturbed each SIGNALING pathway is given both the over-representation of DE genes and the propagation of their fold-changes through the signed wiring.

**Approach:** SPIA combines pNDE (the classical over-representation evidence) with pPERT (the probability of the observed total accumulated perturbation tA, computed by propagating log2 fold-changes through KGML activation/inhibition edges) into a single global pG, then FDR-corrects it. It needs a NAMED vector of DE fold-changes plus the universe, and is defined only for signaling maps. graphite is the modern route: it harmonizes node IDs, resolves complexes/families, removes compounds, and can run SPIA over Reactome topology too.

```r
library(SPIA)
sig <- de[de$padj < 0.05, ]   # DE genes only
map <- bitr(sig$gene, 'SYMBOL', 'ENTREZID', org.Hs.eg.db)   # bitr drops/many-to-one: MERGE, never assign as names
de_vec <- setNames(sig$log2FoldChange[match(map$SYMBOL, sig$gene)], map$ENTREZID)
de_vec <- de_vec[!duplicated(names(de_vec))]
res <- spia(de=de_vec, all=universe, organism='hsa', nB=2000, plots=FALSE)   # nB=2000 bootstraps for pPERT
# output cols: Name, ID, pSize, NDE, pNDE, tA, pPERT, pG, pGFdr, pGFWER, Status, KEGGLINK
# Status reports inferred Activated / Inhibited from the sign of tA

# graphite route (decouples from KEGG's bundled data; works on Reactome too)
library(graphite)
db <- pathways('hsapiens', 'kegg')
db <- convertIdentifiers(db, 'ENTREZID')
prepareSPIA(db, 'kegg_hsa_spia')              # writes the pathway dataset file
gr <- runSPIA(de=de_vec, all=universe, 'kegg_hsa_spia')
```

SPIA aborts if more than ~1% of the DE IDs are absent from `all`, so build the universe from the same ID space. The standalone SPIA package also ships a frozen `hsaSPIA` data object that is an OLDER snapshot than a live enrichKEGG query - do not mix the two in one comparison.

## Pin the KEGG Release for Reproducibility

**Goal:** Freeze the KEGG data a result depends on so the analysis is reproducible and runs offline.

**Approach:** Snapshot the current KEGG sets into a GSON object, persist it, and run enrichment against the snapshot with the generic enricher/GSEA (which accept a `gson` argument); record the access date. Do NOT use use_internal_data=TRUE for this.

```r
library(gson)                                      # GSON class + write.gson/read.gson
k <- gson_KEGG('hsa')                              # gson_KEGG is exported by clusterProfiler; downloads current KEGG sets
k@accessed_date <- as.character(Sys.Date())        # the accessed_date slot survives write/read; a base attr() does not
write.gson(k, file.path(tempdir(), 'kegg_hsa.gson'))
k <- read.gson(file.path(tempdir(), 'kegg_hsa.gson'))

kk_pinned  <- enricher(sig_entrez, gson=k, universe=universe)   # frozen ORA, offline, reproducible
gsea_pinned <- GSEA(geneList, gson=k)                            # frozen GSEA against the snapshot
```

## Compare Multiple Conditions

**Goal:** See shared and condition-specific KEGG pathways across groups in one faceted figure.

**Approach:** Pass named gene lists to compareCluster with fun='enrichKEGG'; it fits one model and produces a faceted dotplot. Compare pathway-ID SETS across conditions, never raw p-values (they depend on sample size, DE gene count, and the KEGG release).

```r
clusters <- list(up=up_entrez, down=down_entrez)
ck <- compareCluster(geneClusters=clusters, fun='enrichKEGG', organism='hsa', keyType='ncbi-geneid')
ck <- setReadable(ck, OrgDb=org.Hs.eg.db, keyType='ENTREZID')
# dotplot(ck) -> enrichment-visualization for the plot grammar
```

## Overlay Data on the KEGG Map (pathview)

pathview downloads a KEGG pathway's KGML and image, joins per-gene values to the nodes, and writes a colored map PNG/PDF (a KEGG-specific operation owned here; generic dot/cnet/emap plots route to enrichment-visualization). It writes files to the working directory and queries KEGG live.

```r
library(pathview)
vals <- setNames(de$log2FoldChange, de$entrez)
pathview(gene.data=vals, pathway.id='hsa04110', species='hsa', gene.idtype='entrez')   # writes hsa04110.pathview.png
```

## Per-Method Failure Modes

### ENSEMBL/SYMBOL passed to enrichKEGG
**Trigger:** feeding OrgDb-style ENSEMBL or SYMBOL IDs to enrichKEGG/gseKEGG. **Mechanism:** KEGG's keyType is 'kegg'/'ncbi-geneid'/'ncbi-proteinid'/'uniprot', not an OrgDb keytype, so no IDs join. **Symptom:** zero enriched pathways, no error. **Fix:** bitr to Entrez and set keyType='ncbi-geneid' (eukaryotes), or pass locus tags with keyType='kegg' (prokaryotes).

### Live-query result treated as reproducible
**Trigger:** reporting an enrichKEGG/gseKEGG/SPIA result without pinning the release. **Mechanism:** the REST query returns the CURRENT KEGG, which changes over time. **Symptom:** a rerun months later yields a different pathway list. **Fix:** snapshot with gson_KEGG, run enricher/GSEA against the gson, and record the access date.

### use_internal_data=TRUE believed to pin current KEGG
**Trigger:** setting use_internal_data=TRUE for reproducibility. **Mechanism:** it loads the deprecated 2012 KEGG.db, not a current pin (and may simply fail). **Symptom:** stale or absent pathways unlike the live result. **Fix:** use a gson snapshot instead; treat KEGG.db as legacy-only.

### SPIA on metabolic maps
**Trigger:** running SPIA/graphite topology on glycolysis or other metabolic maps. **Mechanism:** metabolic maps are compound-mediated and give no clean signed gene->gene graph. **Symptom:** meaningless perturbation scores. **Fix:** restrict SPIA to signaling maps; use enrichKEGG/gseKEGG for metabolism.

### Whole-database universe in ORA
**Trigger:** omitting `universe`. **Mechanism:** the default background is all KEGG-annotated genes, biased toward well-studied, metabolically central genes. **Symptom:** inflated significance for pathways enriched in measured/expressed genes (the tissue-specificity artifact). **Fix:** set universe to the genes that could have been called DE, in the same ID type.

### Locus-tag / strain mismatch in prokaryotes
**Trigger:** locus tags from a re-annotated genome or a different strain than KEGG's reference. **Mechanism:** the gene-ID join is exact; drifted locus tags do not match KEGG's `pae`/`eco` genome. **Symptom:** many genes silently dropped, weak or empty enrichment. **Fix:** confirm the organism code and reference genome with search_kegg_organism; align locus tags to KEGG's annotation, or route through KO.

### bitr/OrgDb forced onto bacteria
**Trigger:** running bitr() or setReadable() on a prokaryote. **Mechanism:** no org.*.eg.db exists for most bacteria and there is no Entrez==KEGG identity. **Symptom:** bitr fails or empties the gene list; setReadable errors. **Fix:** pass locus tags directly with keyType='kegg'; keep raw IDs (no setReadable).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| pvalueCutoff=0.05 | enrichKEGG/gseKEGG default | filters on p.adjust by default; standard FDR gate |
| qvalueCutoff=0.2 | clusterProfiler default | secondary q-value gate on enrichResult |
| pAdjustMethod='BH' | clusterProfiler default | Benjamini-Hochberg FDR; less conservative than Bonferroni for discovery |
| minGSSize=10 | enrichKEGG default | drop tiny sets that overfit and give unstable p-values |
| maxGSSize=500 | enrichKEGG default | drop very broad sets that always 'enrich' |
| nB=2000 | SPIA default | bootstrap replicates for the pPERT null; raise for stable small p-values |
| SPIA aborts if >1% of DE IDs absent from `all` | Tarca 2009 *Bioinformatics* 25:75 | the perturbation null requires the DE genes live in the universe |
| set.seed before gseKEGG/SPIA | reproducibility | gseKEGG seed=FALSE and SPIA bootstrap are stochastic; fix the seed |
| ID-conversion loss > ~15% | practice heuristic | report the bitr conversion rate; heavy loss makes the result unreliable |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| enrichKEGG returns 0 pathways | ENSEMBL/SYMBOL passed, or wrong organism code, or KEGG API unreachable | bitr to Entrez + keyType='ncbi-geneid'; verify code with search_kegg_organism; check network |
| `setReadable` errors | no OrgDb for the organism (prokaryote) | skip setReadable; keep raw KEGG IDs |
| `gson=` rejected by enrichKEGG | enrichKEGG/gseKEGG have no gson argument | pass the gson to the generic enricher()/GSEA() instead |
| Different pathways on rerun | live KEGG changed between runs | pin with a gson snapshot and record the access date |
| SPIA: "more than 1% of de IDs not in all" | DE IDs not a subset of the universe | build de and all from the same ID space |
| SPIA gives nonsense on glycolysis | topology on a metabolic map | use enrichKEGG/gseKEGG; SPIA is signaling-only |
| Bacterial list gives 0 hits | Entrez/bitr forced onto a prokaryote | pass locus tags with keyType='kegg', no OrgDb |

## References

- Kanehisa M, Goto S. 2000. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Res* 28:27-30.
- Kanehisa M, Furumichi M, Sato Y, et al. 2023. KEGG for taxonomy-based analysis of pathways and genomes. *Nucleic Acids Res* 51:D587-D592.
- Wu T, Hu E, Xu S, et al. 2021. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation* 2:100141.
- Tarca AL, Draghici S, Khatri P, et al. 2009. A novel signaling pathway impact analysis (SPIA). *Bioinformatics* 25:75-82.
- Draghici S, Khatri P, Tarca AL, et al. 2007. A systems biology approach for pathway level analysis. *Genome Res* 17:1537-1545.
- Sales G, Calura E, Cavalieri D, Romualdi C. 2012. graphite - a Bioconductor package to convert pathway topology to gene network. *BMC Bioinformatics* 13:20.
- Luo W, Brouwer C. 2013. Pathview: an R/Bioconductor package for pathway-based data integration and visualization. *Bioinformatics* 29:1830-1831.
- Khatri P, Sirota M, Butte AJ. 2012. Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Comput Biol* 8:e1002375.

## Related Skills

- go-enrichment - Hypergeometric ORA and the background-universe problem
- gsea - GSEA running-sum engine and ranking-metric choice (gseKEGG)
- reactome-pathways - Reactome curated-pathway enrichment (reproducible local DB)
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - Dot/bar/cnet/emap/ridge plots of enrichment results
- differential-expression/de-results - Source of the gene list and the fold-changes
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
