---
name: bio-pathway-enrichment-foundations
description: Chooses the enrichment generation before any tool runs, mapping the input shape to a method class - a pre-selected gene list plus a background to over-representation analysis (ORA, hypergeometric), a ranked statistic for all genes to gene set enrichment (GSEA), a signed signaling topology to pathway-topology (SPIA) - then making the null explicit (competitive vs self-contained, gene vs subject sampling) and running a trustworthiness checklist (testable-gene universe, FDR, redundancy collapse, leading-edge check, version reporting). Covers why every clusterProfiler GSEA is the inter-gene-correlation-uncorrected competitive null, why the background not the gene list decides ORA significance, and why no method is universally best. Use when deciding ORA vs GSEA vs topology, which gene-set DB, whether a result is trustworthy, or which null a tool computes. For ORA see go-enrichment, GSEA see gsea, databases kegg-pathways/reactome-pathways/wikipathways; the ranking comes from differential-expression/de-results.
tool_type: r
primary_tool: clusterProfiler
---

## Version Compatibility

Reference examples tested with: clusterProfiler 4.18+, org.Hs.eg.db 3.22+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

This skill is the conceptual spine for the category; it routes the per-tool calls to the sibling skills rather than running them, so the binding versions are the clusterProfiler stack the siblings share. Note that KEGG and WikiPathways query LIVE databases (internet-dependent, not reproducible across data releases) while GO and Reactome use local version-pinned annotation - record the database release and access date for any result.

# Enrichment Foundations -- Choosing the Generation and Defending the Result

**"Which pathways are enriched here?"** -> Pick the enrichment GENERATION from the input shape, make the null explicit, and run the trustworthiness checklist - because an enrichment result is not a discovery, it is a claim about a gene list versus a specific database version under a null that was probably not chosen on purpose.
- R: assemble the input, decide ORA (`enrichGO`/`enrichKEGG`) vs GSEA (`gseGO`/`GSEA`) by whether a ranking exists, then defend the background, FDR, and leading edge

Scope: the meta-decision - which generation (ORA/FCS/topology) is valid, which null each tool computes, which database fits the question, and whether the result survives the checklist. ORA mechanics -> go-enrichment. GSEA mechanics, ranking metric, CAMERA/ROAST -> gsea. Per-database content, IDs, topology -> kegg-pathways, reactome-pathways, wikipathways. Redundancy collapse and plots -> enrichment-visualization. The gene list / ranking statistic itself -> differential-expression/de-results.

## The Single Most Important Modern Insight -- An Enrichment Result Is a Claim Conditioned on a Database, a Universe, and a Null, Not a Discovery

The same gene list, run against a different database version, with a different background universe, under a competitive instead of a self-contained null, with gene-permutation instead of phenotype-permutation, gives a different list of "enriched pathways" - and all of them can be technically valid p-values, because they answer DIFFERENT questions (Goeman & Buhlmann 2007 *Bioinformatics* 23:980). The decision-grade posture treats the deliverable as a defended claim, not a printout, and three facts drive every choice:

1. **The first decision is the GENERATION, and the input sets it, not taste.** A pre-selected list with no ranking -> ORA (1st generation, hypergeometric). A ranked statistic for nearly all genes -> functional class scoring / GSEA (2nd generation), preferred over ORA for that input shape (do not binarize a ranking into a list to force ORA). A ranked list plus a curated signed signaling graph -> pathway topology / SPIA (3rd generation) (Khatri 2012 *PLoS Comput Biol* 8:e1002375).
2. **The null is usually wrong by default and unannounced.** ORA and every clusterProfiler/fgsea preranked GSEA test a COMPETITIVE null by GENE sampling, which assumes genes are independent. They are not - co-regulated pathway genes are positively correlated, so the set-statistic variance is underestimated and the p-values are anti-conservative (Goeman & Buhlmann 2007). Treat these as a discovery screen; CAMERA (Wu & Smyth 2012 *Nucleic Acids Res* 40:e133) is the correlation-aware fix.
3. **The background, not the gene list, decides ORA significance.** The hypergeometric p-value is fully determined by the universe (the denominator). Using the whole genome when only the expressed/testable genes were eligible to be DE is the single most common published error (Wijesooriya 2022 *PLoS Comput Biol* 18:e1009935), and clusterProfiler's default universe is silently all-annotated-genes.

## The Three Generations (Khatri Taxonomy)

The classification every method falls into, set by the available input (Khatri 2012 *PLoS Comput Biol* 8:e1002375).

| Generation | Input required | Question / mechanism | Owns / route to |
|------------|----------------|----------------------|-----------------|
| 1st - Over-Representation (ORA) | a gene LIST + a background universe | is set S over-represented among the hits more than chance? hypergeometric / Fisher on a 2x2 table | enrichGO/enrichKEGG/enrichPathway/enrichWP -> go-enrichment, kegg-pathways, reactome-pathways, wikipathways |
| 2nd - Functional Class Scoring (FCS) | a per-gene STATISTIC for (nearly) all genes, no cutoff | are set S genes systematically shifted along the ranking? running-sum ES + permutation | gseGO/gseKEGG/GSEA/fgsea -> gsea |
| 3rd - Pathway Topology (PT) | ranked list + a curated signed/directed graph | is the pathway IMPACTED given where in the cascade the change lands? perturbation propagation | SPIA/graphite (KEGG/Reactome, model organisms) -> kegg-pathways |

ORA throws away magnitude (binarizes at a cutoff), assumes gene independence, and is acutely sensitive to the background. FCS fixes the first two: it uses the full ranking and detects coordinated WEAK signals ORA misses. PT is the most biologically faithful and the least general (needs a curated directed graph). The agent's first move is to name the generation from the input, never from preference.

## Tool Taxonomy

| Source / method | Citation | Mechanism / role | When |
|-----------------|----------|------------------|------|
| GO ORA (enrichGO) | Ashburner 2000 *Nat Genet* 25:25; Wu 2021 *The Innovation* 2:100141 | hypergeometric vs the GO DAG; local OrgDb | function annotation of a gene LIST; broadest coverage |
| Preranked GSEA (gseGO/GSEA) | Subramanian 2005 *PNAS* 102:15545 | running-sum ES over a ranked vector; gene permutation | all genes carry a statistic; no cutoff; subtle coordinated shifts |
| fgsea engine | (engine behind gseGO/gseKEGG) | fast preranked permutation | the default engine; preranked = gene-sampling |
| CAMERA | Wu & Smyth 2012 *Nucleic Acids Res* 40:e133 | competitive test that estimates and corrects inter-gene correlation (VIF) | a CALIBRATED competitive test when gene-sampling anti-conservatism matters -> gsea |
| ROAST / fry | (limma) | self-contained rotation test; reports directional vs mixed | "is this set involved at all" on the real matrix -> gsea |
| GSVA / ssGSEA | (transform, not a test) | per-sample pathway activity score | features for clustering/modeling; never a p-value -> gsea |
| GOseq | Young 2010 *Genome Biol* 11:R14 | length-bias-corrected ORA (Wallenius null) | RNA-seq with gene-length selection bias -> go-enrichment |
| SPIA / graphite | Tarca 2009 *Bioinformatics* 25:75 | over-representation + perturbation on the signed KEGG graph | signaling directionality / causality -> kegg-pathways |
| MSigDB (H/C2/C5) | Liberzon 2015 *Cell Syst* 1:417; Subramanian 2005 *PNAS* 102:15545 | curated gene-set collections fed via TERM2GENE | Hallmark H is the low-redundancy default; C2:CP/C5 double-count KEGG/GO -> gsea |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| A ranked statistic (t, signed -log10 p, shrunken log2FC) for nearly all genes | GSEA (gseGO/gseKEGG/GSEA) -> gsea | uses the full ranking; no arbitrary cutoff; catches coordinated weak signals |
| A pre-selected list (co-expression module, GWAS loci, screen hits) + a defensible background | ORA (enrichGO/enrichKEGG) -> go-enrichment, kegg-pathways | no ranking exists; binarized membership is all that is available |
| Broad function annotation | enrichGO / gseGO -> go-enrichment, gsea | GO is the broadest local resource |
| Metabolic / signaling pathway maps | enrichKEGG / gseKEGG -> kegg-pathways | KEGG pathway maps (live DB) |
| Reaction-level, peer-reviewed, reproducible offline | enrichPathway -> reactome-pathways | local version-pinned reactome.db |
| Community-curated, broad species coverage | enrichWP -> wikipathways | versioned GMT; variable curation depth |
| A signed signaling topology + fold changes, human/mouse | SPIA / graphite -> kegg-pathways | directionality and cascade position, not just membership |
| RNA-seq with strong gene-length bias | GOseq -> go-enrichment | length-aware null |
| Per-sample pathway activity SCORE (a model feature, not a test) | ssGSEA / GSVA -> gsea | a transform; its scores are not enrichment p-values |
| A calibrated competitive test (correlation-aware) on the expression matrix | CAMERA -> gsea | corrects the gene-sampling anti-conservatism |
| The DE list / ranking statistic itself | -> differential-expression/de-results | that is upstream, not enrichment |
| FDR / p.adjust method theory | -> experimental-design/multiple-testing | the rationale for BH vs Bonferroni lives there |

Default when uncertain: if a ranking exists for almost all genes, run GSEA; otherwise run ORA with the testable-gene universe, FDR-correct, collapse redundancy before interpreting, and inspect the leading-edge genes.

## The Null-Hypothesis Layer

The conceptual core a tutorial skips. Two orthogonal distinctions (Goeman & Buhlmann 2007 *Bioinformatics* 23:980) classify essentially every method and expose why the common ones are mis-calibrated.

**What is being tested (competitive vs self-contained).** A SELF-CONTAINED null is "no gene in set S is associated with the phenotype" - S against zero effect, genes outside S irrelevant; rejecting it means "something in S is differential." A COMPETITIVE null is "the genes in S are no more associated than the genes not in S" - S against its complement; rejecting it means "S is MORE enriched than background." These are different scientific questions: a competitive test can be non-significant simply because everything is changing, and a self-contained test can be significant for a set no more interesting than average. State which one the result answers.

**Where the randomness comes from (gene-sampling vs subject-sampling).** SUBJECT sampling permutes the phenotype labels across samples; each permuted dataset preserves the gene-gene covariance, so this null respects the real correlation structure but needs a sample-level design with adequate replication. GENE sampling treats genes as exchangeable (draw random same-size sets, or permute gene labels); it is the only option with a bare ranked list and no sample data, and it assumes genes are INDEPENDENT.

**The load-bearing consequence.** A competitive test calibrated by gene sampling has the WRONG variance: the variance of the mean of n positively-correlated scores is (sigma^2/n)[1 + (n-1)*rho-bar], and gene-sampling implicitly sets rho-bar = 0, so it underestimates the variance and the p-values are anti-conservative (too many false positives) exactly on the correlated pathway sets of interest. **Every clusterProfiler GSEA (gseGO/gseKEGG/GSEA) and fgsea is preranked = gene-permutation = this uncorrected competitive null.** clusterProfiler does NOT expose phenotype permutation, so in this ecosystem the calibrated subject-sampling GSEA is not available - which is acceptable for a discovery screen but must be stated. CAMERA (Wu & Smyth 2012 *Nucleic Acids Res* 40:e133) estimates the inter-gene correlation and corrects the statistic; ROAST/fry give a self-contained subject-sampling test on the real matrix. Phenotype permutation preserves correlation but needs adequate n.

## Is My Result Trustworthy -- The Checklist

Run before reporting any enrichment, in this order. Each item maps to a documented, published failure (Wijesooriya 2022; Reimand 2019 *Nat Protoc* 14:482).

1. **Correct testable-gene universe.** For ORA, the universe is the genes that passed the SAME filter that produced the list (e.g. the rownames of the filtered DESeq2/edgeR results), NOT the genome. Leaving `universe=` at clusterProfiler's default measures expression bias, not enrichment.
2. **FDR correction, adjusted p reported.** Report `p.adjust`/`qvalue`, never nominal `pvalue`. `pAdjustMethod='BH'` is the clusterProfiler default and is valid under positive dependence.
3. **Redundancy != replication.** A list of 40 significant GO terms is usually ~3 biological stories told 40 times (GO's true-path rule and pathway overlap). Collapse with semantic similarity (`simplify`) or EnrichmentMap BEFORE interpreting, and report clusters of terms, not the raw count -> enrichment-visualization.
4. **Leading-edge / multifunctionality check.** Inspect the core genes (GSEA `core_enrichment`, ORA `geneID`). If the same 3-5 multifunctional hub genes explain the top 20 sets, that is ONE finding wearing twenty pathway costumes, not twenty (Gillis & Pavlidis 2011 *PLoS One* 6:e17258).
5. **Version and parameter reporting.** Record package versions, the database/ontology release + access date (GO/KEGG annotations drift, so interpretation drifts - Tomczak 2018 *Sci Rep* 8:5115), the ranking metric, `pAdjustMethod`, and the universe. Without this the result is unreproducible.
6. **Prefer consensus over single-method certainty.** No method is universally best, so a robust finding survives a second method class (EGSEA logic, Alhamdoosh 2017 *Bioinformatics* 33:414).

## The Gene-Set Database Landscape

GO and MSigDB are the cross-cutting collections everything else samples; KEGG/Reactome/WikiPathways are pathway sources owned by their own skills. Route the per-database mechanics OUT; do not re-teach them here.

- **GO** (Ashburner 2000 *Nat Genet* 25:25) is a controlled VOCABULARY of function (BP/MF/CC), each a directed acyclic graph - not mechanistic pathways. The true-path rule (a gene annotated to a term is annotated to all ancestors) creates the redundancy go-enrichment owns. -> go-enrichment
- **MSigDB** (Liberzon 2015 *Cell Syst* 1:417) is the curated meta-collection GSEA samples. Hallmark (H) is 50 low-redundancy sets and the sensible first-pass default; C2:CP and C5:GO double-count KEGG/Reactome/GO if those are also run natively. MSigDB is human-centric; access via `msigdbr` and feed a TERM2GENE frame to `GSEA`. -> gsea
- **KEGG** - manually drawn metabolic + signaling maps, topology-capable, live REST API. -> kegg-pathways
- **Reactome** - expert-curated, externally reviewed reactions, local reactome.db, reproducible. -> reactome-pathways
- **WikiPathways** - community-curated, broad species, versioned GMT, variable curation depth. -> wikipathways

Never treat "GO enrichment" and "KEGG enrichment" as the same kind of object: GO has thousands of redundant terms (a heavy multiple-testing burden), pathway DBs have hundreds of mechanistic sets.

## Per-Method Failure Modes

### Whole-genome background on an ORA
**Trigger:** running enrichGO/enrichKEGG with `universe=` at the default while only ~12k genes were expressed/testable. **Mechanism:** the hypergeometric p-value is fully set by the denominator; an over-large universe inflates over-representation of any set whose members happen to be expressed in the tissue. **Symptom:** a long list of tissue-specific terms with tiny p-values that runs without warning. **Fix:** set `universe` to the genes that passed the same filter as the DE list.

### Preranked-GSEA p-value reported as a calibrated competitive test
**Trigger:** reporting a gseGO/fgsea p.adjust as evidence the set is specially involved. **Mechanism:** preranked GSEA is competitive + gene-sampling, the inter-gene-correlation-uncorrected null with the wrong variance. **Symptom:** an inflated false-positive rate under permuted phenotypes (Geistlinger 2021 *Brief Bioinform* 22:545). **Fix:** treat it as a discovery screen; run CAMERA for a calibrated competitive test, or ROAST/fry for a self-contained one.

### Redundancy mistaken for replication
**Trigger:** reading 40 overlapping significant terms as 40 independent findings. **Mechanism:** GO's true-path rule and pathway overlap mean the same genes drive many sets; BH is valid but controls the FDR of tests, not the redundancy of findings. **Symptom:** the top terms are obvious parent/child variants of one theme. **Fix:** collapse with `simplify`/EnrichmentMap and report term clusters -> enrichment-visualization.

### Multifunctionality artifact
**Trigger:** a few highly-expressed hub genes (easily called DE) belong to dozens of sets. **Mechanism:** multifunctionality alone predicts apparent enrichment (Gillis & Pavlidis 2011). **Symptom:** the same 3-5 genes appear in the `core_enrichment`/`geneID` of most top sets. **Fix:** inspect the leading edge; if a handful of genes explain the top sets, report one finding, not many.

### Study/annotation bias read as biology
**Trigger:** "cancer / immune / apoptosis pathways enriched" across unrelated experiments. **Mechanism:** genes are annotated in proportion to how studied they are, not how important (Haynes 2018 *Sci Rep* 8:1362), so well-annotated sets light up everywhere. **Symptom:** the same canonical pathways enrich regardless of the biology. **Fix:** treat densely-annotated-set enrichment as the field's null result; weight specific, less-studied hits and confirm with the leading edge.

### Unreported database version
**Trigger:** no record of the GO/KEGG/Reactome release or access date. **Mechanism:** annotations evolve between releases and the interpretation drifts (Tomczak 2018). **Symptom:** the same code returns different pathways next year and the paper cannot be reproduced. **Fix:** record package + database release + access date in a provenance block.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Universe = the testable genes (same filter as the DE list), not the genome | Wijesooriya 2022 *PLoS Comput Biol* 18:e1009935 | the ORA p-value is fully determined by the denominator; the genome inflates it |
| pvalueCutoff=0.05 (filters on p.adjust by default) | clusterProfiler default | standard FDR gate; the enrichResult cutoff applies to the adjusted p, not nominal |
| qvalueCutoff=0.2 | clusterProfiler default | secondary q-value gate on top of p.adjust |
| pAdjustMethod='BH' | clusterProfiler default; Benjamini-Hochberg | valid under positive dependence (overlapping sets); Bonferroni is over-conservative here |
| minGSSize=10, maxGSSize=500 | enrichGO/gseGO defaults | drop tiny sets that overfit and huge sets that always enrich |
| GSEA permutation seed fixed (e.g. set.seed(123)) | reproducibility | permutation p-values drift run-to-run without a fixed seed |
| Leading-edge concentration: if 3-5 genes explain the top ~20 sets | Gillis & Pavlidis 2011 *PLoS One* 6:e17258 | a multifunctionality artifact, not independent findings |
| No single method is best | Tarca 2013 *PLoS One* 8:e79217; Geistlinger 2021 *Brief Bioinform* 22:545 | benchmarks find no method dominates sensitivity + prioritization + specificity; do not claim one beats another by a number |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Implausibly many tissue-specific terms enriched | whole-genome background | set `universe` to the testable genes |
| GSEA hit fails to replicate | preranked gene-sampling null (discovery screen) | confirm with CAMERA/ROAST; report which null was run |
| 40 significant GO terms, all variants of one theme | GO redundancy read as replication | `simplify`/EnrichmentMap; report clusters |
| Top 20 sets all driven by the same few genes | multifunctional hub genes | inspect leading edge; report one finding |
| Same canonical pathways enrich regardless of biology | annotation/study bias | weight specific hits; treat well-studied sets as the null |
| Result not reproducible next year | live DB changed / version unrecorded | pin and report database release + access date |
| Reporting a ssGSEA/GSVA score as a p-value | a transform, not a hypothesis test | use the scores as model features only |

## References

- Khatri P, Sirota M, Butte AJ. 2012. Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Comput Biol* 8:e1002375.
- Goeman JJ, Buhlmann P. 2007. Analyzing gene expression data in terms of gene sets: methodological issues. *Bioinformatics* 23:980-987.
- Subramanian A, Tamayo P, Mootha VK, et al. 2005. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102:15545-15550.
- Wu D, Smyth GK. 2012. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res* 40:e133.
- Young MD, Wakefield MJ, Smyth GK, Oshlack A. 2010. Gene ontology analysis for RNA-seq: accounting for selection bias. *Genome Biol* 11:R14.
- Tarca AL, Draghici S, Khatri P, et al. 2009. A novel signaling pathway impact analysis. *Bioinformatics* 25:75-82.
- Tarca AL, Bhatti G, Romero R. 2013. A comparison of gene set analysis methods in terms of sensitivity, prioritization and specificity. *PLoS One* 8:e79217.
- Geistlinger L, Csaba G, Santarelli M, et al. 2021. Toward a gold standard for benchmarking gene set enrichment analysis. *Brief Bioinform* 22:545-556.
- Wijesooriya K, Jadaan SA, Perera KL, Kaur T, Ziemann M. 2022. Urgent need for consistent standards in functional enrichment analysis. *PLoS Comput Biol* 18:e1009935.
- Reimand J, Isserlin R, Voisin V, et al. 2019. Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap. *Nat Protoc* 14:482-517.
- Gillis J, Pavlidis P. 2011. The impact of multifunctional genes on "guilt by association" analysis. *PLoS One* 6:e17258.
- Haynes WA, Tomczak A, Khatri P. 2018. Gene annotation bias impedes biomedical research. *Sci Rep* 8:1362.
- Tomczak A, Mortensen JM, Winnenburg R, et al. 2018. Interpretation of biological experiments changes with evolution of the Gene Ontology and its annotations. *Sci Rep* 8:5115.
- Alhamdoosh M, Ng M, Wilson NJ, et al. 2017. Combining multiple tools outperforms individual methods in gene set enrichment analyses. *Bioinformatics* 33:414-424.
- Ashburner M, Ball CA, Blake JA, et al. 2000. Gene Ontology: tool for the unification of biology. *Nat Genet* 25:25-29.
- Liberzon A, Birger C, Thorvaldsdottir H, et al. 2015. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Syst* 1:417-425.
- Wu T, Hu E, Xu S, et al. 2021. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The Innovation* 2:100141.

## Related Skills

- go-enrichment - ORA mechanics: enrichGO, the GO DAG, the background universe, GOseq length bias
- gsea - FCS mechanics: ranking metric, running-sum ES, leading edge, CAMERA/ROAST, GSVA/ssGSEA
- kegg-pathways - KEGG pathway/module enrichment and SPIA topology
- reactome-pathways - Reactome curated-pathway ORA and GSEA
- wikipathways - WikiPathways community-pathway enrichment
- enrichment-visualization - Dot/cnet/emap plots and redundancy collapse
- differential-expression/de-results - Source of the gene list and the ranking statistic
- experimental-design/multiple-testing - FDR / p.adjust method theory
- workflows/expression-to-pathways - End-to-end DE-to-enrichment pipeline
