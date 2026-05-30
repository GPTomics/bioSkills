---
name: bio-genome-intervals-overlap-significance
description: Tests whether two genomic interval sets overlap (colocalize) more than expected by chance using a permutation test against a structured-genome null model. Covers bedtools fisher (analytic 2x2 screen), bedtools shuffle + jaccard permutation, GAT (isochore/GC-conditioned simulation with FDR), regioneR (flexible permutation, randomizeRegions vs circularRandomizeRegions, localZScore), LOLA (universe-relative Fisher against a region database), and GREAT/rGREAT (regulatory-domain binomial + hypergeometric for ontology-from-regions). Stresses the universe/background choice, matched background, blacklist exclusion, and multiple-testing control. Use when asking whether peaks/regions are enriched at enhancers/TFBS/features, scoring region-set colocalization or region-set enrichment, comparing CNV/SV concordance, or turning an overlap count into a defensible p-value.
tool_type: mixed
primary_tool: regioneR
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, pybedtools 0.10+, regioneR 1.36+ (Bioconductor 3.18+), GAT 1.3+, LOLA 1.30+, rGREAT 2.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('regioneR')` then `?permTest` to verify parameters

Bioconductor packages (regioneR, LOLA, rGREAT) are version-pinned to the Bioconductor release, not just the package version - record the Bioconductor release with results. rGREAT 2.x runs GREAT locally with general background handling; rGREAT 1.x only proxied the (whole-genome-default) web server. If code throws an error, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Overlap Significance

**"My peaks overlap enhancers a lot - is that more than chance?"** -> Compare the observed overlap to a null distribution from a structured-genome model, not to a uniform-random expectation, and report a permutation p-value/z-score, not a raw count.
- CLI: `bedtools fisher -a A -b B -g genome.txt` (fast screen); `gat-run.py --segments=A --annotations=B --workspace=accessible.bed --isochores=gc.bed --num-samples=10000`
- Python: `BedTool(A).shuffle(g='genome.txt', excl='blacklist.bed').jaccard(B)` looped for a null (pybedtools)
- R: `permTest(A=peaks, B=enhancers, randomize.function=circularRandomizeRegions, evaluate.function=numOverlaps, genome='hg38', mask=blacklist, ntimes=1000)` (regioneR)

## The Single Most Important Modern Insight -- A Raw Overlap Count Is an Observation in Search of a Null, and the Universe Choice Dominates

"847 of 1,000 peaks overlap enhancers" means nothing until the analysis can say what that number would have been *by chance* - and "by chance" is almost never "place the regions uniformly at random on the genome." The genome is structured: genes cluster, GC varies in megabase isochores, mappability is uneven, half the genome is repeat/gap, and query regions are drawn from a *biased universe* (open chromatin, callable space, exons). Two tracks that share nothing but a gene-rich, high-GC, high-mappability habitat overlap far more than uniform-random expectation, and a naive test returns p < 1e-300. The co-localization is real; the *interpretation* ("functional association") is false. Three load-bearing moves:

1. **The universe/background is the lever that moves the whole answer - bigger than the test choice.** Across LOLA (`userUniverse`), GREAT (background regions), GAT (`--workspace`), and regioneR (the `mask`/`resampleRegions` universe), the most consequential choice is *the set of regions the query could have come from*. ATAC/ChIP peaks can only be called in accessible chromatin; their honest universe is "all accessible regions," not the genome. Testing against the whole genome merely rediscovers that open chromatin is gene-rich - every gene-associated annotation lights up, none of it specific. The difference between LOLA/GAT/regioneR on the *same correct universe* is second-order; the difference between a correct universe and a whole-genome universe on the *same tool* is often p~1 vs p~1e-200. Care belongs on the background, not on tool selection.

2. **A correct null preserves three things, or it manufactures significance:** (a) the regions' **size distribution** (a 50 kb domain hits anything; a 200 bp peak rarely does - relocate intervals of the observed sizes, do not sprinkle points); (b) an **accessible workspace** excluding assembly gaps, centromeres, the ENCODE blacklist (Amemiya 2019), and unmappable bins; (c) **local structure** - GC/isochore, gene density, and clustering (GAT `--isochores`; regioneR `circularRandomizeRegions` for autocorrelated regions). The "right" answer is usually *less* significant than the naive test - that deflation is the methodology working.

3. **One ten-minute sanity move catches most false claims:** shuffle the query within the same workspace and re-run the same overlap pipeline. If shuffled regions also overlap the annotation a lot, the "enrichment" is workspace geography, not biology.

## Method Taxonomy

| Tool | Citation | Null model | When |
|------|----------|-----------|------|
| bedtools fisher | Quinlan 2010 *Bioinformatics* | analytic 2x2; estimates the unobserved "in-neither" cell from mean interval size + genome size; ignores genome structure | fast triage screen only - never the reported result |
| bedtools shuffle + jaccard/numOverlaps | Quinlan 2010 *Bioinformatics* | DIY size-preserving permutation; `-incl`/`-excl` make it matched | entry-level permutation in a shell/Python pipeline; full control, more code |
| GAT | Heger 2013 *Bioinformatics* | per-isochore size-preserving simulation in a workspace; GC/composition conditioned; built-in FDR across annotations | composition-aware enrichment vs many tracks with multiple-testing control |
| regioneR | Gel 2016 *Bioinformatics* | flexible permutation; `randomizeRegions` (mask) vs `circularRandomizeRegions` (preserves clustering) vs `resampleRegions` (real universe); any evaluator; `localZScore` | publication-grade, R-native, autocorrelated regions, where-in-the-region probing |
| LOLA | Sheffield & Bock 2016 *Bioinformatics* | NOT a shuffle - Fisher's exact of query vs a region database, relative to a `userUniverse` | ranking enrichment of a peak set against ENCODE/Roadmap region collections |
| GREAT / rGREAT | McLean 2010 *Nat Biotechnol*; Gu & Hubschmann 2023 *Bioinformatics* | gene regulatory-domain (basal 5 kb up/1 kb down, extend to 1 Mb); binomial-over-regions AND hypergeometric-over-genes | GO/ontology enrichment *from regions* (cis-regulatory), not from a gene list |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Quick "is this even worth permuting?" | `bedtools fisher` | one command, analytic; p~1 means stop; low means permute |
| Publication-grade colocalization of two region sets | regioneR `permTest` with a `mask` | flexible, reports p + z-score; circular randomization for clustered regions |
| Need GC/isochore/composition control + FDR over many tracks | GAT with `--workspace` + `--isochores` | per-isochore sampling removes the GC confounder; built-in qvalue |
| Enrich a peak set against a region database (TFBS/chromatin) | LOLA `runLOLA` with a correct `userUniverse` | universe-relative Fisher; ranked odds ratios + q-values |
| GO/ontology terms FROM regions (cis-regulatory) | rGREAT with a real background | regulatory-domain model; require binomial AND hypergeometric |
| GWAS/eQTL statistical colocalization (shared causal variant) | -> causal-genomics/colocalization-analysis | DISTINCT problem: coloc/SuSiE on summary stats, not interval overlap |
| Gene-list (not region) ontology enrichment | -> pathway-analysis/go-enrichment | start from genes; GREAT is the region-based analog |
| Query regions are autocorrelated/clustered | regioneR `circularRandomizeRegions` | rotating the set preserves inter-region spacing; uniform randomization inflates significance |
| CNV/SV concordance between call sets | `bedtools intersect -f 0.5 -r` (50% reciprocal) | "same event" convention; one-sided fractions let a giant call swallow a tiny one |
| Peaks not yet called | -> chip-seq/peak-calling, atac-seq/atac-peak-calling | this skill operates on existing interval sets |

## bedtools fisher - The Analytic Screen (weakest null)

```bash
bedtools fisher -a peaks.bed -b enhancers.bed -g genome.txt   # both inputs sorted; genome file required
```

`fisher` builds a 2x2 table (in-A/not x in-B/not) and runs Fisher's exact test. The trap: it cannot observe the "in-neither" cell (there is no negative class of intervals that do not exist), so it **estimates the table totals from a heuristic on mean interval size and genome size**, assuming intervals are independent points uniformly placeable across the genome - exactly the assumption a structured genome violates. The bedtools docs warn it is prone to **inflation** and advise validating any low p-value by simulation. Treat it as triage: p~1 -> stop; p~1e-50 -> run GAT or regioneR, because the structure-corrected p could be anywhere from 1e-30 to 0.3.

## Permutation Null with bedtools shuffle (entry-level)

**Goal:** Decide whether two interval sets overlap more than chance, controlling for feature size and the accessible workspace.

**Approach:** Compute the observed overlap (jaccard or count), then shuffle one set N times within an include-list / outside a blacklist, recompute each time, and locate the observed value in the resulting null distribution.

```python
import pybedtools

N_PERMUTATIONS = 1000   # >=1000 gives a stable empirical p down to ~0.001; fewer cannot resolve small p
a = pybedtools.BedTool('peaks.bed')
b = pybedtools.BedTool('enhancers.bed').sort()
observed = a.sort().jaccard(b)['jaccard']
null = [a.shuffle(g='genome.txt', incl='accessible.bed', excl='blacklist.bed', chrom=True).sort().jaccard(b)['jaccard'] for _ in range(N_PERMUTATIONS)]
p = (sum(x >= observed for x in null) + 1) / (N_PERMUTATIONS + 1)   # +1 avoids a p of exactly 0 (Phipson & Smyth 2010)
```

`-incl` restricts placement to the accessible workspace (the universe); `-excl` avoids gaps/blacklist; `-chrom` keeps each region on its own chromosome (preserves per-chromosome density). Without `-incl`/`-excl` this collapses to a uniform-random null - the wrong one.

## GAT - Isochore/GC-Conditioned Simulation

**Goal:** Test enrichment against one or many annotation tracks while conditioning on GC/composition and controlling FDR.

**Approach:** Provide the query segments, the annotations, an accessible workspace, and an isochore segmentation; GAT samples size-matched segments *per isochore* and compares observed to sampled overlap, reporting fold, empirical p, and FDR-adjusted q across annotations.

```bash
gat-run.py \
  --segments=peaks.bed \
  --annotations=features.bed \
  --workspace=accessible.bed \
  --isochores=gc_bins.bed \
  --num-samples=10000 \
  --counter=nucleotide-overlap \
  --log=gat.log > gat_results.tsv
```

`--isochores` subdivides the workspace (GC bins, chromatin state, or mappability) so sampling happens *per isochore*, preserving GC/composition confounding rather than averaging it away - GAT's signature over a plain shuffle. Set `--num-samples` (>=10000 for stable small q) and `--counter` explicitly; do not rely on defaults.

## regioneR - Flexible Permutation in R

**Goal:** Get a publication-grade colocalization p-value and z-score, with a null that preserves the structure the query trivially has.

**Approach:** Mask the genome to the workspace, choose a randomizer that concedes the right structure (uniform vs clustering-preserving vs real-universe), permute N times scoring overlaps, then probe *where* the association lives with localZScore.

```r
# Reference: regioneR 1.36+ (Bioconductor 3.18+) | Verify API if version differs
library(regioneR)

N_TIMES <- 1000   # permutation count; >=1000 for a stable empirical p (Gel 2016)
peaks <- toGRanges('peaks.bed')
enhancers <- toGRanges('enhancers.bed')
gam <- getGenomeAndMask(genome = 'hg38', mask = toGRanges('blacklist.bed'))

pt <- permTest(A = peaks, B = enhancers,
               randomize.function = circularRandomizeRegions,   # preserves clustering; use randomizeRegions for non-autocorrelated query
               evaluate.function = numOverlaps,
               genome = gam$genome, mask = gam$mask,
               ntimes = N_TIMES, count.once = TRUE)
pt$numOverlaps$pval; pt$numOverlaps$zscore

lz <- localZScore(A = peaks, B = enhancers, pt = pt, window = 10000, step = 500)   # sharp vs diffuse positional association
```

`circularRandomizeRegions` rotates the whole set around the genome, preserving inter-region spacing/clustering - the honest null when the query is autocorrelated (CpG islands, TAD-restricted peaks); plain `randomizeRegions` breaks clustering and inflates significance. `resampleRegions` draws from a supplied real universe. The `mask` is the workspace control.

## LOLA - Universe-Relative Region-Set Enrichment

```r
# Reference: LOLA 1.30+ | Verify API if version differs
library(LOLA)
regionDB <- loadRegionDB('LOLACore/hg38')
userSets <- readBed('peaks.bed')
userUniverse <- readBed('all_called_regions.bed')   # the candidate pool the query was drawn from - NOT the whole genome
res <- runLOLA(userSets, userUniverse, regionDB, cores = 4)   # odds ratio + p + q per reference set, rankable
```

LOLA is not a shuffle: it tests the query against many reference region sets (ENCODE TFBS, Roadmap chromatin) with a Fisher's exact test **relative to `userUniverse`**. The universe is everything - the LOLA vignette recommends either the *union of all query sets across the experiment* ("regions that were in play") or the assay's full candidate pool (e.g. all tested DHS/called peaks); pick the one that honestly bounds where the query could have come from. A whole-genome universe inflates every enrichment.

## rGREAT - Ontology Enrichment From Regions

```r
# Reference: rGREAT 2.4+ | Verify API if version differs
library(rGREAT)
res <- great(toGRanges('peaks.bed'), gene_sets = 'GO:BP', tss_source = 'txdb:hg38',
             background = toGRanges('accessible.bed'))   # supply a real background, not the whole-genome default
tb <- getEnrichmentTable(res)   # has Binom + Hyper p/adjp columns
```

GREAT assigns each gene a regulatory domain (basal 5 kb upstream / 1 kb downstream, extended up to 1 Mb to the next gene's basal domain), maps regions to those domains, then runs a **binomial-over-regions** test AND a **hypergeometric-over-genes** test - by design, with opposite biases. **Trust a term only if both fire:** a single gene with a huge regulatory domain attracts regions by target size and lights up the binomial; the hypergeometric (counting that gene once) calls the bluff. Supply a real background - the binomial assumes regions are independent and uniformly placeable, which clustered ChIP/ATAC peaks violate (Fulcher 2021 demonstrates orders-of-magnitude false-positive inflation from spatial autocorrelation in genomic enrichment analysis - a transferable critique of region-to-gene-category tests).

## Per-Method Failure Modes

### Whole-genome universe for a biased query
**Trigger:** testing ATAC/ChIP peaks (or DMRs, or capture-panel regions) against the whole genome. **Mechanism:** the query could only have come from accessible/assayed space, which is gene-rich; the genome universe credits that geography as enrichment. **Symptom:** every gene-associated annotation is "significant," none specific, p absurdly small. **Fix:** set the universe to the callable/candidate pool (LOLA `userUniverse`, GAT `--workspace`, regioneR `mask`/`resampleRegions`).

### bedtools fisher reported as the result
**Trigger:** putting a `bedtools fisher` p-value in a figure or reviewer reply. **Mechanism:** its analytic 2x2 estimates the unobserved cell from a uniform-placement heuristic, ignoring genome structure; prone to inflation. **Symptom:** spuriously tiny p that evaporates under a competent permutation null. **Fix:** use fisher only to triage; validate any low p with GAT/regioneR.

### Uniform null on clustered query regions
**Trigger:** `randomizeRegions` / plain `shuffle` on autocorrelated regions (CpG islands, tandem families, TAD-restricted peaks). **Mechanism:** uniform placement destroys the clustering the observed data has, so permuted overlap is too low. **Symptom:** inflated significance vs a structure-preserving null. **Fix:** `circularRandomizeRegions` (regioneR) or per-chromosome/per-class shuffling.

### No blacklist / no workspace exclusion
**Trigger:** shuffling across the full genome including gaps, centromeres, and the ENCODE blacklist. **Mechanism:** the null places regions where reads/peaks could never occur, lowering expected overlap. **Symptom:** enrichment that is really mappability artifact. **Fix:** exclude the ENCODE blacklist (Amemiya 2019) and assembly gaps via `-excl`/`mask` before any test.

### Ignoring multiple testing across tracks/terms
**Trigger:** testing a query against many annotation tracks or many GO terms and reading raw p-values. **Mechanism:** dozens-to-thousands of tests inflate the family-wise false-positive rate. **Symptom:** a long list of "significant" hits dominated by chance. **Fix:** use GAT's built-in FDR, FDR-adjust LOLA ranks, and read GREAT's binomial+hypergeometric q-values; require both GREAT tests.

### GREAT term significant by only one test
**Trigger:** reporting a GO term significant by binomial OR hypergeometric alone. **Mechanism:** the two tests have opposite biases (domain size vs gene count). **Symptom:** a term driven by one large-domain gene, or by gene-counting alone. **Fix:** require significance by both; distrust single-test hits.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Permutation N >= 1000 | empirical p resolution (Gel 2016; Phipson & Smyth 2010) | stable empirical p down to ~0.001; the `(hits+1)/(N+1)` estimator avoids p=0 |
| GAT --num-samples >= 10000 | GAT practice | needed for stable small q across many annotations |
| 50% reciprocal overlap (`-f 0.5 -r`) for CNV/SV concordance | field convention | "same event"; one-sided fractions let a giant call swallow a tiny one |
| Universe = candidate/callable pool, not the genome | LOLA/GAT/regioneR design | the dominant lever; whole-genome universe inflates every enrichment |
| GREAT basal 5 kb up / 1 kb down, extend to 1 Mb | McLean 2010 default | regulatory-domain model; user-configurable, report the values used |
| ENCODE blacklist excluded before any test | Amemiya 2019 | high-signal artifact regions otherwise manufacture overlap |
| FDR control across tracks/terms | multiple-testing | many-track / many-term tests inflate false positives |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Every annotation "significant," p absurdly small | whole-genome universe for a biased query | set the universe to the callable/accessible pool |
| `bedtools fisher` error | inputs not sorted, or missing genome file | `sort -k1,1 -k2,2n`; pass `-g genome.txt` |
| Empty/zero overlap in shuffle null | chrom naming mismatch (`chr1` vs `1`) across query/genome/blacklist | harmonize chromosome naming across all files |
| permTest much more significant than expected | uniform randomizer on clustered regions | use `circularRandomizeRegions` |
| GAT reports no GC effect | no `--isochores` supplied | pass an isochore/GC segmentation of the workspace |
| GREAT term looks real but is fragile | significant by only one of the two tests | require binomial AND hypergeometric |
| regioneR `toGRanges`/`getGenomeAndMask` error | genome name not recognized / mask chrom mismatch | use a supported genome id or supply explicit GRanges; match chrom naming |

## References

- Quinlan AR, Hall IM. 2010. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26:841-842.
- Heger A, Webber C, Goodson M, Ponting CP, Lunter G. 2013. GAT: a simulation framework for testing the association of genomic intervals. *Bioinformatics* 29:2046-2048.
- Gel B, Diez-Villanueva A, Serra E, Buschbeck M, Peinado MA, Malinverni R. 2016. regioneR: an R/Bioconductor package for the association analysis of genomic regions based on permutation tests. *Bioinformatics* 32:289-291.
- Sheffield NC, Bock C. 2016. LOLA: enrichment analysis for genomic region sets and regulatory elements in R and Bioconductor. *Bioinformatics* 32:587-589.
- McLean CY, Bristor D, Hiller M, Clarke SL, Schaar BT, Lowe CB, Wenger AM, Bejerano G. 2010. GREAT improves functional interpretation of cis-regulatory regions. *Nat Biotechnol* 28:495-501.
- Gu Z, Hubschmann D. 2023. rGREAT: an R/Bioconductor package for functional enrichment on genomic regions. *Bioinformatics* 39:btac745.
- Amemiya HM, Kundaje A, Boyle AP. 2019. The ENCODE blacklist: identification of problematic regions of the genome. *Sci Rep* 9:9354.
- Fulcher BD, Arnatkeviciute A, Fornito A. 2021. Overcoming false-positive gene-category enrichment in the analysis of spatially resolved transcriptomic brain atlas data. *Nat Commun* 12:2669.
- Phipson B, Smyth GK. 2010. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Stat Appl Genet Mol Biol* 9:Article 39.

## Related Skills

- interval-arithmetic - The intersect/shuffle/jaccard/fisher mechanics this skill turns into a test
- bed-file-basics - BED format, coordinate systems, and sorting the inputs every test requires
- proximity-operations - Nearest-feature assignment when the question is distance, not overlap enrichment
- chip-seq/peak-calling - Source of the peak query sets tested for enrichment
- chip-seq/peak-annotation - Assign enriched peaks to genes/features
- atac-seq/atac-peak-calling - Source of ATAC peak sets and the accessible-region universe
- pathway-analysis/go-enrichment - Gene-list ontology enrichment; GREAT is the region-based analog
- causal-genomics/colocalization-analysis - GWAS/eQTL statistical colocalization (shared causal variant) - a distinct problem from interval overlap
- data-visualization/genome-tracks - Render the query and annotation tracks behind an enrichment claim
