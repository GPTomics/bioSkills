---
name: bio-methylation-methylkit
description: Imports Bismark coverage or cytosine-report files into the methylKit object model, then runs the import-to-results spine - filterByCoverage, normalizeCoverage, unite/destrand, calculateDiffMeth, getMethylDiff - for both per-CpG (DMC) and fixed-tile (DMR) differential methylation, plus tileMethylCounts, PCA/correlation/clustering QC, and assocComp/removeComp batch handling. Covers the silent default traps that shape the false-positive rate: overdispersion='none' does no correction while 'MN' forces the F-test (ignoring test='Chisq'), adjust defaults to SLIM not BH, getMethylDiff defaults difference=25/qvalue=0.01, cov.bases=0 admits single-CpG tiles, and pool destroys biological replication. Use when importing bisulfite count tables, filtering/normalizing/uniting methylation samples, running methylKit differential testing, or QC-ing methylomes. For per-site test-choice (count vs continuous) see differential-cpg-testing; for selection-aware region FDR (dmrseq/DSS) see dmr-detection.
tool_type: r
primary_tool: methylKit
---

## Version Compatibility

Reference examples tested with: methylKit 1.28+, GenomicRanges 1.54+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The `assembly=` string (e.g. `hg38`) is metadata only - methylKit never checks it. The genome build is real elsewhere: coordinates must match the alignment genome, and annotation packages (`TxDb.Hsapiens.UCSC.hg38.knownGene`, annotatr `build_annotations(genome='hg38')`) are genome-build-specific. methylKit's `overdispersion`, `test`, and `adjust` defaults have shifted across Bioconductor releases - run `?calculateDiffMeth` on the installed build before trusting any default.

# methylKit Analysis

**"Analyze methylation across my samples"** -> Import per-cytosine counts into a methylRawList, then filter, normalize, unite, and test - because each of those steps is a modeling decision that sets which CpGs survive and how many false positives the test emits, not boilerplate.
- R: `methRead(pipeline='bismarkCoverage')` -> `filterByCoverage()` -> `normalizeCoverage()` -> `unite(destrand=)` -> `calculateDiffMeth(overdispersion='MN')` -> `getMethylDiff()`

Scope: the methylKit OBJECT MODEL and the short-read bisulfite import-to-results workflow, for BOTH per-CpG (DMC) and fixed-tile (DMR) results. Which per-site test to use (count vs continuous, beta vs M) -> differential-cpg-testing. Selection-aware region callers (dmrseq/DSS/metilene) and region FDR -> dmr-detection. Long-read MM/ML modBAM input -> long-read-sequencing/nanopore-methylation (its counts pipe back into this object model).

## The Single Most Important Modern Insight -- The Object Model Is the Analysis

Coverage filtering, normalization, destranding, and the overdispersion model are not setup before the "real" test - they ARE the test. Each silently changes which CpGs exist and what the p-value means, and methylKit's defaults are tuned for nothing in particular. Three corollaries every misuse violates:

1. **The defaults offer no protection.** `calculateDiffMeth` defaults to `overdispersion='none'` - a plain logistic LRT that assumes binomial-only variance and over-calls under biological replication. Two healthy replicates differ at a CpG far more than coin-flip sampling predicts; only `overdispersion='MN'` adds the between-replicate (beta-binomial) layer. The paper discusses overdispersion; the function does not apply it unless told.
2. **The knobs interact and several are silent.** `overdispersion='MN'` automatically switches to the F-test, so a passed `test='Chisq'` is ignored with no warning. `adjust` defaults to SLIM (methylKit's own q-method), not BH, so counts are not comparable to a DSS/limma BH analysis. `tileMethylCounts(cov.bases=0)` lets a one-CpG window become a "region." None of these throw an error; the result just quietly changes.
3. **Coverage is the substrate, not the answer.** A single CpG's methylation percentage is a count ratio; below ~10x it is a coin flip. Filtering the low tail (noise) and the high tail (PCR/repeat artifacts) before testing decides the result more than the test does.

Organize the workflow around defending these, not around calling functions in order.

## The Import-to-Results Spine

Run these in order. Skipping or reordering them changes the result silently.

### 1. Import: methRead with the pipeline matching the input

**Goal:** Load per-cytosine counts into a methylRawList, choosing the parser that matches the Bismark output format.

**Approach:** `pipeline='bismarkCoverage'` reads `.cov`/`.cov.gz` (chr/start/end/%meth/numC/numT - NO strand, so destranding is limited); `pipeline='bismarkCytosineReport'` reads the CX/CpG report (carries strand + context, enables proper destranding). `treatment` is an integer vector (0/1, or 0/1/2 for multi-group). `context='CpG'` only - never destrand CHG/CHH downstream.

```r
library(methylKit)
file_list <- list('ctrl1.cov.gz', 'ctrl2.cov.gz', 'treat1.cov.gz', 'treat2.cov.gz')
sample_ids <- list('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2')
meth_obj <- methRead(file_list, sample.id=sample_ids, treatment=c(0,0,1,1),
                     assembly='hg38', context='CpG', pipeline='bismarkCoverage')
# dbtype='tabix', save.db=TRUE gives disk-backed methylRawDB objects for large WGBS
```

### 2. Filter and normalize BEFORE uniting (and before tiling)

**Goal:** Drop unreliable and artifactual CpGs per sample, then remove library-size-driven coverage differences so a deeper sample does not look more "confident."

**Approach:** `filterByCoverage(lo.count, hi.perc)` removes the noisy low tail and the artifactual high tail; `normalizeCoverage` scales coverage between samples. Both are per-sample and must precede `unite` and `tileMethylCounts`.

```r
meth_filt <- filterByCoverage(meth_obj, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
meth_norm <- normalizeCoverage(meth_filt, method='median')
```

### 3. Unite: destrand only for CpG, only with strand info

**Goal:** Build the per-base table of CpGs covered across samples for testing.

**Approach:** `unite` keeps CpGs covered in ALL samples; `min.per.group=2L` relaxes that to >=2 per group (keeps more sites, allows missingness). `destrand=TRUE` merges the + and - strand counts of a CpG dyad - valid ONLY for symmetric CpG context AND only meaningful when strand is present (cytosine report). On `.cov` (bismarkCoverage, no strand) destranding is limited; on CHG/CHH it is wrong.

```r
meth_united <- unite(meth_norm, destrand=TRUE)            # destrand only if strand info present
meth_united <- unite(meth_norm, min.per.group=2L)         # allow missingness across replicates
```

### 4. QC the united object before testing

**Goal:** Confirm samples cluster by biology, not by batch, before believing any DMC.

**Approach:** Run correlation, PCA, and clustering on the united (% methylation) object. A control clustering with the treated group, or PC1 tracking sequencing batch, means the contrast is confounded.

```r
getCorrelation(meth_united, plot=TRUE)
PCASamples(meth_united)
clusterSamples(meth_united, dist='correlation', method='ward.D', plot=TRUE)
```

### 5. Test: calculateDiffMeth with overdispersion correction

**Goal:** Test each CpG for a group difference using a model that accounts for between-replicate overdispersion.

**Approach:** With replicates, set `overdispersion='MN'`, which automatically uses the F-test (the `test=` argument is then ignored - passing `test='Chisq'` alongside MN does not produce a chi-square test). Set `adjust='BH'` if the q-values must be comparable to other tools; the default SLIM is methylKit-specific. `getMethylDiff` filters by effect size AND q.

```r
diff_meth <- calculateDiffMeth(meth_united, overdispersion='MN', adjust='BH', mc.cores=4)
dmcs <- getMethylDiff(diff_meth, difference=25, qvalue=0.01)             # all DMCs
dmcs_hyper <- getMethylDiff(diff_meth, difference=25, qvalue=0.01, type='hyper')
# positive meth.diff = hyper in the higher-treatment group
```

## Tile-Based Regions (a fast screen, not selection-corrected inference)

**Goal:** Aggregate CpGs into fixed windows for a quick region-level scan.

**Approach:** Tile AFTER filter/normalize, raise `cov.bases` so a window needs real CpG support, then flow through the same unite -> calculateDiffMeth -> getMethylDiff path. The same `getMethylDiff` returns DMCs on a per-base object and (window) DMRs on a tiled object.

```r
tiles <- tileMethylCounts(meth_norm, win.size=1000, step.size=1000, cov.bases=3)  # cov.bases>=3
tiles_united <- unite(tiles, destrand=FALSE)              # tiles are not strand objects
diff_tiles <- calculateDiffMeth(tiles_united, overdispersion='MN', adjust='BH', mc.cores=4)
dmrs <- getMethylDiff(diff_tiles, difference=25, qvalue=0.01)
```

The per-tile q is a per-test SLIM/BH value: it does NOT model correlation between tiles and is NOT corrected for the region-selection step. Fixed windows also split or merge true DMRs at arbitrary boundaries. Treat methylKit tiles as a defensible screen; for rigorous region FDR (a permutation null that survives region selection) go to dmr-detection (dmrseq).

## Batch, Multi-Group, and the Single-Factor Limit

methylKit's `calculateDiffMeth` is a single-factor 2-group test. For known batch, remove the associated principal components before testing; for >2 groups, subset to pairwise contrasts. Complex designs (covariates, multi-factor) exceed what methylKit models - move to dmr-detection (DSS multiFactor / dmrseq covariates) or a continuous limma-on-M path (differential-cpg-testing).

```r
sample_anno <- data.frame(batch=c('a','a','b','b'))
as_comp <- assocComp(meth_united, sample_anno)            # which PCs track the covariate
meth_corrected <- removeComp(meth_united, comp=1)         # drop the batch PC, then test
meth_AB <- reorganize(meth_united, sample.ids=c('ctrl_1','ctrl_2','treat_1','treat_2'),
                      treatment=c(0,0,1,1))               # subset/relabel for a pairwise contrast
```

`pool(meth_united, sample.ids=...)` sums replicate counts into one pseudo-sample per group. This DESTROYS biological replication - the test then has no within-group variance estimate and its p-values are meaningless for inference. Use it only for no-replicate exploratory visualization, never for the reported test.

## Per-Method Failure Modes

### Overdispersion left at the default
**Trigger:** `calculateDiffMeth` with replicates and no `overdispersion=` argument. **Mechanism:** the default `'none'` is a binomial-only logistic LRT that ignores between-replicate variance. **Symptom:** implausibly many significant CpGs; q-values far smaller than a beta-binomial tool gives on the same data. **Fix:** `overdispersion='MN'` (which uses the F-test) whenever replicates exist.

### MN plus test='Chisq'
**Trigger:** passing both `overdispersion='MN'` and `test='Chisq'`. **Mechanism:** MN forces the F-test; `test=` is silently ignored. **Symptom:** the reported "chi-square test" was never run. **Fix:** drop `test=` when using MN; `test=` is honored only with `overdispersion='none'`.

### SLIM read as BH
**Trigger:** comparing methylKit q-value counts to a DSS/limma BH analysis. **Mechanism:** `adjust` defaults to SLIM, methylKit's own sliding-linear-model q-method. **Symptom:** DMC counts disagree with another tool's BH results at the "same" q. **Fix:** set `adjust='BH'` for cross-tool comparability and state the method used.

### cov.bases=0 tiles
**Trigger:** `tileMethylCounts` at the default `cov.bases=0`. **Mechanism:** a window with one covered CpG becomes a "region." **Symptom:** thousands of single-CpG "DMRs," many noisy. **Fix:** raise `cov.bases` to >=3.

### Tiling or testing the raw object
**Trigger:** `tileMethylCounts(meth_obj, ...)` or uniting before filtering/normalizing. **Mechanism:** low-coverage and library-size artifacts propagate into the tiles and the test. **Symptom:** artifactual regions; deeper samples look hyper-confident. **Fix:** `filterByCoverage` then `normalizeCoverage` BEFORE tiling/uniting.

### Destranding the wrong context or input
**Trigger:** `unite(destrand=TRUE)` on `.cov` input or on CHG/CHH context. **Mechanism:** `.cov` carries no strand; non-CpG dyads are not symmetric. **Symptom:** double-counting or wrong merges. **Fix:** destrand CpG only, ideally from the cytosine report.

### pool() then test
**Trigger:** `pool()` followed by `calculateDiffMeth`. **Mechanism:** pooling removes within-group variance. **Symptom:** tiny p-values with no biological meaning. **Fix:** never pool for the reported test; keep replicates separate.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `lo.count` = 10 | convention (methylKit tutorial) | below ~10x a single-CpG percentage is a coin flip; not a derived value |
| `hi.perc` = 99.9 | convention | drops the top 0.1% coverage (PCR/repeat artifacts) |
| `overdispersion` = 'MN' with replicates | Akalin 2012 *Genome Biol* 13:R87 | adds the beta-binomial between-replicate layer; default 'none' over-calls |
| `adjust` = 'BH' (default SLIM) | Akalin 2012 *Genome Biol* 13:R87 | BH for comparability; SLIM is methylKit-specific |
| `getMethylDiff` difference=25, qvalue=0.01 | methylKit defaults | 25% is tutorial convention, NOT derived; justify per feature/coverage/purity and report it |
| `tileMethylCounts` cov.bases >= 3 | nuance | default 0 admits single-CpG tiles; require real CpG support |
| `win.size`=`step.size`=1000 (default) | methylKit defaults | `step < win` gives overlapping (sliding) tiles |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Implausibly many DMCs | `overdispersion='none'` under replication | set `overdispersion='MN'` |
| Reported chi-square never ran | `test='Chisq'` with MN | MN forces F; drop `test=` |
| Counts disagree with another tool at same q | default `adjust='SLIM'` | set `adjust='BH'`, state the method |
| Thousands of single-CpG "regions" | `cov.bases=0` | raise `cov.bases` to >=3 |
| Destrand error / double counts | destrand on `.cov` or non-CpG | destrand CpG only, from cytosine report |
| Meaningless tiny p-values | `pool()` before testing | keep replicates; never pool for inference |
| Deeper sample looks more confident | no `normalizeCoverage` | normalize before unite/test |

## References

- Akalin A, Kormaksson M, Li S, Garrett-Bakelman FE, Figueroa ME, Melnick A, Mason CE. 2012. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. *Genome Biol* 13:R87.
- Krueger F, Andrews SR. 2011. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. *Bioinformatics* 27:1571-1572.
- Robinson MD, Kahraman A, Law CW, Lindsay H, Nowicka M, Weber LM, Zhou X. 2014. Statistical methods for detecting differentially methylated loci and regions. *Front Genet* 5:324.

## Related Skills

- methylation-calling - Produces the coverage/cytosine reports read here
- differential-cpg-testing - Per-site statistical model choice (count vs continuous)
- dmr-detection - Selection-aware region callers (dmrseq/DSS) beyond methylKit tiles
- pathway-analysis/go-enrichment - Functional annotation of differentially methylated genes
- long-read-sequencing/nanopore-methylation - Long-read MM/ML calling; pipe counts into this object model
- workflows/methylation-pipeline - End-to-end bisulfite pipeline
