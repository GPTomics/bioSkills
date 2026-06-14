---
name: bio-population-genetics-population-structure
description: Infers and describes population structure with PCA (plink2 --pca, smartpca/EIGENSOFT, FlashPCA2), model-based clustering (ADMIXTURE, fastSTRUCTURE), FST estimators (Weir-Cockerham vs Hudson), and f-statistics (f3/f4/D via AdmixTools/admixr), plus Python plotting of PCs and Q barplots. Every output is a model-conditioned description of variance, not truth: PCs conflate ancestry with LD/inversions/relatedness/batch, ADMIXTURE Q-values are panel- and K-dependent artifacts, and CV-minimum K is a guide not the true population count. FST must combine SNPs as a ratio of averages (sum numerators / sum denominators), never an average of per-SNP FST; negative per-SNP FST is normal and must not be clamped. f3/f4/D need a block jackknife or the significance is fake. Use when running PCA, ADMIXTURE, FST, or f-statistics on QC'd genotypes. For QC and KING relatedness see plink-basics; for LD pruning see linkage-disequilibrium; for array-based Python pipelines see scikit-allel-analysis.
tool_type: mixed
primary_tool: plink2
---

## Version Compatibility

Reference examples tested with: PLINK 2.0 (alpha 6+), ADMIXTURE 1.3+, EIGENSOFT 7.2+, scikit-allel 1.3+, numpy 1.26+, pandas 2.2+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Version traps that change results, not just syntax: ADMIXTURE is a standalone CLI (`admixture --cv input.bed K`), never an R or Python package, and `--cv` defaults to 5-fold. plink2 `--pca` builds on the variance-standardized relationship matrix from `--make-rel`/`--make-grm` and has no Tracy-Widom test; smartpca does. plink2 `--make-king` gives kinship (cutoff 0.0884 = second-degree), not the deprecated PI_HAT from PLINK 1.9 `--genome`. f-statistics live in AdmixTools (CLI) wrapped by admixr (R), not in plink. The single source of truth for versions is this block, not headings.

# Population Structure

**"Analyze the population structure in my genotypes"** -> Project genotype covariance into continuous axes or discrete clusters, after pruning the artifacts the model would otherwise mistake for ancestry, and attach the uncertainty machinery any interpreted statistic requires.
- CLI: `plink2 --pca 20 approx` (top eigenvectors of the relationship matrix; LD-pruned, relatives removed first)
- CLI: `admixture --cv data_pruned.bed 3` (maximum-likelihood mixing weights for K abstract clusters)
- Python: `allel.hudson_fst(ac1, ac2)` then `num.sum()/den.sum()` (FST as a ratio of averages)

Scope: PCA, model-based clustering (ADMIXTURE/fastSTRUCTURE), FST estimators, and f-/D-statistics, with Python plotting. QC and KING relatedness route to plink-basics; LD pruning to linkage-disequilibrium; array-scale Python diversity/FST windows to scikit-allel-analysis; phased haplotype work to phasing-imputation/haplotype-phasing; introgression detection to comparative-genomics/introgression-detection.

## The Single Most Important Insight -- every structure method returns a model-conditioned description of variance, not truth

1. The output is a deterministic function of three silent choices: which samples are in the panel, which SNPs survive ascertainment/QC/LD-pruning, and which model is imposed (continuous PCs vs K discrete clusters vs a tree-with-admixture); change any one and the "answer" changes.
2. PCA does not find ancestry, it finds the directions of greatest genotype covariance, which conflate ancestry with LD blocks, inversions, relatedness, batch, and differential missingness, so the mandatory work is pruning those out before reading axes.
3. ADMIXTURE Q-values are not ancestry fractions but maximum-likelihood weights on K abstract allele-frequency vectors that exist only because K of them were requested, and the CV-minimum K is a prediction-accuracy guide, not the true number of populations.
4. Any interpreted statistic needs its own uncertainty: FST combines across SNPs as a ratio of averages (never an average of per-SNP FST), negative per-SNP FST is kept not clamped, and every f3/f4/D needs a block-jackknife standard error or the significance is fabricated.

## Tool Taxonomy

| Method | Citation | Mechanism / role | When |
|--------|----------|------------------|------|
| plink2 `--pca` | Patterson 2006; Price 2006 | Eigenvectors of the variance-standardized relationship matrix; `approx` for large N | Stratification covariates, gross structure, QC outliers |
| smartpca (EIGENSOFT) | Patterson 2006 | PCA plus Tracy-Widom significance, outlier removal, `lsqproject` projection | Rigorous PCA, aDNA projection, per-PC p-values |
| FlashPCA2 | Abraham 2017 | Randomized PCA for biobank N (>100k) | Very large cohorts |
| ADMIXTURE | Alexander 2009 | Fast ML point estimate of Q (ancestry weights) and P (cluster frequencies) | Genome-wide discrete ancestry proportions |
| fastSTRUCTURE | Raj 2014 | Variational Bayes clustering with `chooseK.py` K guidance | Fast K exploration |
| Hudson FST | Hudson 1992; Bhatia 2013 | Per-SNP heterozygosity estimator; ratio of averages | Pairwise differentiation, unequal sample sizes, rare variants / SNP-array ascertainment |
| Weir-Cockerham FST | Weir & Cockerham 1984 | ANOVA estimator (a/(a+b+c)); ratio of averages | Classic variance-partition framing, balanced n |
| f3 / f4 / D | Patterson 2012; Durand 2011 | Drift-distance tests of admixture and tree-ness; block jackknife | Admixture detection, gene-flow tests |
| TreeMix | Pickrell 2012 | ML tree plus migration edges from frequency covariance | Tree + migration hypotheses |
| admixr / AdmixTools | Petr 2019; Patterson 2012 | Reproducible R wrappers for qp3Pop/qpDstat/qpAdm/qpGraph | f-statistics and admixture graphs |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Stratification covariates for GWAS | plink2 `--pca` (`approx` >5000) | assumption-light, fast, feeds `--glm` directly |
| Per-PC significance, outlier removal, aDNA projection | smartpca | Tracy-Widom test plus `lsqproject` for shrinkage-robust projection |
| Biobank-scale PCA (>100k) | plink2 `--pca approx` or FlashPCA2 | randomized algorithms scale; exact PCA does not |
| Discrete ancestry proportions | ADMIXTURE over a K span (`--cv` as guide) | present a span, never the CV argmin alone |
| Fast K exploration | fastSTRUCTURE + `chooseK.py` | variational, very fast |
| Pairwise differentiation, unequal n or rare variants | Hudson FST, ratio of averages | Bhatia 2013; WC is sensitive to n, population count, and rare variants |
| "Is population C admixed?" | f3(C; A,B) (qp3Pop / admixr) | f3 < 0 proves admixture; block-jackknife Z |
| "Is there gene flow / introgression?" | D / f4 (qpDstat, ABBA-BABA) | `|Z|` > 3; cannot separate from ancient structure alone |
| Clinal / spatial structure | EEMS, Mantel | clines are isolation-by-distance, not discrete demes |
| Relatedness QC before everything | plink2 `--king-cutoff` (see plink-basics) | a relative cluster grabs a spurious PC |

## Mandatory Preprocessing Before PCA

**Goal:** Compute PCs that track ancestry rather than inversions, relatedness, or batch.

**Approach:** LD-prune, exclude long-range-LD/inversion regions by coordinate, remove relatives before computing axes, drop very-low-MAF variants, then run PCA on the survivors.

```bash
# LD-prune so a single dense block cannot dominate a PC (route detail to linkage-disequilibrium).
plink2 --bfile data --indep-pairwise 50 5 0.1 --out prune

# Exclude long-range-LD regions and inversions that survive pruning and create karyotype PCs.
# range_lrld.txt holds MHC chr6:25-35 Mb, 8p23.1, 17q21.31, and LCT/2q21 in plink --exclude range format.
plink2 --bfile data --extract prune.prune.in --exclude range range_lrld.txt --maf 0.01 \
    --make-bed --out data_for_pca

# PCA on the LD-pruned, inversion-stripped, relatedness-pruned set. approx is near-required above ~50k.
plink2 --bfile data_for_pca --pca 20 approx --out pca
# Outputs: pca.eigenvec (FID IID PC1..PCn), pca.eigenval (variance per PC).
```

Relatives must be removed BEFORE computing axes (a cluster of cousins forms its own high-covariance PC); compute the KING cutoff with `plink2 --king-cutoff 0.0884` from plink-basics, build PCs on the unrelated set, then project relatives back. plink2 has no Tracy-Widom test: feed the eigenvalues to smartpca `twstats` or read a scree elbow to decide which PCs are real.

## Rigorous PCA and Projection (smartpca)

**Goal:** Attach per-PC significance and project new/ancient samples without shrinkage artifacts.

**Approach:** Run smartpca with outlier iterations and Tracy-Widom output for the reference build, and use `lsqproject: YES` to place additional samples robustly to missingness.

```bash
# smartpca parameter file (key params verified against the EIGENSOFT POPGEN README):
#   numoutevec: 20            # PCs to output (default 10)
#   numoutlieriter: 5         # outlier-removal iterations (default 5; 0 disables)
#   outliersigmathresh: 6.0   # SD threshold for outlier removal (default 6.0)
#   lsqproject: YES           # least-squares projection, robust to missing data (aDNA standard)
#   poplistname: ref_pops.txt # which populations build the axes (others are projected)
#   altnormstyle: NO          # NO = Price 2006 EIGENSTRAT normalization; YES = Patterson 2006
smartpca -p smartpca.par
# Tracy-Widom test on the eigenvalues: only PCs with p < ~0.05 plus a scree elbow are interpretable.
twstats -t twtable -i out.eval -o out.tw
```

Projected scores shrink toward the origin, worse with a small reference panel (<~5000) and more missing data, so an ancient sample plotting "between" two clusters may be shrunk, not admixed; `lsqproject` is the standard fix. plink2 projects via `--pca allele-wts` then `--score` on the `.eigenvec.allele` weights, but does not correct shrinkage.

## ADMIXTURE Over a K Span

**Goal:** Estimate discrete ancestry proportions while treating K as a model-selection choice, not a discovery.

**Approach:** Run ADMIXTURE on LD-pruned data across a span of K with cross-validation, plot CV error as a guide, and check Q stability across seeds before interpreting any single K.

```bash
# admixture is a standalone CLI; --cv defaults to 5-fold. -jN threads, -BN bootstrap SEs.
for K in $(seq 2 8); do
    admixture --cv -j4 data_pruned.bed "$K" 2>&1 | tee "log_K${K}.out"
done
# Outputs per K: data_pruned.K.Q (N x K ancestry weights) and data_pruned.K.P (cluster frequencies).
# CV error prints as: CV error (K=3): 0.512 -- a guide, never "the true number of populations".
grep -h "CV error" log_K*.out
```

Supervised mode (`admixture --supervised data_pruned.bed K`, reading `data_pruned.pop` with one label per individual, blank for unknowns) fixes labeled individuals to their population but assumes the reference populations are themselves unadmixed. Replicate Q at the same K can land in different local optima; align cluster labels across runs with CLUMPP or pong before averaging or plotting, and treat unstable Q as a sign of mis-specified K, not noise to smooth away.

## FST as a Ratio of Averages

**Goal:** Estimate pairwise differentiation without the average-of-ratios bias.

**Approach:** Compute per-SNP Hudson numerators and denominators, keep negative numerators, then divide summed numerators by summed denominators across SNPs.

```python
import allel
import numpy as np

# ac1, ac2 are AlleleCountsArrays for the two populations at the same SNPs (see scikit-allel-analysis).
num, den = allel.hudson_fst(ac1, ac2)   # per-SNP Hudson numerator and denominator (Bhatia 2013)
fst = num.sum() / den.sum()             # RATIO OF AVERAGES across SNPs; never np.mean(num/den)
# Negative per-SNP numerators are normal sampling behavior near FST=0 and stay in the sum.
# Use Weir-Cockerham only with balanced n: allel.weir_cockerham_fst returns a, b, c variance components.
a, b, c = allel.weir_cockerham_fst(genotype_array, subpops)
fst_wc = np.nansum(a) / np.nansum(a + b + c)   # still a ratio of averages
```

The Hudson estimator is preferred under sample-size asymmetry because Weir-Cockerham's finite-sample correction makes it sensitive to n and to the number of populations; the two can disagree enough to cross a Wright differentiation band. SNP-array ascertainment compresses and warps FST relative to whole-genome sequencing, so cross-study comparisons require matched ascertainment.

## f-statistics with a Block Jackknife (admixr)

**Goal:** Test admixture and gene flow with honest standard errors.

**Approach:** Run f3/f4/D through AdmixTools (wrapped by admixr) so each statistic carries a block-jackknife SE, and read Z, never a raw point estimate.

```r
library(admixr)
# admixr wraps AdmixTools (Petr 2019); each call returns the statistic with a block-jackknife SE and Z.
data <- eigenstrat('prefix')
res_f3 <- f3(A = 'PopA', B = 'PopB', C = 'PopC', data = data)   # f3(C; A,B) < 0 with Z < -3 proves C admixed
res_d  <- d(W = 'PopW', X = 'PopX', Y = 'PopY', Z = 'PopZ', data = data)  # |Z| > 3 indicates gene flow
```

A significantly negative f3 proves the target is admixed (no tree produces a negative f3), but a non-negative f3 is inconclusive, not proof of a clean tree: post-admixture drift in the target itself (a bottleneck after the admixture event), or heavily drifted sources, adds a positive term that can mask the negative cross-product even when admixture is real. A nonzero D or f4 is evidence of a tree violation, not specifically recent introgression: symmetric ancient structure mimics the same ABBA/BABA asymmetry (Durand 2011), so separating them needs admixture-LD decay or explicit modeling.

## Per-Method Failure Modes

### PCA tracks an inversion, not a deme
**Trigger:** PCA on LD-pruned data with MHC/8p23/17q21.31/LCT regions still in. **Mechanism:** megabase-long LD in inversions survives `--indep-pairwise` and dominates a PC. **Symptom:** a PC loads almost entirely on one chromosome arm and splits samples by karyotype. **Fix:** `--exclude range` the long-range-LD and inversion regions by coordinate before PCA.

### Relatives grab a principal component
**Trigger:** computing PCs before relatedness pruning. **Mechanism:** a cluster of relatives forms a high-covariance bundle. **Symptom:** a tight outlier cluster on a top PC that is not a real population. **Fix:** `--king-cutoff 0.0884` first, build PCs on the unrelated set, project relatives back.

### Projection shrinkage misread as admixture
**Trigger:** projecting new/ancient/low-coverage samples onto reference PCs naively. **Mechanism:** projected scores are biased toward the origin, worse with small panels and missing data. **Symptom:** a sample plots "between" two clusters and is narrated as admixed. **Fix:** smartpca `lsqproject: YES`; do not interpret shrunk scores as intermediacy.

### CV-minimum K over-splits
**Trigger:** picking K at the CV argmin when the curve plateaus or keeps falling. **Mechanism:** CV error is prediction accuracy, not a population count. **Symptom:** uninterpretable extra clusters at high K. **Fix:** run K across a span with at least 10 seeds each (`-s`), align replicates with pong/CLUMPP, and report the K where CV error plateaus AND Q is seed-stable; if CV and stability disagree, present both and let sampling design plus orthogonal evidence pick the interpreted K.

### Label switching corrupts averaged barplots
**Trigger:** averaging Q across runs or seeds without alignment. **Mechanism:** cluster labels permute arbitrarily between runs. **Symptom:** a smeared, meaningless mean barplot. **Fix:** align labels with CLUMPP or pong before averaging; treat unstable Q as a mis-specification warning.

### Average-of-ratios FST
**Trigger:** combining per-SNP FST as `mean(num/den)`. **Mechanism:** low-MAF SNPs have tiny denominators and dominate the mean. **Symptom:** badly biased genome-wide FST. **Fix:** ratio of averages, `num.sum()/den.sum()`; never clamp negative per-SNP values first.

### f-statistic significance without a jackknife
**Trigger:** naive SNP-level standard errors for f3/f4/D. **Mechanism:** LD correlates neighboring SNPs, so per-SNP SEs are far too small. **Symptom:** everything looks significant. **Fix:** block-jackknife SE (drop ~5 cM blocks); report Z with `|Z|` > 3.

### Clinal sample forced into discrete clusters
**Trigger:** running K-cluster ADMIXTURE on an isolation-by-distance continuum. **Mechanism:** smooth clines have no discrete demes. **Symptom:** phantom populations and "admixed" intermediates that are really IBD. **Fix:** describe clinal structure with EEMS / Mantel, not a STRUCTURE barplot.

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| LD pruning for PCA/ADMIXTURE | `--indep-pairwise 50 5 0.1` (range 0.05-0.2) | near-independence so structure is not double-counted (see linkage-disequilibrium) |
| MAF floor before PCA | drop MAF < ~0.01 (often < 0.05) | plink2 docs: very-low-MAF variants destabilize PCA |
| Relatedness removal | KING `--king-cutoff 0.0884` (2nd-degree) | KING boundaries 0.354/0.177/0.0884/0.0442 = MZ/1st/2nd/3rd |
| PC significance | Tracy-Widom p < 0.05 plus scree elbow | Patterson 2006; TW null for the largest eigenvalue |
| smartpca outlier removal | `outliersigmathresh 6.0`, `numoutlieriter 5` | EIGENSOFT defaults |
| ADMIXTURE CV | default 5-fold; `--cv=10` for stability | Alexander/Shringarpure manual |
| f3/f4/D significance | `|Z|` > 3 (block jackknife) | Patterson 2012 convention (~3 SE) |
| FST bands (Wright guideposts) | 0-0.05 little, 0.05-0.15 moderate, 0.15-0.25 great, >0.25 very great | heuristic only; estimator-, MAF-, ascertainment-dependent |
| Negative per-SNP FST | keep (do not clamp) | unbiased estimators yield negatives near FST=0 by sampling |

Thresholds are conventions, not laws; the FST bands predate SNP arrays, and the estimator, MAF spectrum, ascertainment, and sample size each move FST by a whole band. Verify current best practice before applying numbers blindly.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `import admixture` fails | ADMIXTURE is a CLI, not a package | run `admixture --cv data.bed K` on the shell |
| CV argmin reported as the true K | reading CV error as a population count | present a K span; CV is a guide, check Q stability |
| A PC tracks one chromosome arm | inversion/long-range-LD region left in | `--exclude range` MHC/8p23/17q21.31/LCT before PCA |
| Outlier cluster is "a new population" | relatives not removed before PCA | `--king-cutoff 0.0884` first, project relatives back |
| Genome-wide FST is biased high | average-of-ratios and/or clamped negatives | ratio of averages; keep negative per-SNP values |
| WC and Hudson FST disagree | unequal sample sizes between populations | prefer Hudson under sample-size asymmetry (Bhatia 2013) |
| Everything is f3/D-significant | naive (non-jackknife) standard errors | use block-jackknife SE; report Z, `|Z|` > 3 |
| Non-negative f3 read as "no admixture" | post-admixture drift in the target (or drifted sources) masks the negative term | non-negative f3 is inconclusive, not a clean tree |
| Smeared mean Q barplot | label switching across replicates | align with CLUMPP/pong before averaging |

## References

1. Patterson N, Price AL, Reich D. Population structure and eigenanalysis. PLoS Genetics 2006; 2(12):e190. DOI:10.1371/journal.pgen.0020190.
2. Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D. Principal components analysis corrects for stratification in genome-wide association studies. Nature Genetics 2006; 38(8):904-909. DOI:10.1038/ng1847.
3. Alexander DH, Novembre J, Lange K. Fast model-based estimation of ancestry in unrelated individuals. Genome Research 2009; 19(9):1655-1664. DOI:10.1101/gr.094052.109.
4. Raj A, Stephens M, Pritchard JK. fastSTRUCTURE: variational inference of population structure in large SNP data sets. Genetics 2014; 197(2):573-589. DOI:10.1534/genetics.114.164350.
5. Weir BS, Cockerham CC. Estimating F-statistics for the analysis of population structure. Evolution 1984; 38(6):1358-1370. DOI:10.1111/j.1558-5646.1984.tb05657.x.
6. Hudson RR, Slatkin M, Maddison WP. Estimation of levels of gene flow from DNA sequence data. Genetics 1992; 132(2):583-589. DOI:10.1093/genetics/132.2.583.
7. Bhatia G, Patterson N, Sankararaman S, Price AL. Estimating and interpreting FST: the impact of rare variants. Genome Research 2013; 23(9):1514-1521. DOI:10.1101/gr.154831.113.
8. Patterson N, Moorjani P, Luo Y, Mallick S, Rohland N, Zhan Y, Genschoreck T, Webster T, Reich D. Ancient admixture in human history. Genetics 2012; 192(3):1065-1093. DOI:10.1534/genetics.112.145037.
9. Durand EY, Patterson N, Reich D, Slatkin M. Testing for ancient admixture between closely related populations. Molecular Biology and Evolution 2011; 28(8):2239-2252. DOI:10.1093/molbev/msr048.
10. Pickrell JK, Pritchard JK. Inference of population splits and mixtures from genome-wide allele frequency data. PLoS Genetics 2012; 8(11):e1002967. DOI:10.1371/journal.pgen.1002967.
11. Petr M, Vernot B, Kelso J. admixr - R package for reproducible analyses using ADMIXTOOLS. Bioinformatics 2019; 35(17):3194-3195. DOI:10.1093/bioinformatics/btz030.
12. Abraham G, Qiu Y, Inouye M. FlashPCA2: principal component analysis of Biobank-scale genotype datasets. Bioinformatics 2017; 33(17):2776-2778. DOI:10.1093/bioinformatics/btx299.

## Related Skills

- plink-basics - QC, KING relatedness pruning, and fileset preparation before structure analysis
- linkage-disequilibrium - LD pruning the SNP set that PCA and ADMIXTURE require
- scikit-allel-analysis - array-scale FST and diversity windows in Python
- phasing-imputation/haplotype-phasing - phased haplotypes for haplotype-based structure methods
- comparative-genomics/introgression-detection - D/f4 introgression scans beyond pairwise tests
