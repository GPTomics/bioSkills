---
name: bio-population-genetics-scikit-allel-analysis
description: In-memory Python population genetics with scikit-allel - GenotypeArray/HaplotypeArray/AlleleCountsArray, diversity (pi, theta, Tajima's D), SFS, FST (Weir-Cockerham, Hudson, Patterson), f3/D admixture stats, LD pruning, PCA, and selection scans (iHS, XP-EHH, nSL, Garud H). Nearly every statistic is a ratio or density with one silent denominator bug in two faces: omit is_accessible= and per-base pi/theta divide by total span not accessible bp (deflated 2-5x); average per-SNP FST instead of sum(a)/(sum(a)+sum(b)+sum(c)) and the estimate is rare-variant-biased - scikit-allel returns the (a,b,c) and (num,den) components on purpose to force ratio-of-sums. to_n_alt default fill=0 imputes missing to reference; sfs() is unfolded and wants derived not alt counts; iHS/XP-EHH need phased data and standardization. Use when computing population-genetics statistics in Python, scanning for selection, or building array pipelines. For PLINK QC see plink-basics; for VCF input see variant-calling/vcf-basics.
tool_type: python
primary_tool: scikit-allel
---

## Version Compatibility

Reference examples tested with: scikit-allel 1.3.13+, numpy 1.26+, zarr 2.18+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Version traps that change results, not just syntax: scikit-allel is in MAINTENANCE mode (latest line v1.3.x, e.g. v1.3.13 Sep 2024); the README names sgkit (xarray+dask) as the successor but states it is "not yet at feature parity", so scikit-allel remains the pragmatic choice for established stat workflows. `average_patterson_f3`/`average_patterson_d` are the current names; pre-2020 code used `blockwise_patterson_*` (gone). `to_n_alt()` default is `fill=0`, not `fill=-1`. `read_vcf` is eager and loads the whole file into RAM. The single source of truth for versions is this block, not headings.

# scikit-allel Analysis

**"Analyze population genetics in Python"** -> Read a VCF into array structures, then compute frequency-, diversity-, differentiation-, and haplotype-based statistics with correct denominators.
- Python: `allel.GenotypeArray`, `allel.AlleleCountsArray`, `allel.windowed_diversity(..., is_accessible=)`, `allel.average_hudson_fst`, `allel.pca`

Scope: in-memory and dask/zarr scikit-allel analysis - the array-API mechanics of the data model, diversity/SFS, FST and f/D admixture statistics, LD pruning, PCA, and selection-scan computation. PLINK-format QC routes to plink-basics; PCA/ADMIXTURE via CLI tools to population-structure; selection-scan DESIGN (standardization, outlier calling, demographic confounding) to selection-statistics; phased input to phasing-imputation/haplotype-phasing; VCF generation to variant-calling/vcf-basics.

## The Single Most Important Insight -- one denominator bug, two faces

1. Almost every statistic scikit-allel computes is a RATIO or a DENSITY, and the two ways the denominator is silently wrong both return a finite, plausible, unflagged number.
2. Face (a) - the per-base denominator: pi, Watterson's theta, and Dxy divide a segregating-site sum by the number of bases; omit `is_accessible=` and they divide by total span (`stop-start+1`) instead of callable bp, deflating values 2-5x AND distorting the genome-wide landscape because callability varies per window.
3. Face (b) - FST is a ratio of variance components: the correct multi-locus estimate is `sum(a)/(sum(a)+sum(b)+sum(c))`, NOT `mean(per_snp_fst)` (which is rare-variant-dominated and biased, Bhatia 2013); scikit-allel returns the `(a,b,c)` components from `weir_cockerham_fst` and `(num,den)` from `hudson_fst`/`patterson_fst` precisely to force ratio-of-sums.
4. Both bugs run silently with no exception or warning, so correctness is a design decision the analyst makes, not something the library enforces.

## Tool Taxonomy -- the data model and the scaling path

| Object / path | Shape / form | Role | When |
|---------------|--------------|------|------|
| `GenotypeArray` | `(n_variants, n_samples, ploidy)` int8, -1 = missing | the fundamental call array | diploid genotypes from `calldata/GT` |
| `HaplotypeArray` | `(n_variants, n_haplotypes)` | phased chromosomes | iHS/XP-EHH/nSL/Garud H (REQUIRE phasing) |
| `AlleleCountsArray` | `(n_variants, n_alleles)` int32 | currency of all frequency stats | `gt.count_alleles()`; ignores -1 |
| `to_n_alt` 012 matrix | `(n_variants, n_samples)` | input to PCA/LD | `gt.to_n_alt(fill=...)` (default fill=0 imputes to ref) |
| in-memory numpy | dense, all in RAM | fast, simple | fits-in-memory regions/chromosomes |
| `GenotypeDaskArray` + zarr | chunked, on-disk, lazy | out-of-core / parallel | biobank-scale; `vcf_to_zarr` once then dask |
| scikit-allel | maintenance mode, v1.3.x | established stat workflows | the pragmatic default today |
| sgkit | xarray+dask, active | successor, NOT yet feature-parity | greenfield biobank-scale infrastructure |

## Decision Tree by Scenario

| Scenario | Use | Why |
|----------|-----|-----|
| Per-base pi/theta/Dxy | `windowed_diversity(..., is_accessible=mask)` | without the mask the per-base denominator is total span, not callable bp |
| Genome-wide FST point estimate + SE | `average_hudson_fst(ac1, ac2, blen)` | ratio-of-sums + block-jackknife done correctly; Hudson is robust to unequal n (Bhatia 2013) |
| FST landscape across the genome | `moving_hudson_fst` / `windowed_weir_cockerham_fst` | per-window ratio aggregation, not `mean(per_snp_fst)` |
| SFS without a confident ancestral allele | `sfs_folded(ac)` | folds on minor-allele count; `sfs()` is unfolded and treats ALT as derived |
| Test if pop C is admixed | `average_patterson_f3(acc, aca, acb, blen)` | significantly negative f3 (z < ~-3) is the formal admixture test; C goes FIRST |
| LD-prune before PCA | `locate_unlinked(gn)` iterated ~3 rounds | one pass leaves residual LD; PCs otherwise track LD blocks/inversions |
| Selection scan on phased data | `ihs`/`nsl` then `standardize_by_allele_count` (DAF bins); `xpehh` then genome-wide `standardize` | raw scores are uninterpretable; the standardization differs by statistic |
| Whole-genome callset (tens of M SNPs) | `vcf_to_zarr` + `GenotypeDaskArray` | `read_vcf` is eager and OOMs; dask materializes per chunk |

## Reading Genotypes and Counting Alleles

**Goal:** Load a VCF region into a GenotypeArray and derive the allele-count currency, missing-aware.

**Approach:** Read only the needed fields (read_vcf is eager), wrap GT, count alleles per site (missing ignored), and get per-population counts in one pass with count_alleles_subpops.

```python
import allel
import numpy as np

callset = allel.read_vcf('data.vcf.gz', fields=['samples', 'calldata/GT', 'variants/POS', 'variants/CHROM'], region='2L:1-5000000')
gt = allel.GenotypeArray(callset['calldata/GT'])   # (n_variants, n_samples, 2); -1 = missing
pos = callset['variants/POS']

ac = gt.count_alleles()                            # ignores -1, so per-site allele number varies
subpops = {'pop1': [0, 1, 2, 3, 4], 'pop2': [5, 6, 7, 8, 9]}
ac_subpops = gt.count_alleles_subpops(subpops)     # one pass, consistent variant axis
ac1, ac2 = ac_subpops['pop1'], ac_subpops['pop2']
```

## Per-base Diversity with the Accessibility Mask

**Goal:** Compute pi and Watterson's theta as honest per-base quantities, not span-deflated ones.

**Approach:** Pass a boolean callability mask (one entry per base, from coverage/mappability, NOT from variant positions) as is_accessible; inspect the returned n_bases per window to confirm the denominator.

```python
# is_accessible: bool array over genomic positions (a callable-loci mask), NOT the VCF variant sites.
pi = allel.sequence_diversity(pos, ac, is_accessible=is_accessible)
theta_w = allel.watterson_theta(pos, ac, is_accessible=is_accessible)

# windowed_diversity returns 4 values; n_bases is the accessible-bp denominator PER window.
pi_w, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size=100000, is_accessible=is_accessible)
# windowed_tajima_d returns 3 values (no n_bases): Tajima's D is dimensionless, no is_accessible.
D, td_windows, td_counts = allel.windowed_tajima_d(pos, ac, size=100000)
```

## FST as a Ratio of Sums

**Goal:** Get a genome-wide FST point estimate with a jackknife SE, and a per-window landscape, without the mean-of-ratios bias.

**Approach:** Let average_hudson_fst do the ratio-of-sums plus block-jackknife; if hand-aggregating, sum the components THEN divide; size blocks (blen) to exceed the LD scale.

```python
# Genome-wide estimate + standard error (ratio-of-sums + delete-one-block jackknife):
fst, se, vb, vj = allel.average_hudson_fst(ac1, ac2, blen=2000)   # blen must exceed the LD decay length

# Hand-aggregating Hudson correctly (NEVER mean of per-SNP fst):
num, den = allel.hudson_fst(ac1, ac2)
fst_manual = np.sum(num) / np.sum(den)

# Weir-Cockerham returns per-allele components (a, b, c); aggregate over BOTH axes:
a, b, c = allel.weir_cockerham_fst(gt, subpops=[[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]])
fst_wc = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))

# Landscape: per-window FST is already ratio-aggregated within each window.
fst_windows = allel.moving_hudson_fst(ac1, ac2, size=1000)
```

## Admixture: f3 and D (block-jackknifed)

**Goal:** Formally test whether a population is admixed (f3) or whether gene flow violates a tree (D / ABBA-BABA).

**Approach:** Call the average_* form for the jackknife z-score; put the test population FIRST in f3; treat a significantly negative f3 (z < ~-3) as admixture and |z| > ~3 for D as treeness violation.

```python
# f3(C; A, B): TEST population C is the FIRST argument. Returns (f3, se, z, vb, vj).
f3, se3, z3, vb3, vj3 = allel.average_patterson_f3(acc, aca, acb, blen=2000)
# Significantly negative f3 (z3 < ~-3) => C is admixed between A and B.

# D-statistic (ABBA-BABA): returns (d, se, z, vb, vj). |z| > ~3 flags gene flow.
d, sed, zd, vbd, vjd = allel.average_patterson_d(aca, acb, acc, acd, blen=2000)
```

## LD Pruning and PCA

**Goal:** Project samples onto ancestry axes that reflect drift, not LD blocks or inversions.

**Approach:** Convert to a missing-free 012 matrix, LD-prune iteratively with locate_unlinked, mask known inversions by position, then run Patterson-scaled PCA (randomized at scale).

```python
gn = gt.to_n_alt(fill=-1)                          # default fill=0 imputes missing to REF; use -1 then handle
gn = np.where(gn < 0, 0, gn)                        # impute-to-reference is a deliberate choice here

# Iterate LD pruning ~3 rounds; one pass leaves residual LD. Returns a KEEP mask (True = unlinked).
for _ in range(3):
    keep = allel.locate_unlinked(gn, size=100, step=20, threshold=0.1)
    gn = gn[keep]

coords, model = allel.randomized_pca(gn, n_components=10, scaler='patterson', random_state=0)
explained = model.explained_variance_ratio_        # scree; coords is (n_samples, n_components)
```

## Selection Scans on Phased Haplotypes

**Goal:** Score the genome for recent selection from haplotype structure.

**Approach:** Reshape phased genotypes to a HaplotypeArray, compute the raw scan, then standardize - iHS/nSL binned by derived-allele frequency (`standardize_by_allele_count`), XP-EHH genome-wide (`standardize`); the raw scores are not directly interpretable.

```python
h = gt.to_haplotypes()                              # VALID only if data are PHASED
ihs_raw = allel.ihs(h, pos, min_maf=0.05)           # unstandardized
ihs_std, bins = allel.standardize_by_allele_count(ihs_raw, ac[:, 1])   # bin by DERIVED count; ac[:,1] is derived only if REF is ancestral (polarize first); |z| > 2 flags candidates
h1, h12, h123, h2_h1 = allel.garud_h(h)             # soft-vs-hard-sweep haplotype-homozygosity stats
```

## Per-Function Failure Modes

### sequence_diversity / watterson_theta without is_accessible
**Trigger:** calling per-base diversity with no callability mask. **Mechanism:** divides the numerator by `stop-start+1` (total span) instead of callable bp. **Symptom:** pi/theta deflated 2-5x and the genome-wide landscape distorted because callability varies per window. **Fix:** pass `is_accessible=` from a coverage/mappability callable-loci mask and inspect the returned `n_bases`.

### mean(per_snp_fst) instead of ratio-of-sums
**Trigger:** averaging per-SNP `a/(a+b+c)` for the genome-wide FST. **Mechanism:** mean-of-ratios is dominated by low-frequency SNPs with tiny noisy denominators. **Symptom:** a biased FST that differs from published estimates of the same comparison (Bhatia 2013). **Fix:** `sum(a)/(sum(a)+sum(b)+sum(c))`, or `average_hudson_fst`/`average_weir_cockerham_fst` which do it plus a jackknife SE.

### to_n_alt default fill=0
**Trigger:** `gt.to_n_alt()` with no `fill`. **Mechanism:** missing calls become 0 alt alleles = homozygous reference. **Symptom:** PCA/LD silently biased toward the reference allele. **Fix:** `to_n_alt(fill=-1)` then impute deliberately, or pre-filter for high call rate; state the imputation choice.

### sfs() fed alt counts and run unfolded
**Trigger:** `allel.sfs(ac[:, 1])` without a confident ancestral allele. **Mechanism:** `sfs()` is unfolded and treats the ALT count as the DERIVED count; ALT != DERIVED. **Symptom:** mis-polarized spectrum biasing demographic/DFE inference. **Fix:** use `sfs_folded(ac)` when polarization is uncertain; use `sfs(dac)` only with a verified ancestral allele.

### iHS/XP-EHH/nSL on unphased or unstandardized data
**Trigger:** running selection scans on unphased genotypes or reporting raw scores. **Mechanism:** these stats need phased haplotype structure, and raw output is on an unstandardized scale. **Symptom:** meaningless scans; un-binned scores not comparable across the genome. **Fix:** require PHASED input and standardize - `standardize_by_allele_count` (DAF bins) for iHS/nSL, genome-wide `standardize` for XP-EHH.

### blen smaller than the LD scale
**Trigger:** a small `blen` in any `average_*` FST or `average_patterson_f3/_d`. **Mechanism:** blocks within an LD region are correlated, so the delete-one-block jackknife under-estimates the SE. **Symptom:** spurious-significant f3 admixture / D-statistics. **Fix:** size `blen` to exceed the LD decay length (multi-Mb / >~1 cM for humans).

### read_vcf on a whole-genome callset
**Trigger:** `allel.read_vcf('genome.vcf.gz')` with no region/fields limits. **Mechanism:** read_vcf is eager and materializes the entire file in RAM. **Symptom:** out-of-memory crash on biobank-scale data. **Fix:** `vcf_to_zarr` once, then `GenotypeDaskArray` for out-of-core counting/filtering; limit `read_vcf(fields=, region=)`.

## Quantitative / correctness notes

| Item | Value / rule | Rationale |
|------|--------------|-----------|
| Accessibility deflation | multiplicative AND per-window | a 40%-accessible window deflates pi ~2.5x, a 90% one ~1.1x - the relative landscape is wrong |
| FST estimator default | Hudson for unequal n / rare variants | Bhatia 2013 recommends the ratio estimator robust to sample-size imbalance |
| Jackknife block size | `blen` > LD decay length | too-small blocks are correlated -> anticonservative SE -> false significance |
| f3 admixture | z < ~-3 (negative) | a significantly negative f3(C; A, B) is the formal admixture test for C |
| D / ABBA-BABA | \|z\| > ~3 | conventional treeness-violation / gene-flow threshold |
| iHS/XP-EHH/nSL | standardized \|z\| > 2, in CLUSTERS | sweeps show clusters of extreme binned z-scores, not isolated SNPs |
| LD pruning rounds | ~3 iterations of `locate_unlinked` | one pass leaves residual LD; expect to discard most SNPs |
| PCA scaler | `'patterson'` (default) | centers then divides each SNP by sqrt(p(1-p)); equal expected variance under drift |

Thresholds are conventions, not laws - inspect distributions and verify current best practice before applying numbers blindly.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| pi/theta look 2-5x too small | `is_accessible=` omitted | pass a callable-loci mask; check the returned `n_bases` |
| FST disagrees with published value | `mean(per_snp_fst)` aggregation | `sum(num)/sum(den)` or `average_hudson_fst(ac1, ac2, blen)` |
| `hudson_fst`/`patterson_fst` "FST" out of range | treating the first return as FST | they return `(num, den)`; aggregate `np.sum(num)/np.sum(den)` |
| `ValueError` unpacking `windowed_tajima_d` | expecting 4 values | it returns 3 `(D, windows, counts)`; `windowed_diversity` returns 4 |
| PCA skewed toward reference allele | `to_n_alt()` default `fill=0` | `to_n_alt(fill=-1)` then impute deliberately, or pre-filter |
| `pca` raises on -1/NaN | missing values in the 012 matrix | impute or filter; the patterson scaler cannot handle missing |
| f3 admixture test makes no sense | wrong argument order | `average_patterson_f3(acc, aca, acb, blen)` - test pop C is FIRST |
| `blockwise_patterson_f3` AttributeError | old name | use `average_patterson_f3` / `average_patterson_d` |
| raw iHS values uninterpretable | not standardized | `standardize_by_allele_count(score, aac)` binned by DAF |
| MemoryError on `read_vcf` | eager whole-genome read | `vcf_to_zarr` + `GenotypeDaskArray`; limit `fields=`/`region=` |

## References

1. Miles A, Harding N, et al. scikit-allel: explore and analyse genetic variation. Zenodo; cite-all-versions DOI:10.5281/zenodo.597309. Docs: scikit-allel.readthedocs.io; successor sgkit: github.com/sgkit-dev/sgkit.
2. Weir BS, Cockerham CC. Estimating F-statistics for the analysis of population structure. Evolution 1984; 38(6):1358-1370. DOI:10.1111/j.1558-5646.1984.tb05657.x.
3. Hudson RR, Slatkin M, Maddison WP. Estimation of levels of gene flow from DNA sequence data. Genetics 1992; 132(2):583-589. PMID:1427045.
4. Bhatia G, Patterson N, Sankararaman S, Price AL. Estimating and interpreting FST: the impact of rare variants. Genome Research 2013; 23(9):1514-1521. PMID:23861382.
5. Patterson N, Moorjani P, Luo Y, Mallick S, Rohland N, Zhan Y, Genschoreck T, Webster T, Reich D. Ancient admixture in human history. Genetics 2012; 192(3):1065-1093. DOI:10.1534/genetics.112.145037.
6. Patterson N, Price AL, Reich D. Population structure and eigenanalysis. PLoS Genetics 2006; 2(12):e190. DOI:10.1371/journal.pgen.0020190.
7. Tajima F. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics 1989; 123(3):585-595. PMID:2513255.
8. Voight BF, Kudaravalli S, Wen X, Pritchard JK. A map of recent positive selection in the human genome. PLoS Biology 2006; 4(3):e72. DOI:10.1371/journal.pbio.0040072.
9. Sabeti PC, et al. Genome-wide detection and characterization of positive selection in human populations. Nature 2007; 449:913-918. DOI:10.1038/nature06250.
10. Ferrer-Admetlla A, Liang M, Korneliussen T, Nielsen R. On detecting incomplete soft or hard selective sweeps using haplotype structure. Molecular Biology and Evolution 2014; 31(5):1275-1291. DOI:10.1093/molbev/msu077.
11. Garud NR, Messer PW, Buzbas EO, Petrov DA. Recent selective sweeps in North American Drosophila melanogaster show signatures of soft sweeps. PLoS Genetics 2015; 11(2):e1005004. DOI:10.1371/journal.pgen.1005004.

## Related Skills

- selection-statistics - selection-scan design, standardization, and demography-aware outlier interpretation (this skill owns the array mechanics)
- population-structure - PCA and ADMIXTURE via PLINK2/FlashPCA2
- linkage-disequilibrium - LD pruning and clumping
- plink-basics - PLINK-format QC before array-based analysis
- variant-calling/vcf-basics - VCF generation and manipulation before loading
- phasing-imputation/haplotype-phasing - phased haplotypes required for iHS/XP-EHH/nSL
