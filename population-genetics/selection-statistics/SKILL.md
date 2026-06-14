---
name: bio-population-genetics-selection-statistics
description: Scans genomes for natural selection with SFS tests (Tajima's D, Fay & Wu H, Zeng E, SweepFinder2 CLR), haplotype tests (iHS, nSL, XP-EHH, Rsb, H12), and differentiation (FST, PBS) using scikit-allel, selscan, and SweepFinder2. No single statistic separates selection from demography at one locus, so the deliverable is empirical genome-wide outliers plus multiple orthogonal signals, not an absolute cutoff. iHS detects incomplete sweeps and collapses to zero at fixation while XP-EHH catches fixed sweeps; iHS/nSL standardize within derived-allele-frequency bins but XP-EHH gets a genome-wide z-score; derived-allele tests need substitution-model polarization; background selection mimics FST and CLR. Use when computing selection statistics like FST, Tajima's D, iHS, or XP-EHH, or scanning for selective sweeps. For phasing inputs see phasing-imputation/haplotype-phasing; for dN/dS see comparative-genomics/positive-selection.
tool_type: mixed
primary_tool: scikit-allel
---

## Version Compatibility

Reference examples tested with: scikit-allel 1.3+, numpy 1.26+, selscan 2.0+, SweepFinder2 1.0+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Version traps that change results, not just syntax: `allel.standardize_by_allele_count(score, aac, ...)` is for iHS and nSL (it bins by derived-allele count), while XP-EHH and iHH12 use plain genome-wide `allel.standardize(score)`; conflating them is a real bug. `allel.nsl(h)` takes no `pos`/`map_pos` (nSL is map-free by construction) whereas `allel.ihs(h, pos, map_pos=...)` distorts without a genetic map. selscan 2.0 adds `--unphased` for multilocus genotypes; selscan 1.x requires phased haplotypes. The single source of truth for versions is this block, not headings.

# Selection Statistics

**"Scan my population for signatures of natural selection"** -> Contrast each locus against the genome-wide neutral expectation that already carries the demographic history, using statistics that read complementary features of a sweep.
- Python: `allel.ihs()`, `allel.xpehh()`, `allel.nsl()`, `allel.garud_h()`, `allel.windowed_tajima_d()`, `allel.hudson_fst()` (scikit-allel)
- CLI: `selscan --ihs|--xpehh|--nsl` then `norm` for standardization; `SweepFinder2 -lrb` for CLR with a background-selection map

Scope: selection-scan DESIGN - SFS neutrality tests, haplotype/EHH sweep statistics, differentiation outliers (FST, PBS), CLR scans, and polygenic Qx, with standardization and demography-aware outlier interpretation. PCA/ADMIXTURE population assignment routes to population-structure; LD/EHH mechanics to linkage-disequilibrium; the scikit-allel array-API mechanics (signatures, accessibility masks, data loading) to scikit-allel-analysis; phasing to phasing-imputation/haplotype-phasing; dN/dS and McDonald-Kreitman to comparative-genomics/positive-selection.

## The Single Most Important Insight -- no single statistic separates selection from demography at one locus

1. A bottleneck removes variation genome-wide and skews the SFS toward rare variants exactly as a sweep does, so a single Tajima's D cannot tell them apart; recent expansion gives genome-wide negative D that mimics a post-sweep signal, and population structure gives positive D that mimics balancing selection.
2. The field's response is structural, not statistical: rank windows and flag empirical genome-wide OUTLIERS (the bulk absorbs the shared demography), intersect MULTIPLE ORTHOGONAL signals (a sweep distorts haplotype, SFS, diversity, and differentiation together while each confound acts differently), and calibrate against an explicit demographic-model SIMULATION (dadi/momi/fastsimcoal2 + msprime/SLiM).
3. iHS and nSL detect INCOMPLETE/ongoing sweeps and collapse to ~0 once the allele fixes; XP-EHH and Rsb catch FIXED/near-fixed sweeps by borrowing a second population, so they are COMPLEMENTARY not redundant (high XP-EHH with near-zero iHS is the textbook completed-sweep signature).
4. An absolute cutoff such as "Tajima's D < -2 means a sweep" is wrong by construction; -2 means nothing without the genome-wide context, which is why the entire skill is outlier-based.

## Tool Taxonomy

| Family | Statistic / tool | Citation | Reads / role | Phasing |
|--------|------------------|----------|--------------|---------|
| SFS | Tajima's D (`allel.tajima_d`) | Tajima 1989 | theta_pi - theta_W; rare-vs-intermediate variant skew | no |
| SFS | Fay & Wu H, Zeng E | Fay & Wu 2000; Zeng 2006 | high-frequency-DERIVED excess; needs polarization | no |
| SFS-CLR | SweepFinder2 (`-f/-l/-lr/-lrb`) | DeGiorgio 2016 | local SFS vs empirical background; localizes target; `-lrb` deconfounds BGS | no |
| SFS-CLR | SweeD, OmegaPlus | Pavlidis 2013; Alachiotis 2012 | parallel CLR; LD-shoulder omega | no |
| Haplotype | iHS (`allel.ihs`, selscan `--ihs`) | Voight 2006 | INCOMPLETE sweeps; standardize by DAF bin | YES |
| Haplotype | nSL (`allel.nsl`, selscan `--nsl`) | Ferrer-Admetlla 2014 | incomplete hard+soft; map-free | YES |
| Haplotype | XP-EHH (`allel.xpehh`, selscan `--xpehh`) | Sabeti 2007 | FIXED/near-fixed sweeps; genome-wide z-score | YES |
| Haplotype | Rsb | Tang 2007 | XP-EHH niche; iES ratio | YES |
| Haplotype | H12 / H2H1 (`allel.garud_h`) | Garud 2015 | SOFT sweeps; hard-vs-soft tendency | YES |
| Differentiation | FST (`allel.hudson_fst`, `allel.weir_cockerham_fst`) | Weir & Cockerham 1984; Bhatia 2013 | local-adaptation divergence | no |
| Differentiation | PBS | Yi 2010 | which branch the divergence is on (3 pops) | no |
| Polygenic | Qx | Berg & Coop 2014 | over-dispersed polygenic scores across pops | no |

## Decision Tree by Scenario

| Scenario | Use | Why |
|----------|-----|-----|
| Localized recent sweep, ongoing/incomplete, one population, genetic map available | iHS | reads the frequency contrast between long derived and short ancestral haplotypes |
| Same but no reliable genetic map | nSL | integrates over segregating-site count; map-free, robust to recombination-rate variation |
| Sweep complete/fixed in one of two populations | XP-EHH or Rsb | maximum power exactly where iHS is blind (no within-population contrast remains) |
| Selection on standing variation (soft sweep) | H12 / H2H1 | collapses co-rising haplotypes; hard-sweep stats lose power |
| Localize the target with a demographic/BGS model | SweepFinder2 CLR (`-lrb`) or SweeD; OmegaPlus | empirical background SFS as null; B-value map deconfounds background selection |
| Differentiation outliers, unequal sample sizes | FST Hudson estimator + PBS for 3 pops | Hudson is robust to unequal n; PBS localizes the branch |
| Quick SFS reconnaissance | windowed Tajima's D + pi + Fay & Wu H (substitution-model polarization) | cheap, but never an absolute cutoff |
| Polygenic trait shift across populations | Qx with within-family GWAS effect sizes only | sweep scans are blind to coordinated tiny shifts; standard GWAS betas carry stratification bias |
| Ancient/recurrent coding selection across species | hand off to comparative-genomics/positive-selection | dN/dS and MK are a different timescale, not a within-population sweep |

## FST and PBS - Differentiation

**Goal:** Quantify allele-frequency divergence between populations and flag local-adaptation candidates without being fooled by rare variants or unequal samples.

**Approach:** Count alleles per subpopulation, compute the Hudson estimator per SNP, and report mean FST as a ratio-of-averages over windows after an MAF filter (rare variants deflate FST, Bhatia 2013).

```python
import allel
import numpy as np

callset = allel.read_vcf('data.vcf.gz')
gt = allel.GenotypeArray(callset['calldata/GT'])
pos = callset['variants/POS']

subpops = {'pop1': [0, 1, 2, 3, 4], 'pop2': [5, 6, 7, 8, 9]}
ac_subpops = gt.count_alleles_subpops(subpops)
ac1, ac2 = ac_subpops['pop1'], ac_subpops['pop2']

# MAF filter first: rare variants systematically deflate FST (Bhatia 2013).
maf = np.minimum(ac1.to_frequencies()[:, 1], ac2.to_frequencies()[:, 1])
keep = (maf > 0.05) | (1 - maf > 0.05)

# Hudson estimator: robust to unequal sample sizes. Mean FST = ratio-of-averages, never mean of per-SNP ratios.
num, den = allel.hudson_fst(ac1[keep], ac2[keep])
fst_mean = np.nansum(num) / np.nansum(den)

# Windowed scan for outlier localization.
fst_win, windows, n_snps = allel.windowed_hudson_fst(pos[keep], ac1[keep], ac2[keep], size=100000, step=50000)
outlier = fst_win > np.nanpercentile(fst_win, 99)
```

PBS (Population Branch Statistic, Yi 2010) needs three populations: convert three pairwise FST values to branch lengths via T = -log(1 - FST) and isolate the focal branch as PBS = (T_12 + T_13 - T_23) / 2. A long focal branch in a small/bottlenecked population can be drift, not selection.

## iHS and nSL - Within-Population Incomplete Sweeps

**Goal:** Detect ongoing sweeps from extended haplotype homozygosity around derived core alleles.

**Approach:** Filter to segregating biallelic SNPs, compute the raw integrated-EHH score, then standardize WITHIN derived-allele-frequency bins (the raw score depends strongly on derived frequency and is uninterpretable unbinned).

```python
import allel
import numpy as np

h = gt.to_haplotypes()
ac = h.count_alleles()
flt = (ac[:, 0] > 1) & (ac[:, 1] > 1)
h_flt, pos_flt, ac_flt = h.compress(flt, axis=0), pos[flt], ac.compress(flt, axis=0)

# iHS needs a genetic map (map_pos) where recombination rate varies; physical distance distorts iHH.
ihs_raw = allel.ihs(h_flt, pos_flt, min_maf=0.05, include_edges=True)
# Standardize WITHIN derived-allele-count bins (NOT genome-wide) - the most common iHS bug.
ihs_std, _ = allel.standardize_by_allele_count(ihs_raw, ac_flt[:, 1])
ihs_outlier = np.abs(ihs_std) > 2

# nSL: map-free (no pos/map_pos), robust to recombination-rate variation; still bin-standardized.
nsl_raw = allel.nsl(h_flt)
nsl_std, _ = allel.standardize_by_allele_count(nsl_raw, ac_flt[:, 1])
```

A null iHS does NOT mean no selection: once the derived allele fixes, the ancestral haplotype is gone and iHS collapses to zero. Switch to XP-EHH/Rsb for completed sweeps. Report the windowed PROPORTION of |iHS|>2 SNPs, not single noisy hits.

## XP-EHH - Cross-Population Completed Sweeps

**Goal:** Catch fixed or near-fixed sweeps in one of two populations, where within-population haplotype tests are blind.

**Approach:** Compute the cross-population integrated-EHH ratio on shared positions, then apply a plain GENOME-WIDE z-score (XP-EHH is not strongly correlated with derived frequency, so frequency-bin standardization is the wrong rule).

```python
import allel
import numpy as np

h1 = h_flt.take(pop1_hap_idx, axis=1)
h2 = h_flt.take(pop2_hap_idx, axis=1)

xpehh_raw = allel.xpehh(h1, h2, pos_flt, include_edges=True)
# Genome-wide z-score for XP-EHH/iHH12 - NOT standardize_by_allele_count (that is for iHS/nSL).
xpehh_std = allel.standardize(xpehh_raw)
completed_sweep = np.abs(xpehh_std) > 2
```

A shared ancestral sweep at the same locus in BOTH populations cancels in the contrast (false negative). Differential phasing quality between the two populations systematically biases the difference statistic.

## H12 - Soft Sweeps

**Goal:** Recover sweeps on standing variation that iHS and CLR miss because no single long haplotype dominates.

**Approach:** Compute Garud's H over haplotype-frequency spectra; H12 collapses the two most common haplotypes into one, and H2/H1 tends (not classifies) toward hard-vs-soft.

```python
h1, h12, h123, h2_h1 = allel.garud_h(h_flt)
h12_win = allel.moving_garud_h(h_flt, size=100)  # SNP-count windows, not bp
```

## SFS-CLR with selscan and SweepFinder2

```bash
# selscan haplotype scans (phased VCF + genetic map), then norm for standardization.
selscan --ihs --vcf phased.vcf --map genetic.map --out scan        # within-pop incomplete sweeps
selscan --xpehh --vcf pop1.vcf --vcf-ref pop2.vcf --map genetic.map --out xp  # completed sweeps
selscan --nsl --vcf phased.vcf --out nsl_scan                       # no map; map-free
selscan --ihs --unphased --vcf unphased.vcf --map genetic.map --out scan_unphased  # selscan 2.0 multilocus-genotype

# norm bins iHS/nSL by derived-allele frequency (--bins default 100); XP-EHH/iHH12 get a genome-wide z-score.
norm --ihs --files scan.ihs.out --bins 100
norm --xpehh --files xp.xpehh.out

# SweepFinder2 CLR: build a genome-WIDE background SFS, then scan with recombination + B-value (BGS) map.
SweepFinder2 -f genome.freq genome.spect                            # empirical background spectrum
# -lrb takes N1 (current ingroup Ne), N2 (ancestral Ne), T (divergence time in generations) before OutFile.
SweepFinder2 -lrb 2000 region.freq genome.spect region.rec bvalue.map "$N1" "$N2" "$T" out.clr  # G=2000 grid; -lrb deconfounds BGS
```

## Polarization and Polygenic Selection

Derived-allele tests (Fay & Wu H, Zeng E, unfolded SFS, the iHS sign) are POLARIZATION traps: one mispolarized CpG site masquerades as a high-frequency-derived variant and fakes a strongly negative H. Polarize with a probabilistic substitution model and two or more outgroups (Hernandez 2007), never a single chimp allele. Conservation scores (phyloP, GERP, phastCons) measure CONSTRAINT (purifying selection), the opposite sign of a recent positive sweep, and must not be read as "under selection." Polygenic Qx (Berg & Coop 2014) inherits every stratification bias in the GWAS effect sizes it uses: the celebrated European height-selection signal was shown to be a stratification artifact (Sohail 2019; Berg 2019), so use within-family/sibling effect estimates and treat standard-GWAS-beta Qx as provisional.

## Per-Method Failure Modes

### Absolute Tajima's D cutoff
**Trigger:** flagging windows on |D|>2 as selection. **Mechanism:** growth/bottleneck/structure produce the same sign as a sweep/balancing selection (non-identifiable at one locus). **Symptom:** genome-wide false positives that are pure demography. **Fix:** use the genome-wide empirical distribution and intersect orthogonal statistics; report D relative to the genomic background.

### Unstandardized or mis-standardized iHS
**Trigger:** reporting raw iHS, or applying a genome-wide z-score to iHS. **Mechanism:** raw iHH ratio depends strongly on derived-allele frequency. **Symptom:** uninterpretable scores, false tails at particular frequencies. **Fix:** `standardize_by_allele_count` (DAF bins) for iHS/nSL; reserve `standardize` (genome-wide) for XP-EHH/iHH12.

### iHS null read as no-selection
**Trigger:** concluding neutrality from absent iHS signal. **Mechanism:** a completed sweep fixes the derived allele, erasing the within-population contrast iHS needs. **Symptom:** real fixed sweeps missed. **Fix:** add XP-EHH/Rsb against a second population for completed sweeps.

### Single-outgroup polarization
**Trigger:** deriving ancestral state from one chimp allele for H/E/unfolded SFS. **Mechanism:** recurrent mutation at CpG/hypermutable sites biases the single-allele ancestral estimator. **Symptom:** spurious strongly-negative Fay & Wu H (fake sweep). **Fix:** substitution-model polarization with >=2 outgroups (Hernandez 2007; est-sfs).

### Background selection as sweep
**Trigger:** reading an FST or CLR peak in a low-recombination region as positive selection. **Mechanism:** purifying selection against linked deleterious variants reduces local diversity. **Symptom:** FST/PBS/CLR peaks with no positive selection. **Fix:** supply SweepFinder2 a B-value map (`-lrb`); control for recombination rate; treat low-recombination outliers with suspicion.

### Switch error in haplotype stats
**Trigger:** EHH statistics on poorly phased data. **Mechanism:** a switch error truncates long haplotypes and deflates EHH. **Symptom:** iHS/nSL/XP-EHH biased toward null; XP-EHH biased by differential phasing between populations. **Fix:** read-backed/high-quality phasing (see phasing-imputation/haplotype-phasing), or selscan 2.0 `--unphased` multilocus-genotype mode.

### Polygenic Qx on stratified GWAS betas
**Trigger:** running Qx on standard meta-analysis effect sizes. **Mechanism:** residual population stratification correlates effect sizes with ancestry axes. **Symptom:** spurious polygenic-adaptation clines (the height story). **Fix:** within-family/sibling effect sizes; treat PC-corrected meta-analysis betas as insufficient.

## Quantitative Thresholds

| Quantity | Typical value | Rationale |
|----------|---------------|-----------|
| iHS / nSL flag | \|standardized score\| > 2 | ~top 2.5% two-sided of N(0,1); an empirical convention, NOT a calibrated p-value (Voight 2006) |
| Empirical outlier tail | top 1% (top 0.1% for stringency) | genome-wide bulk absorbs demography; percentile is a sensitivity/specificity tradeoff, not an alpha |
| iHS/nSL standardization | within derived-allele-frequency bins (selscan `--bins 100`) | raw score is frequency-dependent; binning makes scores comparable |
| XP-EHH/iHH12 standardization | genome-wide z-score | not frequency-correlated; bin-standardizing inverts the rule |
| Window size | 10-100 kb, 50% step (scale to LD decay) | wide enough for stable SFS estimates, narrow enough to localize |
| Haplotype-stat MAF | drop core alleles MAF < ~0.05 | EHH on singletons is undefined/uninformative (`min_maf=0.05`) |
| SweepFinder2 grid G | finer than expected sweep width (~1-2 kb in humans) | coarse grids step over narrow sweeps; significance from demographic-model simulation |
| FST mean | ratio-of-averages, MAF-filtered, Hudson estimator | rare variants deflate FST; arithmetic mean of per-SNP ratios is biased (Bhatia 2013) |

Thresholds are conventions, not laws - inspect distributions and calibrate against a demographic-model simulation or the genome-wide empirical tail before quoting any number.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "D < -2 = sweep" | absolute cutoff ignoring demography | use the genome-wide empirical tail; D sign never identifies cause |
| iHS scores incomparable across SNPs | unstandardized or genome-wide-standardized iHS | `standardize_by_allele_count` by DAF bin |
| XP-EHH never standardized | applying the iHS rule or no standardization | `allel.standardize(xpehh)` genome-wide z-score |
| "no iHS so no selection" | completed sweep has no within-pop contrast | add XP-EHH/Rsb against a second population |
| spurious negative Fay & Wu H | single-outgroup polarization at CpG sites | substitution-model polarization, >=2 outgroups (Hernandez 2007) |
| FST/CLR peak misread as positive selection | background selection in low-recombination region | B-value map (`-lrb`); recombination control |
| `allel.nsl(h, pos)` TypeError | nSL is map-free; takes no pos/map_pos | call `allel.nsl(h)` |
| FST tail driven by rare variants | no MAF filter; mean of per-SNP ratios | MAF filter + Hudson ratio-of-averages |
| Qx height-style false signal | stratified GWAS effect sizes | within-family/sibling effect estimates |
| "high phyloP = under selection" | conflating constraint with positive selection | phyloP/GERP measure purifying constraint (opposite sign) |

## References

1. Tajima F. Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics 1989; 123(3):585-595.
2. Fay JC, Wu CI. Hitchhiking under positive Darwinian selection. Genetics 2000; 155(3):1405-1413.
3. Zeng K, Fu YX, Shi S, Wu CI. Statistical tests for detecting positive selection by utilizing high-frequency variants. Genetics 2006; 174(3):1431-1439.
4. DeGiorgio M, Huber CD, Hubisz MJ, Hellmann I, Nielsen R. SweepFinder2: increased sensitivity, robustness and flexibility. Bioinformatics 2016; 32(12):1895-1897.
5. Voight BF, Kudaravalli S, Wen X, Pritchard JK. A map of recent positive selection in the human genome. PLoS Biology 2006; 4(3):e72.
6. Sabeti PC, Varilly P, Fry B, et al. Genome-wide detection and characterization of positive selection in human populations. Nature 2007; 449(7164):913-918.
7. Ferrer-Admetlla A, Liang M, Korneliussen T, Nielsen R. On detecting incomplete soft or hard selective sweeps using haplotype structure. Molecular Biology and Evolution 2014; 31(5):1275-1291.
8. Garud NR, Messer PW, Buzbas EO, Petrov DA. Recent selective sweeps in North American Drosophila melanogaster show signatures of soft sweeps. PLoS Genetics 2015; 11(2):e1005004.
9. Yi X, Liang Y, Huerta-Sanchez E, et al. Sequencing of 50 human exomes reveals adaptation to high altitude. Science 2010; 329(5987):75-78.
10. Weir BS, Cockerham CC. Estimating F-statistics for the analysis of population structure. Evolution 1984; 38(6):1358-1370.
11. Bhatia G, Patterson N, Sankararaman S, Price AL. Estimating and interpreting FST: the impact of rare variants. Genome Research 2013; 23(9):1514-1521.
12. Hernandez RD, Williamson SH, Bustamante CD. Context dependence, ancestral misidentification, and spurious signatures of natural selection. Molecular Biology and Evolution 2007; 24(8):1792-1800.
13. Berg JJ, Coop G. A population genetic signal of polygenic adaptation. PLoS Genetics 2014; 10(8):e1004412.
14. Sohail M, Maier RM, Ganna A, et al. Polygenic adaptation on height is overestimated due to uncorrected stratification in genome-wide association studies. eLife 2019; 8:e39702.
15. Berg JJ, Harpak A, Sinnott-Armstrong N, et al. Reduced signal for polygenic adaptation of height in UK Biobank. eLife 2019; 8:e39725.

## Related Skills

- scikit-allel-analysis - genotype/haplotype array loading and allele-count basics
- population-structure - PCA and ADMIXTURE for population assignment before FST/PBS
- linkage-disequilibrium - EHH and LD mechanics underlying haplotype statistics
- phasing-imputation/haplotype-phasing - phased haplotypes that all EHH statistics require
- comparative-genomics/positive-selection - dN/dS and McDonald-Kreitman for recurrent coding selection
- comparative-genomics/introgression-detection - archaic/admixture signals that confound differentiation outliers
