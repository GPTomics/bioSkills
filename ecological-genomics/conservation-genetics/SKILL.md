---
name: bio-ecological-genomics-conservation-genetics
description: Assesses genetic health of populations for conservation with Ne estimation across time horizons (LDNe NeEstimator V2 option-file API + SNeP physical-linkage correction; recent trajectory via GONE/GONE2; deep history via Stairway Plot 2 / dadi / fastsimcoal2 / PSMC), F-statistics, runs of homozygosity binned by length class to date inbreeding, genetic-load decomposition (Bertorelle 2022 realized vs masked), the modern 100/1000 Ne rule (Frankham 2014), Ne/Nc 2-6 orders of magnitude in marine fish (Hauser & Carvalho 2008), tree-sequence forward simulations (SLiM 4 + pyslim + tskit), and the Sukumaran-Knowles caveat against MSC methods for management-unit definition. Use when estimating Ne by time horizon, detecting inbreeding via F_ROH, decomposing genetic load, justifying conservation thresholds, distinguishing ESU/MU/DPS, configuring NeEstimator V2, or correcting LDNe physical linkage.
tool_type: mixed
primary_tool: hierfstat
---

## Version Compatibility

Reference examples tested with: hierfstat 0.5+, adegenet 2.1+, detectRUNS 2.0+, poppr 2.9+, NeEstimator V2.1+, GONE2, SNeP 1.1+, SLiM 4+, msprime 1.3+, bcftools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Conservation Genetics

**"Assess the genetic health of my endangered population"** -> Estimate Ne by time horizon (contemporary LD, recent trajectory, deep demographic history), measure inbreeding via F_ROH with length-class thresholds, decompose realized vs masked genetic load, and interpret against the modern 100/1000 Ne rule.
- R: `hierfstat::basic.stats()` for F-statistics and diversity
- R: `detectRUNS::consecutiveRUNS.run()` or `bcftools roh -G30` for ROH detection
- CLI: `Ne2.exe option_file.ne2` for NeEstimator V2 contemporary Ne (option-file driven, NOT CLI flags)
- CLI: `gone2` for recent Ne trajectory from genome-wide LD (requires populated cM map)

## The Single Most Important Modern Insight -- There is No Single Best Ne Estimator

The most common methodological failure in conservation genetics is treating "Ne" as one quantity. **It is not.** Each Ne estimator captures a different time horizon with different biases:
- **LD-based (NeEstimator V2, Waples 2006)**: Last 1-few generations; biased downward at small N
- **GONE / GONE2 (Santiago 2020)**: Recent 100-200 generations trajectory; requires populated cM genetic map
- **Heterozygote-excess**: Last 1 generation only; very low power
- **Temporal (Jorde-Ryman)**: Between two time-sampled cohorts
- **Stairway Plot 2 / dadi / fastsimcoal2**: 10^3 to 10^5 generations from SFS
- **PSMC**: 10^4 to 10^6 generations from single high-coverage genome

A second cornerstone: **Ne/Nc ratio is NOT 0.1 universally.** Frankham 1995 reported a median of 0.1 in vertebrates, but Hauser & Carvalho 2008 *Fish Fish* 9:333-362 documented Ne/Nc spanning 2-6 orders of magnitude smaller than census (i.e., 10^-2 to 10^-6) in marine fish with sex-biased and sweepstakes-recruitment reproduction. Conservation papers that assume Ne/Nc = 0.1 are wrong for many taxa.

A third: **the 50/500 rule was revised to 100/1000** by Frankham et al. 2014 *Biol Conserv* 170:56-63. The 50/500 numbers came from 1980 and underestimate the Ne needed for adaptive maintenance. The modern thresholds are Ne >= 100 for short-term inbreeding-fitness protection and Ne >= 1000 for long-term adaptive potential.

## Algorithmic Taxonomy

| Method | Time horizon | Strength | Bias mode |
|--------|--------------|----------|-----------|
| LDNe (NeEstimator V2 LD method) | Last 1-few generations | Single sample; SNP-friendly | Biased downward at small N; sensitive to physical linkage |
| SNeP (LD with physical-linkage correction) | Last 1-few generations | Corrects LDNe for chromosomal linkage in genomic data | Requires phased genotypes |
| GONE / GONE2 | Recent 100-200 generations | Detects bottleneck/expansion trajectory | Requires populated cM genetic map column (NOT physical position) |
| Heterozygote-excess | Last 1 generation | Detects recent bottleneck | Very low power; requires high heterozygosity |
| Temporal (Jorde-Ryman, Pollak) | Between two cohorts | Direct drift estimate | Requires temporal samples |
| Sibship (Wang COLONY) | Current generation | Pedigree-based from kinship | Computationally expensive |
| dadi (Gutenkunst 2009) | 10^3 - 10^5 generations | Diffusion SFS-based; analytical | Local-optima trap; <=3 populations practical |
| fastsimcoal2 (Excoffier 2013) | 10^3 - 10^5 generations | Simulation-based composite-likelihood | Optimization needs >= 50 replicates |
| moments (Jouganous 2017) | 10^3 - 10^5 generations | Moment-based ODE; faster than dadi | Same dimensionality limits |
| momi2 (Kamm 2020) | 10^3 - 10^5 generations | Analytical SFS via moments | Stiff at very recent times |
| Stairway Plot 2 | 10^3 - 10^5 generations | Folded SFS; no parametric model | Requires explicit mutation rate and generation time |
| PSMC | 10^4 - 10^6 generations | Whole-genome pairwise coalescent | Requires >= 20x WGS coverage |
| MSMC2 | Multi-individual deep coalescent | Better recent resolution than PSMC | Requires phased data; computational |
| msprime (Kelleher 2016) | Simulation (not inference) | Gold-standard neutral simulator | Used to generate data under a fitted model |
| SLiM 3/4 (Haller & Messer 2019) | Forward simulation | Non-WF, age structure, spatial | Eidos scripting; tree-sequence recording essential for speed |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Contemporary Ne from a single SNP sample | NeEstimator V2 LD method with Pcrit = 0.02 | Standard; widely accepted by conservation reviewers |
| LDNe applied to RAD-seq / WGS | Use SNeP (Barbato 2015) OR thin SNPs to >= 1 cM apart | LDNe assumes inter-locus r^2 reflects demography only; physical linkage inflates it |
| Recent Ne trajectory (bottleneck detection) | GONE2 with populated cM genetic-map column | GONE2 fails silently if cM column is zero (PLINK default uses physical position) |
| Deep demographic history (mammal/vertebrate) | PSMC with >= 20x WGS or Stairway Plot 2 from SFS | PSMC requires high coverage; Stairway Plot 2 needs only the SFS |
| Multi-population demographic inference | dadi / fastsimcoal2 with >= 50 independent optimization replicates | Likelihood surfaces have local optima; single-replicate inference is unreliable |
| Inbreeding from genomic data | F_ROH with length-class thresholds (bcftools roh OR detectRUNS) | F_IS is genome-averaged; F_ROH partitions inbreeding by time |
| Genetic-load assessment | Decompose realized vs masked via Bertorelle 2022 framework | Single "load" estimate hides the dynamics |
| Define management unit (MU) or ESU | Moritz 1994 reciprocal monophyly + nuclear divergence | Do NOT use BPP/BFD* species-delimitation methods (Sukumaran-Knowles 2017 oversplitting) |
| Forward simulation of selection + demography | SLiM 4 with tree-sequence recording + pyslim/tskit | Tree-sequence recording is 5-100x faster than mutation tracking |
| Comparing Ne across species | Cite Ne/Nc taxonomic variation (Hauser & Carvalho 2008) | The Ne/Nc = 0.1 default is wrong for many taxa |
| Genetic rescue decision | Frankham 2015 meta-analysis; usually favors rescue | Outbreeding-depression fear is overstated for most cases |

## Genetic Diversity and F-Statistics

**Goal:** Compute standard population-genetic metrics with Weir-Cockerham estimators and bootstrap confidence intervals.

**Approach:** Convert VCF to genind via adegenet, then to hierfstat format. Run `basic.stats` for Fis/Fst/Ho/He and `pairwise.WCfst` with bootstrapping for population differentiation.

```r
library(hierfstat)
library(adegenet)
library(poppr)

# VCF/genepop/PLINK -> genind via adegenet
data_genind <- read.genepop('populations.gen')
data_hf <- genind2hierfstat(data_genind)

# F-statistics with Weir-Cockerham estimators
bstats <- basic.stats(data_hf)
cat('Overall Fis:', bstats$overall['Fis'], '\n')
cat('Overall Fst:', bstats$overall['Fst'], '\n')

# Pairwise FST with bootstrap CIs
pw_fst <- pairwise.WCfst(data_hf)
boot_fst <- boot.ppfst(data_hf, nboot = 1000)

# Rarefied allelic richness (corrects for unequal sample sizes)
ar <- allelic.richness(data_hf)

# Private alleles unique to each population
pa <- private_alleles(data_genind, count.alleles = TRUE)
```

## Runs of Homozygosity with Length-Class Dating

**Goal:** Detect autozygous segments and date the inbreeding events by ROH length class.

**Approach:** Use `bcftools roh -G30` (HMM-based, density-independent) for VCF input OR `detectRUNS::consecutiveRUNS.run()` for PLINK. Bin ROH by length to date inbreeding events: > 16 Mb means parents/grandparents; 4-16 Mb means ~5 generations; 1-4 Mb means ~10-20 generations; < 1 Mb is deep background.

```r
library(detectRUNS)

# detectRUNS: SNP-by-SNP scanning from PLINK files
runs <- consecutiveRUNS.run(
    genotypeFile = 'genotypes.ped',
    mapFile = 'genotypes.map',
    minSNP = 20,        # Minimum SNPs per run; prevents short-stretch false positives
    minLengthBps = 1e6, # 1 Mb minimum; shorter ROH usually not IBD-derived
    maxGap = 1e6,
    maxOppRun = 1,
    maxMissRun = 2
)

# F_ROH per individual = sum(ROH length) / autosomal genome length
# > 16 Mb ROH: parents/grandparents inbreeding (very recent)
# 4-16 Mb: ~5 generations back
# 1-4 Mb: ~10-20 generations (historical bottleneck)
# < 1 Mb: deep background, ancient homozygosity
froh <- Froh_inbreeding(runs, mapFile = 'genotypes.map', genome_wide = TRUE)

# Bin ROH by length to estimate inbreeding timing
runs_summary <- summaryRuns(runs, genotypeFile = 'genotypes.ped',
                            mapFile = 'genotypes.map')
plot_DistributionRuns(runs)
```

```bash
# Alternative: bcftools roh with HMM (density-independent; better for low-density data)
# -G30: per-sample genotype likelihood threshold (30 = high-quality calls)
# --AF-tag AF: use AF tag from VCF for allele frequencies (else --AF-file)
bcftools roh -G30 --AF-tag AF -o roh_results.txt input.vcf.gz
```

## NeEstimator V2 — Option-File API

**Goal:** Estimate contemporary Ne via LD method from a single sample (Do et al. 2014 *Mol Ecol Resour* 14:209-214).

**Approach:** NeEstimator V2 is option-file driven, NOT CLI-flag driven. Create a `.ne2` text option file specifying input format, methods (LD/HetExcess/Coancestry/Temporal), Pcrit thresholds, and output. Run with `Ne2.exe option_file.ne2`.

```text
# Example NeEstimator V2 option file (option_file.ne2)
# Line order matters; defaults can cause silent failures
# See NeEstimator V2 documentation: https://github.com/bunop/NeEstimator2.X

# Input format
1 0                                   # 1=GenePop, 2=FSTAT; second value reserved
input_genotypes.gen

# Methods (1 = run; 0 = skip)
1 0 0 0                               # LD only; others (HetExcess, Coancestry, Temporal) off

# Pcrit cutoffs (allele frequencies below threshold excluded)
3                                     # number of Pcrit values
0.05 0.02 0.01                        # 0.02 is standard; 0.05 conservative; 0.01 sensitive

# Mating system: 0 = random mating; 1 = monogamy
0

# Output
output_results.txt
```

```bash
# Run NeEstimator V2 (Java-based)
java -jar NeEstimator.jar option_file.ne2
# Verify INFO output line: confirms which methods ran and which were skipped

# Common failure: "Ne = infinity" reported
# Diagnosis: insufficient drift signal (population larger than method can detect)
# Try multiple Pcrit; with genomic data add SNeP physical-linkage correction
```

## GONE2 — Recent Ne Trajectory with the cM-Column Trap

**Goal:** Estimate the Ne trajectory over the last ~100-200 generations from genome-wide LD across recombination rates.

**Approach:** GONE2 requires a genetic map with POPULATED cM positions (not just physical position). PLINK MAP files default to cM=0; if this is not corrected, GONE2 will silently produce nonsense. Verify the cM column is populated from a genetic-map file or use Hi-C-derived recombination map.

```bash
# GONE2 input: PLINK BED/BIM/FAM or VCF + populated MAP
# CRITICAL: BIM/MAP cM column must be populated; default cM=0 produces silent failure

# Check cM column populated
head genotypes.bim   # column 3 should NOT be all zeros

# Run GONE2
# -t 4: threads
# -u 0.05: upper recombination rate bound (exclude pairs with r > 0.05)
# Smaller -u focuses on more recent generations
./gone2 -t 4 -u 0.05 genotypes.vcf

# Output: OUTPUT_GONE2 with generation, Ne, CI_low, CI_high
```

```r
# Parse GONE2 output and plot Ne trajectory
gone_out <- read.table('OUTPUT_GONE2', header = TRUE, sep = '\t')

pdf('gone2_ne_trajectory.pdf', width = 8, height = 5)
plot(gone_out$generation, gone_out$Ne,
     type = 'l', lwd = 2, col = 'blue',
     xlab = 'Generations ago', ylab = 'Effective population size (Ne)',
     main = 'Recent Ne Trajectory (GONE2)', log = 'y')
polygon(c(gone_out$generation, rev(gone_out$generation)),
        c(gone_out$CI_low, rev(gone_out$CI_high)),
        col = adjustcolor('blue', alpha = 0.2), border = NA)
dev.off()
```

## SNeP — Physical-Linkage-Corrected LDNe for Genomic Data

**Goal:** Apply physical-linkage correction to LDNe when SNPs come from RAD-seq or WGS (where inter-locus r^2 reflects chromosomal linkage, not just demography).

**Approach:** SNeP (Barbato 2015 *Front Genet* 6:109) implements Waples & Do's physical-linkage correction. Use when SNPs are not pre-thinned to >= 1 cM apart.

```bash
# SNeP input: PLINK PED/MAP or VCF + map file
# Multi-threaded; supports several corrections (sample size, mutation, phasing, recombination)
./SNeP1.1 -ped genotypes.ped -map genotypes.map -threads 4 \
          -mutationrate 1.4e-8 -out snep_results.txt

# Output: Ne estimates per recombination-rate bin
```

## Deep Demographic History

**Goal:** Reconstruct Ne trajectory over thousands of generations from genome-wide variation.

**Approach:** Pick the appropriate tool by data type. PSMC for a single high-coverage diploid genome. Stairway Plot 2 for SFS from population samples. dadi or fastsimcoal2 for multi-population history with composite likelihood.

```bash
# --- PSMC: single high-coverage WGS (>= 20x) ---
bcftools mpileup -C50 -Q 30 -q 30 -f reference.fa sample.bam | \
    bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 > consensus.fq
fq2psmcfa -q20 consensus.fq > consensus.psmcfa
psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o sample.psmc consensus.psmcfa
psmc_plot.pl -u 1.4e-8 -g 5 -p sample_psmc_plot sample.psmc

# --- Stairway Plot 2: SFS-based (works with RAD-seq or WGS) ---
# Build SFS from VCF (use easySFS or vcf2sfs)
# Create blueprint.txt specifying nseq, L, mu, generation time
java -cp stairway_plot_v2.jar Stairbuilder blueprint.txt
bash blueprint.sh

# --- fastsimcoal2: multi-population composite-likelihood SFS ---
# REQUIRES >= 50 independent optimization replicates; single-run is unreliable
for i in {1..50}; do
    fsc27 -t template.tpl -e template.est -m -L 50 -n 100000 -q
done
# Compare likelihoods across replicates; report best AND distribution
```

## Genetic Load — Realized vs Masked Decomposition

**Goal:** Decompose total genetic load into realized (currently expressed) and masked (heterozygous, potential) components per Bertorelle 2022 framework.

**Approach:** Use SnpEff/VEP to annotate variants for predicted functional effect. Compute realized load from homozygous-derived-allele counts at deleterious positions; compute masked load from heterozygous counts. Hedrick & Garcia-Dorado 2016 distinguish purging vs drift dynamics; report both.

```r
# Conceptual workflow (Bertorelle et al. 2022 NRG 23:492-503)
# Realized load: count_homozygous_deleterious / total_deleterious_sites
# Masked load: count_heterozygous_deleterious / total_deleterious_sites
# Total load = realized + 0.5 * masked (assuming partial dominance)

# Annotate VCF with SnpEff or VEP first to classify deleterious vs neutral
# Then per individual:
# realized_load <- sum(genotype == 2 & severity == 'HIGH') / sum(severity == 'HIGH')
# masked_load <- sum(genotype == 1 & severity == 'HIGH') / sum(severity == 'HIGH')

# Cite Hedrick & Garcia-Dorado 2016 TREE 31:940-952 for purging vs drift:
# Strong-s alleles purge under inbreeding (potentially good)
# Weak-s alleles fix by drift (definitely bad)
# Net effect depends on selection-coefficient distribution
# Empirical purging example: Robinson 2018 Curr Biol 28:3487-3494 (Channel Island foxes)
# Simulation-based load assessment: Kyriazis 2021 Evol Lett 5:33-47
```

## Forward Simulation with Tree-Sequence Recording

**Goal:** Simulate selection and non-Wright-Fisher demography efficiently using SLiM 4 with tskit tree-sequence recording.

**Approach:** SLiM 4 with `treeSeqOutput()` records the genealogy; `pyslim` + `tskit` adds neutral mutations a posteriori. This is 5-100x faster than mutation-tracking (Haller, Galloway, Kelleher, Messer, Ralph 2019 *Mol Ecol Resour* 19:552-566). The Eidos scripting language is NOT Python.

```python
# Reference: SLiM 4+, pyslim 1.0+, tskit 0.5+, msprime 1.3+
# Note: SLiM scripts (.slim files) use Eidos, not Python.
# See https://messerlab.org/slim/ for full Eidos syntax.

# After running a SLiM simulation with treeSeqOutput():
import tskit
import pyslim
import msprime

ts = tskit.load('simulation.trees')
ts = pyslim.update(ts)  # update to current SLiM/pyslim conventions

# Recapitate: add coalescence above the SLiM tree root (neutral burn-in)
ts_recap = pyslim.recapitate(ts, recombination_rate=1e-8,
                              ancestral_Ne=10000, random_seed=42)

# Add neutral mutations after-the-fact via msprime
ts_mut = msprime.sim_mutations(ts_recap, rate=1e-8, random_seed=42)

# Export VCF for downstream analyses
with open('simulation.vcf', 'w') as f:
    ts_mut.write_vcf(f)
```

## Per-Method Failure Modes

### LDNe applied to RAD-seq SNPs without physical-linkage correction

**Trigger:** Running NeEstimator V2 LD method on RAD-seq or WGS genotypes thinned only by MAF, with multiple SNPs per chromosome at close physical distance. For RAD-seq biases more broadly, see Andrews 2016 *Nat Rev Genet* 17:81-92.

**Mechanism:** LDNe assumes inter-locus r^2 reflects only demographic LD (drift, random mating). With physically linked SNPs, much r^2 is chromosomal, NOT demographic. The estimator interprets this as "more drift has happened" and reports a smaller Ne.

**Symptom:** LDNe estimate dramatically lower than ecologically reasonable; downward bias compared with temporal or other estimators on the same population.

**Fix:** Thin SNPs to >= 1 cM apart OR use SNeP (Barbato 2015) which corrects for physical linkage explicitly. NeEstimator V2 has a `Chrom` flag for chromosome-aware LDNe; document its use.

### GONE2 silent failure with PLINK MAP cM=0 (the default)

**Trigger:** Running GONE/GONE2 on PLINK output where the MAP file's cM column was never populated (defaults to 0 because PLINK uses physical position).

**Mechanism:** GONE/GONE2 use the cM positions to interpret recombination distances; cM=0 for all loci means the algorithm cannot distinguish linked from unlinked SNPs and produces nonsense.

**Symptom:** GONE2 output shows extreme Ne values or completely flat trajectory; results inconsistent across runs or with other methods.

**Fix:** Populate the cM column from a species-specific genetic map, OR build one from Hi-C data, OR use the linkage-rate approximation from PLINK `--cm-map`. Verify with `head genotypes.bim` — column 3 should not be all zeros.

### NeEstimator silently fails with malformed option file

**Trigger:** Modifying line order in the `.ne2` option file or omitting expected parameters.

**Mechanism:** NeEstimator V2 reads the option file by line order with rigid parsing. Reordered lines or skipped parameters cause the program to misinterpret subsequent lines without raising clear errors.

**Symptom:** "INFO" output line shows fewer methods running than expected; results file is empty or has impossible values.

**Fix:** Strictly follow the documented option-file format. The INFO line confirms which methods actually ran. Cross-reference with the official documentation (https://github.com/bunop/NeEstimator2.X).

### dadi/fastsimcoal2 single-replicate inference at local optimum

**Trigger:** Running dadi or fastsimcoal2 with only 1-5 optimization replicates and reporting the result.

**Mechanism:** Likelihood surfaces for demographic inference have local optima. A single optimization may converge to a local maximum that is much worse than the global maximum.

**Symptom:** Reviewer asks "how many replicates?"; parameter estimates inconsistent across studies; demographic events implausible.

**Fix:** Run >= 50 independent replicates with random starting parameters; report the best-likelihood replicate AND the spread; flag if the top replicates disagree.

### Ne/Nc = 0.1 assumption applied to fish or other high-fecundity taxon

**Trigger:** Converting Ne to Nc via ratio = 0.1 for a marine fish, broadcast spawner, or other species with high reproductive skew.

**Mechanism:** Frankham 1995 median Ne/Nc = 0.1 was derived from primarily terrestrial vertebrates. Hauser & Carvalho 2008 documented Ne/Nc spanning 2-6 orders of magnitude smaller than census (10^-2 to 10^-6) in marine fish with sex-biased and sweepstakes recruitment. The 0.1 default produces wildly wrong Nc estimates.

**Symptom:** Census size estimate disagrees with field observations by orders of magnitude.

**Fix:** Use taxon-specific Ne/Nc ratios from the literature. For marine fish with sweepstakes recruitment, Hauser & Carvalho 2008 documented Ne/Nc spanning 2-6 orders of magnitude smaller than census (10^-2 to 10^-6). For long-lived mammals, 0.1 - 0.3 is typical. Cite Hauser & Carvalho 2008 when discussing Ne/Nc.

## Quantitative Thresholds

| Threshold | Value | Source / rationale |
|-----------|-------|-------------------|
| Modern Ne for inbreeding protection | Ne >= 100 | Frankham 2014 Biol Conserv 170:56-63 (revised from 50) |
| Modern Ne for adaptive maintenance | Ne >= 1000 | Frankham 2014 (revised from 500) |
| F_ROH for very recent inbreeding | ROH > 16 Mb | Parents/grandparents-scale; standard ROH-length-class convention |
| F_ROH for historical bottleneck | 1-4 Mb ROH | ~10-20 generations back |
| LDNe Pcrit (allele frequency cutoff) | 0.02 standard; 0.05 conservative; 0.01 sensitive | NeEstimator V2 convention |
| GONE2 minimum SNPs | 10,000 | Reliable trajectory inference |
| GONE2 minimum N | 50 diploids | Statistical power floor |
| PSMC minimum coverage | >= 20x WGS | Single-genome SMC requires high coverage |
| dadi / fsc2 optimization replicates | >= 50 | Local-optima trap below this |
| Ne/Nc default (use cautiously) | 0.1 vertebrate median | Hauser & Carvalho 2008 documented 2-6 orders of magnitude variation in marine fish (10^-2 to 10^-6) |
| SLiM forward-sim with tree sequences | treeSeqOutput() mandatory for speed | 5-100x faster than mutation tracking |

## Common errors

| Error | Cause | Solution |
|-------|-------|----------|
| LDNe reports Ne = infinity | Population too large for method to detect drift | Try multiple Pcrit; switch to SNeP for genomic data; acknowledge "Ne very large" |
| GONE2 nonsense trajectory | cM column = 0 in BIM/MAP file (PLINK default) | Populate cM with species genetic map or Hi-C inference |
| NeEstimator silent skip | Malformed option file or wrong method index | Check INFO line in output; align with documented format |
| bcftools roh empty output | Missing -G30 flag or --AF-tag | Add `-G30 --AF-tag AF` |
| detectRUNS missing F_ROH | Wrong mapFile path or chromosome naming mismatch | Verify map file format matches PED file |
| dadi numerical error at recent times | Mixing float32 in SFS construction | Use numpy float64 explicitly |
| SLiM script error "Eidos not Python" | Trying to use Python syntax | SLiM uses Eidos; consult SLiM manual |
| msprime PopulationConfiguration deprecation | Old msprime <1.0 API in newer install | Use `msprime.Demography()` constructor |

## References

- Waples RS (2006) A bias correction for estimates of effective population size based on linkage disequilibrium. *Conserv Genet* 7(2):167-184. doi:10.1007/s10592-005-9100-y
- Do C, Waples RS, Peel D, Macbeth GM, Tillett BJ, Ovenden JR (2014) NeEstimator V2. *Mol Ecol Resour* 14(1):209-214. doi:10.1111/1755-0998.12157
- Santiago E, Novo I, Pardiñas AF, Saura M, Wang J, Caballero A (2020) Recent demographic history inferred by high-resolution analysis of linkage disequilibrium (GONE). *Mol Biol Evol* 37(12):3642-3653. doi:10.1093/molbev/msaa169
- Barbato M, Orozco-terWengel P, Tapio M, Bruford MW (2015) SNeP: physical-linkage-corrected LDNe. *Front Genet* 6:109. doi:10.3389/fgene.2015.00109
- Frankham R, Bradshaw CJA, Brook BW (2014) Revised 100/1000 Ne rules. *Biol Conserv* 170:56-63. doi:10.1016/j.biocon.2013.12.036
- Hauser L, Carvalho GR (2008) Paradigm shifts in marine fisheries genetics. *Fish Fish* 9(4):333-362. doi:10.1111/j.1467-2979.2008.00299.x
- Hedrick PW, García-Dorado A (2016) Understanding inbreeding depression, purging, and genetic rescue. *Trends Ecol Evol* 31(12):940-952. doi:10.1016/j.tree.2016.09.005
- Bertorelle G, Raffini F, Bosse M, Bortoluzzi C, Iannucci A, Trucchi E, Morales HE, van Oosterhout C (2022) Genetic load. *Nat Rev Genet* 23(8):492-503. doi:10.1038/s41576-022-00448-x
- Gutenkunst RN, Hernandez RD, Williamson SH, Bustamante CD (2009) Inferring joint demographic history with dadi. *PLoS Genet* 5(10):e1000695. doi:10.1371/journal.pgen.1000695
- Excoffier L, Dupanloup I, Huerta-Sánchez E, Sousa VC, Foll M (2013) Robust demographic inference (fastsimcoal2). *PLoS Genet* 9(10):e1003905. doi:10.1371/journal.pgen.1003905
- Jouganous J, Long W, Ragsdale AP, Gravel S (2017) Inferring the joint demographic history of multiple populations: beyond the diffusion approximation (moments). *Genetics* 206(3):1549-1567. doi:10.1534/genetics.117.200493
- Kamm J, Terhorst J, Durbin R, Song YS (2020) Efficiently inferring the demographic history of many populations with allele count data (momi2). *J Am Stat Assoc* 115(531):1472-1487. doi:10.1080/01621459.2019.1635482
- Kelleher J, Etheridge AM, McVean G (2016) Efficient coalescent simulation (msprime). *PLoS Comput Biol* 12(5):e1004842. doi:10.1371/journal.pcbi.1004842
- Haller BC, Messer PW (2019) SLiM 3: forward genetic simulation beyond Wright-Fisher. *Mol Biol Evol* 36(3):632-637. doi:10.1093/molbev/msy228
- Haller BC, Galloway J, Kelleher J, Messer PW, Ralph PL (2019) Tree-sequence recording in SLiM. *Mol Ecol Resour* 19(2):552-566. doi:10.1111/1755-0998.12968
- Andrews KR, Good JM, Miller MR, Luikart G, Hohenlohe PA (2016) Harnessing RADseq for population genomics. *Nat Rev Genet* 17(2):81-92. doi:10.1038/nrg.2015.28
- Frankham R (2015) Genetic rescue meta-analysis. *Mol Ecol* 24(11):2610-2618. doi:10.1111/mec.13139
- Frankham R (1995) Effective population size / adult population size ratios in wildlife. *Genet Res* 66:95-107. doi:10.1017/S0016672300034455
- Moritz C (1994) Defining 'Evolutionarily Significant Units' for conservation. *Trends Ecol Evol* 9(10):373-375. doi:10.1016/0169-5347(94)90057-4
- Robinson JA, Brown C, Kim BY, Lohmueller KE, Wayne RK (2018) Purging of strongly deleterious mutations explains long-term persistence and absence of inbreeding depression in island foxes. *Curr Biol* 28(21):3487-3494.e4. doi:10.1016/j.cub.2018.08.066
- Kyriazis CC, Wayne RK, Lohmueller KE (2021) Strongly deleterious mutations are a primary determinant of extinction risk due to inbreeding depression. *Evol Lett* 5(1):33-47. doi:10.1002/evl3.209

## Related Skills

- ecological-genomics/landscape-genomics - Adaptive variation and genotype-environment associations
- ecological-genomics/species-delimitation - Taxonomic unit definition (cite Sukumaran-Knowles caveat for ESU/MU vs species)
- population-genetics/population-structure - Population stratification and STRUCTURE/ADMIXTURE
- population-genetics/selection-statistics - Genome-wide selection signatures
- variant-calling/vcf-basics - VCF preparation from RAD-seq or WGS
