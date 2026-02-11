---
name: bio-ecological-genomics-conservation-genetics
description: Assesses genetic health of populations for conservation using effective population size estimation (GONE2 for recent Ne trajectory, NeEstimator for contemporary Ne, Stairway Plot 2 and PSMC for historical Ne), F-statistics (hierfstat), runs of homozygosity (detectRUNS), and genetic diversity metrics. Use when estimating effective population size, detecting inbreeding or bottlenecks, or assessing genetic diversity in threatened species from microsatellite or SNP data.
tool_type: mixed
primary_tool: hierfstat
---

# Conservation Genetics

Assesses genetic health of populations through diversity metrics, effective population size estimation, inbreeding detection, and demographic history reconstruction.

## Genetic Diversity with hierfstat

Basic population genetics statistics using Weir & Cockerham estimators:

```r
library(hierfstat)
library(adegenet)

# genind object from adegenet (VCF, PLINK, or genepop input)
# read.genepop() or vcfR::read.vcfR() + vcfR2genind()
data_genind <- read.genepop('populations.gen')

# Convert to hierfstat format
data_hf <- genind2hierfstat(data_genind)

# Basic F-statistics: Fis (inbreeding), Fst (differentiation), Ho, He
bstats <- basic.stats(data_hf)
cat('Overall Fis:', bstats$overall['Fis'], '\n')
cat('Overall Fst:', bstats$overall['Fst'], '\n')

# Per-population statistics
pop_stats <- data.frame(
    Ho = colMeans(bstats$Ho, na.rm = TRUE),
    He = colMeans(bstats$Hs, na.rm = TRUE),
    Fis = colMeans(bstats$Fis, na.rm = TRUE)
)
pop_stats
```

## Pairwise Fst

```r
# Bootstrap confidence intervals for pairwise Fst
# nboot=1000: standard for publication-quality CIs
pw_fst <- pairwise.WCfst(data_hf)
pw_fst

# Bootstrap p-values
boot_fst <- boot.ppfst(data_hf, nboot = 1000)
boot_fst$ll  # lower 95% CI
boot_fst$ul  # upper 95% CI
```

## Allelic Richness

Rarefied allelic richness corrects for unequal sample sizes:

```r
# min.n: rarefaction to smallest population sample size
ar <- allelic.richness(data_hf)
cat('Rarefied allelic richness per population:\n')
colMeans(ar$Ar, na.rm = TRUE)
```

## Private Alleles

```r
library(poppr)

# Private alleles unique to each population
pa <- private_alleles(data_genind, count.alleles = TRUE)
private_counts <- rowSums(pa > 0)
cat('Private alleles per population:\n')
print(private_counts)
```

## Runs of Homozygosity (ROH)

Detects autozygous segments indicating recent inbreeding:

```r
library(detectRUNS)

# Input: PLINK .ped/.map files
# consecutiveRuns: SNP-by-SNP scanning (preferred for SNP arrays/WGS)
runs <- consecutiveRuns(
    genotypeFile = 'genotypes.ped',
    mapFile = 'genotypes.map',
    minSNP = 20,        # min SNPs in a run; 20 prevents false positives from short homozygous stretches
    minLengthBps = 1e6, # 1 Mb minimum; shorter ROH are often not due to IBD
    maxGap = 1e6,       # max gap between consecutive SNPs within a run
    maxOppWindow = 1,   # max opposing homozygotes allowed (genotyping errors)
    maxMissWindow = 2   # max missing genotypes allowed per window
)

# Summary statistics
summary_runs <- summaryRuns(runs, genotypeFile = 'genotypes.ped',
                            mapFile = 'genotypes.map')

# F_ROH: inbreeding coefficient from ROH
# F_ROH = total ROH length / autosomal genome length
# <0.0625: background; 0.0625-0.125: offspring of half-cousins; >0.125: moderate inbreeding
genome_length_bp <- 2.5e9  # adjust to species genome size
froh <- Froh_inbreeding(runs, mapFile = 'genotypes.map',
                         genome_wide = TRUE, genomeLength = genome_length_bp)

cat('F_ROH per individual:\n')
print(froh)
```

### ROH Length Classes

| ROH Length | Approximate Ancestor | Interpretation |
|-----------|---------------------|----------------|
| > 16 Mb | Parents or grandparents | Very recent inbreeding |
| 4-16 Mb | ~5 generations back | Recent inbreeding |
| 1-4 Mb | ~10-20 generations back | Historical inbreeding/bottleneck |
| < 1 Mb | Deep background | Ancient homozygosity |

```r
# Plot ROH length distribution by population
plot_ROHclasses(runs, class_breaks = c(0, 1e6, 4e6, 16e6, Inf))
```

## Effective Population Size (Ne)

### GONE2: Recent Ne Trajectory (Linkage Disequilibrium)

Estimates Ne changes over the last ~200 generations from LD patterns:

```r
library(GONE2)

# Input: PLINK bed/bim/fam files
# Requires at least 10,000 SNPs and 50 individuals for reliable estimates
gone_result <- gone('genotypes', hc = 0.05, num_threads = 4)

# hc=0.05: maximum recombination distance in Morgans (~5 cM); default inter-SNP distance cutoff
# Smaller hc focuses on more recent generations

# Plot Ne trajectory
pdf('gone2_ne_trajectory.pdf', width = 8, height = 5)
plot(gone_result$generation, gone_result$Ne,
     type = 'l', lwd = 2, col = 'blue',
     xlab = 'Generations ago', ylab = 'Effective population size (Ne)',
     main = 'Recent Ne Trajectory (GONE2)', log = 'y')
dev.off()

cat('Current Ne estimate:', gone_result$Ne[1], '\n')
```

### NeEstimator: Contemporary Ne (LD Method)

```bash
# NeEstimator v2 command-line usage
# --method LD: linkage disequilibrium method (single time point)
# --pcrit 0.02: allele frequency threshold; excludes rare alleles
# 0.02 balances bias (lower pcrit) vs precision (higher pcrit)
NeEstimator --input genotypes.gen --method LD --pcrit 0.02 \
    --output ne_results.txt

# Parse results
# Reports point estimate and 95% CI (jackknife or parametric)
```

### Stairway Plot 2: Demographic History from SFS

Reconstructs Ne over thousands of generations using the site frequency spectrum:

```bash
# Step 1: Generate folded SFS from VCF
# Blueprint file specifies SFS, mutation rate, generation time
# nseq: number of haploid sequences (2 * n_diploid_individuals)
# L: total number of monomorphic + polymorphic sites analyzed
# mutation_rate: per-generation per-site (e.g., 1.4e-8 for vertebrates)
# year_per_generation: species generation time

# Step 2: Run Stairway Plot 2
java -cp stairway_plot_v2.jar Stairbuilder blueprint.txt

# Step 3: Execute the generated bash script
bash blueprint.sh

# Output: Ne trajectory plot and table
```

### PSMC: Whole-Genome Pairwise Coalescent

For a single diploid genome (whole-genome sequencing):

```bash
# Consensus FASTQ from BAM
bcftools mpileup -C50 -Q 30 -q 30 -f reference.fa sample.bam | \
    bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 > consensus.fq

# PSMC input format
fq2psmcfa -q20 consensus.fq > consensus.psmcfa

# Run PSMC
# -N25: number of iterations (25 is standard)
# -t15: maximum 2N0 coalescent time
# -r5: initial theta/rho ratio
# -p: time interval pattern (4+25*2+4+6: standard for most organisms)
psmc -N25 -t15 -r5 -p '4+25*2+4+6' -o sample.psmc consensus.psmcfa

# Bootstrap (100 replicates for confidence intervals)
seq 100 | xargs -P 4 -I {} \
    psmc -N25 -t15 -r5 -p '4+25*2+4+6' -b \
    -o sample_bootstrap_{}.psmc consensus.psmcfa

# Plot with species-specific mutation rate and generation time
# -u: mutation rate per generation per site
# -g: generation time in years
psmc_plot.pl -u 1.4e-8 -g 5 -p sample_psmc_plot sample.psmc
```

## Ne Interpretation Thresholds

| Ne Value | Conservation Status | Rationale |
|----------|-------------------|-----------|
| Ne < 50 | Critical | Inbreeding depression risk (50/500 rule, Franklin 1980) |
| Ne = 50-500 | Vulnerable | Sufficient to avoid inbreeding but limited adaptive potential |
| Ne > 500 | Viable | Adequate for long-term adaptive evolution |
| Ne/N ratio | ~0.1-0.3 typical | Ne is usually 10-30% of census size |

## Bottleneck Detection

```r
# Heterozygosity excess test (Cornuet & Luikart 1996)
# Under mutation-drift equilibrium, He ~ expected from allele number
# After bottleneck, alleles lost faster than He, creating excess
library(hierfstat)

# Compare observed He against equilibrium He expected from allele number
# After bottleneck, alleles lost faster than He, creating transient excess
# IAM (infinite alleles model) for SNPs; TPM preferred for microsatellites
bstats <- basic.stats(data_hf)
he_obs <- bstats$perloc$Hs

# Heq under IAM: expected He given observed allele count k per locus
# Heq_IAM = 1 - 1/k (for diploid data)
k_alleles <- sapply(data_hf[, -1], function(x) length(unique(na.omit(x))))
heq_iam <- 1 - 1 / k_alleles

# Wilcoxon signed-rank test: He > Heq indicates bottleneck
wilcox.test(he_obs, heq_iam, paired = TRUE, alternative = 'greater')
```

## Related Skills

- landscape-genomics - Adaptive variation and genotype-environment associations
- species-delimitation - Taxonomic unit definition for conservation
- population-genetics/population-structure - Population stratification
- population-genetics/selection-statistics - Selection signatures
- variant-calling/vcf-basics - VCF input from RADseq/WGS
