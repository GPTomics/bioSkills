---
name: bio-ecological-genomics-landscape-genomics
description: Tests genotype-environment associations and identifies loci under local adaptation using LFMM2 (LEA), pcadapt outlier detection, OutFLANK Fst-based selection scans, and redundancy analysis. Detects adaptive genetic variation correlated with environmental variables while controlling for population structure. Use when identifying adaptive loci across environmental gradients, testing for signatures of local adaptation, or predicting genetic vulnerability to climate change with gradientForest.
tool_type: r
primary_tool: LEA
---

# Landscape Genomics

Identifies loci under local adaptation by testing genotype-environment associations while controlling for population structure.

## Population Structure Estimation with LEA

Determine K (number of ancestral populations) before running GEA:

```r
library(LEA)

# Convert VCF to lfmm/geno format
vcf2lfmm('variants.vcf', 'genotypes.lfmm')
vcf2geno('variants.vcf', 'genotypes.geno')

# sNMF for K estimation (faster than STRUCTURE)
snmf_result <- snmf('genotypes.geno', K = 1:10, repetitions = 5,
                     entropy = TRUE, project = 'new')

# Select K with minimum cross-entropy across K values
# cross.entropy(obj, K) returns per-run values for that K; take min per K
ce_values <- sapply(1:10, function(k) min(cross.entropy(snmf_result, K = k)))
plot(1:10, ce_values, xlab = 'K', ylab = 'Cross-entropy', pch = 19, col = 'blue')
best_K <- which.min(ce_values)
```

## LFMM2 Genotype-Environment Association

Latent factor mixed model tests association between each locus and environmental variables while correcting for population structure via latent factors:

```r
library(LEA)

genotypes <- read.lfmm('genotypes.lfmm')
env_vars <- read.env('environment.env')

# K from sNMM cross-entropy (latent factors = population structure proxy)
# K=3: determined from sNMF analysis above
lfmm_result <- lfmm2(input = genotypes, env = env_vars, K = 3)

# Extract p-values for each environmental variable
# Genomic inflation factor (GIF) calibration corrects for residual confounding
# Target lambda ~ 1.0; >1.5 suggests insufficient structure correction
pvalues <- lfmm2.test(lfmm_result, input = genotypes, env = env_vars,
                       full = TRUE, genomic.control = TRUE)

gif <- median(qchisq(1 - pvalues$pvalues[, 1], df = 1)) / qchisq(0.5, df = 1)
cat('Genomic inflation factor (lambda):', gif, '\n')
# lambda should be ~1.0; if >1.5, increase K; if <0.5, decrease K

# Adjust p-values for multiple testing
# q-value < 0.05: Storey FDR control, more powerful than BH for genomic data
library(qvalue)
qvals <- qvalue(pvalues$pvalues[, 1])$qvalues

# Candidate loci under selection
# q < 0.05: 5% false discovery rate
candidates <- which(qvals < 0.05)
cat('Candidate adaptive loci:', length(candidates), '\n')
```

### Multiple Environmental Variables

```r
# env_vars matrix: samples x variables (e.g., temperature, precipitation, altitude)
# LFMM2 tests each variable separately but fits latent factors jointly
for (i in 1:ncol(env_vars)) {
    pvals <- lfmm2.test(lfmm_result, input = genotypes, env = env_vars,
                         full = TRUE, genomic.control = TRUE)
    qvals_i <- qvalue(pvals$pvalues[, i])$qvalues
    cat('Variable', i, '- candidates (q<0.05):', sum(qvals_i < 0.05), '\n')
}
```

## pcadapt Outlier Detection

Detects loci under selection based on Mahalanobis distances from principal components, without requiring environmental data:

```r
library(pcadapt)

# Read genotypes (bed/bim/fam or lfmm format)
geno <- read.pcadapt('genotypes.bed', type = 'bed')

# Select K from screeplot (number of significant PCs)
# Look for elbow where eigenvalues flatten
x_scree <- pcadapt(geno, K = 20)
plot(x_scree, option = 'screeplot')

# Run pcadapt with selected K
# K=3: number of principal components capturing population structure
x <- pcadapt(geno, K = 3)

# QQ-plot to check calibration (should follow diagonal except tail)
plot(x, option = 'qqplot')

# Manhattan-style plot of -log10(p-values)
plot(x, option = 'manhattan')

# Identify outliers using q-values
# q < 0.05: 5% FDR threshold
library(qvalue)
qvals <- qvalue(x$pvalues)$qvalues
outliers <- which(qvals < 0.05)
cat('pcadapt outlier loci:', length(outliers), '\n')
```

## OutFLANK Fst-Based Selection Scan

Identifies Fst outliers by fitting a chi-squared distribution to the trimmed Fst distribution:

```r
library(OutFLANK)

# fst_data: matrix with columns FSTNoCorr, T1, T2, indexOrder, LocusName, He
# From hierfstat or custom calculation
fst_mat <- MakeDiploidFSTMat(SNPmat = genotype_matrix,
                              locusNames = snp_names,
                              popNames = pop_assignments)

# Trim extremes to estimate neutral Fst distribution
# LeftTrimFraction=0.05, RightTrimFraction=0.05: standard 5% trim on each tail
# Hmin=0.1: exclude low-heterozygosity loci (unreliable Fst estimates)
outflank_result <- OutFLANK(fst_mat, LeftTrimFraction = 0.05,
                             RightTrimFraction = 0.05, Hmin = 0.1,
                             NumberOfSamples = length(unique(pop_assignments)),
                             qthreshold = 0.05)

# Visualize Fst distribution and outliers
OutFLANKResultsPlotter(outflank_result, withOutliers = TRUE,
                        NoCorr = TRUE, Hmin = 0.1)

# Extract outlier loci
outlier_idx <- which(outflank_result$results$OutlierFlag == TRUE)
cat('OutFLANK Fst outliers:', length(outlier_idx), '\n')
```

## RDA-Based GEA

Redundancy analysis detects multilocus adaptive signatures correlated with environment:

```r
library(vegan)

# allele_freq: matrix of allele frequencies (individuals x loci)
# Impute missing values with population means
allele_freq[is.na(allele_freq)] <- apply(allele_freq, 2,
    function(x) mean(x, na.rm = TRUE))[col(allele_freq)[is.na(allele_freq)]]

# Condition on population structure (Q-matrix from sNMF)
# ~ env_vars + Condition(structure): partial RDA controlling for structure
rda_result <- rda(allele_freq ~ temperature + precipitation + altitude +
                  Condition(as.matrix(q_matrix)), data = env_data)

# Permutation test
anova(rda_result, permutations = 999)

# Identify outlier loci from RDA loadings
# z-score > 3: conservative threshold for candidate adaptive loci
loadings <- scores(rda_result, choices = 1:3, display = 'species')
zscores <- apply(loadings, 2, function(x) (x - mean(x)) / sd(x))
rda_candidates <- which(apply(abs(zscores), 1, max) > 3)
cat('RDA candidate loci:', length(rda_candidates), '\n')
```

## gradientForest: Allele Turnover Prediction

Predicts how allele frequencies shift along environmental gradients using random forests:

```r
library(gradientForest)
library(terra)

# allele_data: allele frequencies at candidate loci
# env_predictors: environmental variables at sampling sites
gf <- gradientForest(cbind(env_predictors, allele_data),
                     predictor.vars = colnames(env_predictors),
                     response.vars = colnames(allele_data),
                     ntree = 500, trace = FALSE)

# Variable importance: which environment drives most allele turnover
plot(gf, plot.type = 'Overall.Importance')

# Cumulative importance curves for each predictor
plot(gf, plot.type = 'Cumulative.Importance')

# Predict genetic offset under climate change
# current_rasters: current bioclimatic layers
# future_rasters: projected climate layers
current <- rast('bioclim_current.tif')
future <- rast('bioclim_2070_ssp585.tif')

current_vals <- extract(current, sampling_coords)
future_vals <- extract(future, sampling_coords)

# Genetic offset = Euclidean distance in transformed environmental space
current_transformed <- predict(gf, current_vals)
future_transformed <- predict(gf, future_vals)
genetic_offset <- sqrt(rowSums((current_transformed - future_transformed)^2))
```

## Environmental Data Extraction with terra

```r
library(terra)

# WorldClim bioclimatic variables (1 km resolution)
bioclim <- rast('wc2.1_30s_bio.tif')

# Extract values at sampling coordinates
# coords: data.frame with longitude and latitude columns
coords_vect <- vect(coords, geom = c('longitude', 'latitude'), crs = 'EPSG:4326')
env_values <- extract(bioclim, coords_vect)

# Common bioclimatic variables for GEA:
# BIO1: mean annual temperature
# BIO4: temperature seasonality
# BIO12: annual precipitation
# BIO15: precipitation seasonality
```

## Related Skills

- conservation-genetics - Population genetic health assessment
- community-ecology - Environmental gradient analysis for species
- population-genetics/selection-statistics - Selection scans in human data
- population-genetics/population-structure - Population structure inference
- variant-calling/vcf-basics - VCF input preparation
