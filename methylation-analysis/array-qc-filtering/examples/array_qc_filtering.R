# Reference: minfi 1.48+, sesame 1.20+, maxprobes 0.0.2+, ChAMP 2.32+ | Verify API if version differs
#
# Probe filtering + sample-level QC + EPICv2 replicate collapse + cross-version harmonization
# for Illumina Infinium methylation arrays. Inputs come from array-preprocessing (a corrected,
# detection-masked GenomicRatioSet / RGChannelSet / beta matrix); this script decides which
# probes and samples to trust and how to merge across array versions. The minfi/sesame blocks
# assume those objects exist; the standalone simulation at the bottom runs end-to-end to
# demonstrate the rs-SNP identity-clustering and sex-mismatch logic without array packages.

# The minfi/sesame/maxprobes functions below require those Bioconductor packages; the standalone
# simulation at the bottom is base-R and runs without them, so the loads are attached lazily.
suppressWarnings(suppressMessages({require(minfi); require(maxprobes); require(sesame)}))

lo_detp <- 0.01      # detection-p above this = failed probe (masking applied upstream in array-preprocessing)
maf_drop <- 0        # 0 drops any annotated SNP at CpG/SBE; raise to keep rare variants

# --- Probe filtering on a corrected GenomicRatioSet (gset from array-preprocessing) ---
filter_probes <- function(gset, drop_sex = TRUE) {
    start_n <- nrow(gset)
    gset <- dropLociWithSnps(gset, snps = c('CpG', 'SBE'), maf = maf_drop)
    gset <- maxprobes::dropXreactiveLoci(gset)
    if (drop_sex) {
        anno <- getAnnotation(gset)
        gset <- gset[!(anno$chr %in% c('chrX', 'chrY')), ]
    }
    attr(gset, 'attrition') <- c(start = start_n, retained = nrow(gset))
    gset
}

# --- EPICv2 replicate collapse (sesame) ---
# betasCollapseToPfx averages replicate designs to one value per cg core ID (it takes betas only).
# To keep the best-detection replicate instead, collapse at the SigDF stage from IDATs:
# openSesame(idats, func = getBetas, collapseToPfx = TRUE, collapseMethod = 'minPval').
collapse_epicv2 <- function(betas) {
    betasCollapseToPfx(betas)
}

# --- Cross-version harmonization: collapse, intersect, liftover before merging ---
harmonize_versions <- function(betas_450k, betas_epic, betas_epicv2) {
    v2 <- collapse_epicv2(betas_epicv2)
    shared <- Reduce(intersect, list(rownames(betas_450k), rownames(betas_epic), rownames(v2)))
    # mLiftOver harmonizes EPICv2 (hg38) onto the 450K/EPIC (hg19) build before a coordinate merge:
    # v2 <- mLiftOver(v2, target_platform = 'HM450')
    cbind(betas_450k[shared, ], betas_epic[shared, ], v2[shared, ])
}

# --- Sample-level identity QC ---
sex_mismatches <- function(gmset, sheet_sex) {
    predicted <- getSex(gmset)$predictedSex
    which(predicted != sheet_sex)
}

fingerprint_clusters <- function(rgset) {
    snp_betas <- getSnpBeta(rgset)        # rs genotyping probes (65 on 450K, ~59 on EPIC); identity, independent of methylation
    hclust(dist(t(snp_betas)))
}

# --- Standalone runnable demo: identity clustering + sex mismatch on simulated rs-SNP genotypes ---
set.seed(1)
n_snp <- 59
donors <- replicate(6, sample(c(0, 0.5, 1), n_snp, replace = TRUE))   # 6 distinct genotypes
colnames(donors) <- paste0('donor', 1:6)

# Build 8 samples where sample8 is secretly a relabeled replicate of donor3 (a swap to catch).
sample_idx <- c(1, 2, 3, 4, 5, 6, 1, 3)
labels <- c(paste0('s', 1:6), 's7_dupOf_s1', 's8_swap_labeledNew')
snp_matrix <- donors[, sample_idx] + matrix(rnorm(n_snp * 8, sd = 0.02), n_snp, 8)
colnames(snp_matrix) <- labels

clusters <- hclust(dist(t(snp_matrix)))
identity_pairs <- cutree(clusters, h = 0.5)        # samples sharing a genotype fall in one cluster
duplicated_groups <- names(which(table(identity_pairs) > 1))
flagged <- identity_pairs[identity_pairs %in% duplicated_groups]
cat('rs-SNP identity clusters with >1 member (candidate swaps/duplicates):\n')
print(flagged)

predicted_sex <- c('M', 'F', 'M', 'F', 'M', 'F', 'M', 'F')
sheet_sex <- c('M', 'F', 'M', 'F', 'M', 'F', 'M', 'M')   # s8 sheet says M, intensity predicts F
sex_flag <- which(predicted_sex != sheet_sex)
cat('\nsamples whose predicted sex disagrees with the sample sheet:', labels[sex_flag], '\n')
