# Reference: MultiAssayExperiment 1.36+, sva 3.50+ | Verify API if version differs
# Cross-omic harmonization: assemble the container, equalize block contribution, correct batch
# ONCE per omic (only after a confounding check), and triage missingness by mechanism. Assumes
# each block is ALREADY per-omic normalized (VST / log2 / M-values) in its own category.

library(MultiAssayExperiment)
library(sva)

set.seed(1)
n <- 50
samples <- paste0('S', seq_len(n))
sample_info <- DataFrame(row.names=samples,
                         Condition=factor(rep(c('case', 'control'), length.out=n)),
                         Batch=factor(rep(c('b1', 'b2'), each=n / 2)))

vst_rna  <- matrix(rnorm(2000 * n, 8, 2), nrow=2000, dimnames=list(paste0('gene', 1:2000), samples))
norm_prot <- matrix(rnorm(300 * n, 20, 1), nrow=300, dimnames=list(paste0('prot', 1:300), samples))
m_values <- matrix(rnorm(5000 * n), nrow=5000, dimnames=list(paste0('cg', 1:5000), samples))

mae <- MultiAssayExperiment(experiments=ExperimentList(RNA=vst_rna, Protein=norm_prot, Methylation=m_values),
                            colData=sample_info)
cat('subjects with every omic:', sum(complete.cases(mae)), 'of', nrow(colData(mae)), '\n')

drop_constant <- function(mat, min_sd=1e-8) mat[apply(mat, 1, sd, na.rm=TRUE) > min_sd, ]   # per-feature scaling blows up zero-variance features
scale_per_view <- function(mat) mat / sqrt(sum(apply(mat, 1, var, na.rm=TRUE)))             # equalize block contribution, keep within-block feature ratios

blocks <- lapply(experiments(mae), function(x) scale_per_view(drop_constant(as.matrix(x))))
cat('per-block variance after scaling:\n')
print(round(sapply(blocks, function(x) sum(apply(x, 1, var))), 3))

cat('\nbatch x condition (confounding gate; any 0 cell = collinear, do NOT correct):\n')
print(with(as.data.frame(colData(mae)), table(Batch, Condition)))

mod <- model.matrix(~ Condition, data=as.data.frame(colData(mae)))                          # protect biology
rna_bc <- ComBat(dat=vst_rna, batch=colData(mae)$Batch, mod=mod, par.prior=TRUE)            # ONE omic at a time, never a stacked matrix
cat('\nRNA batch-corrected:', nrow(rna_bc), 'x', ncol(rna_bc), '\n')

prot_miss <- norm_prot
prot_miss[sample(length(prot_miss), 0.1 * length(prot_miss))] <- NA                         # simulate sporadic gaps
keep <- rowMeans(is.na(prot_miss)) < 0.30                                                    # drop features missing in >30% of samples
cat('proteins kept after missingness filter:', sum(keep), 'of', nrow(prot_miss), '\n')
