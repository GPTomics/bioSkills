# Reference: EpiDISH 2.18+, minfi 1.48+ | Verify API if version differs
# Reference-based cell-type deconvolution of a bulk methylation beta matrix,
# then using the fractions as EWAS covariates (drop one to avoid collinearity)
# and a cell-type-resolved interaction test with CellDMC.

library(EpiDISH)

data(centDHSbloodDMC.m)    # shipped 7-immune-cell adult whole-blood reference centroid

simulate_bulk <- function(ref, n_samples) {
    n_ct <- ncol(ref)
    fractions <- matrix(rgamma(n_samples * n_ct, shape = 2), nrow = n_samples)
    fractions <- fractions / rowSums(fractions)    # project onto the simplex (rows sum to 1)
    beta <- ref %*% t(fractions)
    noise <- matrix(rnorm(length(beta), sd = 0.02), nrow = nrow(beta))    # mild array noise
    beta <- pmin(pmax(beta + noise, 0), 1)    # betas are bounded to [0, 1]
    colnames(beta) <- paste0('sample', seq_len(n_samples))
    list(beta = beta, truth = fractions)
}

set.seed(1)
n_samples <- 40
sim <- simulate_bulk(centDHSbloodDMC.m, n_samples)

# RPC is the robust option (downweights noisy CpGs); estF is samples x cell types.
out <- epidish(beta.m = sim$beta, ref.m = centDHSbloodDMC.m, method = 'RPC')
fractions <- out$estF
round(head(fractions), 3)

recovery_r <- diag(cor(fractions, sim$truth))    # per-cell-type estimate-vs-truth correlation
round(recovery_r, 3)

# Use fractions as EWAS covariates: drop one cell type because all K sum to ~1 (collinear).
covariates <- fractions[, -ncol(fractions), drop = FALSE]

# Cell-type-resolved test: which cell type carries a phenotype-associated signal.
phenotype <- rbinom(n_samples, 1, 0.5)
celldmc_res <- CellDMC(beta.m = sim$beta, pheno.v = phenotype, frac.m = fractions)
table(celldmc_res$dmct[, 'DMC'])    # binary 0/1: is this CpG a DMC in any cell type (per-cell-type cols are -1/0/1)

cat('mean fraction per cell type (distrust csDM for types below ~0.05):\n')
round(colMeans(fractions), 3)
