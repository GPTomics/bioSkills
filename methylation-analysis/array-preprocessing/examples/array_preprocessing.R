# Reference: sesame 1.20+, minfi 1.48+ | Verify API if version differs
# Preprocess Illumina Infinium methylation IDATs into a corrected beta/M matrix.
# Two interchangeable routes are shown: sesame (EPICv2-safe default) and minfi (450K/EPICv1).
# Run only the route matching the data; both expect a directory of raw _Grn.idat/_Red.idat pairs.

library(sesame)
library(minfi)

idat_dir <- 'idat_dir'

detp_cutoff <- 0.01     # detection p above this = signal indistinguishable from background, mask the probe
funnorm_pcs <- 2        # first 2 control-probe PCs absorb technical variation without erasing biology

logit_m <- function(betas) log2(betas / (1 - betas))

preprocess_sesame <- function(dir) {
    sesameDataCache()
    betas <- openSesame(dir, prep = 'QCDPB', func = getBetas)
    list(beta = betas, mval = logit_m(betas))
}

preprocess_minfi <- function(dir, global_differences = FALSE) {
    rg <- read.metharray.exp(base = dir)
    detP <- detectionP(rg)
    grSet <- if (global_differences) preprocessFunnorm(rg, nPCs = funnorm_pcs) else preprocessQuantile(rg)
    beta <- getBeta(grSet)      # GenomicRatioSet holds precomputed betas (offset already applied upstream)
    mval <- getM(grSet)
    common <- intersect(rownames(beta), rownames(detP))
    beta[common, ][detP[common, colnames(beta)] > detp_cutoff] <- NA
    list(beta = beta, mval = mval)
}

if (interactive()) {
    res <- preprocess_sesame(idat_dir)
    str(res$beta)
}
