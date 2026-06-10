# Reference: sva 3.50+, bacon 1.30+, limma 3.58+ | Verify API if version differs
# EWAS skeleton: build a confounder-aware model, soak up unwanted variation with SVA,
# run the per-CpG test on M-values, then correct genomic inflation with BACON and
# monitor the smoking positive control through the whole pipeline.
# The per-site test mechanics live in differential-cpg-testing; this script is the
# DESIGN-and-inference layer (covariates, SVA, lambda/BACON, threshold, positive control).

library(sva)
library(limma)
library(bacon)

ahrr_probe <- 'cg05575921'   # smoking positive control; must survive correction (Joehanes 2016)
fwer_threshold <- 9e-8       # EPIC genome-wide FWER (Mansell 2019); use 2.4e-7 for 450K (Saffari 2018)

# --- synthetic stand-in: M-values (CpGs x samples) + a sample sheet. Replace with real data. ---
set.seed(1)
n <- 120
n_cpg <- 2000
pheno <- factor(rep(c('control', 'case'), each = n / 2), levels = c('control', 'case'))   # control is the reference -> coef 'phenocase'
age <- round(rnorm(n, 55, 8))
sex <- factor(sample(c('F', 'M'), n, replace = TRUE))
cell_gran <- pmin(pmax(rnorm(n, 0.6, 0.1), 0), 1)   # granulocyte fraction (dominant confounder)
smoker <- rbinom(n, 1, 0.4)

mvals <- matrix(rnorm(n_cpg * n, 0, 1), nrow = n_cpg, dimnames = list(paste0('cg', seq_len(n_cpg)), NULL))
rownames(mvals)[1] <- ahrr_probe
mvals[ahrr_probe, ] <- mvals[ahrr_probe, ] - 2.5 * smoker   # plant the smoking effect at AHRR
mvals <- mvals + matrix(rep(3 * cell_gran, each = n_cpg), nrow = n_cpg)   # composition leaks everywhere

sheet <- data.frame(pheno = pheno, age = age, sex = sex, cell_gran = cell_gran, smoker = smoker)

# --- full and null design: SVs are built orthogonal to the variable of interest ---
mod <- model.matrix(~ pheno + age + sex + cell_gran + smoker, data = sheet)
mod0 <- model.matrix(~ age + sex + cell_gran + smoker, data = sheet)

n_sv <- num.sv(mvals, mod, method = 'be')   # Buja-Eyuboglu estimate of how many SVs to keep
svobj <- sva(mvals, mod, mod0, n.sv = n_sv)
mod_sv <- cbind(mod, svobj$sv)

# --- per-CpG test on M-values with empirical-Bayes moderation ---
fit <- eBayes(lmFit(mvals, mod_sv))
tt <- topTable(fit, coef = 'phenocase', number = Inf, sort.by = 'none')
zscores <- sign(tt$logFC) * qnorm(tt$P.Value / 2, lower.tail = FALSE)

# --- genomic inflation BEFORE correction (report, do not divide by it) ---
chisq <- qchisq(tt$P.Value, df = 1, lower.tail = FALSE)
lambda <- median(chisq) / qchisq(0.5, df = 1)

# --- BACON: empirical-null bias AND inflation, then re-derive p-values ---
bc <- bacon(zscores)
cat(sprintf('lambda=%.3f  BACON bias=%.3f  BACON inflation=%.3f\n', lambda, bias(bc), inflation(bc)))
p_bacon <- pval(bc)

# --- positive-control check: the planted smoking signal must surface at AHRR ---
# AHRR (cg05575921) is a smoking marker, so it must appear in the SMOKER coefficient (where the
# effect was planted), not the case/control contrast. If it does not, the model or the data is broken.
tt_smoke <- topTable(fit, coef = 'smoker', number = Inf, sort.by = 'none')
ahrr_rank <- rank(tt_smoke$P.Value)[match(ahrr_probe, rownames(mvals))]
cat(sprintf('AHRR %s smoking p=%.2e (rank %.0f of %d) -- %s\n', ahrr_probe,
            tt_smoke$P.Value[match(ahrr_probe, rownames(mvals))], ahrr_rank, n_cpg,
            if (ahrr_rank <= 5) 'positive control recovered' else 'BROKEN: smoking signal not recovered at AHRR'))

# --- thresholds: FWER headline + BH-FDR discovery ---
fwer_hits <- sum(p_bacon < fwer_threshold)
fdr_hits <- sum(p.adjust(p_bacon, method = 'BH') < 0.05)
cat(sprintf('FWER hits (p<%.0e)=%d   BH-FDR hits (q<0.05)=%d\n', fwer_threshold, fwer_hits, fdr_hits))
