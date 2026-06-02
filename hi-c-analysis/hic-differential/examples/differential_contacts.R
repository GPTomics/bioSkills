# Differential bin-pair Hi-C contacts: SCC reproducibility gate, then
# distance-stratified between-sample normalization and a replicate-aware edgeR test.
# Reference: multiHiCcompare 1.20+, hicrep (Bioconductor) 1.26+ | Verify API if version differs

library(multiHiCcompare)
library(hicrep)

# --- 0. Inputs ---------------------------------------------------------------
# Each replicate is a sparse upper-triangular table: chr, region1(bp), region2(bp), IF.
# chr coded 1-22, 23=X, 24=Y. Replace these reads with your own per-replicate tables.
c1_r1 <- read.table('cond1_rep1.tsv', header = FALSE, col.names = c('chr', 'region1', 'region2', 'IF'))
c1_r2 <- read.table('cond1_rep2.tsv', header = FALSE, col.names = c('chr', 'region1', 'region2', 'IF'))
c2_r1 <- read.table('cond2_rep1.tsv', header = FALSE, col.names = c('chr', 'region1', 'region2', 'IF'))
c2_r2 <- read.table('cond2_rep2.tsv', header = FALSE, col.names = c('chr', 'region1', 'region2', 'IF'))

RESOL <- 50000          # bin size (bp); must match the input tables
MAX_DIST <- 5000000     # SCC max interaction distance: covers compartment/TAD scale, caps sparse long-range
ZERO_P <- 0.8           # drop bin-pairs >80% zero across samples: sparse pairs break NB and waste FDR
A_MIN <- 5              # mean-IF independent filter: contrast-independent, preserves FDR validity
PADJ_CUT <- 0.1         # genome-wide multiple testing over millions of bin-pairs needs FDR control
LOGFC_CUT <- 1         # report >=2-fold changes

# --- 1. Replicate QC gate: HiCRep SCC ---------------------------------------
# Plain Pearson on Hi-C always looks reproducible (shared P(s) decay inflates r>0.9 even
# between unrelated maps). SCC stratifies by distance and smooths for sparsity.
# Current Bioconductor signature: get.scc(dat, resol, max) on a 4-col (mid1, mid2, IF_A, IF_B)
# table; the older TaoYang-dev signature is get.scc(mat1, mat2, resol, h, lbr, ubr). Verify with ?get.scc.
# The positional cbind below assumes both replicates share identical bin-pair support (true here
# because the demo tables are row-aligned). On real sparse pairs, merge on (region1, region2) first.
within_dat  <- data.frame(c1_r1$region1, c1_r1$region2, c1_r1$IF, c1_r2$IF)
between_dat <- data.frame(c1_r1$region1, c1_r1$region2, c1_r1$IF, c2_r1$IF)
scc_within  <- get.scc(within_dat,  resol = RESOL, max = MAX_DIST)$scc
scc_between <- get.scc(between_dat, resol = RESOL, max = MAX_DIST)$scc
cat(sprintf('SCC within-condition: %.3f  between-condition: %.3f\n', scc_within, scc_between))
stopifnot(scc_within > scc_between)   # if within does not exceed between, differences are replicate noise

# --- 2. Build the experiment and normalize on the M-D plot ------------------
hicexp <- make_hicexp(c1_r1, c1_r2, c2_r1, c2_r2,
                      groups = c(0, 0, 1, 1),
                      zero.p = ZERO_P, A.min = A_MIN, filter = TRUE)
hicexp <- cyclic_loess(hicexp, span = NA)   # span=NA -> GCV chooses the per-stratum loess span

# --- 3. Replicate-aware test and differential bin-pairs ---------------------
hicexp <- hic_exactTest(hicexp)             # 2-group; use hic_glm(hicexp, design) for covariates
res <- results(hicexp)                      # chr, region1, region2, D, logFC, logCPM, p.value, p.adj
sig <- topDirs(hicexp, logfc_cutoff = LOGFC_CUT, logcpm_cutoff = 1,
               p.adj_cutoff = PADJ_CUT, return_df = 'pairedbed')

cat(sprintf('Tested bin-pairs: %d  Significant (p.adj<%.2f, |logFC|>%.0f): %d\n',
            nrow(res), PADJ_CUT, LOGFC_CUT, nrow(sig)))
write.table(sig, 'differential_contacts.bedpe', sep = '\t', quote = FALSE, row.names = FALSE)

# --- 4. Diagnostic: M should center on 0 at every distance D ----------------
# A residual M-trend at large D means the between-sample normalization failed at long range.
pdf('md_diagnostic.pdf'); MD_hicexp(hicexp); dev.off()
