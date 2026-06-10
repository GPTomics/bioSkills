# Reference: dmrseq 1.22+, bsseq 1.38+, DSS 2.50+ | Verify API if version differs
# HEADLINE selection-aware DMR call: dmrseq builds its null by PERMUTING condition
# labels and RE-RUNNING region detection, so its qval is a region FDR that accounts
# for the region-selection step (post-selection inference). DSS callDMR corroborates
# (beta-binomial dispersion shrinkage). Cross-tool OVERLAP is the evidence statement;
# the two q-values control different objects and are not directly comparable.
library(dmrseq)
library(bsseq)
library(DSS)
library(GenomicRanges)

cov_files <- c('ctrl1.cov.gz', 'ctrl2.cov.gz', 'ctrl3.cov.gz',
               'treat1.cov.gz', 'treat2.cov.gz', 'treat3.cov.gz')
condition <- c('ctrl', 'ctrl', 'ctrl', 'treat', 'treat', 'treat')

bs <- read.bismark(cov_files, colData = DataFrame(condition = condition),
                   rmZeroCov = TRUE, strandCollapse = TRUE)

# dmrseq requires non-zero coverage in EVERY sample at every retained locus.
bs <- bs[rowSums(getCoverage(bs) == 0) == 0, ]

# cutoff=0.1 only SEEDS candidate detection (10% smoothed difference); significance is
# the permutation statistic, so dmrseq does NOT hard-threshold delta-beta. Do NOT run
# BSmooth() first - dmrseq smooths the methylation DIFFERENCE internally.
candidate_cutoff <- 0.1
dmrs <- dmrseq(bs, testCovariate = 'condition', cutoff = candidate_cutoff)

qval_floor <- 0.05   # REGION FDR (selection-aware), not a per-CpG FDR
sig_dmrseq <- dmrs[dmrs$qval < qval_floor]
sprintf('dmrseq: %d candidate regions, %d at qval<%.2f (selection-aware region FDR)',
        length(dmrs), length(sig_dmrseq), qval_floor)

build_dss <- function(cov_file) {
    d <- read.table(cov_file, col.names = c('chr', 'start', 'end', 'meth_pct', 'numC', 'numT'))
    data.frame(chr = d$chr, pos = d$start, N = d$numC + d$numT, X = d$numC)
}
dss_obj <- makeBSseqData(lapply(cov_files, build_dss), sample_ids <- sub('.cov.gz', '', cov_files))

dml <- DMLtest(dss_obj, group1 = sample_ids[condition == 'ctrl'],
               group2 = sample_ids[condition == 'treat'], smoothing = TRUE)

# callDMR delta DEFAULTS to 0 (no effect-size floor) - SET it so tiny shifts are not called.
dss_delta <- 0.1
dss_dmrs <- callDMR(dml, delta = dss_delta, p.threshold = 1e-5,
                    minlen = 50, minCG = 3, dis.merge = 100, pct.sig = 0.5)

dss_gr <- GRanges(dss_dmrs$chr, IRanges(dss_dmrs$start, dss_dmrs$end))
overlap <- length(subsetByOverlaps(sig_dmrseq, dss_gr))
sprintf('DSS callDMR: %d DMRs | dmrseq-DSS overlap: %d regions', nrow(dss_dmrs), overlap)

sig_dmrseq
