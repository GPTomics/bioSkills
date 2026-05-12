# Reference: TwoSampleMR 0.5.11+ | Verify API if version differs
## Bidirectional MR with Steiger pre-filter
##
## Forward MR: instrument exposure E -> test effect on outcome Y
## Reverse MR: instrument outcome Y -> test effect on exposure E
## Significant forward + null reverse strengthens the forward causal claim.
## A bidirectionally significant pair flags ambiguity (feedback loop, shared
## confounder, or genuine reciprocal causation).
##
## Steiger pre-filter (Hemani 2017 PLoS Genet 13:e1007081) removes SNPs whose
## per-SNP r^2 with the outcome exceeds r^2 with the exposure, ie SNPs that
## are likely outcome-direction instruments contaminating the forward set.

library(TwoSampleMR)

## Forward: exposure -> outcome
dat_fwd <- harmonise_data(exposure_dat = exposure_instruments,
                          outcome_dat  = outcome_for_forward)
dat_fwd <- steiger_filtering(dat_fwd)
dat_fwd_pass <- dat_fwd[dat_fwd$steiger_dir == TRUE, ]

res_fwd <- mr(dat_fwd_pass,
              method_list = c('mr_ivw', 'mr_egger_regression',
                              'mr_weighted_median', 'mr_weighted_mode'))
dir_fwd <- directionality_test(dat_fwd_pass)

cat('=== Forward MR (E -> Y) ===\n')
print(res_fwd[, c('method', 'nsnp', 'b', 'se', 'pval')])
cat('Steiger correct direction:', dir_fwd$correct_causal_direction, '\n')
cat('SNPs retained after Steiger:', nrow(dat_fwd_pass), '/', nrow(dat_fwd), '\n')

## Reverse: outcome -> exposure
dat_rev <- harmonise_data(exposure_dat = outcome_instruments,
                          outcome_dat  = exposure_for_reverse)
dat_rev <- steiger_filtering(dat_rev)
dat_rev_pass <- dat_rev[dat_rev$steiger_dir == TRUE, ]

res_rev <- mr(dat_rev_pass,
              method_list = c('mr_ivw', 'mr_egger_regression',
                              'mr_weighted_median', 'mr_weighted_mode'))

cat('\n=== Reverse MR (Y -> E) ===\n')
print(res_rev[, c('method', 'nsnp', 'b', 'se', 'pval')])

cat('\n=== Interpretation ===\n')
fwd_sig <- res_fwd$pval[res_fwd$method == 'Inverse variance weighted'] < 0.05
rev_sig <- res_rev$pval[res_rev$method == 'Inverse variance weighted'] < 0.05
if (fwd_sig && !rev_sig) {
    cat('Forward significant, reverse null: forward causal claim strengthened\n')
} else if (fwd_sig && rev_sig) {
    cat('Bidirectional significance: feedback / shared confounder; run LHC-MR\n')
} else if (!fwd_sig && rev_sig) {
    cat('Reverse-only significance: re-examine instrument-exposure direction\n')
} else {
    cat('Neither direction significant\n')
}
