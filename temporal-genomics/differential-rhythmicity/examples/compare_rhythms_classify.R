# Reference: compareRhythms 1.0+ | Verify API if version differs
# compareRhythms classifies each feature's rhythm change between two groups DIRECTLY into
# gain / loss / change / same, built to replace the detect-then-Venn anti-pattern (which
# overestimates reprogramming because two independently thresholded lists differ mostly
# from threshold noise near p~0.05). Simulated two-condition data; method='mod_sel' (BIC
# model selection) needs no DE package, unlike the 'limma'/'voom'/'deseq2'/'edger' methods.

suppressPackageStartupMessages(library(compareRhythms))

set.seed(1)

period <- 24
zt <- rep(seq(0, 44, by = 4), each = 2)   # 12 timepoints x 2 replicates over ~2 full cycles

# exp_design: one row per sample, a numeric 'time' column and a 2-level factor 'group'.
exp_design <- data.frame(
  time = rep(zt, times = 2),
  group = factor(rep(c('WT', 'KO'), each = length(zt)), levels = c('WT', 'KO')))  # WT is the reference level

sim <- function(amp_wt, amp_ko, ph_wt, ph_ko, noise = 0.3) {
  amp <- ifelse(exp_design$group == 'WT', amp_wt, amp_ko)
  ph <- ifelse(exp_design$group == 'WT', ph_wt, ph_ko)
  8 + amp * cos(2 * pi * (exp_design$time - ph) / period) + rnorm(nrow(exp_design), 0, noise)
}

# Rows = features (rownames become the reported ids); columns = samples matching exp_design rows.
data <- rbind(
  same     = sim(1.0, 1.0, 6, 6),    # unchanged rhythm
  loss     = sim(1.0, 0.0, 6, 6),    # loss of rhythm in KO
  gain     = sim(0.0, 1.0, 6, 6),    # gain of rhythm in KO
  change   = sim(1.2, 0.3, 6, 6),    # amplitude change in KO
  arrhythm = sim(0.0, 0.0, 6, 6))    # rhythmic in neither
colnames(data) <- rownames(exp_design) <- paste0('s', seq_len(nrow(exp_design)))

# amp_cutoff = peak-to-trough amplitude floor: a feature must clear it in >=1 group to enter DR results.
# criterion='bic' penalizes model size more than aic, favoring the simpler (parsimonious) DR category.
res <- compareRhythms(data, exp_design = exp_design, period = period,
                      method = 'mod_sel', amp_cutoff = 0.5, criterion = 'bic')

cat('Differential-rhythmicity classification (relative to WT reference):\n')
print(res)   # data.frame: id + category (gain / loss / change / same)
