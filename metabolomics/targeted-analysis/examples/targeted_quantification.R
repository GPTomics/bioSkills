# Reference: R 4.3+, ggplot2 3.5+ | Verify API if version differs
# Demonstrates the central targeted-quant claims on synthetic data:
#   weighted (1/x^2) calibration accepted by back-calculated %RE not R-squared,
#   LLOQ from accuracy, IS normalization, and ion-ratio confirmation.
library(ggplot2)

out_dir <- file.path(tempdir(), 'targeted_quant')
dir.create(out_dir, showWarnings = FALSE)

# ICH M10 acceptance: each calibrator within +-15% back-calc, +-20% at the LLOQ;
# SANTE/2020/12830 ion-ratio tolerance +-30% relative for LC-MS/MS.
re_tol_pct <- 15
re_tol_lloq_pct <- 20
ion_ratio_tol <- 0.30

# Synthetic standards spanning three orders of magnitude. Areas carry constant-CV
# (heteroscedastic) noise so absolute variance grows with concentration.
set.seed(7)
conc <- c(1, 5, 10, 25, 50, 100, 250, 500, 1000)
true_slope <- 0.0049
istd_area <- rnorm(length(conc), 1e5, 2e3)
analyte_area <- true_slope * conc * istd_area * rnorm(length(conc), 1, 0.04)
standards <- data.frame(conc, analyte_area, istd_area)
standards$ratio <- standards$analyte_area / standards$istd_area

# Fit unweighted vs 1/x^2-weighted and compare low-end accuracy, the whole point.
fit_ols <- lm(ratio ~ conc, data = standards)
fit_w <- lm(ratio ~ conc, data = standards, weights = 1 / standards$conc^2)

back_calc <- function(fit) (standards$ratio - coef(fit)[1]) / coef(fit)[2]
re_pct <- function(bc) (bc - standards$conc) / standards$conc * 100

standards$re_ols <- re_pct(back_calc(fit_ols))
standards$re_w <- re_pct(back_calc(fit_w))

cat('R-squared (unweighted):', round(summary(fit_ols)$r.squared, 4), '\n')
cat('Low-end %RE unweighted vs weighted (lowest 3 levels):\n')
print(round(head(standards[, c('conc', 're_ols', 're_w')], 3), 1))

# Accept the weighted curve per-level; LLOQ = lowest passing calibrator.
tol <- ifelse(standards$conc == min(standards$conc), re_tol_lloq_pct, re_tol_pct)
standards$pass <- abs(standards$re_w) <= tol
lloq <- min(standards$conc[standards$pass])
cat('LLOQ (lowest calibrator within tolerance):', lloq, 'nM\n')

cal_line <- data.frame(conc = standards$conc, ratio = predict(fit_w))
cal_plot <- ggplot(standards, aes(conc, ratio)) +
  geom_point(size = 3, colour = 'steelblue') +
  geom_line(data = cal_line, colour = 'firebrick') +
  scale_x_log10() + scale_y_log10() + theme_bw() +
  labs(title = 'Weighted calibration (1/x^2)', x = 'Concentration (nM)', y = 'Analyte/IS ratio')
ggsave(file.path(out_dir, 'calibration_curve.png'), plot = cal_plot, width = 7, height = 5)

# IS normalization on samples: matrix effect cancels in the analyte/IS ratio.
samples <- data.frame(
  sample = paste0('S', 1:6),
  condition = rep(c('control', 'treated'), each = 3),
  analyte_area = c(24500, 23800, 25200, 61000, 58500, 63000),
  istd_area = c(1.00e5, 9.8e4, 1.01e5, 9.9e4, 1.00e5, 1.02e5),
  quantifier_area = c(24500, 23800, 25200, 61000, 58500, 63000),
  qualifier_area = c(9600, 9300, 9900, 23800, 9200, 24600)
)
samples$ratio <- samples$analyte_area / samples$istd_area
samples$conc <- (samples$ratio - coef(fit_w)[1]) / coef(fit_w)[2]
samples$conc[samples$conc < lloq] <- NA

# Ion-ratio confirmation: a drifted qualifier/quantifier ratio means an interference.
# Calibrator ratio ~0.39 here; S5's collapsed qualifier (0.16) falls outside +-30%.
cal_ratio <- 0.39
samples$ion_ratio <- samples$qualifier_area / samples$quantifier_area
samples$id_confirmed <- abs(samples$ion_ratio - cal_ratio) / cal_ratio <= ion_ratio_tol

cat('\nQuantified samples (S5 qualifier collapsed -> id_confirmed FALSE):\n')
print(samples[, c('sample', 'condition', 'conc', 'ion_ratio', 'id_confirmed')])

write.csv(samples[, c('sample', 'condition', 'conc', 'ion_ratio', 'id_confirmed')],
          file.path(out_dir, 'sample_concentrations.csv'), row.names = FALSE)
cat('\nOutputs written to', out_dir, '\n')
