# Reference: mgcv 1.9+ | Verify API if version differs
library(mgcv)

set.seed(42)

# --- Simulate time-course data with two conditions ---
# 10 timepoints over 0-48h, denser early to capture rapid changes.
timepoints <- c(0, 2, 4, 8, 12, 18, 24, 30, 36, 48)
n_replicates <- 3

gene_df <- data.frame()
for (cond in c('control', 'treated')) {
    for (rep in seq_len(n_replicates)) {
        for (t in timepoints) {
            base <- 8 + 3 * sin(pi * t / 48)
            effect <- if (cond == 'treated') 1.5 * (1 - exp(-t / 10)) else 0
            # SD=0.4: moderate biological + technical noise on a LOG-expression scale,
            # which is what makes a Gaussian GAM valid here (see the NB block for raw counts).
            expr <- base + effect + rnorm(1, 0, 0.4)
            gene_df <- rbind(gene_df, data.frame(time = t, expression = expr,
                                                 condition = cond, replicate = rep))
        }
    }
}

# --- Basic GAM: temporal trend ---
# k=6: basis-dimension CEILING (max flexibility), NOT a knot count and < 10 unique timepoints.
# method='REML': better-behaved objective than GCV, resists under/over-smoothing (Wood 2011).
# Independence: with 3 replicates/timepoint here plain gam() is defensible; for a SINGLE series with
# autocorrelated residuals, switch to gamm(correlation = corAR1(form = ~time)) or bam(rho=, AR.start=).
fit_basic <- gam(expression ~ s(time, k = 6, bs = 'tp'), data = gene_df, method = 'REML')

cat('=== Basic GAM (temporal trend) ===\n')
print(summary(fit_basic))
# edf near 1 => penalty shrank to linear; edf near k-1 => nearly full flexibility.
# The smooth-term p-value is APPROXIMATE (Wood 2013): treat as categorical, not a magnitude.

cat('\n=== Model diagnostics (k.check) ===\n')
# k-index < 1 with small p => residual pattern the basis is too rigid to catch.
# Correct response: double k and refit; if edf barely moves, suspect autocorrelation.
print(k.check(fit_basic))

# --- Raw counts: NB family with a library-size offset ---
# Gaussian is wrong on raw counts (overdispersed, mean-variance coupled, can predict < 0).
# nb() estimates theta by REML; offset(log(libsize)) makes s(time) model rate, not depth.
count_df <- gene_df
count_df$library_size <- round(runif(nrow(count_df), 8e6, 1.2e7))
count_df$counts <- rpois(nrow(count_df), lambda = exp(count_df$expression - 8) *
                         count_df$library_size / 1e7)
fit_nb <- gam(counts ~ s(time, k = 6) + offset(log(library_size)),
              data = count_df, family = nb(), method = 'REML')
cat('\n=== NB GAM on raw counts (test column is Chi.sq, not F) ===\n')
print(summary(fit_nb)$s.table)

# --- Condition comparison: ordered-factor difference smooth ---
# Ordered factor => s(time, by=condition) is the difference smooth (treated minus control);
# its single p-value directly tests divergence. The parametric main effect is REQUIRED
# because centered smooths cannot carry the group's overall level.
gene_df$condition <- as.ordered(gene_df$condition)
fit_diff <- gam(expression ~ condition + s(time, k = 6) + s(time, k = 6, by = condition),
                data = gene_df, method = 'REML')
cat('\n=== Condition comparison (ordered-factor difference smooth) ===\n')
print(summary(fit_diff))

# --- Model comparison: linear vs GAM ---
fit_linear <- gam(expression ~ time, data = gene_df, method = 'REML')
# Delta AIC > 2: meaningful improvement; > 10: strong. Decides non-linear vs linear.
cat(sprintf('\nAIC linear: %.1f, AIC GAM: %.1f, Delta AIC: %.1f\n',
            AIC(fit_linear), AIC(fit_basic), AIC(fit_linear) - AIC(fit_basic)))

# --- Prediction with pointwise intervals (within sampled range only) ---
lvls <- levels(gene_df$condition)
grid_ctrl <- data.frame(time = seq(0, 48, length.out = 200),
                        condition = ordered('control', levels = lvls))
grid_treat <- data.frame(time = seq(0, 48, length.out = 200),
                         condition = ordered('treated', levels = lvls))
pred_ctrl <- predict(fit_diff, newdata = grid_ctrl, se.fit = TRUE)
pred_treat <- predict(fit_diff, newdata = grid_treat, se.fit = TRUE)
# 1.96*SE is a POINTWISE band; overlapping bands are NOT a divergence test (use the p-value).
grid_ctrl$fitted <- pred_ctrl$fit; grid_ctrl$lower <- pred_ctrl$fit - 1.96 * pred_ctrl$se.fit
grid_ctrl$upper <- pred_ctrl$fit + 1.96 * pred_ctrl$se.fit
grid_treat$fitted <- pred_treat$fit; grid_treat$lower <- pred_treat$fit - 1.96 * pred_treat$se.fit
grid_treat$upper <- pred_treat$fit + 1.96 * pred_treat$se.fit

# --- Visualization (written to a temp file, then removed: no stray outputs) ---
out_pdf <- tempfile(fileext = '.pdf')
pdf(out_pdf, width = 12, height = 5)
par(mfrow = c(1, 2))
is_ctrl <- gene_df$condition == 'control'
plot(gene_df$time[is_ctrl], gene_df$expression[is_ctrl], pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.5), cex = 0.8, xlab = 'Time (hours)', ylab = 'Expression',
     main = 'GAM trajectory by condition', ylim = c(6, 14))
points(gene_df$time[!is_ctrl], gene_df$expression[!is_ctrl], pch = 17,
       col = rgb(0.8, 0.2, 0.2, 0.5), cex = 0.8)
lines(grid_ctrl$time, grid_ctrl$fitted, col = 'blue', lwd = 2)
polygon(c(grid_ctrl$time, rev(grid_ctrl$time)), c(grid_ctrl$lower, rev(grid_ctrl$upper)),
        col = rgb(0.2, 0.4, 0.8, 0.15), border = NA)
lines(grid_treat$time, grid_treat$fitted, col = 'red', lwd = 2)
polygon(c(grid_treat$time, rev(grid_treat$time)), c(grid_treat$lower, rev(grid_treat$upper)),
        col = rgb(0.8, 0.2, 0.2, 0.15), border = NA)
legend('topleft', c('Control', 'Treated'), col = c('blue', 'red'), lwd = 2,
       pch = c(19, 17), bty = 'n')
plot(fit_basic, shade = TRUE, shade.col = rgb(0, 0, 0, 0.1), xlab = 'Time (hours)',
     ylab = 's(time)', main = 'Smooth term (basic GAM)')
rug(gene_df$time)
dev.off()
unlink(out_pdf)

# --- Genome-wide GAM fitting + BH FDR (5 demo genes) ---
cat('\n=== Genome-wide GAM fitting (5 demo genes) ===\n')
demo_results <- data.frame()
for (g in seq_len(5)) {
    demo_expr <- 8 + rnorm(nrow(gene_df), 0, 0.5)
    if (g <= 3) demo_expr <- demo_expr + 2 * sin(pi * gene_df$time / 48)
    demo_df <- data.frame(expression = demo_expr, time = gene_df$time)
    fit_g <- gam(expression ~ s(time, k = 6), data = demo_df, method = 'REML')
    s_tab <- summary(fit_g)$s.table
    demo_results <- rbind(demo_results, data.frame(gene = paste0('gene_', g),
        edf = s_tab[, 'edf'], p_value = s_tab[, 'p-value']))
}
# q<0.05: standard FDR floor. Per-gene smooth p-values are approximate, so the FDR is too.
demo_results$q_value <- p.adjust(demo_results$p_value, method = 'BH')
print(demo_results)
