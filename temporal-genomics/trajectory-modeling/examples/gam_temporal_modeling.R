library(mgcv)

set.seed(42)

# --- Simulate time-course data with two conditions ---
# 10 timepoints: 0-48h, denser early sampling to capture rapid changes
timepoints <- c(0, 2, 4, 8, 12, 18, 24, 30, 36, 48)
n_replicates <- 3

gene_df <- data.frame()
for (cond in c('control', 'treated')) {
    for (rep in seq_len(n_replicates)) {
        for (t in timepoints) {
            base <- 8 + 3 * sin(pi * t / 48)
            if (cond == 'treated') {
                effect <- 1.5 * (1 - exp(-t / 10))
            } else {
                effect <- 0
            }
            # SD = 0.4: moderate biological + technical noise for log-expression
            expr <- base + effect + rnorm(1, 0, 0.4)
            gene_df <- rbind(gene_df, data.frame(
                time = t, expression = expr,
                condition = cond, replicate = rep
            ))
        }
    }
}

gene_df$condition <- as.factor(gene_df$condition)
gene_df$is_treated <- as.numeric(gene_df$condition == 'treated')

# --- Basic GAM: temporal trend ---
# k=6: basis dimension; appropriate for 10 unique timepoints (k < n_unique)
# method='REML': restricted maximum likelihood; more robust than GCV for smooth estimation
fit_basic <- gam(expression ~ s(time, k = 6), data = gene_df, method = 'REML')

cat('=== Basic GAM (temporal trend) ===\n')
summary(fit_basic)

# gam.check: residual diagnostics and basis dimension adequacy
# k-index >= 1.0 indicates sufficient basis dimension
cat('\n=== Model diagnostics ===\n')
gam.check(fit_basic)

# --- Condition comparison with difference smooth ---
# s(time, k=6): shared baseline smooth
# s(time, k=6, by=is_treated): deviation smooth for treatment effect over time
# The p-value of the second smooth tests whether conditions differ
fit_diff <- gam(
    expression ~ is_treated + s(time, k = 6) + s(time, k = 6, by = is_treated),
    data = gene_df, method = 'REML'
)

cat('\n=== Condition comparison (difference smooth) ===\n')
summary(fit_diff)

# --- Model comparison: linear vs GAM ---
fit_linear <- gam(expression ~ time, data = gene_df, method = 'REML')

# AIC comparison: lower AIC = better model
# Difference > 2: meaningful improvement; > 10: strong improvement
aic_linear <- AIC(fit_linear)
aic_gam <- AIC(fit_basic)
cat(sprintf('\nAIC linear: %.1f, AIC GAM: %.1f, Delta AIC: %.1f\n',
            aic_linear, aic_gam, aic_linear - aic_gam))

# --- Prediction with confidence intervals ---
new_ctrl <- data.frame(time = seq(0, 48, length.out = 200), is_treated = 0)
new_treat <- data.frame(time = seq(0, 48, length.out = 200), is_treated = 1)

pred_ctrl <- predict(fit_diff, newdata = new_ctrl, se.fit = TRUE)
pred_treat <- predict(fit_diff, newdata = new_treat, se.fit = TRUE)

# 1.96 * SE: 95% pointwise confidence interval
new_ctrl$fitted <- pred_ctrl$fit
new_ctrl$lower <- pred_ctrl$fit - 1.96 * pred_ctrl$se.fit
new_ctrl$upper <- pred_ctrl$fit + 1.96 * pred_ctrl$se.fit

new_treat$fitted <- pred_treat$fit
new_treat$lower <- pred_treat$fit - 1.96 * pred_treat$se.fit
new_treat$upper <- pred_treat$fit + 1.96 * pred_treat$se.fit

# --- Visualization ---
pdf('gam_trajectory_results.pdf', width = 12, height = 5)
par(mfrow = c(1, 2))

plot(gene_df$time[gene_df$condition == 'control'],
     gene_df$expression[gene_df$condition == 'control'],
     pch = 19, col = rgb(0.2, 0.4, 0.8, 0.5), cex = 0.8,
     xlab = 'Time (hours)', ylab = 'Expression',
     main = 'GAM trajectory by condition', ylim = c(6, 14))
points(gene_df$time[gene_df$condition == 'treated'],
       gene_df$expression[gene_df$condition == 'treated'],
       pch = 17, col = rgb(0.8, 0.2, 0.2, 0.5), cex = 0.8)

lines(new_ctrl$time, new_ctrl$fitted, col = 'blue', lwd = 2)
polygon(c(new_ctrl$time, rev(new_ctrl$time)),
        c(new_ctrl$lower, rev(new_ctrl$upper)),
        col = rgb(0.2, 0.4, 0.8, 0.15), border = NA)

lines(new_treat$time, new_treat$fitted, col = 'red', lwd = 2)
polygon(c(new_treat$time, rev(new_treat$time)),
        c(new_treat$lower, rev(new_treat$upper)),
        col = rgb(0.8, 0.2, 0.2, 0.15), border = NA)

legend('topleft', c('Control', 'Treated'), col = c('blue', 'red'),
       lwd = 2, pch = c(19, 17), bty = 'n')

plot(fit_basic, shade = TRUE, shade.col = rgb(0, 0, 0, 0.1),
     xlab = 'Time (hours)', ylab = 's(time)',
     main = 'Smooth term (basic GAM)')
rug(gene_df$time)

dev.off()
cat('\nPlot saved to gam_trajectory_results.pdf\n')

# --- Genome-wide GAM fitting example ---
cat('\n=== Genome-wide GAM fitting (5 demo genes) ===\n')
n_demo <- 5
demo_results <- data.frame()
for (g in seq_len(n_demo)) {
    demo_expr <- 8 + rnorm(nrow(gene_df), 0, 0.5)
    if (g <= 3) {
        demo_expr <- demo_expr + 2 * sin(pi * gene_df$time / 48)
    }
    demo_df <- data.frame(expression = demo_expr, time = gene_df$time)
    fit_g <- gam(expression ~ s(time, k = 6), data = demo_df, method = 'REML')
    s_tab <- summary(fit_g)$s.table
    demo_results <- rbind(demo_results, data.frame(
        gene = paste0('gene_', g), edf = s_tab[, 'edf'],
        F_stat = s_tab[, 'F'], p_value = s_tab[, 'p-value']
    ))
}
demo_results$q_value <- p.adjust(demo_results$p_value, method = 'BH')
print(demo_results)
