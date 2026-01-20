library(ggplot2)
library(dplyr)

# === CONFIGURATION ===
output_dir <- 'results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. STANDARD CURVE DATA ===
cat('Setting up standard curves...\n')

standards <- data.frame(
    concentration = c(0, 1, 5, 10, 25, 50, 100, 250, 500, 1000),
    analyte_area = c(50, 4800, 24500, 49000, 122000, 245000, 490000, 1220000, 2450000, 4900000),
    istd_area = c(100000, 98000, 102000, 99000, 101000, 100000, 98000, 101000, 99000, 100000)
)

# Response ratio
standards$response_ratio <- standards$analyte_area / standards$istd_area

# === 2. FIT STANDARD CURVE ===
cat('Fitting calibration curve...\n')

# Weighted linear regression (1/x^2 weighting)
weights <- 1 / (standards$concentration + 1)^2
fit <- lm(response_ratio ~ concentration, data = standards, weights = weights)

cat('Slope:', coef(fit)[2], '\n')
cat('Intercept:', coef(fit)[1], '\n')
cat('R-squared:', summary(fit)$r.squared, '\n')

# Back-calculated concentrations
standards$back_calc <- (standards$response_ratio - coef(fit)[1]) / coef(fit)[2]
standards$accuracy <- ifelse(standards$concentration > 0,
                             standards$back_calc / standards$concentration * 100, NA)

# === 3. STANDARD CURVE PLOT ===
ggplot(standards, aes(x = concentration, y = response_ratio)) +
    geom_point(size = 3, color = 'blue') +
    geom_smooth(method = 'lm', se = TRUE, color = 'red') +
    theme_bw() +
    labs(title = sprintf('Standard Curve (R² = %.4f)', summary(fit)$r.squared),
         x = 'Concentration (nM)', y = 'Response Ratio')
ggsave(file.path(output_dir, 'standard_curve.png'), width = 8, height = 6)

# === 4. SAMPLE QUANTIFICATION ===
cat('Quantifying samples...\n')

samples <- data.frame(
    sample = paste0('Sample', 1:12),
    condition = rep(c('Control', 'Treatment'), each = 6),
    analyte_area = c(52000, 48000, 55000, 51000, 49000, 53000,
                     125000, 118000, 132000, 121000, 128000, 119000),
    istd_area = c(100000, 98000, 101000, 99000, 100000, 102000,
                  99000, 101000, 100000, 98000, 101000, 99000)
)

# Calculate response ratio and concentration
samples$response_ratio <- samples$analyte_area / samples$istd_area
samples$concentration <- (samples$response_ratio - coef(fit)[1]) / coef(fit)[2]

# Account for dilution
dilution_factor <- 10
samples$concentration_original <- samples$concentration * dilution_factor

cat('Sample concentrations calculated\n')

# === 5. QC ASSESSMENT ===
cat('Assessing QC samples...\n')

qc_samples <- data.frame(
    level = rep(c('Low', 'Medium', 'High'), each = 3),
    nominal = rep(c(25, 250, 750), each = 3),
    measured = c(23.5, 26.2, 24.8, 242, 258, 251, 738, 762, 745)
)

qc_summary <- qc_samples %>%
    group_by(level, nominal) %>%
    summarise(
        mean = mean(measured),
        sd = sd(measured),
        cv = sd(measured) / mean(measured) * 100,
        accuracy = mean(measured) / nominal * 100,
        .groups = 'drop'
    )

cat('\nQC Summary:\n')
print(qc_summary)

# QC acceptance
qc_summary$cv_pass <- qc_summary$cv < 15
qc_summary$accuracy_pass <- qc_summary$accuracy > 85 & qc_summary$accuracy < 115

# === 6. LOD/LOQ CALCULATION ===
residuals_sd <- sd(residuals(fit))
slope <- coef(fit)[2]

LOD <- 3.3 * residuals_sd / slope
LOQ <- 10 * residuals_sd / slope

cat('\nLOD:', round(LOD, 2), 'nM\n')
cat('LOQ:', round(LOQ, 2), 'nM\n')

# === 7. GROUP COMPARISON ===
cat('\nComparing groups...\n')

group_summary <- samples %>%
    group_by(condition) %>%
    summarise(
        mean = mean(concentration_original),
        sd = sd(concentration_original),
        n = n(),
        .groups = 'drop'
    )

print(group_summary)

# t-test
ttest <- t.test(concentration_original ~ condition, data = samples)
cat('p-value:', format(ttest$p.value, digits = 3), '\n')

# === 8. RESULTS PLOT ===
ggplot(samples, aes(x = condition, y = concentration_original, fill = condition)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 2) +
    theme_bw() +
    labs(title = sprintf('Metabolite Concentration (p = %.3f)', ttest$p.value),
         x = 'Condition', y = 'Concentration (nM)')
ggsave(file.path(output_dir, 'group_comparison.png'), width = 6, height = 5)

# === 9. EXPORT RESULTS ===
write.csv(samples, file.path(output_dir, 'sample_concentrations.csv'), row.names = FALSE)
write.csv(qc_summary, file.path(output_dir, 'qc_summary.csv'), row.names = FALSE)

cat('\nResults saved to', output_dir, '\n')
