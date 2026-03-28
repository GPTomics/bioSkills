# Reference: limma 3.58+ | Verify API if version differs
# Per-CpG differential methylation with limma on M-values
# Recommended when sample sizes are small (n=3-5 per group)

library(limma)

beta_matrix <- read.csv('beta_values.csv', row.names = 1)

# M-value: logit transform (base 2) for statistical testing
# 1e-6 offset: prevents log(0) and division by zero at beta=0 or beta=1
offset <- 1e-6
m_values <- log2((beta_matrix + offset) / (1 - beta_matrix + offset))

group <- factor(c(rep('case', 6), rep('ctrl', 6)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast_matrix <- makeContrasts(case - ctrl, levels = design)

fit <- lmFit(m_values, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
# trend=TRUE: models intensity-dependent prior variance
#   (methylation variance differs across M-value range)
# robust=TRUE: protects against outlier CpGs inflating variance estimates
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

results <- topTable(fit2, number = Inf, adjust.method = 'BH', sort.by = 'none')

# Delta-beta from original beta values (not from M-value logFC)
# logFC on M-value scale does not map linearly to beta differences
delta_beta <- rowMeans(beta_matrix[, group == 'case']) -
              rowMeans(beta_matrix[, group == 'ctrl'])
results$delta_beta <- delta_beta
results$significant <- ifelse(results$adj.P.Val < 0.05, 'TRUE', 'FALSE')

write.csv(results, 'dmc_limma_results.csv', row.names = TRUE)

n_sig <- sum(results$significant == 'TRUE')
cat(sprintf('CpGs tested: %d, significant (adj.P.Val < 0.05): %d\n', nrow(results), n_sig))
