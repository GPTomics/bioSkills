# Reference: limma 3.58+ | Verify API if version differs
# Per-CpG differential methylation with limma moderated-t on M-values.
# This is the ARRAY / CONTINUOUS path (450K/EPIC, or high uniform-coverage sequencing).
# For bisulfite sequencing COUNTS prefer a beta-binomial count model (DSS / methylKit MN)
# that uses coverage as precision instead of collapsing to a continuous beta.

library(limma)

beta_matrix <- read.csv('beta_values.csv', row.names = 1)

# M-value: logit transform (base 2); test on M, report effect on beta.
# 1e-3 offset: boundary-safe so beta = 0 or 1 does not blow up the logit
offset <- 1e-3
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

out_file <- tempfile(fileext = '.csv')   # write to a temp path, not the working directory
write.csv(results, out_file, row.names = TRUE)

n_sig <- sum(results$significant == 'TRUE')
cat(sprintf('CpGs tested: %d, significant (adj.P.Val < 0.05): %d\n', nrow(results), n_sig))
