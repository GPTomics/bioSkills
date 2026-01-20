library(MOFA2)
library(ggplot2)

# === CONFIGURATION ===
output_dir <- 'results/'
dir.create(output_dir, showWarnings = FALSE)

# === 1. LOAD DATA ===
cat('Loading data...\n')
rna <- read.csv('rnaseq_normalized.csv', row.names = 1)
protein <- read.csv('proteomics_normalized.csv', row.names = 1)
metab <- read.csv('metabolomics_normalized.csv', row.names = 1)

# Find common samples
common_samples <- Reduce(intersect, list(rownames(rna), rownames(protein), rownames(metab)))
cat('Common samples:', length(common_samples), '\n')

rna <- rna[common_samples, ]
protein <- protein[common_samples, ]
metab <- metab[common_samples, ]

# === 2. FEATURE SELECTION ===
cat('Selecting variable features...\n')
select_top_var <- function(data, n) {
    vars <- apply(data, 2, var, na.rm = TRUE)
    data[, names(sort(vars, decreasing = TRUE))[1:min(n, ncol(data))]]
}

rna_var <- select_top_var(rna, 2000)
protein_var <- select_top_var(protein, 1000)
metab_var <- select_top_var(metab, 500)

# === 3. CREATE MOFA OBJECT ===
data_list <- list(
    RNA = t(as.matrix(rna_var)),
    Protein = t(as.matrix(protein_var)),
    Metabolome = t(as.matrix(metab_var))
)

mofa <- create_mofa(data_list)

# === 4. CONFIGURE MODEL ===
data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- TRUE

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 10

train_opts <- get_default_training_options(mofa)
train_opts$maxiter <- 1000
train_opts$seed <- 42

mofa <- prepare_mofa(mofa, data_options = data_opts,
                     model_options = model_opts, training_options = train_opts)

# === 5. TRAIN MODEL ===
cat('Training MOFA model...\n')
mofa <- run_mofa(mofa, outfile = file.path(output_dir, 'mofa_model.hdf5'), use_basilisk = TRUE)

# === 6. ANALYZE RESULTS ===
cat('Analyzing results...\n')

# Variance explained
plot_variance_explained(mofa, max_r2 = 15)
ggsave(file.path(output_dir, 'variance_explained.png'), width = 10, height = 6)

# Factor values
factors <- get_factors(mofa)[[1]]

# Top weights
plot_top_weights(mofa, view = 'RNA', factors = 1:3, nfeatures = 10)
ggsave(file.path(output_dir, 'top_weights_rna.png'), width = 10, height = 8)

plot_top_weights(mofa, view = 'Protein', factors = 1:3, nfeatures = 10)
ggsave(file.path(output_dir, 'top_weights_protein.png'), width = 10, height = 8)

# === 7. EXPORT ===
write.csv(factors, file.path(output_dir, 'factor_values.csv'))
weights <- get_weights(mofa, as.data.frame = TRUE)
write.csv(weights, file.path(output_dir, 'all_weights.csv'), row.names = FALSE)

cat('Analysis complete! Results in', output_dir, '\n')
