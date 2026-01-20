library(MOFA2)

# Simulate multi-omics data (replace with real data)
set.seed(42)
n_samples <- 100
rna <- matrix(rnorm(n_samples * 500), nrow = 500, dimnames = list(paste0('Gene', 1:500), paste0('Sample', 1:n_samples)))
protein <- matrix(rnorm(n_samples * 200), nrow = 200, dimnames = list(paste0('Protein', 1:200), paste0('Sample', 1:n_samples)))

data_list <- list(RNA = rna, Protein = protein)

# Create and configure MOFA
mofa <- create_mofa(data_list)
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 10

train_opts <- get_default_training_options(mofa)
train_opts$seed <- 42

mofa <- prepare_mofa(mofa, model_options = model_opts, training_options = train_opts)
mofa <- run_mofa(mofa)

# Results
cat('Factors learned:', mofa@dimensions$K, '\n')
var_exp <- get_variance_explained(mofa)
print(var_exp$r2_total)

# Export factor values
factors <- get_factors(mofa, as.data.frame = TRUE)
write.csv(factors, 'mofa_factors.csv', row.names = FALSE)
