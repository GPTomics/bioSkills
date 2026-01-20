library(MOFA2)

# --- ALTERNATIVE: Use the CLL dataset from MOFA2 package ---
# The CLL (Chronic Lymphocytic Leukemia) dataset is the canonical MOFA2 example
# with matched RNA-seq, drug response, methylation, and mutation data:
#
# data('CLL_data')  # Loads list with mRNA, Drugs, Methylation, Mutations
# mofa <- create_mofa(CLL_data)
#
# For more details see: https://biofam.github.io/MOFA2/CLL.html

# Simulated multi-omics data with shared latent structure
set.seed(42)
n_samples <- 100
n_factors <- 5

# Create shared latent factors (what MOFA will try to recover)
latent_factors <- matrix(rnorm(n_samples * n_factors), nrow = n_samples)

# RNA: 500 genes, ~60% influenced by shared factors
rna_loadings <- matrix(rnorm(500 * n_factors, 0, 0.5), nrow = 500)
rna_loadings[1:300, ] <- rna_loadings[1:300, ] * 3  # Stronger signal in subset
rna <- rna_loadings %*% t(latent_factors) + matrix(rnorm(500 * n_samples, 0, 0.5), nrow = 500)
dimnames(rna) <- list(paste0('Gene', 1:500), paste0('Sample', 1:n_samples))

# Protein: 200 proteins, correlated with RNA via shared factors
protein_loadings <- matrix(rnorm(200 * n_factors, 0, 0.5), nrow = 200)
protein_loadings[1:100, ] <- protein_loadings[1:100, ] * 2.5
protein <- protein_loadings %*% t(latent_factors) + matrix(rnorm(200 * n_samples, 0, 0.5), nrow = 200)
dimnames(protein) <- list(paste0('Protein', 1:200), paste0('Sample', 1:n_samples))

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
