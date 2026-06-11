# Reference: MOFA2 1.12+ | Verify API if version differs
# Unsupervised cross-omic factor discovery. The deliverable is the per-view variance
# decomposition plus factors that are annotated, checked against technical covariates, and
# robust across seeds - not the factors themselves. Each view is features-by-samples and
# already per-omic normalized / variance-stabilized (data-harmonization owns that).

library(MOFA2)

set.seed(42)
n <- 60
samples <- paste0('S', seq_len(n))
metadata <- data.frame(row.names=samples,
                       condition=factor(rep(c('case', 'control'), length.out=n)),
                       batch=factor(rep(c('b1', 'b2'), each=n / 2)))

rna  <- matrix(rnorm(2000 * n), nrow=2000, dimnames=list(paste0('gene', 1:2000), samples))   # features x samples
prot <- matrix(rnorm(300 * n), nrow=300, dimnames=list(paste0('prot', 1:300), samples))
meth <- matrix(rnorm(800 * n), nrow=800, dimnames=list(paste0('cg', 1:800), samples))

mofa <- create_mofa(list(RNA=rna, Protein=prot, Methylation=meth))

data_opts  <- get_default_data_options(mofa)
model_opts <- get_default_model_options(mofa)
train_opts <- get_default_training_options(mofa)

model_opts$num_factors <- 15                  # over-specify; ARD prunes inactive factors
model_opts$likelihoods <- c(RNA='gaussian', Protein='gaussian', Methylation='gaussian')   # counts transformed upstream -> gaussian
train_opts$drop_factor_threshold <- 0.01      # drop factors explaining <1% variance in ALL views
train_opts$seed <- 42                          # set for reproducibility; confirm headline factors recur on retrain
data_opts$scale_views <- TRUE                  # equalize per-view variance when feature counts cannot be balanced

mofa <- prepare_mofa(mofa, data_options=data_opts, model_options=model_opts, training_options=train_opts)
mofa <- run_mofa(mofa, outfile=file.path(tempdir(), 'model.hdf5'), use_basilisk=TRUE)

var_exp <- get_variance_explained(mofa)        # $r2_per_factor is the central output: shared = 2+ views, view-specific = 1
print(var_exp$r2_total)

md <- metadata[unlist(samples_names(mofa)), ]
md$sample <- rownames(md)                       # MOFA2's samples_metadata<- requires a literal 'sample' column
samples_metadata(mofa) <- md
correlate_factors_with_covariates(mofa, covariates=c('condition', 'batch'))   # a factor that correlates with batch IS a batch factor
