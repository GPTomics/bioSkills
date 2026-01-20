library(flowCore)
library(flowWorkspace)
library(ggcyto)

# Load preprocessed data
fs <- read.flowSet(list.files('preprocessed', pattern = '\\.fcs$', full.names = TRUE))
gs <- GatingSet(fs)

# 1. Lymphocyte gate (FSC vs SSC)
lymph_gate <- polygonGate(filterId = 'Lymphocytes',
    .gate = matrix(c(50000, 20000, 150000, 200000, 250000, 150000, 200000, 50000),
                   ncol = 2, dimnames = list(NULL, c('FSC-A', 'SSC-A'))))
gs_pop_add(gs, lymph_gate, parent = 'root')

# 2. Singlet gate (FSC-A vs FSC-H)
singlet_gate <- polygonGate(filterId = 'Singlets',
    .gate = matrix(c(50000, 50000, 250000, 250000, 250000, 200000, 50000, 20000),
                   ncol = 2, dimnames = list(NULL, c('FSC-A', 'FSC-H'))))
gs_pop_add(gs, singlet_gate, parent = 'Lymphocytes')

# 3. CD3+ gate
cd3_gate <- rectangleGate(filterId = 'CD3+', 'CD3' = c(2, Inf))
gs_pop_add(gs, cd3_gate, parent = 'Singlets')

# 4. CD4/CD8 quadrant gate
quad_gate <- quadGate(filterId = 'CD4_CD8', 'CD4' = 2, 'CD8' = 2)
gs_pop_add(gs, quad_gate, parent = 'CD3+')

# Recompute
recompute(gs)

# Statistics
stats <- gs_pop_get_stats(gs)
print(stats)

# Plot hierarchy
plot(gs)

# Save
save_gs(gs, 'gating_set')
cat('Saved gating set\n')
