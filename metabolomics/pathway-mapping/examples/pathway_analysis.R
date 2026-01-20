library(MetaboAnalystR)

# Significant metabolites (HMDB IDs)
metabolites <- c(
    'HMDB0000122',  # Glucose
    'HMDB0000190',  # Lactate
    'HMDB0000243',  # Pyruvate
    'HMDB0001206',  # Succinate
    'HMDB0000695',  # Fumarate
    'HMDB0000254'   # Citrate
)

cat('Input metabolites:', length(metabolites), '\n')

# Initialize
mSet <- InitDataObjects('conc', 'pathora', FALSE)
mSet <- SetOrganism(mSet, 'hsa')

# Map data
mSet <- Setup.MapData(mSet, metabolites)
mSet <- CrossReferencing(mSet, 'hmdb')

# Check mapping
cat('Mapped metabolites:', mSet$dataSet$cmpd.count, '\n')

# Pathway analysis
mSet <- SetKEGG.PathLib(mSet, 'hsa', 'current')
mSet <- SetMetabolomeFilter(mSet, FALSE)
mSet <- CalculateOraScore(mSet, 'rbc', 'hyperg')

# Results
results <- mSet$analSet$ora.mat
results <- as.data.frame(results)
results$pathway <- rownames(results)

# Filter significant
sig_pathways <- results[results[, 'Raw p'] < 0.05, ]
sig_pathways <- sig_pathways[order(sig_pathways[, 'Raw p']), ]

cat('\nSignificant pathways (p<0.05):\n')
print(sig_pathways[, c('pathway', 'Total', 'Hits', 'Raw p')])

# Save
write.csv(results, 'pathway_results.csv', row.names = FALSE)
cat('\nResults saved to pathway_results.csv\n')
