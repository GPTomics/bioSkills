library(flowCore)

# Load FCS files
fcs_files <- list.files('data', pattern = '\\.fcs$', full.names = TRUE)
cat('Found', length(fcs_files), 'FCS files\n')

# Read as flowSet
fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)
cat('Loaded', length(fs), 'samples\n')

# Sample info
cat('\nSample names:\n')
print(sampleNames(fs))

# Channel info
cat('\nChannels:\n')
print(colnames(fs[[1]]))

# Events per sample
cat('\nEvents per sample:\n')
print(fsApply(fs, nrow))

# Parameter details
params <- pData(parameters(fs[[1]]))
print(params[, c('name', 'desc', 'range')])
