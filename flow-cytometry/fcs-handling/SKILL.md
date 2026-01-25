---
name: bio-flow-cytometry-fcs-handling
description: Read and manipulate Flow Cytometry Standard (FCS) files. Covers loading data, accessing parameters, and basic data exploration. Use when loading and inspecting flow or mass cytometry data before preprocessing.
tool_type: r
primary_tool: flowCore
---

# FCS File Handling

## Load FCS Files

```r
library(flowCore)

# Read single FCS file
fcs <- read.FCS('sample.fcs', transformation = FALSE, truncate_max_range = FALSE)

# File info
print(fcs)

# Parameter names
colnames(fcs)  # Short names
pData(parameters(fcs))  # Full metadata including descriptions
```

## Load Multiple Files

```r
# Read multiple files into flowSet
files <- list.files('data', pattern = '\\.fcs$', full.names = TRUE)
fs <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)

# Sample names
sampleNames(fs)

# Access individual samples
fcs1 <- fs[[1]]
```

## Access Expression Data

```r
# Get expression matrix
expr <- exprs(fcs)
head(expr)

# Dimensions
dim(expr)  # cells x channels

# Channel statistics
summary(expr)

# Get specific channels
cd4_expr <- expr[, 'CD4']
```

## Channel Metadata

```r
# Parameter information
params <- pData(parameters(fcs))
print(params)

# Parameter columns:
# - name: short name (e.g., "FL1-A")
# - desc: description (e.g., "CD4")
# - range: max value
# - minRange: min value

# Get channel mapping
channel_map <- setNames(params$desc, params$name)
```

## Rename Channels

```r
# Rename using descriptions
rename_channels <- function(fcs) {
    params <- pData(parameters(fcs))
    new_names <- ifelse(is.na(params$desc) | params$desc == '',
                        params$name, params$desc)
    colnames(fcs) <- new_names
    return(fcs)
}

fcs_renamed <- rename_channels(fcs)
```

## Subsetting Data

```r
# Subset by cells (rows)
fcs_subset <- fcs[1:1000, ]

# Subset by channels (columns)
fcs_markers <- fcs[, c('CD4', 'CD8', 'CD3')]

# Subset by expression values
high_cd4 <- fcs[exprs(fcs)[, 'CD4'] > 1000, ]
```

## Merge flowSets

```r
# Combine multiple flowSets
fs_combined <- rbind2(fs1, fs2)

# Or concatenate into single flowFrame
all_data <- fsApply(fs, exprs)
all_data <- do.call(rbind, all_data)
```

## Write FCS Files

```r
# Write single file
write.FCS(fcs, 'output.fcs')

# Write flowSet
write.flowSet(fs, outdir = 'output_dir')
```

## Sample Metadata

```r
# Add sample annotations
pData(fs) <- data.frame(
    name = sampleNames(fs),
    condition = c('Control', 'Control', 'Treatment', 'Treatment'),
    patient = c('P1', 'P2', 'P1', 'P2')
)

# Access
pData(fs)
```

## Basic Visualization

```r
library(ggcyto)

# Density plot
autoplot(fcs, 'FSC-A')

# Bivariate plot
autoplot(fcs, 'CD4', 'CD8')

# Multiple samples
autoplot(fs, 'CD4', 'CD8')
```

## Check Data Quality

```r
# Time parameter check
if ('Time' %in% colnames(fcs)) {
    time <- exprs(fcs)[, 'Time']
    plot(time, type = 'l', main = 'Acquisition Time')
}

# Event count per file
fsApply(fs, nrow)

# Check for saturated events
saturation <- apply(exprs(fcs), 2, function(x) mean(x == max(x)) * 100)
print(saturation)
```

## Convert to Data Frame

```r
# For use with tidyverse
library(tidyverse)

df <- as.data.frame(exprs(fcs))
df$sample <- 'sample1'

# From flowSet
df_all <- fsApply(fs, function(f) {
    d <- as.data.frame(exprs(f))
    d$sample <- identifier(f)
    d
}, simplify = FALSE)
df_all <- bind_rows(df_all)
```

## Related Skills

- compensation-transformation - Apply compensation and transforms
- gating-analysis - Define cell populations
- clustering-phenotyping - Unsupervised analysis
