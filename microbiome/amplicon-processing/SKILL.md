---
name: bio-microbiome-amplicon-processing
description: Amplicon sequence variant (ASV) inference from 16S rRNA or ITS amplicon sequencing using DADA2. Covers quality filtering, error learning, denoising, and chimera removal. Use when processing demultiplexed amplicon FASTQ files to generate an ASV table for downstream analysis.
tool_type: r
primary_tool: dada2
---

# Amplicon Processing with DADA2

## Complete DADA2 Workflow

```r
library(dada2)

path <- 'raw_reads'
fnFs <- sort(list.files(path, pattern = '_R1_001.fastq.gz', full.names = TRUE))
fnRs <- sort(list.files(path, pattern = '_R2_001.fastq.gz', full.names = TRUE))
sample_names <- sapply(strsplit(basename(fnFs), '_'), `[`, 1)

# Quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

## Quality Filtering and Trimming

```r
filtFs <- file.path('filtered', paste0(sample_names, '_F_filt.fastq.gz'))
filtRs <- file.path('filtered', paste0(sample_names, '_R_filt.fastq.gz'))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Filter parameters depend on amplicon region and read length
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(240, 160),      # Trim to quality scores
                     maxN = 0,                     # No ambiguous bases
                     maxEE = c(2, 2),              # Max expected errors
                     truncQ = 2,                   # Truncate at first Q <= 2
                     rm.phix = TRUE,               # Remove PhiX
                     compress = TRUE,
                     multithread = TRUE)
```

## Error Rate Learning

```r
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Visualize error rates
plotErrors(errF, nominalQ = TRUE)
```

## Sample Inference (Denoising)

```r
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Check results
dadaFs[[1]]
```

## Merge Paired Reads

```r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Check merge success
head(mergers[[1]])
```

## Construct Sequence Table

```r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Check length distribution
table(nchar(getSequences(seqtab)))
```

## Remove Chimeras

```r
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus',
                                     multithread = TRUE, verbose = TRUE)

# Percentage retained
sum(seqtab_nochim) / sum(seqtab)
```

## Track Reads Through Pipeline

```r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab_nochim))
colnames(track) <- c('input', 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim')
rownames(track) <- sample_names
track
```

## ITS-Specific Processing

```r
# For ITS, use cutadapt to remove primers first (variable length amplicons)
# Then skip truncLen (don't truncate ITS to fixed length)

out_its <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                         maxN = 0, maxEE = c(2, 2), truncQ = 2,
                         minLen = 50,  # Minimum length
                         rm.phix = TRUE, compress = TRUE, multithread = TRUE)
```

## Related Skills

- taxonomy-assignment - Assign taxonomy to ASVs
- read-qc/quality-reports - Pre-DADA2 quality assessment
- diversity-analysis - Analyze ASV table
