# tximport Workflow - Usage Guide

## Overview
tximport imports transcript-level quantifications and summarizes them to the gene level for differential expression analysis with DESeq2 or edgeR.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('tximport')
BiocManager::install('tximeta')  # Optional: metadata-aware import
BiocManager::install('GenomicFeatures')  # For creating tx2gene from GTF
```

## Quick Start
Tell your AI agent what you want to do:
- "Import my Salmon quantification results into R"
- "Summarize transcript counts to gene level for DESeq2"
- "Create a tx2gene mapping from my GTF file"

## Example Prompts
### Basic Import
> "Import Salmon quant.sf files from my samples and create a gene-level count matrix"

> "Load my kallisto results into R for differential expression analysis"

### tx2gene Creation
> "Create a transcript-to-gene mapping from my Ensembl GTF file"

> "Build tx2gene using biomaRt for human transcripts"

### Downstream Analysis
> "Import Salmon results and create a DESeqDataSet ready for analysis"

> "Use tximeta to automatically link my Salmon results to annotations"

## What the Agent Will Do
1. Locate all quantification output files (quant.sf or abundance.tsv)
2. Create or load the tx2gene mapping
3. Import transcript-level counts with tximport
4. Summarize to gene level with proper length scaling
5. Create a DESeqDataSet or DGEList for downstream analysis

## Why tximport?

1. **Corrects for transcript length bias** - Longer genes get more reads
2. **Aggregates transcripts to genes** - Most DE tools work at gene level
3. **Preserves uncertainty** - Passes length information to DESeq2/edgeR
4. **Fast** - No re-counting required

## Creating tx2gene

The tx2gene data frame must have exactly two columns:
- `TXNAME` - Transcript identifiers matching your quant files
- `GENEID` - Gene identifiers for summarization

### Method 1: From Ensembl GTF
```r
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('Homo_sapiens.GRCh38.110.gtf.gz')
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- select(txdb, k, 'GENEID', 'TXNAME')
write.csv(tx2gene, 'tx2gene.csv', row.names = FALSE)
```

### Method 2: From biomaRt
```r
library(biomaRt)
mart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
tx2gene <- getBM(
    attributes = c('ensembl_transcript_id_version', 'ensembl_gene_id_version'),
    mart = mart
)
colnames(tx2gene) <- c('TXNAME', 'GENEID')
```

### Method 3: Quick Parse from IDs
```r
# If transcript IDs are like ENST00000456328.2|ENSG00000223972.5
quant <- read.table('sample1_quant/quant.sf', header = TRUE)
tx2gene <- data.frame(
    TXNAME = quant$Name,
    GENEID = sapply(strsplit(quant$Name, '\\|'), '[', 2)
)
```

## Complete Workflow Example

```r
library(tximport)
library(DESeq2)

# 1. Set up files
sample_names <- c('ctrl1', 'ctrl2', 'treat1', 'treat2')
files <- file.path(paste0(sample_names, '_quant'), 'quant.sf')
names(files) <- sample_names

# 2. Load tx2gene
tx2gene <- read.csv('tx2gene.csv')

# 3. Import
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)

# 4. Create metadata
coldata <- data.frame(
    condition = factor(c('control', 'control', 'treated', 'treated')),
    row.names = sample_names
)

# 5. Create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# 6. Filter and analyze
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds)
```

## tximeta Alternative

tximeta automatically links to annotation databases:

```r
library(tximeta)

coldata <- data.frame(
    files = files,
    names = sample_names,
    condition = c('control', 'control', 'treated', 'treated')
)

# Import with automatic annotation
se <- tximeta(coldata)

# Summarize to gene level
gse <- summarizeToGene(se)

# Get gene symbols
library(org.Hs.eg.db)
gse <- addIds(gse, 'SYMBOL', gene = TRUE)
```

## Tips
- Match versions - tx2gene IDs must exactly match quant file IDs
- Use `ignoreTxVersion = TRUE` if transcript versions don't match
- Save tx2gene once and reuse for all analyses with the same annotation
- Check import results - verify row counts match expected gene count
- Use tximeta for automatic annotation linking when working with standard references
