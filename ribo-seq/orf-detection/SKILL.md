---
name: bio-ribo-seq-orf-detection
description: Detect and quantify translated ORFs from Ribo-seq data including uORFs and novel ORFs using RiboCode and ORFquant. Use when identifying translated regions beyond annotated coding sequences or quantifying ORF-level translation.
tool_type: mixed
primary_tool: RiboCode
---

# ORF Detection

## RiboCode Workflow

```bash
# Step 1: Prepare annotation
prepare_transcripts \
    -g annotation.gtf \
    -f genome.fa \
    -o ribocode_annot

# Step 2: Run RiboCode
RiboCode \
    -a ribocode_annot \
    -c config.txt \
    -l 27,28,29,30 \
    -o output_prefix

# config.txt format:
# SampleName  AlignmentFile  Stranded
# sample1     sample1.bam    yes
```

## One-Step RiboCode

```bash
# All-in-one command
RiboCode_onestep \
    -g annotation.gtf \
    -r riboseq.bam \
    -f genome.fa \
    -l 27,28,29,30 \
    -o output_dir
```

## RiboCode Output

| File | Description |
|------|-------------|
| *_ORF_result.txt | Detected ORFs with coordinates |
| *_ORF_result.html | Interactive visualization |
| *_binomial_test.txt | Statistical test results |

## Parse RiboCode Results

```python
import pandas as pd

def load_ribocode_orfs(filepath):
    '''Load RiboCode ORF predictions'''
    df = pd.read_csv(filepath, sep='\t')

    # ORF categories
    categories = {
        'annotated': df[df['ORF_type'] == 'annotated'],
        'uORF': df[df['ORF_type'] == 'uORF'],
        'dORF': df[df['ORF_type'] == 'dORF'],
        'novel': df[df['ORF_type'].isin(['novel', 'noncoding'])]
    }

    return df, categories
```

## Alternative: RibORF

```bash
# RibORF uses random forest classifier
RibORF.py \
    -f genome.fa \
    -r riboseq.bam \
    -g annotation.gtf \
    -o output_dir
```

## Manual ORF Detection

```python
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence, min_length=30):
    '''Find all ORFs in a sequence'''
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']

    orfs = []
    seq = str(sequence).upper()

    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon == start_codon:
                # Find next stop codon
                for j in range(i + 3, len(seq) - 2, 3):
                    if seq[j:j+3] in stop_codons:
                        orf_length = j - i + 3
                        if orf_length >= min_length:
                            orfs.append({
                                'start': i,
                                'end': j + 3,
                                'frame': frame,
                                'length': orf_length,
                                'sequence': seq[i:j+3]
                            })
                        break

    return orfs

def detect_translated_orfs(orfs, coverage_data, min_coverage=10):
    '''Filter ORFs by Ribo-seq coverage'''
    translated = []
    for orf in orfs:
        cov = coverage_data[orf['start']:orf['end']]
        if sum(cov) >= min_coverage:
            translated.append(orf)
    return translated
```

## uORF Analysis

```python
def find_uorfs(transcript, cds_start):
    '''Find upstream ORFs before main CDS'''
    utr5 = transcript[:cds_start]
    uorfs = find_orfs(utr5)

    # Classify uORFs
    for uorf in uorfs:
        if uorf['end'] <= cds_start:
            uorf['type'] = 'contained'  # Fully in 5' UTR
        else:
            uorf['type'] = 'overlapping'  # Overlaps main CDS

    return uorfs
```

## ORF Categories

| Type | Description |
|------|-------------|
| annotated | Known CDS in annotation |
| uORF | Upstream of main CDS |
| dORF | Downstream of main CDS |
| internal | Within CDS, different frame |
| noncoding | In annotated non-coding RNA |
| novel | Unannotated region |

## ORFquant for ORF Quantification

ORFquant provides transcript-level and ORF-level quantification from Ribo-seq data.

### Installation

```r
# Install from Bioconductor
BiocManager::install('ORFik')
# ORFquant is part of the ORFik ecosystem
```

### Basic ORF Quantification

```r
library(ORFik)
library(GenomicFeatures)

# Load annotation
txdb <- makeTxDbFromGFF('annotation.gtf')

# Load Ribo-seq data
riboseq <- fimport('riboseq.bam')

# Get CDS regions
cds <- cdsBy(txdb, by = 'tx', use.names = TRUE)

# Calculate ORF-level RPKM
# fpkm: Fragments Per Kilobase per Million mapped reads
orf_counts <- countOverlaps(cds, riboseq)
orf_lengths <- sum(width(cds))
total_reads <- length(riboseq)
orf_fpkm <- (orf_counts * 1e9) / (orf_lengths * total_reads)
```

### P-site Corrected Quantification

```r
library(ORFik)

# Load with P-site offset correction
# p_offsets=c(12,12,12): P-site offset for 28-30nt reads. Determine from metagene.
riboseq <- fimport('riboseq.bam', p_offsets = c(12, 12, 12), lengths = 28:30)

# Count P-sites per ORF
psite_counts <- countOverlaps(cds, riboseq)
```

### Detect and Quantify Novel ORFs

```r
library(ORFik)

# Find candidate ORFs in 5' UTRs
utr5 <- fiveUTRsByTranscript(txdb, use.names = TRUE)
uorf_candidates <- findORFs(utr5, startCodon = 'ATG', longestORF = FALSE,
                            minimumLength = 9)  # 9 codons minimum

# Quantify uORFs
uorf_counts <- countOverlaps(uorf_candidates, riboseq)

# Filter by coverage
# min_count=10: Minimum reads for confident detection.
active_uorfs <- uorf_candidates[uorf_counts >= 10]
```

### ORFquant Output Interpretation

```r
# Create ORF summary table
orf_summary <- data.frame(
    orf_id = names(cds),
    length = sum(width(cds)),
    counts = orf_counts,
    fpkm = orf_fpkm
)

# Classify by expression
# fpkm>1: Low expression threshold. Adjust based on library depth.
orf_summary$expressed <- orf_summary$fpkm > 1

write.csv(orf_summary, 'orf_quantification.csv', row.names = FALSE)
```

### Compare ORF Expression Across Conditions

```r
library(DESeq2)

# Build count matrix for multiple samples
orf_count_matrix <- cbind(
    sample1 = countOverlaps(cds, riboseq1),
    sample2 = countOverlaps(cds, riboseq2),
    sample3 = countOverlaps(cds, riboseq3),
    sample4 = countOverlaps(cds, riboseq4)
)

# Run DESeq2 for differential translation
coldata <- data.frame(condition = c('control', 'control', 'treatment', 'treatment'))
dds <- DESeqDataSetFromMatrix(orf_count_matrix, coldata, ~ condition)
dds <- DESeq(dds)
results <- results(dds)
```

## Related Skills

- ribosome-periodicity - Validate ORF calling
- translation-efficiency - Quantify ORF translation
- differential-expression - Compare ORF expression
