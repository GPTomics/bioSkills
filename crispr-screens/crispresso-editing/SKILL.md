---
name: bio-crispr-screens-crispresso-editing
description: CRISPResso2 for analyzing CRISPR gene editing outcomes. Quantifies indels, HDR efficiency, and generates comprehensive editing reports. Use when analyzing amplicon sequencing data from CRISPR editing experiments to assess editing efficiency.
tool_type: cli
primary_tool: CRISPResso2
---

# CRISPResso2 Editing Analysis

## Basic Analysis

```bash
# Analyze single amplicon
CRISPResso \
    --fastq_r1 sample_R1.fastq.gz \
    --fastq_r2 sample_R2.fastq.gz \
    --amplicon_seq AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT \
    --guide_seq CTGCCCACAGGTGAGGAGGT \
    --output_folder crispresso_output \
    --name sample1

# Output includes:
# - Editing efficiency statistics
# - Indel distribution
# - Allele frequency plots
```

## With HDR Template

```bash
# Analyze HDR editing
CRISPResso \
    --fastq_r1 hdr_sample_R1.fastq.gz \
    --fastq_r2 hdr_sample_R2.fastq.gz \
    --amplicon_seq AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGAGGCTCGGCTGTGGTT \
    --guide_seq CTGCCCACAGGTGAGGAGGT \
    --expected_hdr_amplicon_seq AATGTCCCCCAATGGGAAGTTCATCTGGCACTGCCCACAGGTGAGGAGGTCATGATCCCCTTCTGGAGCTCCCAACGGGCCGTGGTCTGGTTCATCATCTGTAAGAATGGCTTCAAGATGCTCGGCTGTGGTT \
    --output_folder hdr_output \
    --name hdr_sample
```

## Batch Analysis

```bash
# Create batch file (tab-separated)
# batch.txt:
# name    fastq_r1    fastq_r2    amplicon_seq    guide_seq
# sample1 s1_R1.fq.gz s1_R2.fq.gz AMPLICON1       GUIDE1
# sample2 s2_R1.fq.gz s2_R2.fq.gz AMPLICON2       GUIDE2

CRISPRessoBatch \
    --batch_settings batch.txt \
    --output_folder batch_output \
    --n_processes 8
```

## Pool Analysis (Multiple Guides)

```bash
# Analyze pooled amplicons
CRISPRessoPooled \
    --fastq_r1 pooled_R1.fastq.gz \
    --fastq_r2 pooled_R2.fastq.gz \
    --amplicon_file amplicons.txt \
    --output_folder pooled_output \
    --n_processes 8

# amplicons.txt format:
# amplicon_name    amplicon_seq    guide_seq
```

## WGS Analysis

```bash
# Analyze off-target editing from WGS
CRISPRessoWGS \
    --bam aligned.bam \
    --reference genome.fa \
    --regions_file targets.bed \
    --output_folder wgs_output
```

## Parse Results in Python

```python
import pandas as pd
import json

# Load mapping statistics
with open('crispresso_output/CRISPResso_mapping_statistics.txt') as f:
    stats = {}
    for line in f:
        key, value = line.strip().split('\t')
        stats[key] = value

print(f"Reads aligned: {stats['READS_ALIGNED']}")
print(f"Reads aligned %: {stats['READS_ALIGNED_PERCENTAGE']}")

# Load quantification
quant = pd.read_csv('crispresso_output/CRISPResso_quantification_of_editing_frequency.txt', sep='\t')
print(quant)

# Load allele frequency
alleles = pd.read_csv('crispresso_output/Alleles_frequency_table.zip', compression='zip', sep='\t')
print(f"Unique alleles: {len(alleles)}")
print(alleles.head(10))
```

## Key Output Files

```
CRISPResso_output/
├── CRISPResso_mapping_statistics.txt    # Read mapping stats
├── CRISPResso_quantification_of_editing_frequency.txt  # Summary
├── Alleles_frequency_table.zip          # All allele sequences
├── CRISPResso_RUNNING_LOG.txt           # Analysis log
├── Indel_histogram.png                  # Indel size distribution
├── Insertion_deletion_substitution.png  # Edit type pie chart
├── Alleles_frequency_table.png          # Top allele bar plot
└── CRISPResso2_info.json               # Machine-readable summary
```

## Quantify Specific Outcomes

```bash
# Define expected outcomes
CRISPResso \
    --fastq_r1 sample_R1.fastq.gz \
    --amplicon_seq AMPLICON \
    --guide_seq GUIDE \
    --coding_seq CODING_REGION \
    --quantification_window_size 5 \
    --quantification_window_center -3 \
    --output_folder output
```

## Base Editing Analysis

```bash
# For base editors (CBE/ABE)
CRISPResso \
    --fastq_r1 base_edit_R1.fastq.gz \
    --amplicon_seq AMPLICON \
    --guide_seq GUIDE \
    --base_editor_output \
    --conversion_nuc_from C \
    --conversion_nuc_to T \
    --output_folder base_edit_output
```

## Prime Editing Analysis

```bash
# For prime editing
CRISPResso \
    --fastq_r1 prime_edit_R1.fastq.gz \
    --amplicon_seq AMPLICON \
    --guide_seq GUIDE \
    --prime_editing_pegRNA_spacer_seq SPACER \
    --prime_editing_pegRNA_extension_seq EXTENSION \
    --prime_editing_pegRNA_scaffold_seq SCAFFOLD \
    --output_folder prime_edit_output
```

## Compare Samples

```bash
# Compare two CRISPResso runs
CRISPRessoCompare \
    --crispresso_output_folder_1 sample1_output \
    --crispresso_output_folder_2 sample2_output \
    --output_folder comparison_output
```

## Related Skills

- screen-qc - QC for editing experiments
- read-alignment/bwa-alignment - Align reads for WGS analysis
- variant-calling/variant-calling - Detect editing-induced variants
