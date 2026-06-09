---
name: bio-workflows-metagenomics-pipeline
description: End-to-end shotgun metagenomics workflow from FASTQ to taxonomic and functional profiles, orchestrating controls/host depletion, Kraken2+Bracken classification, MetaPhlAn marker profiling, and HUMAnN functional profiling. Covers the controls-first ordering, why Kraken2 read counts are not abundances and MetaPhlAn cell fractions do not equal Bracken read fractions, and the consistent-pipeline framing. Use when profiling shotgun metagenomic samples end to end, or chaining classification, abundance, and function. For resistome see metagenomics/amr-detection; for strains see metagenomics/strain-tracking; for assembly see genome-assembly/metagenome-assembly.
tool_type: cli
primary_tool: Kraken2
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - metagenomics/contamination-controls
  - metagenomics/kraken-classification
  - metagenomics/metaphlan-profiling
  - metagenomics/abundance-estimation
  - metagenomics/functional-profiling
  - metagenomics/metagenome-visualization
qc_checkpoints:
  - after_qc: "Q30 >80%, host reads removed"
  - after_classification: "Classification rate >60%, known taxa dominant"
  - after_functional: "Pathway coverage reasonable, unmapped <50%"
---

## Version Compatibility

Reference examples tested with: Bowtie2 2.5.3+, Bracken 2.9+, HUMAnN 3.8+, Kraken2 2.1+, MetaPhlAn 4.1+, fastp 0.23+, matplotlib 3.8+, pandas 2.2+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metagenomics Pipeline

**"Analyze my metagenomic samples from FASTQ to taxonomic and functional profiles"** -> Orchestrate controls and host depletion, Kraken2/Bracken taxonomic classification, MetaPhlAn profiling, and HUMAnN3 functional analysis - reporting results relative to a consistent pipeline, never as a direct observation of the community.

Complete workflow from metagenomic FASTQ to taxonomic and functional profiles. Every result is a position in a choice-chain (extraction -> depletion -> depth -> classifier -> database -> normalization); hold the chain constant across a study and report each link.

## Workflow Overview

```
FASTQ files (+ extraction blanks, mock)
    |
    v
[0. QC, Host Removal & Controls] --> fastp + Hostile/Bowtie2(T2T) + blanks/decontam + Nonpareil depth check
    |
    v
[1. Taxonomic Classification]
    |
    +---> Kraken2 (+confidence, +hit-groups) + Bracken -> read fraction
    |
    +---> MetaPhlAn 4 (marker-based, pinned --index) -> cell fraction (NOT comparable to Bracken %)
    |
    v
[2. Functional Profiling] --> HUMAnN (potential, not activity; keep UNMAPPED)
    |
    v
Taxonomic profiles + Pathway abundances (+ AMR/strain via their own skills)
```

## Primary Path: Kraken2 + Bracken + HUMAnN

### Step 0: Quality Control, Host Removal, and Controls

Carry extraction blanks and a mock through the whole workflow; host-deplete against T2T-CHM13; confirm depth with Nonpareil. See metagenomics/contamination-controls for the controls/decontam detail.

```bash
# QC with fastp (trimming mechanics: read-qc/fastp-workflow)
for sample in sample1 sample2 sample3; do
    fastp -i ${sample}_R1.fastq.gz -I ${sample}_R2.fastq.gz \
        -o trimmed/${sample}_R1.fq.gz -O trimmed/${sample}_R2.fq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --html qc/${sample}_fastp.html
done

# Remove host reads - Hostile with a T2T-CHM13 index removes >99.5% host with low microbial loss.
# Report the reads removed; host depletion can halve usable depth.
for sample in sample1 sample2 sample3; do
    hostile clean --fastq1 trimmed/${sample}_R1.fq.gz --fastq2 trimmed/${sample}_R2.fq.gz \
        --index human-t2t-hla --aligner bowtie2 --output host_removed/
done
# Then run decontam on the classifier output table using the blanks (contamination-controls),
# and confirm depth adequacy with Nonpareil before interpreting any non-detection.
```

### Step 2A: Kraken2 Classification

```bash
# Classify reads. Raise --confidence above the default 0 to suppress single-k-mer false positives,
# and require >=2 hit groups. The database defines what can be detected.
for sample in sample1 sample2 sample3; do
    kraken2 --db kraken2_db \
        --threads 8 \
        --paired \
        --confidence 0.1 \
        --minimum-hit-groups 2 \
        --report kraken/${sample}.report \
        --output kraken/${sample}.output \
        host_removed/${sample}_R1.fq.gz \
        host_removed/${sample}_R2.fq.gz
done
```

### Step 2B: Bracken Abundance Estimation

```bash
# Estimate species abundance
for sample in sample1 sample2 sample3; do
    bracken -d kraken2_db \
        -i kraken/${sample}.report \
        -o bracken/${sample}.species.txt \
        -r 150 \
        -l S \
        -t 10
done

# Combine samples into abundance matrix
combine_bracken_outputs.py \
    --files bracken/*.species.txt \
    -o bracken/combined_species.txt
```

### Step 2C: Alternative - MetaPhlAn Profiling

```bash
# Profile with MetaPhlAn 4. Pin --index (DB version is a batch variable). MetaPhlAn % is a cell
# fraction - do NOT merge it with Bracken read fractions. In 4.2 --bowtie2out is renamed --mapout.
for sample in sample1 sample2 sample3; do
    metaphlan host_removed/${sample}_R1.fq.gz,host_removed/${sample}_R2.fq.gz \
        --bowtie2out metaphlan/${sample}.bowtie2.bz2 \
        --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
        --input_type fastq \
        --nproc 8 \
        -o metaphlan/${sample}_profile.txt
done

# Merge profiles
merge_metaphlan_tables.py metaphlan/*_profile.txt > metaphlan/merged_abundance.txt
```

### Step 3: Functional Profiling with HUMAnN

```bash
# Run HUMAnN
for sample in sample1 sample2 sample3; do
    # Concatenate paired reads
    cat host_removed/${sample}_R1.fq.gz host_removed/${sample}_R2.fq.gz > \
        host_removed/${sample}_concat.fq.gz

    humann --input host_removed/${sample}_concat.fq.gz \
        --output humann/${sample} \
        --threads 8 \
        --metaphlan-options "--bowtie2db metaphlan_db"
done

# Normalize and join tables
humann_renorm_table --input humann/sample1/sample1_pathabundance.tsv \
    --output humann/sample1/sample1_pathabundance_cpm.tsv \
    --units cpm

humann_join_tables --input humann \
    --output humann/merged_pathabundance.tsv \
    --file_name pathabundance
```

### Visualization

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load Bracken species table
species = pd.read_csv('bracken/combined_species.txt', sep='\t', index_col=0)

# Top 20 species heatmap
top20 = species.sum(axis=1).nlargest(20).index
plt.figure(figsize=(12, 8))
sns.heatmap(species.loc[top20], cmap='viridis', annot=False)
plt.title('Top 20 Species Abundance')
plt.tight_layout()
plt.savefig('top20_species_heatmap.pdf')

# Stacked bar plot
species_norm = species.div(species.sum()) * 100
top10 = species_norm.sum(axis=1).nlargest(10).index
other = species_norm.loc[~species_norm.index.isin(top10)].sum()

plot_data = species_norm.loc[top10].T
plot_data['Other'] = other
plot_data.plot(kind='bar', stacked=True, figsize=(10, 6))
plt.ylabel('Relative Abundance (%)')
plt.legend(bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.savefig('species_barplot.pdf')
```

## Parameter Recommendations

| Step | Parameter | Value |
|------|-----------|-------|
| fastp | --length_required | 50 (metagenomic reads) |
| Kraken2 | --confidence | 0.1-0.4 (default 0.0 over-classifies; see metagenomics/kraken-classification) |
| Kraken2 | --minimum-hit-groups | 2 (cut single-region false positives) |
| Bracken | -r | Read length (e.g., 150; must match the DB build) |
| Bracken | -l | S (species) or G (genus) |
| Bracken | -t | 10 (min reads threshold) |
| MetaPhlAn | --min_cu_len | 2000 (default) |
| HUMAnN | --threads | 8+ |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| Low classification rate | Database mismatch, novel organisms | Try different database, check sample type |
| High unclassified | Novel microbes, host contamination | Remove host, use larger database |
| High host reads | Incomplete host removal | Use multiple host reference genomes |
| HUMAnN slow | Large files | Increase threads, pre-filter reads |

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

THREADS=8
KRAKEN_DB="kraken2_standard_db"
HOST_INDEX="human_bt2_index"
SAMPLES="sample1 sample2 sample3"
OUTDIR="metagenomics_results"

mkdir -p ${OUTDIR}/{trimmed,host_removed,kraken,bracken,metaphlan,humann,qc}

# Step 1: QC
echo "=== QC ==="
for sample in $SAMPLES; do
    fastp -i ${sample}_R1.fastq.gz -I ${sample}_R2.fastq.gz \
        -o ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -O ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --length_required 50 \
        --html ${OUTDIR}/qc/${sample}_fastp.html -w ${THREADS}
done

# Host removal
echo "=== Host Removal ==="
for sample in $SAMPLES; do
    bowtie2 -p ${THREADS} -x ${HOST_INDEX} \
        -1 ${OUTDIR}/trimmed/${sample}_R1.fq.gz \
        -2 ${OUTDIR}/trimmed/${sample}_R2.fq.gz \
        --un-conc-gz ${OUTDIR}/host_removed/${sample}_R%.fq.gz \
        > /dev/null 2> ${OUTDIR}/qc/${sample}_host.log
done

# Step 2: Kraken2
echo "=== Kraken2 ==="
for sample in $SAMPLES; do
    kraken2 --db ${KRAKEN_DB} --threads ${THREADS} --paired \
        --confidence 0.1 --minimum-hit-groups 2 \
        --report ${OUTDIR}/kraken/${sample}.report \
        --output ${OUTDIR}/kraken/${sample}.output \
        ${OUTDIR}/host_removed/${sample}_R1.fq.gz \
        ${OUTDIR}/host_removed/${sample}_R2.fq.gz
done

# Bracken
echo "=== Bracken ==="
for sample in $SAMPLES; do
    bracken -d ${KRAKEN_DB} \
        -i ${OUTDIR}/kraken/${sample}.report \
        -o ${OUTDIR}/bracken/${sample}.species.txt \
        -r 150 -l S -t 10
done

echo "=== Pipeline Complete ==="
echo "Kraken reports: ${OUTDIR}/kraken/"
echo "Bracken abundances: ${OUTDIR}/bracken/"
```

## Related Skills

- database-access/sra-data - Pull metagenomic FASTQ from SRA / ENA (16S amplicon or shotgun)
- database-access/ncbi-datasets-cli - Bulk-pull reference genomes for read mapping
- database-access/remote-homology - DIAMOND --ultra-sensitive for predicted-ORF annotation
- metagenomics/contamination-controls - Host depletion, blanks/decontam, depth checks up front
- metagenomics/kraken-classification - Kraken2 details
- metagenomics/metaphlan-profiling - MetaPhlAn parameters
- metagenomics/abundance-estimation - Bracken options and compositional handling
- metagenomics/functional-profiling - HUMAnN workflow
- metagenomics/amr-detection - Community resistome from the same reads
- metagenomics/strain-tracking - Strain resolution from the same reads
- metagenomics/metagenome-visualization - Plotting and community statistics
