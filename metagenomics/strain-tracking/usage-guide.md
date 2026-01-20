# Strain Tracking Usage Guide

Track bacterial strains at sub-species resolution.

## Prerequisites

```bash
# MASH
conda install -c bioconda mash

# sourmash
pip install sourmash

# fastANI
conda install -c bioconda fastani

# inStrain
pip install instrain
```

## Quick Start

### MASH - Fast Distance Estimation

```bash
# Create sketches
mash sketch -o reference genome1.fasta genome2.fasta

# Calculate distances
mash dist reference.msh query.fasta > distances.tsv

# Screen against RefSeq
mash screen refseq.msh reads.fastq.gz > screen_results.tsv
```

### sourmash - Metagenome Comparisons

```bash
# Compute signatures
sourmash sketch dna -p k=31,scaled=1000 *.fasta

# Compare all-vs-all
sourmash compare *.sig -o comparison.numpy

# Search against database
sourmash gather sample.sig gtdb-rs214.k31.zip
```

## Distance Interpretation

| MASH Distance | ANI | Interpretation |
|---------------|-----|----------------|
| 0.00 | 100% | Same strain |
| < 0.05 | > 95% | Same species |
| 0.05-0.10 | 90-95% | Related species |
| > 0.10 | < 90% | Different species |

## Strain-Level Metagenomics

### inStrain Profiling

```bash
# Profile strain variation
inStrain profile reads.bam reference.fasta -o profile_output

# Compare samples
inStrain compare profile1 profile2 -o comparison_output
```

### Key Metrics

- **popANI**: Population ANI across reads
- **conANI**: Consensus ANI
- **SNV density**: Variation within sample

## Tracking Workflows

### Outbreak Investigation

```bash
# Sketch all isolates
mash sketch -s 10000 -o outbreak *.fasta

# Pairwise distances
mash dist outbreak.msh outbreak.msh > pairwise.tsv

# Cluster close strains (< 0.001)
```

### Longitudinal Strain Tracking

```bash
# For each timepoint sample
for sample in timepoint*.bam; do
    inStrain profile $sample reference.fasta -o ${sample%.bam}_profile
done

# Compare across timepoints
inStrain compare *_profile -o longitudinal_comparison
```

## Large-Scale Comparisons

```bash
# fastANI for many genomes
fastANI \
    --ql genome_list.txt \
    --rl genome_list.txt \
    -o ani_results.tsv \
    -t 16

# Parse results
awk '$3 > 99' ani_results.tsv  # Near-identical strains
```

## See Also

- [MASH documentation](https://mash.readthedocs.io/)
- [sourmash documentation](https://sourmash.readthedocs.io/)
- [inStrain tutorial](https://instrain.readthedocs.io/)
