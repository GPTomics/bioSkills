# Functional Profiling Usage Guide

HUMAnN3 (HMP Unified Metabolic Analysis Network) profiles the functional potential of metagenomic communities.

## What It Does

1. **Maps reads** to species-specific pangenomes (via MetaPhlAn)
2. **Translates unmapped reads** to protein database
3. **Quantifies gene families** (UniRef90 clusters)
4. **Infers pathway abundance** (MetaCyc pathways)

## Installation

```bash
conda create -n humann -c bioconda humann
conda activate humann

# Download databases (~16 GB total)
humann_databases --download chocophlan full /db/humann
humann_databases --download uniref uniref90_diamond /db/humann
```

## Quick Start

```bash
# Single sample
humann --input sample.fq.gz --output sample_out --threads 8

# Multiple samples
for fq in *.fq.gz; do
    humann --input $fq --output $(basename $fq .fq.gz)_out -t 8
done

# Join and normalize
humann_join_tables -i . -o merged_pathways.tsv --file_name pathabundance
humann_renorm_table -i merged_pathways.tsv -o pathways_relab.tsv -u relab
```

## Output Interpretation

### Gene Families (UniRef90)
- Abundance in RPKs (reads per kilobase)
- Can be regrouped to EC, KO, GO, Pfam

### Pathway Abundance
- MetaCyc pathway abundances
- Stratified by contributing species
- UNMAPPED = reads not mapping to any gene
- UNINTEGRATED = genes not in any pathway

### Pathway Coverage
- Fraction of pathway genes detected (0-1)
- Low coverage = incomplete pathway

## Key Parameters

| Parameter | Description |
|-----------|-------------|
| `--threads` | Number of CPUs |
| `--memory-use` | Memory limit (minimum/maximum) |
| `--taxonomic-profile` | Pre-computed MetaPhlAn profile |
| `--bypass-nucleotide-search` | Skip pangenome search |

## Database Options

| Database | Size | Speed | Sensitivity |
|----------|------|-------|-------------|
| uniref90_diamond | 16GB | Fast | Standard |
| uniref50_diamond | 5GB | Faster | Lower |
| uniref90_ec_filtered | 0.8GB | Fastest | EC only |

## Tips

1. **Use MetaPhlAn pre-profile** for speed if you have it
2. **Concatenate paired-end reads** - HUMAnN handles both orientations
3. **Check UNMAPPED rate** - high rates indicate missing database coverage
4. **Normalize before comparison** - use relab or CPM
5. **Consider stratification** - species contributions add biological insight
