# AMR Detection Usage Guide

Identify antimicrobial resistance genes in bacterial genomes and metagenomes.

## Prerequisites

```bash
# AMRFinderPlus
conda install -c bioconda ncbi-amrfinderplus
amrfinder -u  # Update database

# ResFinder
pip install resfinder

# ABRicate
conda install -c bioconda abricate
abricate --setupdb
```

## Quick Start

### AMRFinderPlus (Recommended)

```bash
# From assembled contigs
amrfinder \
    -n contigs.fasta \
    -o amr_results.tsv \
    --plus  # Include stress/virulence genes

# From protein sequences
amrfinder \
    -p proteins.faa \
    -o amr_results.tsv
```

### ABRicate (Quick Screening)

```bash
# Screen against multiple databases
abricate contigs.fasta --db resfinder > resfinder_results.tsv
abricate contigs.fasta --db card > card_results.tsv
abricate contigs.fasta --db ncbi > ncbi_results.tsv

# Summary across samples
abricate --summary *_results.tsv > summary.tsv
```

## Database Comparison

| Database | Focus | Best For |
|----------|-------|----------|
| NCBI (AMRFinderPlus) | Curated, comprehensive | General use |
| ResFinder | Acquired resistance | Clinical isolates |
| CARD | Mechanism-focused | Research |
| ARG-ANNOT | Acquired genes | Surveillance |

## Interpreting Results

### AMRFinderPlus Output

```
Gene symbol | Sequence name | Scope | Element type | Element subtype | Class | Subclass
```

Key fields:
- **Gene symbol**: Standardized resistance gene name
- **Class**: Drug class (e.g., aminoglycoside, beta-lactam)
- **Subclass**: Specific drug (e.g., gentamicin)
- **% Coverage/Identity**: Match quality

### Resistance Mechanisms

| Type | Example | Implication |
|------|---------|-------------|
| Acquired | blaCTX-M | Horizontal transfer |
| Mutational | gyrA | Chromosomal |
| Efflux | mexAB-oprM | Broad spectrum |

## Metagenomic Analysis

```bash
# Assemble first, then screen
megahit -1 reads_1.fq.gz -2 reads_2.fq.gz -o assembly
amrfinder -n assembly/final.contigs.fa -o metagenome_amr.tsv

# Or use read-based approach
shortbred_identify.py \
    --goi card_proteins.faa \
    --ref uniref90.fasta \
    --markers markers.faa
```

## Batch Processing

```bash
for fasta in *.fasta; do
    sample=$(basename $fasta .fasta)
    amrfinder -n $fasta -o ${sample}_amr.tsv --plus
done

# Combine results
cat *_amr.tsv | head -1 > all_amr.tsv
cat *_amr.tsv | grep -v "^Protein" >> all_amr.tsv
```

## See Also

- [AMRFinderPlus wiki](https://github.com/ncbi/amr/wiki)
- [CARD database](https://card.mcmaster.ca/)
