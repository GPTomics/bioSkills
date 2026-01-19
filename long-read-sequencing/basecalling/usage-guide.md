# Basecalling Usage Guide

This guide covers converting raw Nanopore signal data to nucleotide sequences.

## Dorado vs Guppy

- **Dorado**: Current production basecaller, recommended
- **Guppy**: Legacy, being deprecated

## Workflow Position

```
Raw Signal (FAST5/POD5)
    |
    v
Basecalling (Dorado/Guppy)
    |
    v
FASTQ/BAM
    |
    v
Quality Filtering (Chopper)
    |
    v
Alignment (minimap2)
```

## Model Selection

Choose based on speed vs accuracy tradeoff:

- **fast**: Quick preview, lowest accuracy
- **hac**: Balanced, good for most uses
- **sup**: Highest accuracy, slowest

## Chemistry Matching

Match the model to your flow cell chemistry:
- R10.4.1: Current chemistry
- R9.4.1: Legacy chemistry

## Requirements

```bash
# Dorado (from ONT community)
# Download from https://github.com/nanoporetech/dorado

# POD5 tools
pip install pod5

# Quality filtering
conda install -c bioconda chopper nanoplot

# Guppy (requires ONT account)
# Download from ONT community site
```

## Common Issues

### GPU Memory
Use smaller batch size: `--batchsize 32`

### Wrong Model
Check chemistry: R10.4.1 vs R9.4.1

### Slow Performance
Use GPU if available; use fast model for preview
