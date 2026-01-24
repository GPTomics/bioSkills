# Basecalling - Usage Guide

## Overview
Basecalling converts raw Nanopore signal data (FAST5/POD5) to nucleotide sequences. Dorado is the current production basecaller and should be used for all new analyses. Guppy is deprecated and no longer receiving updates.

## Prerequisites
```bash
# Dorado (from ONT community)
# Download from https://github.com/nanoporetech/dorado

# POD5 tools
pip install pod5

# Quality filtering
conda install -c bioconda chopper nanoplot
```

## Quick Start
Tell your AI agent what you want to do:
- "Basecall my POD5 files with Dorado using SUP accuracy"
- "Basecall my FAST5 files with Dorado"

## Example Prompts

### Basic Basecalling
> "Basecall my POD5 files in the raw_data/ folder using Dorado with the SUP model"

> "Run Dorado on my R10.4.1 data with HAC accuracy for faster processing"

### Model Selection
> "What Dorado model should I use for R10.4.1 chemistry?"

> "Basecall with the fast model for a quick preview of my data"

### With Modifications
> "Basecall and detect 5mC methylation using Dorado"

### Quality Filtering After
> "Basecall my data and then filter for Q10+ reads"

## What the Agent Will Do
1. Identify your flow cell chemistry (R10.4.1 vs R9.4.1)
2. Select appropriate model based on speed/accuracy tradeoff
3. Run basecalling with Dorado
4. Optionally filter output for quality

## Model Selection

Choose based on speed vs accuracy tradeoff:
- **fast**: Quick preview, lowest accuracy
- **hac**: Balanced, good for most uses
- **sup**: Highest accuracy, slowest

## Chemistry Matching

Match the model to your flow cell chemistry:
- R10.4.1: Current chemistry
- R9.4.1: Legacy chemistry

## Workflow Position

```
Raw Signal (FAST5/POD5)
    |
    v
Basecalling (Dorado)
    |
    v
FASTQ/BAM
    |
    v
Quality Filtering (chopper)
    |
    v
Alignment (minimap2)
```

## Tips
- Use SUP model for final analysis, HAC for exploratory work, fast for quick checks
- Dorado can output BAM directly with alignment if reference provided
- Use `--batchsize 32` if running out of GPU memory
- Check chemistry (R10.4.1 vs R9.4.1) before selecting model
- POD5 is the new format; FAST5 is legacy but still supported

## Resources
- [Dorado GitHub](https://github.com/nanoporetech/dorado)
- [ONT Community](https://community.nanoporetech.com/)
