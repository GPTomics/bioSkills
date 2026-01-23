# long-read-sequencing

## Overview

Analysis of long-read sequencing data from Oxford Nanopore and PacBio. Covers alignment with minimap2, polishing, variant calling with medaka and Clair3, and structural variant detection with Sniffles.

**Tool type:** cli | **Primary tools:** minimap2, medaka, Clair3, Sniffles, NanoPlot

## Skills

| Skill | Description |
|-------|-------------|
| basecalling | Convert raw signal to sequences with Dorado/Guppy |
| long-read-alignment | Align long reads with minimap2 |
| medaka-polishing | Polish assemblies and call variants with medaka |
| clair3-variants | Deep learning variant calling with Clair3 |
| structural-variants | Detect SVs from long reads |
| long-read-qc | Quality control for long reads |
| isoseq-analysis | PacBio Iso-Seq isoform discovery and quantification |

## Example Prompts

- "Basecall my POD5 files with Dorado"
- "Convert FAST5 to sequences"
- "Run super-accuracy basecalling"
- "Align my Nanopore reads with minimap2"
- "Map PacBio HiFi reads to the reference genome"
- "Call variants with Clair3"
- "Polish my assembly with medaka"
- "Call variants from Nanopore reads with medaka"
- "Find structural variants from my long reads"
- "Detect deletions and insertions with Sniffles"
- "Check the quality of my Nanopore reads"
- "Process my Iso-Seq data"
- "Discover novel isoforms from PacBio"
- "QC my transcript assembly with SQANTI"

## Requirements

```bash
# Dorado (from ONT)
# Download from https://github.com/nanoporetech/dorado

# POD5 tools
pip install pod5

# Alignment and variant calling
conda install -c bioconda minimap2 medaka clair3

# SV callers
conda install -c bioconda sniffles cutesv

# QC tools
conda install -c bioconda nanoplot chopper

# Iso-Seq tools
conda install -c bioconda pbccs lima isoseq3 sqanti3
```

## Related Skills

- **alignment-files** - BAM manipulation
- **variant-calling** - Short-read variant calling
- **genome-assembly** - Long-read assembly
