# Long-Read Sequencing

Analysis of long-read sequencing data from Oxford Nanopore (ONT) and PacBio using minimap2, medaka, and associated tools.

## Overview

This category covers long-read sequencing analysis: minimap2 for fast alignment, medaka for Oxford Nanopore consensus and variant calling, and tools for structural variant detection. Long reads enable assembly of complex regions, detection of structural variants, and direct detection of base modifications.

**Tool type:** `cli`
**Primary tools:** minimap2, medaka, samtools

## Skills

| Skill | Description |
|-------|-------------|
| [long-read-alignment](long-read-alignment/) | Align long reads with minimap2 |
| [medaka-polishing](medaka-polishing/) | Polish assemblies and call variants with medaka |
| [structural-variants](structural-variants/) | Detect SVs from long reads |
| [long-read-qc](long-read-qc/) | Quality control for long reads |

## Workflow

```
Raw Long Reads (FASTQ/FAST5/POD5)
    |
    v
[long-read-qc] -----------> Quality metrics, filtering
    |
    v
[long-read-alignment] ----> Align to reference (minimap2)
    |
    +-------------------+
    |                   |
    v                   v
[medaka-polishing]     [structural-variants]
    |                   |
    v                   v
Consensus/Variants     SV calls (DEL, INS, INV, DUP)
```

## Platform Comparison

| Feature | Oxford Nanopore | PacBio HiFi |
|---------|-----------------|-------------|
| Read length | 10kb - 2Mb | 10-25kb |
| Accuracy | ~95-99% (Q20+) | >99% (Q30+) |
| Error type | Systematic | Random |
| Throughput | High | Moderate |
| Base mods | Direct detection | Kinetics |
| minimap2 preset | -ax map-ont | -ax map-hifi |

## Example Prompts

### Alignment
- "Align my Nanopore reads with minimap2"
- "Map PacBio HiFi reads to the reference genome"
- "Create a sorted BAM from long-read alignment"

### Polishing
- "Polish my assembly with medaka"
- "Call variants from Nanopore reads with medaka"
- "Generate consensus sequence for my region"

### Structural Variants
- "Find structural variants from my long reads"
- "Detect deletions and insertions with sniffles"
- "Call SVs from PacBio alignments"

### Quality Control
- "Check the quality of my Nanopore reads"
- "Generate read length distribution"
- "Filter reads by quality score"

## Requirements

### minimap2

```bash
# Conda installation
conda install -c bioconda minimap2

# Verify
minimap2 --version
```

### medaka

```bash
# Conda installation (recommended)
conda create -n medaka -c conda-forge -c bioconda medaka

# Or pip (may need dependencies)
pip install medaka

# Verify
medaka --version
```

### Additional Tools

```bash
# Structural variants
conda install -c bioconda sniffles cutesv

# Quality control
conda install -c bioconda nanoplot nanostat chopper

# PacBio tools
conda install -c bioconda pbmm2 pbsv
```

## Key Functions

| Tool | Command | Purpose |
|------|---------|---------|
| minimap2 | minimap2 -ax | Align reads |
| medaka | medaka_consensus | Polish assembly |
| medaka | medaka_variant | Haploid variant calling (v2.0+) |
| sniffles | sniffles | SV detection |
| cuteSV | cuteSV | SV detection |
| NanoPlot | NanoPlot | QC visualization |
| chopper | chopper | Read filtering |

## minimap2 Presets

| Preset | Use Case |
|--------|----------|
| map-ont | ONT reads to reference |
| map-hifi | PacBio HiFi to reference |
| map-pb | PacBio CLR to reference |
| splice | Long reads to transcriptome |
| asm5 | Assembly to assembly (~0.1% divergence) |
| asm20 | Assembly to assembly (~5% divergence) |

## Notes

- **Memory requirements** - minimap2 loads index into RAM (~3GB for human)
- **medaka models** - Use correct model for your basecaller version
- **SV calling** - Requires good coverage (>10x recommended)
- **BAM sorting** - Always sort and index BAM files after alignment
- **Reference indexing** - minimap2 can index on-the-fly but pre-indexing is faster

## Related Skills

- **alignment-files** - BAM manipulation
- **variant-calling** - Short-read variant calling
- **sequence-io** - FASTQ handling

## References

- [minimap2 GitHub](https://github.com/lh3/minimap2)
- [medaka GitHub](https://github.com/nanoporetech/medaka)
- [Sniffles GitHub](https://github.com/fritzsedlazeck/Sniffles)
- [NanoPlot GitHub](https://github.com/wdecoster/NanoPlot)
