# long-read-sequencing

## Overview

Analysis of Oxford Nanopore and PacBio long-read data: basecalling, QC, alignment, polishing, small-variant and structural-variant calling, read-backed phasing, direct methylation detection, and full-length isoform analysis.

**Tool type:** cli | **Primary tools:** Dorado, minimap2, medaka, Clair3, Sniffles2, WhatsHap, modkit, SQANTI3

## Skills

| Skill | Description |
|-------|-------------|
| basecalling | Basecall raw signal (POD5/FAST5) with Dorado, requesting modifications at basecall time |
| long-read-qc | Assess read length/quality/identity and run health; filter for the downstream goal |
| long-read-alignment | Map reads with the error-rate-matched minimap2 preset, preserving SV/methylation tags |
| medaka-polishing | Polish ONT assemblies with the basecaller-matched medaka consensus model |
| clair3-variants | Call germline small variants with Clair3 using the chemistry-matched model |
| haplotype-phasing | Read-backed phasing (WhatsHap/LongPhase) and haplotagging the BAM (HP/PS) |
| structural-variants | Detect SVs with Sniffles2/cuteSV; joint-genotype cohorts; Truvari benchmarking |
| nanopore-methylation | Pile up 5mC/5hmC/6mA from MM/ML tags into bedMethyl with modkit |
| isoseq-analysis | Discover, classify, and filter full-length isoforms (PacBio + ONT) with SQANTI3/pigeon |

## Example Prompts

- "Basecall my POD5 files with Dorado and call 5mC methylation"
- "Check the quality and real identity of my Nanopore reads"
- "Align my R10 Nanopore reads with the right minimap2 preset"
- "Polish my ONT assembly with medaka"
- "Call germline variants with Clair3 using the matching model"
- "Phase my variants and haplotag the BAM"
- "Find structural variants and joint-genotype my cohort"
- "Call methylation from my nanopore BAM file"
- "Discover and curate full-length isoforms from my Iso-Seq data"

## Requirements

```bash
# Dorado (from ONT): https://github.com/nanoporetech/dorado
pip install pod5

# QC, alignment, polishing, variant calling
conda install -c bioconda nanoplot cramino chopper filtlong seqkit \
    minimap2 samtools medaka clair3 whatshap longphase

# SV calling and benchmarking
conda install -c bioconda sniffles cutesv truvari

# Methylation
# modkit: https://github.com/nanoporetech/modkit

# Isoform analysis (PacBio + ONT)
conda install -c bioconda isoseq pbmm2 pigeon lima sqanti3 isoquant
```

## Related Skills

- **genome-assembly** - De novo assembly, polishing strategy, and assembly QC for long reads
- **variant-calling** - Short-read and cross-platform variant calling, VCF manipulation
- **alignment-files** - Sort, index, and inspect the BAMs these skills produce
- **methylation-analysis** - DMR statistics downstream of bedMethyl
- **alternative-splicing** - Differential isoform usage downstream of isoform discovery
- **phasing-imputation** - Statistical/reference-panel phasing (vs read-backed here)
