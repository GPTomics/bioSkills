# genome-assembly

## Overview

Assemble genomes and transcriptomes from sequencing reads using short-read, long-read, and hybrid approaches. Includes assembly quality assessment, polishing, scaffolding, and contamination detection.

**Tool type:** cli | **Primary tools:** SPAdes, Flye, hifiasm, QUAST, BUSCO

## Skills

| Skill | Description |
|-------|-------------|
| short-read-assembly | De novo assembly from Illumina reads with SPAdes |
| long-read-assembly | Long-read assembly with Flye and Canu |
| hifi-assembly | High-quality assembly from PacBio HiFi with hifiasm |
| metagenome-assembly | Metagenome assembly with metaFlye and metaSPAdes |
| assembly-polishing | Polish assemblies with Pilon, Racon, and medaka |
| assembly-qc | Assess assembly quality with QUAST and BUSCO |
| scaffolding | Hi-C and optical map scaffolding with YaHS |
| contamination-detection | Detect contamination with CheckM2 and GUNC |

## Example Prompts

- "Assemble my bacterial genome from Illumina reads"
- "Run SPAdes on my paired-end data"
- "Assemble my Nanopore reads with Flye"
- "Assemble HiFi reads with hifiasm"
- "Create a phased assembly with Hi-C data"
- "Assemble a metagenome with metaFlye"
- "Polish my assembly with Pilon"
- "Run QUAST to assess my assembly"
- "Check completeness with BUSCO"
- "Scaffold my assembly with Hi-C"
- "Check for contamination with CheckM2"

## Requirements

```bash
# Assemblers
conda install -c bioconda spades flye canu hifiasm

# Polishing
conda install -c bioconda pilon racon medaka

# QC
conda install -c bioconda quast busco checkm2

# Scaffolding
conda install -c bioconda yahs
```

## Related Skills

- **long-read-sequencing** - Long-read alignment and polishing
- **read-qc** - Preprocess reads before assembly
- **sequence-io** - Work with assembled FASTA files
