# genome-assembly

## Overview

Assemble genomes and transcriptomes from sequencing reads using short-read, long-read, and hybrid approaches. Includes assembly quality assessment and polishing workflows.

**Tool type:** cli | **Primary tools:** SPAdes, Flye, QUAST, BUSCO

## Skills

| Skill | Description |
|-------|-------------|
| short-read-assembly | De novo assembly from Illumina reads with SPAdes |
| long-read-assembly | Long-read assembly with Flye and Canu |
| assembly-polishing | Polish assemblies with Pilon, Racon, and medaka |
| assembly-qc | Assess assembly quality with QUAST and BUSCO |

## Example Prompts

- "Assemble my bacterial genome from Illumina reads"
- "Run SPAdes on my paired-end data"
- "Assemble my Nanopore reads with Flye"
- "Polish my assembly with Pilon"
- "Run QUAST to assess my assembly"
- "Check completeness with BUSCO"
- "Do hybrid assembly with short and long reads"
- "Assemble a metagenome with metaSPAdes"
- "Polish my ONT assembly with medaka then Pilon"
- "Compare assembly statistics across samples"

## Requirements

```bash
# Assemblers
conda install -c bioconda spades flye canu

# Polishing
conda install -c bioconda pilon racon

# QC
conda install -c bioconda quast busco
```

## Related Skills

- **long-read-sequencing** - Long-read alignment and polishing
- **read-qc** - Preprocess reads before assembly
- **sequence-io** - Work with assembled FASTA files
