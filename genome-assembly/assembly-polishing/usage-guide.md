# Assembly Polishing - Usage Guide

## Overview

Polishing improves assembly accuracy by using additional sequencing data to correct errors. Essential for long-read assemblies which have higher raw error rates.

## When to Use This Skill

- After long-read assembly (ONT or PacBio CLR)
- When Illumina data is available for validation
- To improve consensus accuracy
- Before gene annotation

## Installation

```bash
conda install -c bioconda pilon medaka racon
```

## Quick Start

### Illumina Polishing (Pilon)

```bash
bwa index assembly.fa
bwa mem -t 16 assembly.fa R1.fq.gz R2.fq.gz | samtools sort -o aligned.bam
samtools index aligned.bam
pilon --genome assembly.fa --frags aligned.bam --output polished
```

### ONT Polishing (medaka)

```bash
medaka_consensus -i reads.fq.gz -d assembly.fa -o medaka_out -t 8
```

## Recommended Workflow

| Assembly Type | Polishing Steps |
|--------------|-----------------|
| ONT | Racon x2 → medaka → Pilon x2 |
| PacBio CLR | Racon x2 → Pilon x2 |
| PacBio HiFi | Usually none needed |

## When to Stop

- No more changes in Pilon output
- Error rate stabilizes
- BUSCO completeness stops improving

## Resources

- [Pilon](https://github.com/broadinstitute/pilon)
- [medaka](https://github.com/nanoporetech/medaka)
- [Racon](https://github.com/isovic/racon)
