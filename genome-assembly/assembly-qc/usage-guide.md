# Assembly QC - Usage Guide

## Overview

Assembly QC evaluates the quality of genome assemblies using contiguity metrics (QUAST) and gene completeness (BUSCO).

## When to Use This Skill

- After any genome assembly
- Comparing different assemblers
- Evaluating polishing improvement
- Before downstream annotation
- Publishing assembly statistics

## Installation

```bash
conda install -c bioconda quast busco
```

## Quick Start

```bash
# Basic contiguity stats
quast.py assembly.fasta -o quast_output

# Gene completeness
busco -i assembly.fasta -m genome -l bacteria_odb10 -o busco_output
```

## Key Metrics

| Metric | Meaning | Target |
|--------|---------|--------|
| N50 | Contig length at 50% of assembly | Higher is better |
| L50 | Number of contigs in N50 | Lower is better |
| BUSCO Complete | Gene completeness | >95% |
| Misassemblies | Structural errors | 0 |

## BUSCO Lineage Selection

| Organism | Lineage |
|----------|---------|
| Bacteria | bacteria_odb10 |
| Fungi | fungi_odb10 |
| Animals | metazoa_odb10 |
| Plants | viridiplantae_odb10 |

Use `busco --list-datasets` for full list.

## Resources

- [QUAST Manual](http://quast.sourceforge.net/docs/manual.html)
- [BUSCO User Guide](https://busco.ezlab.org/busco_userguide.html)
