# Long-Read Assembly - Usage Guide

## Overview

Long-read assemblers produce highly contiguous assemblies from ONT or PacBio reads, often achieving complete bacterial chromosomes in a single contig.

## When to Use This Skill

- Bacterial genome closure
- Resolving repetitive regions
- Structural variant detection
- Complete chromosome assembly
- When short-read assembly is too fragmented

## Installation

```bash
conda install -c bioconda flye canu hifiasm
```

## Quick Start

```bash
# ONT with Flye (recommended)
flye --nano-raw reads.fq.gz --out-dir output --genome-size 5m

# PacBio HiFi with Hifiasm
hifiasm -o asm -t 16 reads.fq.gz
awk '/^S/{print ">"$2;print $3}' asm.bp.p_ctg.gfa > assembly.fasta
```

## Tool Selection

| Data Type | Recommended |
|-----------|-------------|
| ONT R9/R10 | Flye |
| PacBio CLR | Flye or Canu |
| PacBio HiFi | Hifiasm |
| Metagenome | Flye --meta |

## Coverage Requirements

| Goal | Minimum | Recommended |
|------|---------|-------------|
| Draft | 20x | 30x |
| Complete | 50x | 100x |

## Post-Assembly

Long-read assemblies typically need polishing:
1. medaka (ONT)
2. Pilon (short reads)

See assembly-polishing skill.

## Resources

- [Flye](https://github.com/fenderglass/Flye)
- [Canu](https://github.com/marbl/canu)
- [Hifiasm](https://github.com/chhylp123/hifiasm)
