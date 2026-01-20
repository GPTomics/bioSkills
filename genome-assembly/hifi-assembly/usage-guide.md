# HiFi Assembly Usage Guide

## Overview

hifiasm is the leading assembler for PacBio HiFi reads, producing highly contiguous and accurate genome assemblies with built-in phasing support.

## When to Use

| Scenario | Recommended Approach |
|----------|---------------------|
| Diploid, no phasing data | Basic hifiasm |
| Diploid with Hi-C | hifiasm --h1/--h2 |
| Trio sample available | hifiasm with yak |
| Highly repetitive genome | Add ONT ultra-long |
| Polyploid | Adjust --n-hap |

## Quick Start Prompts

- "Assemble this HiFi dataset with hifiasm"
- "Create a phased assembly using Hi-C data"
- "Run trio-binned assembly with parental reads"
- "Assess my hifiasm assembly quality"

## Requirements

```bash
# Install hifiasm
conda install -c bioconda hifiasm

# For trio binning
conda install -c bioconda yak

# For quality assessment
conda install -c bioconda quast busco seqkit
```

## Input Recommendations

| Parameter | Recommendation |
|-----------|---------------|
| HiFi coverage | 30-60x for mammals |
| Read quality | Q20+ (typically Q30+) |
| Read N50 | >15 kb preferred |
| Hi-C coverage | 30-50x for phasing |

## Output Files

| File | Description |
|------|-------------|
| *.bp.p_ctg.gfa | Primary contigs (consensus) |
| *.bp.a_ctg.gfa | Alternate contigs (haplotigs) |
| *.bp.hap1.p_ctg.gfa | Haplotype 1 (phased) |
| *.bp.hap2.p_ctg.gfa | Haplotype 2 (phased) |
| *.ovlp.* | Overlap information |
| *.ec.* | Error-corrected reads |

## Related Skills

- **genome-assembly/quality-assessment** - Evaluate assemblies
- **genome-assembly/scaffolding** - Chromosome-scale scaffolding
- **long-read-sequencing/read-qc** - Input QC
