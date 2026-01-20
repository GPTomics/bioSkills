# Assembly Polishing - Usage Guide

## Overview

Polishing improves assembly accuracy by using additional sequencing data to correct errors. Essential for long-read assemblies which have higher raw error rates.

## Prerequisites

```bash
conda install -c bioconda pilon medaka racon bwa samtools
```

## Quick Start

Tell your AI agent what you want to do:
- "Polish my ONT assembly with medaka"
- "Use Illumina reads to polish my long-read assembly with Pilon"
- "Run multiple rounds of polishing to improve accuracy"

## Example Prompts

### ONT Polishing
> "Polish my Flye assembly with medaka using the original ONT reads"
> "Run two rounds of Racon followed by medaka on my nanopore assembly"

### Illumina Polishing
> "Polish my assembly with Pilon using paired-end Illumina reads"
> "Run two iterations of Pilon polishing to maximize accuracy"

### Combined Approach
> "Polish my ONT assembly first with medaka then with Illumina using Pilon"

## What the Agent Will Do

1. Assess input assembly and available polishing reads
2. Determine optimal polishing strategy based on data types
3. Align reads to assembly (BWA for Illumina, minimap2 for long reads)
4. Run polishing tool with appropriate parameters
5. Iterate polishing if multiple rounds requested
6. Compare before/after statistics to verify improvement

## Tips

- For ONT assemblies: Racon x2 -> medaka -> Pilon x2 is a recommended workflow
- PacBio HiFi assemblies typically do not need polishing due to high accuracy
- Stop polishing when changes between rounds become minimal
- Monitor BUSCO completeness to track improvement
- Pilon requires significant memory; use `--fix bases` to limit corrections
- medaka models should match your flowcell/basecaller version
