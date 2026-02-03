# PCR Primer Design - Usage Guide

## Overview

This skill covers designing PCR primers using primer3-py, the Python binding for the Primer3 primer design software. Use this for basic PCR, cloning, sequencing primers, and amplicon design.

## Prerequisites

```bash
pip install primer3-py biopython
```

## Quick Start

Ask your AI agent:

- "Design primers for this sequence"
- "Find primers that amplify a 200bp region around position 500"
- "Design primers with Tm around 60C"

## Example Prompts

### Basic Design
> "Design PCR primers for my gene sequence"

> "Find primers for a 150-250bp amplicon"

### Target Specific Regions
> "Design primers flanking positions 100-200"

> "Find primers that amplify exon 3"

### With Constraints
> "Design primers with Tm 58-62C and 45-55% GC"

> "Find primers avoiding the SNP at position 150"

> "Design primers with no more than 3 consecutive identical bases"

### For Cloning
> "Design primers for Gibson assembly with 20bp overlaps"

> "Find primers for restriction cloning with EcoRI sites"

## What the Agent Will Do

1. Load your template sequence
2. Configure primer3 parameters based on your requirements
3. Run primer design
4. Return ranked primer pairs with Tm, GC%, and product size
5. Optionally check for secondary structures

## Tips

- **Tm matching** - Aim for <2C difference between forward and reverse primers
- **GC content** - 40-60% is ideal; avoid GC-rich 3' ends
- **Product size** - For standard PCR, 100-500bp works well
- **3' end** - Should have 1-2 G/C bases (GC clamp) but avoid >3
- **Poly-X runs** - Avoid >4 consecutive identical bases
- **Secondary structure** - Check hairpin and dimer Tm are well below annealing temp
