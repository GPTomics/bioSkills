# Contact Pairs Processing - Usage Guide

## Overview

This skill covers processing Hi-C read pairs using pairtools, including parsing, sorting, deduplication, and filtering to generate valid contact pairs.

## Prerequisites

```bash
conda install -c bioconda pairtools cooler
```

## Quick Start

Tell your AI agent what you want to do:
- "Process my Hi-C BAM file with pairtools"
- "Filter valid pairs from my pairs file"

## Example Prompts

### Basic Processing
> "Parse Hi-C pairs from this BAM"

> "Remove duplicates from my pairs file"

### Filtering
> "Select only UU pair types"

> "Remove self-ligation artifacts"

### Statistics
> "Get statistics from my pairs file"

> "How many valid pairs do I have?"

## What the Agent Will Do

1. Parse alignments to pairs format
2. Sort pairs by position
3. Remove PCR duplicates
4. Filter by pair type
5. Generate statistics

## Tips

- **Pair types** - UU (unique-unique) are highest confidence
- **MAPQ** - Filter for MAPQ >= 30 for quality
- **Self-ligation** - Remove pairs < 1kb apart
- **Pipeline** - Use pipes to chain operations
