# Annotation Transfer - Usage Guide

## Overview

Transfer gene annotations between genome assemblies using Liftoff for same-species annotation liftover and MiniProt for cross-species protein-to-genome alignment. Enables rapid annotation of new assemblies when a high-quality reference annotation already exists.

## Prerequisites

```bash
# Liftoff (same-species liftover)
pip install liftoff

# LiftOn (improved successor)
pip install lifton

# MiniProt (cross-species protein alignment)
conda install -c bioconda miniprot

# Python utilities
pip install gffutils biopython pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Transfer annotations from the reference genome to my new assembly"
- "Map gene models from a related species to my genome"
- "Lift over the GRCh38 annotation to my custom human assembly"

## Example Prompts

### Same-Species Transfer

> "Use Liftoff to transfer the reference annotation to my new assembly of the same species"

> "Lift over the GENCODE annotation from GRCh38 to my T2T-based assembly"

### Cross-Species Transfer

> "Align mouse proteins to my rat genome with MiniProt to transfer annotations"

> "Map proteins from a closely related species to my new genome assembly"

### Quality Assessment

> "Compare the transferred annotation to the original - how many genes transferred successfully?"

> "Check how many transferred gene models have valid open reading frames"

## What the Agent Will Do

1. Determine appropriate tool (Liftoff for same-species, MiniProt for cross-species)
2. Set up chromosome name mapping if needed
3. Run annotation transfer with appropriate parameters
4. Report transfer rate and unmapped features
5. Validate transferred gene models (ORF integrity check)
6. Suggest de novo prediction for regions that failed transfer

## Tips

- **Same species** - Liftoff with strict parameters (sc 0.95, s 0.90) gives best results
- **Close species** - LiftOn combines Liftoff + MiniProt for improved accuracy
- **Distant species** - MiniProt protein alignment with appropriate max intron size (-G)
- **Intron size** - Set -G appropriately: vertebrates 500000, insects 50000, fungi 5000
- **Assembly quality** - Transfer rates directly correlate with target assembly contiguity
- **Chromosome mapping** - Provide explicit mapping if reference and target use different naming
- **Complement with de novo** - Transfer catches known genes; de novo prediction finds novel ones

## Related Skills

- genome-annotation/eukaryotic-gene-prediction - De novo prediction alternative
- comparative-genomics/ortholog-inference - Orthology-based functional transfer
- comparative-genomics/synteny-analysis - Synteny context for annotation transfer
- genome-intervals/gtf-gff-handling - Parse transferred annotation files
