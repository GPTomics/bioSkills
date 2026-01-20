# Library Design - Usage Guide

## Overview
CRISPR library design involves selecting optimal sgRNAs for genetic screens, balancing efficiency, specificity, and practical cloning considerations.

## Prerequisites
```bash
pip install biopython pandas numpy
# For sgRNA scoring
pip install crispor  # or use web interface
```

## Quick Start
Tell your AI agent what you want to do:
- "Design sgRNAs targeting my gene list"
- "Create a focused CRISPR library for pathway genes"
- "Select optimal guides with low off-target potential"

## Example Prompts

### sgRNA Selection
> "Design 4 sgRNAs per gene for my list of 50 kinases. Filter for GC content 40-60% and no poly-T runs."

> "Find the best guides targeting TP53 exons. Prioritize high efficiency scores and low off-targets."

> "Select sgRNAs for a CRISPRa activation library targeting my gene set."

### Library Design
> "Design a focused knockout library for 200 DNA repair genes with 6 guides per gene."

> "Create a genome-wide library design. Include essential gene controls and 500 non-targeting controls."

> "Design a sublibrary from the Brunello library targeting only my genes of interest."

### Validation and QC
> "Check my designed guides for off-target sites in the human genome."

> "Validate that my library has adequate controls: essentials, non-essentials, and non-targeting."

> "Calculate the optimal cell numbers needed for 500x library coverage."

## What the Agent Will Do
1. Parse gene list and retrieve sequences
2. Identify PAM sites and candidate guides
3. Score guides for efficiency and specificity
4. Filter based on design criteria
5. Add control guides (essential, non-targeting)
6. Generate oligo sequences for synthesis

## Tips
- Include 4-6 guides per gene for dropout screens, 6-10 for focused screens
- Optimal GC content is 40-60%, acceptable is 30-70%
- Avoid poly-T runs (>3 Ts) which terminate Pol III transcription
- 5' G preferred for U6 promoter transcription
- Include 5-10% of library as controls

## Design Criteria

| Criterion | Optimal | Acceptable |
|-----------|---------|------------|
| GC content | 40-60% | 30-70% |
| Guide length | 20 nt | 18-22 nt |
| Poly-T runs | None | < 4 Ts |
| 5' nucleotide | G | A, C |
| Off-targets | 0 in exons | < 3 in exons |

## Library Size Guidelines

| Screen Type | Guides/Gene | Cell Coverage |
|-------------|-------------|---------------|
| Dropout | 4-6 | 500x |
| Activation | 3-4 | 300x |
| Focused | 6-10 | 1000x |

## Control Design
- **Essential genes**: 20-50 guides, validates dropout
- **Non-targeting**: 100-500 guides, background estimation
- **Safe harbor**: AAVS1/ROSA26, validates Cas9 activity

## References
- CRISPOR: doi:10.1093/nar/gkw441
- Azimuth 2.0: doi:10.1038/nbt.3437
- Library design review: doi:10.1038/nrg.2017.97
