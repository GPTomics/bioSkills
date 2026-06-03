# hi-c-analysis

## Overview

Analyze Hi-C and chromosome-conformation-capture data: read-pair processing and library QC, cooler matrices, ICE/KR balancing and distance-decay normalization, A/B compartments, TADs, loops, visualization, differential comparison, and protein-directed 3C (HiChIP/PLAC-seq/Capture Hi-C).

**Tool type:** mixed | **Primary tools:** cooler, cooltools, pairtools, HiCExplorer, chromosight, FitHiChIP

## Skills

| Skill | Description |
|-------|-------------|
| hic-data-io | Load, convert, and manipulate Hi-C matrices in cooler format (.cool/.mcool/.hic) |
| contact-pairs | Process Hi-C read pairs with pairtools; decide library quality from cis/trans and orientation |
| matrix-operations | Balance (ICE/KR), compute distance-decay expected/P(s), and choose a resolution for the depth |
| compartment-analysis | Detect and sign-phase A/B compartments; quantify strength with saddle plots |
| tad-detection | Call TAD boundaries from insulation score across a window sweep |
| loop-calling | Detect chromatin loops (cooltools dots, chromosight, Mustache) and validate with APA |
| hic-visualization | Visualize contact matrices with the correct normalization and color scale |
| hic-differential | Compare Hi-C between conditions at the matched scale (bin-pair, compartment, boundary, loop) |
| hichip-plac-loops | Call protein-directed loops from HiChIP/PLAC-seq/Capture Hi-C (FitHiChIP/MAPS/CHiCAGO) |

## Example Prompts

- "Process my Hi-C read pairs and tell me if the library worked"
- "Balance my contact matrix and compute observed/expected"
- "Call A/B compartments and phase them by GC content"
- "Detect TAD boundaries from my Hi-C data"
- "Find chromatin loops and check them with aggregate peak analysis"
- "Plot a contact matrix the right way for showing compartments"
- "Compare Hi-C between treatment and control"
- "Call loops from my H3K27ac HiChIP"

## Requirements

```bash
pip install cooler cooltools bioframe pairtools matplotlib
conda install -c bioconda hicexplorer chromosight pairix bwa-mem2
# Differential (R/Bioconductor): multiHiCcompare, diffHic, HiCRep; dcHiC (standalone)
# Protein-directed 3C: FitHiChIP, MAPS, hichipper (HiChIP/PLAC-seq); CHiCAGO (Capture Hi-C)
```

## Related Skills

- **read-alignment** - Produces the BAM/alignments feeding pairtools
- **chip-seq** - TF/CTCF peaks to anchor loops and boundaries; peak calling for HiChIP anchors
- **atac-seq** - Enhancer-gene linking and co-accessibility complement loop calls
- **genome-intervals** - Coordinate operations and overlap significance for boundaries/anchors
- **genome-assembly** - Hi-C scaffolding of draft genomes (3D-DNA/SALSA2/YaHS)
- **copy-number** - CNV context that confounds balancing and differential contacts
- **data-visualization** - Genome-track rendering around Hi-C features
