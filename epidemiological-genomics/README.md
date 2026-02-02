# epidemiological-genomics

## Overview

Pathogen surveillance and outbreak genomics including strain typing, time-scaled phylogenies, transmission inference, AMR tracking, and variant surveillance.

**Tool type:** mixed | **Primary tools:** mlst, TreeTime, TransPhylo, AMRFinderPlus, nextclade

## Skills

| Skill | Description |
|-------|-------------|
| pathogen-typing | MLST, cgMLST, and SNP-based strain typing |
| phylodynamics | Time-scaled phylogenies with TreeTime, BEAST2 |
| transmission-inference | Transmission networks with TransPhylo |
| amr-surveillance | AMR surveillance with epidemiological context |
| variant-surveillance | Lineage assignment with Nextclade, pangolin |

## Example Prompts

- "Type my Salmonella isolates with MLST"
- "Build a time-scaled tree for this outbreak"
- "Identify transmission clusters in my sequences"
- "Who is the likely source of this outbreak?"
- "Track AMR trends across my surveillance samples"
- "Assign lineages to my new SARS-CoV-2 sequences"

## Requirements

```bash
pip install treetime
conda install -c bioconda mlst amrfinderplus nextclade
# R packages for TransPhylo
install.packages('TransPhylo')
```

## Related Skills

- **phylogenetics** - Tree building methods
- **metagenomics** - Community analysis, AMR detection
- **variant-calling** - Variant identification
