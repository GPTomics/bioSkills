# rna-structure

## Overview
Predict and analyze RNA secondary structures, search for non-coding RNA families, and interpret experimental structure probing data.

**Tool type:** cli | **Primary tools:** ViennaRNA, Infernal, ShapeMapper2

## Skills
| Skill | Description |
|-------|-------------|
| secondary-structure-prediction | Predict RNA secondary structures with ViennaRNA (MFE, partition function, consensus) |
| ncrna-search | Search for ncRNA homologs and classify RNA families with Infernal/Rfam |
| structure-probing | Analyze SHAPE-MaP and DMS-MaPseq experimental structure probing data |

## Example Prompts
- "Predict the secondary structure and MFE for this RNA sequence"
- "Search my transcript against Rfam to identify ncRNA families"
- "Process my SHAPE-MaP data and use reactivities to constrain folding"
- "Build a consensus structure from an alignment of homologous RNAs"
- "Build a custom covariance model for a novel RNA family"

## Requirements
```bash
# ViennaRNA (includes RNAfold, RNAalifold, RNAcofold, Python API)
conda install -c bioconda viennarna

# Infernal + Rfam database
conda install -c bioconda infernal
# Download Rfam CMs (~500 MB)
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz && cmpress Rfam.cm

# ShapeMapper2
conda install -c bioconda shapemapper2

# Python dependencies
pip install biopython pandas matplotlib
```

## Related Skills
- **sequence-manipulation** - Sequence property calculations and reverse complement
- **genome-annotation** - Genome-wide ncRNA annotation pipelines
- **small-rna-seq** - Small RNA sequencing analysis (miRNA, piRNA)
- **epitranscriptomics** - RNA modification detection (m6A, pseudouridine)
