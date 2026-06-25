# rna-structure

## Overview
Predict and analyze RNA secondary structures, search for non-coding RNA families, interpret experimental structure probing data, and validate conserved structures by evolutionary covariation.

**Tool type:** cli | **Primary tools:** ViennaRNA, ShapeMapper2, Infernal, R-scape

## Skills
| Skill | Description |
|-------|-------------|
| secondary-structure-prediction | Predict RNA structure with ViennaRNA, reporting the Boltzmann ensemble (MFE, partition function, centroid/MEA, confidence), with consensus, long-RNA, and pseudoknot routes |
| structure-probing | Process SHAPE-MaP / DMS-MaPseq data into reactivities and use them as soft folding restraints; handle MaP vs RT-stop, normalization, in-cell vs in-vitro, and multiple conformations |
| ncrna-search | Search for ncRNA homologs and classify families with Infernal/Rfam covariance models; choose CM vs BLAST, gathering thresholds, and clan resolution |
| covariation-analysis | Test whether a conserved RNA structure is supported by evolutionary covariation with R-scape, returning a power-aware verdict |

## Example Prompts
- "Fold my RNA and report the ensemble, not just the MFE, with per-base confidence"
- "My RNA is 8 kb or has a pseudoknot; fold it sensibly"
- "Process my SHAPE-MaP data and use the reactivities to constrain folding"
- "Is this protected region base-paired or protein-bound?"
- "Search my transcripts against Rfam to classify ncRNA families"
- "Should I use a covariance model or BLAST for this RNA?"
- "Is my conserved RNA structure actually supported by covariation, or does the alignment lack power?"

## Requirements
```bash
# ViennaRNA (RNAfold, RNAalifold, RNAcofold, RNAplfold, Python API)
conda install -c bioconda viennarna

# Infernal + Rfam database (pre-calibrated: press, do not recalibrate)
conda install -c bioconda infernal
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz && cmpress Rfam.cm

# ShapeMapper2 (Linux only; Docker on macOS), SEISMIC-RNA, R-scape
conda install -c bioconda shapemapper2 rscape
pip install seismic-rna

# Python dependencies
pip install biopython pandas matplotlib numpy
```

## Related Skills
- **genome-annotation** - Genome-wide ncRNA annotation pipelines
- **small-rna-seq** - Small RNA sequencing analysis (miRNA, piRNA) and target prediction
- **epitranscriptomics** - RNA modification detection (m6A, pseudouridine) that confounds probing
- **clip-seq** - Protein-RNA binding site detection (in-cell protection as an RBP footprint)
- **sequence-manipulation** - Sequence property calculations and reverse complement
