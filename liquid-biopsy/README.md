# liquid-biopsy

## Overview
Cell-free DNA and circulating tumor DNA analysis for non-invasive cancer detection, tumor fraction estimation, mutation detection, and treatment monitoring from plasma samples.

**Tool type:** mixed | **Primary tools:** ichorCNA, FinaleToolkit, fgbio, VarDict, MethylDackel

## Skills
| Skill | Description |
|-------|-------------|
| cfdna-preprocessing | Preprocess cfDNA reads with UMI/duplex consensus error suppression |
| analytical-validation | Set and report limits of detection from molecule-counting and error statistics |
| ctdna-mutation-detection | Detect somatic mutations at low VAF with CHIP subtraction |
| fragment-analysis | Analyze fragmentomics patterns for cancer detection |
| tumor-fraction-estimation | Estimate ctDNA fraction from shallow WGS |
| methylation-based-detection | Analyze cfDNA methylation for detection and tissue-of-origin |
| longitudinal-monitoring | Track ctDNA dynamics over treatment for MRD and response |

## Example Prompts
- "Preprocess my plasma cfDNA FASTQ files with UMI consensus calling"
- "How many genome-equivalents do I need to detect a 0.1% variant?"
- "Detect mutations at 0.5% VAF from my targeted panel and filter CHIP"
- "Analyze fragment size distribution for tumor signal"
- "Estimate tumor fraction from my shallow WGS data"
- "Infer tissue-of-origin from cfDNA methylation"
- "Track ctDNA levels across my serial samples for molecular relapse"

## Requirements
```bash
# Python
pip install finaletoolkit pysam pandas numpy scipy statsmodels matplotlib

# R
# ichorCNA is on GitHub, not Bioconductor
devtools::install_github('GavinHaLab/ichorCNA')

# CLI
conda install -c bioconda fgbio vardict-java methyldackel
```

## Related Skills

- **variant-calling** - Somatic variant calling principles
- **copy-number** - CNV detection concepts
- **methylation-analysis** - Methylation processing
