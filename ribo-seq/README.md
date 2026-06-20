# ribo-seq

## Overview

Analyze ribosome profiling (Ribo-seq) data to study translation at codon resolution: preprocessing with UMI handling, periodicity QC and P-site calibration, ORF detection, translation efficiency, ribosome stalling, and initiation-site mapping.

**Tool type:** mixed | **Primary tools:** riboWaltz, RiboCode, ORFquant, riborex, Ribo-TISH, STAR

## Skills

| Skill | Description |
|-------|-------------|
| riboseq-preprocessing | UMI handling, trimming, rRNA depletion, and footprint-aware alignment |
| ribosome-periodicity | Validate 3-nt periodicity and calibrate read-length-specific P-site offsets |
| orf-detection | Detect and quantify translated ORFs (uORFs, novel ORFs) with RiboCode and ORFquant |
| translation-efficiency | Differential TE separating translational control from buffering |
| ribosome-stalling | Detect pausing with local-relative metrics, judging drug artifacts |
| initiation-site-mapping | Map start codons, including non-AUG, from harringtonine/LTM/Ribo-RET data |

## Method Selection

| Question | Skill | Key decision |
|----------|-------|--------------|
| Is my library usable at codon level? | ribosome-periodicity | Frame-0 fraction; per-length P-site offsets |
| Should I deduplicate? | riboseq-preprocessing | Only with UMIs |
| What ORFs are translated? | orf-detection | Periodicity (not coverage); near-cognate starts |
| Where does translation initiate? | initiation-site-mapping | Needs an initiation-drug (TI-seq) library |
| Translational control or buffering? | translation-efficiency | anota2seq mode-of-regulation |
| Is this pause real or a drug artifact? | ribosome-stalling | Flash-frozen no-drug data; A-site, local metrics |

## Example Prompts

- "Preprocess my Ribo-seq with UMI extraction and deduplication"
- "Check 3-nt periodicity and calibrate P-site offsets"
- "Find actively translated uORFs at near-cognate start codons"
- "Tell translational control from buffering"
- "Map translation initiation sites from my harringtonine data"
- "Detect ribosome stalling, accounting for cycloheximide artifacts"
- "Quantify ORF-level translation isoform-aware across conditions"
- "Run a differential translation efficiency analysis"

## Requirements

```bash
# CLI
conda install -c bioconda cutadapt umi_tools star bowtie2 sortmerna samtools ribocode ribotish

# R (periodicity, quantification, differential TE)
# BiocManager::install(c('riboWaltz', 'ORFquant', 'ORFik', 'riborex', 'xtail', 'anota2seq', 'DESeq2'))

# Python alternative paths
pip install plastid
```

## Related Skills

- **read-alignment** - Spliced and unspliced read alignment
- **rna-quantification** - Matched RNA-seq quantification for translation efficiency
- **differential-expression** - Count-based differential testing for ORFs and TE
- **workflows/riboseq-pipeline** - End-to-end orchestration of these skills
