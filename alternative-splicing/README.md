# alternative-splicing

## Overview
Alternative splicing analysis for short-read and long-read RNA-seq, covering event quantification, differential splicing, isoform switching with functional consequences, visualization, sequencing QC, single-cell splicing, deep-learning splice variant prediction, outlier detection in rare disease, and full-isoform analysis from PacBio/ONT long reads.

**Tool type:** mixed | **Primary tools:** rMATS-turbo, leafcutter, MAJIQ V3, SUPPA2, IsoformSwitchAnalyzeR, SpliceAI, FRASER 2.0, FLAIR

## Skills
| Skill | Description |
|-------|-------------|
| splicing-quantification | Quantify PSI for SE/A5SS/A3SS/MXE/RI events; microexons, exitrons, IR vs detained introns |
| differential-splicing | Detect splicing changes between conditions (rMATS-turbo, leafcutter, MAJIQ V3, SUPPA2, Shiba) |
| isoform-switching | DTU framework with NMD/ORF/domain/IDR consequences; IsoformSwitchAnalyzeR v2 + DRIMSeq/DEXSeq/satuRn/stageR |
| sashimi-plots | Visualize splicing with ggsashimi, MAJIQ-VOILA, leafviz, Jutils, pyGenomeTracks |
| splicing-qc | RNA-seq QC for splicing: STAR 2-pass, library prep, depth, MaxEntScan/SpliceAI, novel-junction rate |
| single-cell-splicing | Chemistry-first decision (10X 3' is limited); MARVEL, BRIE2, scQuint, SpliZ, Sierra, Psix |
| splice-variant-prediction | SpliceAI/Pangolin/MMSplice/SpliceVault for variant impact; ClinGen SVI 2023 framework |
| outlier-splicing-detection | FRASER 2.0, OUTRIDER, LeafcutterMD, DROP for rare-disease aberrant splicing detection |
| long-read-splicing | Full-isoform splicing from PacBio Iso-Seq / ONT (FLAIR, IsoQuant, Bambu, SQANTI3, rMATS-long) |

## Example Prompts
- "Quantify exon skipping events from my RNA-seq BAMs with rMATS-turbo"
- "Find differential splicing between tumor and normal - run rMATS and leafcutter, report concordant hits"
- "Identify isoform switches with NMD or domain consequences"
- "Create publication-quality sashimi plots for the top 25 differential events"
- "Audit my RNA-seq design for splicing analysis suitability"
- "Will my 10X 3' single-cell data support splicing analysis?"
- "Predict SpliceAI delta scores for clinical variants and apply ClinGen SVI 2023 thresholds"
- "Run FRASER 2.0 on my rare-disease patient vs control panel"
- "Discover full-length isoforms from PacBio Iso-Seq with IsoQuant and SQANTI3"

## Requirements
```bash
# Python tools
pip install rmats-turbo suppa spliceai pangolin mmsplice brie2 \
 rseqc maxentpy pysam pandas numpy scanpy anndata flair-brookslab isoquant

# R / Bioconductor
BiocManager::install(c(
 'IsoformSwitchAnalyzeR', 'leafcutter', 'FRASER', 'OUTRIDER',
 'DRIMSeq', 'DEXSeq', 'satuRn', 'stageR', 'fishpond', 'tximeta',
 'bambu', 'Sierra', 'MARVEL'
))

# CLI tools
conda install -c bioconda regtools ggsashimi pygenometracks rmats2sashimiplot \
 minimap2 samtools snakemake star

# MAJIQ V3 from majiq.biociphers.org (academic license)
# SQANTI3 from github.com/ConesaLab/SQANTI3
# rMATS-long from github.com/Xinglab/rMATS-long
# DROP via bioconda: mamba create -n drop_env -c conda-forge -c bioconda drop --override-channels
```

## Related Skills
- **read-alignment** - STAR 2-pass alignment is required upstream of splicing analysis
- **rna-quantification** - Salmon/kallisto for SUPPA2 and DTU pipelines
- **differential-expression** - Compare gene-level vs transcript-level vs splicing changes
- **variant-calling** - Clinical variant interpretation with SpliceAI integration
- **long-read-sequencing** - PacBio Iso-Seq pipeline upstream of long-read-splicing analysis
- **single-cell** - Cell typing prerequisite for single-cell splicing
- **clinical-databases** - ClinVar, gnomAD splice constraint annotations
- **read-qc** - Sequencing quality assessment
- **workflows** - splicing-pipeline workflow orchestration
