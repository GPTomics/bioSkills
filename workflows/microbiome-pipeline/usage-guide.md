# Microbiome Pipeline - Usage Guide

## Overview

This workflow takes demultiplexed 16S/ITS amplicon FASTQ through to a confidence-graded differential-abundance result, orchestrating cutadapt primer removal, per-run DADA2 ASV inference, region-matched taxonomy, a placed phylogenetic tree, declared-depth alpha/beta diversity, a consensus of compositionally-aware differential-abundance tools, and optional PICRUSt2 functional prediction. It is an orchestration skill: it sequences the stages and defers every per-step scientific decision to the six microbiome category skills. The central discipline is that the pipeline is a chain of modeling choices, not a conveyor belt - each stage's parameters silently set what the next stage can find, so each choice is declared and defended in order. Shotgun (whole-genome) sequencing is a different assay and lives in workflows/metagenomics-pipeline.

## Prerequisites

```bash
# CLI tools: primer removal, optional SEPP/QIIME2 tree, optional PICRUSt2 prediction
conda install -c bioconda cutadapt
```

```r
# R / Bioconductor analysis stack
BiocManager::install(c('dada2', 'phyloseq', 'ALDEx2', 'ANCOMBC', 'DECIPHER', 'decontam'))
install.packages('vegan')
```

Conceptual prerequisites and notes:
- Reads are DEMULTIPLEXED, paired-end (or CCS), and the PRIMER sequences and amplified REGION (V4, V3-V4, ITS) are known - primer removal needs the primers, and the region sets the truncation budget and the classifier.
- The DADA2 error model is fit PER SEQUENCING RUN; multi-run studies process each run separately, then merge sequence tables, then remove chimeras once.
- The reference database (SILVA 138.1, GTDB r220, UNITE) is a large download and a versioned dependency; record the release and the classifier training region.
- A SEPP/Greengenes2 reference package (or a QIIME2 install) is needed for a placed tree; a de novo tree from short reads is a last resort.
- PICRUSt2 and QIIME2 each install as their own conda environment; functional prediction is optional and hypothesis-generating.
- For low-biomass samples, sequence negative/positive controls: the pipeline removes host mitochondria/chloroplast features after taxonomy and reagent contaminants with decontam before diversity/DA.

## Quick Start

Tell your AI agent what you want to do:
- "Run the full 16S amplicon pipeline on my demultiplexed FASTQ files"
- "Remove primers, denoise per run with DADA2, and build me an ASV table"
- "Assign taxonomy, then compute alpha and beta diversity at a sensible depth"
- "Find differentially abundant taxa as a consensus of two compositional tools"
- "Predict functional potential from my 16S data and report NSTI"

## Example Prompts

### End-to-end pipeline
> "I have demultiplexed paired-end 16S V4 reads from two sequencing runs. Remove the primers, learn the error model per run, denoise, merge the run tables, remove chimeras, assign taxonomy against SILVA, and give me an ASV table with a placed tree."

### Diversity
> "Build a phyloseq object with a SEPP tree, pick a rarefaction depth from the alpha-rarefaction plateau, report which samples get dropped, and test group differences with weighted UniFrac, pairing PERMANOVA with betadisper."

### Differential abundance
> "On the unrarefied counts, run ALDEx2 and ANCOM-BC2 and report the consensus differentially abundant taxa with effect sizes - tell me which tool found what."

### Functional prediction
> "Predict KO and MetaCyc potential from my ASVs with PICRUSt2, gate on NSTI, and report the read fraction dropped - frame it as potential, not activity."

## What the Agent Will Do

1. Confirms the reads are demultiplexed and that the primers and region are known, and runs read QC (read-qc/quality-reports).
2. Removes primers with cutadapt BEFORE any truncation, discarding primerless pairs (read-qc/adapter-trimming, microbiome/amplicon-processing).
3. Runs the DADA2 block (filterAndTrim -> learnErrors -> dada -> mergePairs) ONCE PER sequencing run, choosing truncation lengths within the merge-overlap budget (microbiome/amplicon-processing).
4. Merges the per-run sequence tables by exact sequence string, then removes chimeras once on the combined table.
5. Assigns taxonomy with a region-matched classifier, reporting genus rather than species for 16S (microbiome/taxonomy-assignment).
6. Filters host mitochondria/chloroplast features after taxonomy, and for low-biomass studies removes reagent contaminants with decontam using the negative controls (microbiome/taxonomy-assignment, microbiome/amplicon-processing).
7. Assembles the phyloseq object and attaches a SEPP/Greengenes2 placed tree rather than a de novo tree from short reads (microbiome/diversity-analysis, phylogenetics/tree-io).
8. Rarefies ONLY into the diversity branch at a declared depth, reports the dropped samples, computes alpha/beta diversity, and pairs adonis2 with betadisper (microbiome/diversity-analysis).
9. Runs differential abundance on the UNrarefied counts using two or more compositionally-aware tools and reports the intersection as high-confidence (microbiome/differential-abundance).
10. Optionally runs PICRUSt2, reports the NSTI distribution and the dropped read fraction, and frames every claim as predicted potential (microbiome/functional-prediction).
11. Records the per-stage choices (primers, truncation, error model per run, database release, organelle/decontam filtering, tree method, rarefaction depth, DA tool panel, PICRUSt2 reference) alongside the results.

## QC Checkpoints

| Stage | Check | Action if it fails |
|-------|-------|--------------------|
| Primer removal | Most pairs retain a primer match | Confirm the primer sequences and orientation; do not skip to filtering |
| Denoising | Reads survive filter -> denoise -> merge with no cliff | A merge cliff is a truncLen budget problem, not bad data; keep length, loosen reverse maxEE |
| Chimeras | Small READ fraction removed | A large read loss flags leftover primers, not a chimera storm |
| Taxonomy | Genus assigned for most ASVs | Use a region-matched classifier; do not force species calls for 16S |
| Contamination | Organelle filtered; low-biomass decontaminated | filter Mitochondria/Chloroplast after taxonomy; run decontam with negative controls before diversity/DA |
| Diversity | Sampling depth declared, dropped samples listed | Pick the depth from the alpha-rarefaction plateau, not min(sample_sums) |
| Differential abundance | Consensus of >=2 CoDA tools on unrarefied counts | Never rarefy for DA; never report one tool's hit list alone |

## Tips

- Primer removal comes first: leftover primers corrupt the error model and look chimeric. Remove them with cutadapt before filterAndTrim.
- Fit the error model per sequencing run. Concatenating runs before learnErrors denoises wrong; merge run-level sequence tables instead, and carry run as a batch covariate into DA.
- Set truncation lengths to keep the merge-overlap budget (truncLen_F + truncLen_R >= amplicon length + ~12), not just to chase quality.
- Filter host mitochondria/chloroplast features after taxonomy, and decontaminate low-biomass samples with decontam using negative controls, before diversity and DA (microbiome/taxonomy-assignment, microbiome/amplicon-processing).
- Prefer a SEPP/Greengenes2 placed tree over a de novo tree from short reads for any UniFrac or Faith PD metric.
- Rarefy for diversity, never for differential abundance. Keep the raw counts and rarefy only into the diversity branch.
- Report the rarefaction depth AND which samples were dropped below it; the dropped ones are rarely random.
- Pair every PERMANOVA (adonis2) with betadisper so a dispersion difference is not read as a community shift.
- Treat differential abundance as a consensus of two or more compositional tools (ALDEx2, ANCOM-BC2, LinDA), not one tool's volcano plot - the hit list depends more on the tool than the biology.
- PICRUSt2 predicts potential, not activity, and is circular with taxonomy. Report NSTI and frame functional results as hypothesis-generating; for measured function use shotgun (metagenomics/functional-profiling).
- For shotgun (whole-genome) data, this pipeline does not apply - use workflows/metagenomics-pipeline.

## Related Skills

- microbiome/amplicon-processing - Primer removal, per-run error model, truncLen budget, chimeras, ITS
- microbiome/taxonomy-assignment - Region-matched classifier and reference-database choice
- microbiome/diversity-analysis - Sampling depth, tree choice, metric choice, adonis2 + betadisper
- microbiome/differential-abundance - Compositional DA tools and the consensus deliverable
- microbiome/functional-prediction - PICRUSt2 predicted potential gated on NSTI
- microbiome/qiime2-workflow - The QIIME2 artifact/provenance route for the whole chain
- read-qc/adapter-trimming - cutadapt primer removal mechanics before DADA2
- reporting/automated-qc-reports - Aggregate FastQC/MultiQC across samples
- metagenomics/abundance-estimation - Shared compositional/normalization/rarefaction theory
- workflows/metagenomics-pipeline - The shotgun (WGS) equivalent of this amplicon pipeline
