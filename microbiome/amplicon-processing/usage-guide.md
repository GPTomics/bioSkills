# Amplicon Processing - Usage Guide

## Overview

This skill turns demultiplexed marker-gene FASTQ (16S rRNA, ITS) into an ASV (amplicon sequence variant) feature table plus representative sequences, using DADA2's per-run error model. An ASV is a model-inferred exact sequence conditioned on one sequencing run - not a clustered OTU and not an organism. The decisions that matter (which are invisible in a recipe) are: remove primers before truncating, learn the error model separately per run, choose truncation lengths within the merge-overlap budget, and never fix-truncate variable-length ITS. The output feeds taxonomy assignment, diversity, and differential abundance; the compositional statistics of the resulting table are shared with shotgun metagenomics.

## Prerequisites

```bash
conda install -c bioconda cutadapt itsxpress
```

```r
install.packages('BiocManager')
BiocManager::install(c('dada2', 'decontam'))
```

The QIIME2 path (`qiime dada2 denoise-paired`, `qiime deblur denoise-16S`) installs as its own conda env (`qiime2-amplicon-<release>`); see qiime2-workflow.

Conceptual prerequisites:
- Input is DEMULTIPLEXED paired-end (or single-end) FASTQ, one file set per sample.
- The forward and reverse PCR primer sequences are known (needed for cutadapt).
- The amplicon region and its approximate length are known (V4 ~253 bp, V3-V4 ~460 bp, ITS variable) - this sets the truncation budget.
- The error model is learned PER sequencing run; a multi-run study processes each run separately before merging tables.
- For low-biomass samples (skin, biopsy, BAL, sterile-site swabs), sequence negative controls (extraction blanks, no-template PCR) and a positive mock community so contaminants can be identified with decontam.

## Quick Start

Tell your AI agent what you want to do:
- "Remove primers from my paired-end 16S reads, then run DADA2 to get an ASV table"
- "I have V3-V4 reads on 2x250 and my merge rate is near zero - help"
- "Process two MiSeq runs correctly with a per-run error model and merge the tables"
- "Process my fungal ITS amplicons into ASVs"
- "These are low-biomass skin swabs with extraction blanks - remove reagent contaminants with decontam"

## Example Prompts

### ASV inference
> "I have demultiplexed paired-end 16S V4 reads. Remove the primers with cutadapt, learn the error model per run, denoise, merge pairs, and remove chimeras to give me an ASV table and a read-tracking summary."

### Truncation budget
> "My amplicon is V3-V4 (~460 bp) on 2x250 reads. What truncLen should I use, and how do I keep enough overlap to merge while controlling quality?"

### Multi-run studies
> "My samples came off three sequencing runs. Show me how to learn the error model on each run separately and combine the sequence tables before chimera removal."

### ITS / non-16S
> "I have fungal ITS2 amplicons of variable length. Trim the spacer and infer ASVs without fixed truncation."

### Platform caveats
> "These are NovaSeq reads with binned quality scores. Check whether the DADA2 error model is fit correctly and fix it if not."

### Low-biomass / decontamination
> "These are low-biomass biopsy samples sequenced with extraction-blank and no-template-PCR negative controls. After building the ASV table, run decontam to flag and remove reagent/kit contaminants, and report how many ASVs and reads were removed."

## What the Agent Will Do

1. Confirm the marker, region, primers, read length, and number of sequencing runs.
2. Remove primers with cutadapt (`--discard-untrimmed`) before any quality step.
3. Inspect quality profiles and compute the truncation budget from amplicon and read length.
4. Filter and truncate on expected errors within that budget (per run).
5. Learn the error model from each run separately and inspect `plotErrors` (enforcing monotonicity on binned-quality data).
6. Denoise, merge pairs, and build a per-run sequence table.
7. Merge run-level tables, then run a single chimera removal.
8. For low-biomass studies, classify and remove reagent/kit contaminant ASVs with decontam using the negative controls (and DNA concentration if measured), reporting what was removed.
9. Produce the ASV table, representative sequences, and a read-tracking table; flag whether ASV count overstates richness.

## Tips

- An ASV is a run-conditioned exact sequence, not an organism; do not equate ASV count with species richness without collapsing to a rank.
- Primers off FIRST. A large read fraction removed as "chimeric" usually means primers were not trimmed.
- truncLen is a detection budget: truncLen_F + truncLen_R must exceed the amplicon length by at least 12 bp, or pairs cannot merge.
- A near-zero merge rate is almost always a budget problem, not bad data.
- For ITS, never set a fixed truncLen - trim the spacer with ITSxpress and use truncLen=0.
- On NovaSeq/NextSeq/iSeq, quality scores are binned to ~4 values; inspect plotErrors and enforce a monotonic error fit.
- Use pool='pseudo' when rare or singleton ASVs matter.
- Low-biomass samples can be dominated by the reagent/kit "kitome": never interpret a near-sterile sample without negative controls, and run decontam (prevalence method with controls, frequency with DNA concentration, combined with both). The shotgun analogue is metagenomics/contamination-controls.
- Save the chimera-free sequence table as an RDS for downstream taxonomy and phyloseq work.

## Related Skills

- taxonomy-assignment - Assign taxonomy to the ASVs produced here
- diversity-analysis - Alpha/beta diversity of the resulting community table
- differential-abundance - Compositional DA on the ASV/feature table
- qiime2-workflow - The QIIME2 CLI equivalent of this R workflow
- read-qc/adapter-trimming - cutadapt primer removal before DADA2
- metagenomics/kraken-classification - Shotgun (not amplicon) read classification
- metagenomics/abundance-estimation - Shared compositional/normalization theory
- phylogenetics/tree-io - Phylogenetic tree for UniFrac / Faith PD
- workflows/microbiome-pipeline - End-to-end amplicon pipeline
