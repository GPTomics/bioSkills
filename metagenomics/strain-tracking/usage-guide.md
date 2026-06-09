# Strain Tracking - Usage Guide

## Overview
This skill resolves and compares bacterial strains below the species level from shotgun metagenomes. A strain is not a thing you find - it is a claim your threshold makes: inStrain's popANI >= 99.999% over >= 50% of the genome, or a per-species nGD cutoff, IS the strain definition. ANI tools (MASH/skani/fastANI) answer "are these two genomes the same?" - they saturate and cannot resolve the difference that matters for transmission. Use microdiversity-aware tools (inStrain, StrainPhlAn, MIDAS2, StrainGE) for in-situ strain resolution, and map to your own dRep MAGs, not database genomes.

## Prerequisites
```bash
conda install -c bioconda instrain drep skani bowtie2 samtools
pip install instrain
# StrainPhlAn ships with MetaPhlAn (conda install metaphlan); StrainGE / MIDAS2 separately.
```

Conceptual prerequisites:
- Three distinct tasks: identification, tracking/sharing, deconvolution. Pick the right one.
- Map to dRep-dereplicated MAGs from your own dataset; a distant reference fakes SNVs and corrupts popANI.
- You can only share a strain both samples detect at depth (>= 5x, >= 50% breadth). Absence is not absence.
- Sharing is an undirected edge; direction needs timepoints or contact data, not the genomic comparison.

## Quick Start
Tell your AI agent what you want to do:
- "Detect whether a strain is shared between my paired samples with inStrain popANI"
- "Run a cross-sample transmission survey with StrainPhlAn and derive per-species nGD thresholds"
- "Track a low-abundance pathogen strain down to 0.5x coverage"
- "Compare two assembled MAGs with skani"

## Example Prompts

### Shared-strain detection
> "I have mother and infant gut metagenomes. dRep my MAGs, map reads back, run inStrain profile and compare, and report which genomes are shared at popANI >= 99.999% over at least 50% breadth - and the co-detection rate."

### Transmission survey
> "Run StrainPhlAn across my 200 samples for a target species, build the tree, and derive an nGD threshold that separates same-individual timepoints from unrelated pairs."

### Low-abundance tracking
> "My pathogen is below 1% relative abundance. Use StrainGE to detect and compare its strain across samples."

### Genome comparison (not strains)
> "Compare these two MAGs with skani and tell me if they are the same species - but do not call them the same strain from ANI."

## What the Agent Will Do
1. Identify the task (identification vs tracking vs deconvolution) and pick the matching tool.
2. dRep MAGs from the dataset and map reads to them, not to database genomes.
3. Run inStrain profile and compare on popANI over the co-covered genome fraction, or StrainPhlAn for marker-based surveys.
4. Report the threshold used, the percent genome compared, and co-detection rates alongside sharing.
5. Refuse to infer transmission direction from a single cross-sectional comparison.
6. Use ANI tools only for isolate/MAG comparison and dereplication, never for strain calls.

## Tips
- popANI (microdiversity-aware) detects shared strains conANI/consensus tools miss.
- A sample missing from a StrainPhlAn tree usually means low coverage dropped its markers - not no sharing.
- Prefer skani over fastANI for fragmented MAGs.
- StrainGE pushes detection to ~0.5x; inStrain needs ~5x and 50% breadth.
- SNV tools resolve and compare; they do not separate co-occurring strains - that needs DESMAN or long reads.

## inStrain Key Metrics
- **popANI**: population ANI; a position differs only if the two samples share no alleles (incl. minor) - the shared-strain metric.
- **conANI**: consensus ANI; confounded by microdiversity (conservative sibling).
- **percent_compared** (genome-level `genomeWide_compare.tsv`; the per-scaffold `comparisonsTable.tsv` names it `percent_genome_compared`): genome fraction covered at min_cov in both samples; gate at >= 50%.
- **nucleotide diversity (pi)**: within-population microdiversity.

## Resources
- [inStrain docs](https://instrain.readthedocs.io/)
- [StrainPhlAn (MetaPhlAn) docs](https://github.com/biobakery/MetaPhlAn)
- [skani](https://github.com/bluenote-1577/skani)
- [StrainGE](https://github.com/broadinstitute/StrainGE)

## Related Skills

- metaphlan-profiling - StrainPhlAn builds on MetaPhlAn markers
- kraken-classification - Species presence before strain resolution
- genome-assembly/metagenome-assembly - dRep MAGs to map against; long-read deconvolution
- epidemiological-genomics/amr-surveillance - Isolate outbreak SNP/cgMLST trees
- workflows/metagenomics-pipeline - End-to-end shotgun analysis
