# Functional Profiling - Usage Guide

## Overview
HUMAnN 3 profiles the functional potential of a community by quantifying gene families (UniRef90, in RPK) and MetaCyc pathway abundances, using a tiered search that only translates the reads its fast steps could not place. The output is potential, not activity - it tells you what the community can do, never what it is doing. Two reference databases (ChocoPhlAn, UniRef) and the MetaPhlAn version drive every number, and the UNMAPPED/UNINTEGRATED rows are the denominator, not noise.

## Prerequisites
```bash
conda create -n humann -c bioconda humann
conda activate humann
humann_databases --download chocophlan full /db/humann
humann_databases --download uniref uniref90_diamond /db/humann
humann_databases --download utility_mapping full /db/humann   # for regroup/rename
```

Conceptual prerequisites:
- A metagenome measures capability, not expression. For activity, pair with metatranscriptomics (RNA HUMAnN normalized by matched DNA).
- Host-deplete and quality-trim before HUMAnN; host reads inflate UNMAPPED and waste DIAMOND time.
- Keep UNMAPPED/UNINTEGRATED through normalization; they encode database coverage and often track the phenotype.
- UniRef90 suits gut/host-associated samples; novel/environmental biomes have large genuine UNMAPPED and may need UniRef50 or the assembly route.

## Quick Start
Tell your AI agent what you want to do:
- "Profile gene families and pathways for my metagenome, reusing the MetaPhlAn profile"
- "Regroup gene families to KEGG Orthologs and normalize to CPM"
- "Compare pathways between groups without dropping the unmapped fraction"
- "Decide read-based HUMAnN vs assembly-based functional annotation for my soil samples"

## Example Prompts

### Basic profiling
> "Run HUMAnN 3 on sample.fastq.gz reusing sample_metaphlan.tsv, keep the temp output, then normalize gene families to CPM."

### Differential pathways done right
> "Compare pathway abundance between healthy and disease, keep UNMAPPED and UNINTEGRATED, check whether the unmapped fraction differs by group, and use a compositional method rather than a raw Mann-Whitney."

### Biome-aware database choice
> "My samples are marine and over half the reads are UNMAPPED with UniRef90. Help me decide between dropping to UniRef50 and switching to the assembly route."

### Specialized function
> "I care about carbohydrate-active enzymes. Should I use HUMAnN/UniRef or dbCAN, and do I need to assemble first?"

## What the Agent Will Do
1. Host-deplete and trim, then run HUMAnN reusing a MetaPhlAn profile.
2. Normalize RPK to CPM per sample, join, regroup to KO/EC/GO, and split stratified from unstratified tables.
3. Keep UNMAPPED/UNINTEGRATED and check they do not confound the group comparison.
4. Run differential abundance on the unstratified table with a compositional method (MaAsLin2/ANCOM-BC).
5. State that the result is functional potential and flag the unclassified stratified fraction.
6. Recommend the assembly route or a specialized database (dbCAN, antiSMASH) when context or niche function is the goal.

## Tips
- Concatenate paired-end reads before running; HUMAnN has no native pairing flag.
- Do not `--remove-temp-output` if you want to reuse the MetaPhlAn profile.
- `humann_regroup_table` adds an UNGROUPED row - the analogue of UNINTEGRATED; keep it.
- Run statistics on the unstratified table; stratified features are heavily zero-inflated.
- A large UNMAPPED in a novel biome is the environment, not a QC failure.

## Database Options

| Database | Resolution | When |
|----------|-----------|------|
| uniref90_diamond | finer families, higher specificity | gut/host-associated; well-covered biomes |
| uniref50_diamond | coarser, higher sensitivity to divergent homologs | novel/environmental; high UNMAPPED |
| uniref90_ec_filtered | EC-annotated subset, fastest | when only enzyme-level function is needed |

## Resources
- [HUMAnN 3 User Manual](https://huttenhower.sph.harvard.edu/humann/)
- [bioBakery Tools](https://github.com/biobakery/biobakery)

## Related Skills

- metaphlan-profiling - The taxonomic prescreen HUMAnN reuses
- abundance-estimation - Compositional normalization shared with functional tables
- amr-detection - Dedicated ARG quantification
- metagenome-visualization - Plot and test functional tables
- contamination-controls - Host depletion before HUMAnN
- genome-assembly/metagenome-assembly - Assembly route for contextualized function
- pathway-analysis/kegg-pathways - Organism-centric pathway interpretation
