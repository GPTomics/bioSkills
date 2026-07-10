# Outbreak Pipeline Usage Guide

## Overview

This workflow orchestrates outbreak investigation from whole-genome sequencing data. It combines pathogen typing (MLST, cgMLST, in-silico serotyping; Pangolin / Nextclade lineage assignment for viruses), AMR detection (AMRFinderPlus with species mode, TB-Profiler for Mtb, hAMRonization for cross-tool reporting), recombination-aware phylodynamics (snippy + Gubbins on `core.full.aln`, TreeTime or BEAST2 BDSKY with TempEst R^2 QC), and transmission inference (outbreaker2 / TransPhylo / transcluster with pathogen-specific SNP thresholds). The workflow defers methodological depth to the underlying epidemiological-genomics skills.

## Prerequisites

```bash
conda install -c bioconda mlst chewbbaca ncbi-amrfinderplus hamronization tb-profiler \
    snippy snp-dists gubbins clonalframeml iqtree treetime pangolin nextclade freyja \
    mob_suite plasmidfinder sistr_cmd seqsero2 kleborate

conda install -c bioconda beast2
packagemanager -add BDSKY BEASTLabs feast ORC MASCOT BICEPS

Rscript -e "install.packages(c('TransPhylo', 'outbreaker2', 'BactDating', 'bdskytools', 'coda', 'ape'))"

amrfinder -u
tb-profiler update_tbdb
```

**Required databases:**
- MLST schemes (auto-downloaded by mlst); chewBBACA cgMLST schema (record source + fetch date)
- AMR databases: NCBI ReferenceGeneCatalog (AMRFinderPlus), WHO Mtb 2nd-edition catalogue (TB-Profiler)
- Reference genome for core SNP alignment; for SARS-CoV-2, pin pangolin-data + Nextclade dataset version

## Quick Start

Tell your AI agent what you want to do:
- "Investigate this Salmonella outbreak with genomic data"
- "Build a transmission network from my isolate genomes"
- "Type my bacterial isolates and check for AMR genes"
- "Run phylodynamic analysis on my outbreak samples"

## Example Prompts

### Complete pipeline
> "Run outbreak investigation on my Klebsiella isolate genomes"

> "Build a transmission tree from my hospital outbreak samples"

### Typing and AMR
> "MLST type my E. coli isolates"

> "Screen my isolates for antibiotic resistance genes"

> "Compare AMR profiles across outbreak isolates"

### Phylodynamics
> "Build a time-scaled phylogeny from my outbreak with collection dates"

> "Estimate the clock rate for this outbreak"

### Transmission
> "Infer who infected whom in my outbreak using outbreaker2 with contact-tracing data"

> "Estimate R_e through time for my outbreak using BEAST2 BDSKY with multi-chain convergence diagnostics"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Isolate genomes | FASTA | Assembled genomes for each isolate |
| Collection dates | TSV | Sample names with dates (YYYY-MM-DD) |
| Reference genome | GenBank/FASTA | For core SNP alignment |
| Location data (optional) | TSV | Geographic coordinates for mapping |

## What the Agent Will Do

1. **Typing**: 7-locus MLST + cgMLST via chewBBACA; in-silico serotyping per organism (SISTR/SeqSero2 for Salmonella; Kleborate for Klebsiella; SeroBA for pneumococcus; spa + SCCmec for S. aureus). For SARS-CoV-2, Pangolin `--analysis-mode usher` + Nextclade with pinned dataset version.
2. **AMR**: AMRFinderPlus with `--organism` for species-specific point-mutation panel; TB-Profiler against the WHO 2nd-edition catalogue for *M. tuberculosis* (AMRFinderPlus has no Mtb organism mode); hAMRonization to merge per-tool output to the PHA4GE schema.
3. **Core alignment + recombination masking (bacteria)**: snippy + snippy-core; Gubbins on `core.full.aln` (NOT `core.aln`) to mask recombinant tracts before any clock inference.
4. **Phylogenetics**: IQ-TREE 2 with `+ASC` (ascertainment correction required for SNP-only alignment) on the recombination-masked sites.
5. **Phylodynamics**: TempEst root-to-tip R^2 >= 0.3 minimum first; TreeTime ML for outbreak-scale or BEAST2 BDSKY for posterior R_e; for multi-deme sampling use MASCOT (NOT BEAST DTA which is sampling-biased).
6. **Transmission inference**: outbreaker2 (dense + contact data) OR TransPhylo (sparse / longer time-scale from dated tree); apply pathogen-specific SNP cutoff for cluster definition with the source population caveated.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| snippy --mincov | 10 | Minimum coverage for variant call |
| Gubbins input | core.full.aln | NOT core.aln; variable-only alignment cannot estimate background SNP density |
| IQ-TREE model | GTR+G+ASC | +ASC for SNP-only post-Gubbins input |
| TreeTime --clock-filter | 4 | SD multiplier on root-to-tip residual; TreeTime convention |
| TreeTime R^2 minimum | 0.3 | Below this, temporal signal treated as insufficient (field convention; TempEst itself sets no cutoff) |
| Generation time | pathogen-specific | Bacteria ~7-14 days; TB years; SARS-CoV-2 ~5 days; never hardcode |
| TransPhylo mcmcIterations | 1e5-1e6+ | 10k is a smoke-test only |
| BEAST2 BDSKY origin | > rootHeight | Setting origin == rootHeight biases R_e upward systematically |
| BEAST2 chains | >=3-4 independent | Single-chain ESS necessary but not sufficient |
| Pangolin --analysis-mode | usher | pangoLEARN deprecated mid-2023 |
| Pangolin pangolin-data | pinned + dated | Cross-time-point comparisons require a single pangolin-data version |

## Pathogen-Specific SNP / cgMLST Cluster Thresholds

No universal cutoff. See epidemiological-genomics/transmission-inference for the full table with citations.

| Pathogen | Threshold | Source |
|----------|-----------|--------|
| *M. tuberculosis* (core SNP) | <=12 SNPs (likely); <=5 (recent) | Walker 2013 *Lancet Infect Dis* 13:137 -- UK low-transmission; inflates 2-5x in high-burden |
| *S. aureus* (core SNP) | <=15 SNPs (within hospital) | Coll 2020 *Lancet Microbe* 1:e328 |
| *K. pneumoniae* (KPC) | <=21 SNPs | Field convention; no threshold source |
| *C. difficile* (recombination-masked) | <=2 SNPs (direct) | Eyre 2013 *NEJM* 369:1195 |
| *Salmonella* (cgMLST) | <=5 alleles | EFSA / EnteroBase convention |
| *Listeria* (PulseNet cgMLST) | <=4 alleles | PulseNet protocol |
| SARS-CoV-2 | NOT SNP alone | 0-2 SNPs + epi link + sampling window |
| HIV-1 subtype B | 1.5% TN93 | HIV-TRACE US-CDC default; re-tune for non-B |

## Pathogen-Specific Clock Rates and Generation Times

| Pathogen | Generation Time | Clock Rate | Notes |
|----------|-----------------|------------|-------|
| Salmonella | 7-14 days | ~1-3e-7 subs/site/year | Hospital and foodborne outbreaks |
| E. coli | 7-10 days | ~1-2e-7 | STEC outbreaks |
| Klebsiella | 7-14 days | ~1-3e-7 | Hospital AMR / KPC outbreaks |
| Mtb | years (reactivation possible) | ~0.3-0.5 SNPs/genome/year | Within-host coalescence matters; use TransPhylo |
| S. pneumoniae | days | ~1.5e-6 subs/site/year (recombination-masked) | >5e-6 indicates unmasked recombination |
| SARS-CoV-2 ancestral | ~5 days | ~8e-4 subs/site/year | Rate varies by lineage |

## Tips

- Verify clonality first: MLST should show same or closely related STs; cgMLST distance triages outbreak vs distinct introductions.
- Recombination masking with Gubbins on `core.full.aln` is mandatory for recombining bacteria (S. pneumoniae, N. gonorrhoeae, E. coli, Klebsiella, Campylobacter, H. pylori). Skipping inflates the clock rate 2-5x; the date-randomisation test is NOT a sufficient guard.
- Check temporal signal via TempEst root-to-tip regression (Rambaut 2016). R^2 >= 0.3 is a field convention, not a threshold from that paper. If R^2 < 0.3, time-scaling is not supported.
- AMRFinderPlus has NO Mtb organism mode; for any Mtb resistance call route to TB-Profiler interpreted against the WHO 2nd-edition catalogue. Report Group 3 (Uncertain) mutations explicitly -- never collapse to "S".
- For cross-tool AMR comparison, hAMRonization is mandatory; raw AMRFinderPlus / RGI / ResFinder outputs are not directly comparable (CARD ARO vs NCBI ReferenceGeneCatalog vs ResFinder gene-name registries).
- For SARS-CoV-2 lineage assignment, Pangolin `--analysis-mode usher` is the modern default (pangoLEARN deprecated mid-2023); pin pangolin-data version explicitly for any cross-time-point comparison.
- Pathogen-specific SNP thresholds: never apply Walker 2013 UK TB cutoff in high-prevalence settings without caveat; it inflates apparent recent-transmission rates 2-5x.
- For BEAST2 BDSKY, origin > rootHeight always; setting them equal systematically biases R_e upward (Stadler 2013). Run >=3-4 independent chains and combine only after marginal posteriors overlap.
- R_e (effective reproduction number under current conditions) is NOT R_0 (basic, fully susceptible population). Use R_e (or R_t) consistently in reports.
- Integrate epi data: outbreaker2 with contact-tracing matrix `ctd` materially improves who-infected-whom posterior over genomic-only methods.

## Related Skills

- epidemiological-genomics/pathogen-typing - MLST / cgMLST / serotyping / Pangolin / Nextclade details
- epidemiological-genomics/amr-surveillance - AMRFinderPlus species mode, TB-Profiler, hAMRonization, MGE context
- epidemiological-genomics/phylodynamics - TreeTime / BEAST2 BDSKY / BICEPS / MASCOT / Gubbins recombination masking
- epidemiological-genomics/transmission-inference - outbreaker2 / TransPhylo / pathogen-specific SNP thresholds
- epidemiological-genomics/variant-surveillance - Nextstrain Augur, Pangolin UShER, Freyja wastewater deconvolution
- phylogenetics/divergence-dating - General-purpose dating beyond outbreaks
- phylogenetics/bayesian-inference - BEAST mechanics beyond BDSKY
- comparative-genomics/whole-genome-alignment - Core-genome alignment input
- variant-calling/vcf-basics - Per-isolate variant calls for SNP-typing
- read-alignment/bwa-alignment - Read mapping upstream
- database-access/sra-data - Download outbreak FASTQ from SRA / ENA
- database-access/ncbi-datasets-cli - Bulk-pull pathogen reference genomes
