# epidemiological-genomics

## Overview

Pathogen surveillance and outbreak genomics from raw reads to public-health-actionable output. Covers isolate typing (MLST, cgMLST, in-silico serotyping, MTBC and SARS-CoV-2 lineage assignment), AMR detection with mobile-element context and WHO-catalogue-anchored phenotype prediction, time-scaled phylogenies and R_e estimation under recombination-aware clock models, transmission-tree inference with pathogen-tuned SNP thresholds, and wastewater lineage deconvolution. Distinct from `metagenomics/` (community profiling, not isolate-focused) and `variant-calling/` (general germline / somatic; pathogen-specific surveillance lives here).

**Tool type:** mixed | **Primary tools:** AMRFinderPlus, hAMRonization, mlst, chewBBACA, Pangolin, Nextclade, TreeTime, BEAST2 (BDSKY + MASCOT + BICEPS), TransPhylo, outbreaker2, Freyja, TB-Profiler, Kleborate

## Skills

| Skill | Description |
|-------|-------------|
| amr-surveillance | Acquired and point-mutation AMR detection (AMRFinderPlus with `--organism`, RGI/CARD, ResFinder 4.0, TB-Profiler, Mykrobe), hAMRonization to PHA4GE schema, MGE context (MOB-suite, PlasmidFinder, MobileElementFinder), WHO Mtb 2nd-edition catalogue interpretation, heteroresistance via deep variant calling |
| pathogen-typing | 7-locus MLST, cgMLST / wgMLST (chewBBACA, BIGSdb, EnteroBase HierCC), in-silico serotyping (SISTR / SeqSero2 for Salmonella; SerotypeFinder for E. coli; Kaptive K/O via Kleborate; SeroBA for pneumococcus; spa + SCCmec for S. aureus), MTBC Coll/Napier 90-SNP barcode, SARS-CoV-2 Pangolin UShER + Nextclade with version pinning |
| phylodynamics | Time-scaled trees (TreeTime, BEAST2 with BICEPS / BDSKY / MASCOT / ORC clock); recombination masking via Gubbins / ClonalFrameML mandatory for bacteria; date-randomisation tests; structured-coalescent migration (MASCOT, MASCOT-GLM); UShER + matUtils for pandemic-scale; multi-chain BEAST2 convergence diagnostics |
| transmission-inference | outbreaker2 / TransPhylo / phybreak / BadTrIP / SCOTTI / transcluster; pathogen-specific SNP cluster thresholds (TB <=12 / Walker 2013; MRSA <=15 / Coll 2017; C. diff <=2 / Eyre 2013; Klebs <=21 / Snitkin 2012); within-host diversity and transmission bottleneck (McCrone 2018; Lythgoe 2021); HIV-TRACE subtype-aware clustering; Bayesian source attribution |
| variant-surveillance | Nextstrain Augur + Auspice; Pangolin UShER mode (pangoLEARN deprecated mid-2023); Nextclade clade + QC; wastewater deconvolution (Freyja, COJAC, alcov, lineagespot); ARTIC primer scheme version awareness (V3 / V4.1 / V5.3.2 / Midnight 1200); pangolin-data and Nextclade dataset version pinning; recombinant X-prefix detection |

## Example Prompts

- "Type my Salmonella isolates with 7-locus MLST and chewBBACA cgMLST, then cluster at the EFSA harmonised <=5 allele threshold"
- "Build a time-scaled tree for this Mycobacterium tuberculosis outbreak with TreeTime and estimate R_e with BEAST2 BDSKY"
- "Find acquired AMR determinants in 200 Klebsiella assemblies and tell me which carbapenemases are plasmid-borne via MOB-suite"
- "Run TB-Profiler on these MTB assemblies against the WHO 2nd-edition catalogue and report Group 1 / 2 / 3 / 4 / 5 grading per drug"
- "Infer the most likely transmission tree for this hospital MRSA outbreak using outbreaker2 with the contact-tracing matrix"
- "Assign Pango lineages and Nextclade clades to these 1200 SARS-CoV-2 consensus sequences with full pangolin-data version pinning"
- "Deconvolve these wastewater samples into SARS-CoV-2 lineage frequencies with Freyja; verify barcode date postdates sample collection"
- "Is this Clostridioides difficile cluster a real outbreak or convergent ribotype 027 evolution after Gubbins masking?"
- "Build a Nextstrain Augur build for our influenza H3N2 surveillance from the last flu season with explicit subsampling sensitivity"
- "Mask recombination with Gubbins on this S. pneumoniae core-genome alignment before running BEAST2 BDSKY"
- "Cross-check AMRFinderPlus, ResFinder, and RGI calls for these E. coli isolates and produce a WHO GLASS-formatted report via hAMRonization"
- "Distinguish recent TB transmission from reactivation using TransPhylo with within-host coalescent"

## Requirements

```bash
conda install -c bioconda ncbi-amrfinderplus resfinder rgi abricate hamronization abritamr staramr mlst chewbbaca sistr_cmd seqsero2 serotypefinder kleborate kaptive seroba poppunk pangolin nextclade tb-profiler mykrobe mash skani snippy snp-dists gubbins clonalframeml iqtree treetime freyja cojac usher matutils mob_suite plasmidfinder mefinder hiv-trace augur ivar lofreq
```

```bash
conda install -c bioconda beast2
packagemanager -add BDSKY BEASTLabs feast ORC MASCOT BICEPS
```

```r
install.packages(c('TransPhylo', 'outbreaker2', 'phybreak', 'transcluster',
                   'BactDating', 'bdskytools', 'coda', 'ape', 'igraph', 'ggplot2'))
BiocManager::install('lineagespot')
```

```bash
amrfinder -u
tb-profiler update_tbdb
pangolin --all-versions
nextclade dataset list --tag latest sars-cov-2
freyja --version
```

## Related Skills

- **phylogenetics** - Tree inference (IQ-TREE / RAxML / BEAST) and divergence dating fundamentals that phylodynamics builds on
- **variant-calling** - Per-isolate variant calling that feeds core-SNP alignments for bacterial typing and transmission
- **comparative-genomics** - Pangenome analysis (Panaroo / Roary), core-genome alignment, ANI species delineation
- **metagenomics** - AMR detection in mixed communities (NOT isolates); strain tracking via Kraken2 / StrainPhlAn
- **long-read-sequencing** - Hybrid assembly for plasmid resolution; ONT-only WGS for outbreak triage
- **clinical-databases** - Variant interpretation infrastructure adjacent to pathogen surveillance
- **read-alignment** - Read mapping upstream of variant calling
- **read-qc** - Sequencing QC upstream
- **data-visualization** - Phylodynamic trajectory / lineage frequency / transmission network plotting
- **workflows** - End-to-end outbreak pipelines orchestrating typing + AMR + phylodynamics + transmission
