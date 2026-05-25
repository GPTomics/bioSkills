# eDNA Metabarcoding Usage Guide

## Overview

Processes environmental DNA metabarcoding amplicon data from raw paired-end reads to species occurrence tables, navigating: the ASV vs OTU decision (Callahan 2017 "ASVs replace OTUs" vs Schloss 2021 multi-copy bacterial 16S critique; ASVs are recommended for COI/12S/18S/fungi but should be evaluated cautiously for bacterial 16S), marker and primer choice with taxon-specific bias (mlCOIintF/jgHCO2198 for Leray metazoan COI, MiFish-U/E for fish 12S, 515F/806R V4 for bacteria, ITS2 for fungi), OBITools3 v3 command-name break (`obi stats` plural; DMS-based ingestion; `.tar.gz` taxonomy), tag-jumping quantification with dual-indexing (Schnell 2015; NovaSeq has ~10x higher rates than MiSeq), decontam contamination screening as classifier-not-truth (Davis 2018), the read-counts-are-not-abundance critique (Elbrecht 2015; Lamb 2019 meta-analysis), site-occupancy models for detection-corrected species occurrence (Ficetola 2015), eDNA decay kinetics with temperature dependence (Strickler 2015), and 85% taxonomic-gap reporting honesty.

## Prerequisites

- DADA2 (`BiocManager::install('dada2')`)
- cutadapt (`pip install cutadapt`)
- OBITools3 v3 (`pip install OBITools3`)
- decontam (`BiocManager::install('decontam')`)
- microDecon (`remotes::install_github('donaldtmcknight/microDecon')`)
- occumb and JAGS for occupancy modeling (`install.packages('occumb')`; JAGS from https://mcmc-jags.sourceforge.io)
- metabaR for tag-jumping filtering (`install.packages('metabaR')`)
- Reference database appropriate to marker: MIDORI2 (COI/12S/18S), MitoFish (fish 12S), BOLD (animal COI), SILVA (16S/18S), UNITE (fungal ITS)

## Quick Start

Tell your AI agent what you want to do:

- "Trim primers with cutadapt BEFORE DADA2 (mandatory) and denoise to ASVs"
- "Run the OBITools3 v3 pipeline using `obi stats` (plural) and DMS-based ingestion"
- "Quantify tag-jumping rate before filtering; apply NovaSeq-appropriate thresholds (10x MiSeq)"
- "Use decontam as a SCREENING tool; verify each flagged ASV for biological plausibility"
- "Report presence/absence or relative abundance with explicit mock-community calibration; never read counts as absolute biomass"
- "Apply occumb site-occupancy modeling to correct for false negatives in PCR replicates"

## Example Prompts

### Pipeline Selection

> "I have paired-end FASTQ from eDNA water samples amplified with MiFish 12S primers. Trim primers with cutadapt (linked-adapter mode), denoise with DADA2 to ASVs, remove chimeras, and assign taxonomy against MitoFish with minBoot=80 for genus-level confidence."

> "Run the OBITools3 v3 pipeline on my COI data. NOTE: v3 uses `obi stats` (plural, with space), DMS-based ingestion, and `.tar.gz` taxonomy archive (different from v1)."

### Marker / Primer Decision

> "For metazoan eDNA from coral reefs, use mlCOIintF/jgHCO2198 (Leray 2013) -> 313 bp. For fish eDNA, use MiFish-U (Miya 2015) -> 163-185 bp. For bacterial communities, 515F/806R (Parada modified) V4 region. Acknowledge primer-specific bias in interpretation."

### Tag-Jumping

> "Quantify residual tag-jumping rate from per-ASV cross-sample appearance. For NovaSeq, expect ~10x higher rates than MiSeq; apply heavier filtering thresholds (~0.005-0.01 vs ~0.001-0.005)."

### Contamination Screening (NOT Classification)

> "Run decontam with combined frequency + prevalence method (DNA concentration AND negative controls). Treat the output as SCREENING: manually verify biological plausibility of each flagged ASV before deletion. For low-biomass samples, lower the threshold to 0.05."

### Quantitative Interpretation

> "DO NOT report read counts as abundance. Read counts have weak, non-linear, taxon-specific correlation with biomass (Lamb 2019). Either report presence/absence, relative abundance with caveat, OR calibrate with mock community."

### Site Occupancy

> "I collected 3 PCR replicates per site across 30 sites. Fit a multi-species occupancy model with `occumb` to estimate detection-corrected species occurrence. Cite Ficetola 2015 for the false-presence/absence framework."

## What the Agent Will Do

1. Pick marker and primer pair from the literature (Leray COI, MiFish 12S, 515F/806R 16S, ITS2 fungal, trnL plant) and document the choice with primer-bias caveats
2. Trim primers with cutadapt in linked-adapter mode (`--discard-untrimmed`) BEFORE any DADA2 quality filtering; primer-bearing reads corrupt the DADA2 error model
3. Decide ASV vs OTU: ASVs (DADA2 / UNOISE3) for COI / 12S / 18S / fungi; cautiously evaluate ASVs vs OTUs for bacterial 16S due to multi-copy rRNA (Schloss 2021)
4. Run DADA2 with inspected quality profiles for `truncLen`, `maxEE = c(2, 2)`, learn errors per dataset, denoise, merge pairs with `minOverlap`, build sequence table, remove chimeras via `method='consensus'`
5. Assign taxonomy against marker-appropriate curated reference database (MIDORI2 / MitoFish / BOLD / SILVA / UNITE) with `minBoot = 80` for genus-level confidence
6. For OBITools3: use v3 syntax (`obi stats` plural, DMS-based ingestion, `.tar.gz` taxonomy archive); cite Boyer 2016
7. Quantify tag-jumping rate before filtering; apply per-ASV abundance threshold with `metabaR::tagjumpslayer`; use NovaSeq-appropriate threshold (~0.005-0.01) for patterned-flow-cell libraries
8. Run decontam with combined frequency + prevalence method when DNA concentration and negative controls are available; explicitly flag output as SCREENING; manually verify biological plausibility of flagged ASVs
9. Lower decontam threshold to 0.05 for low-biomass samples (open ocean, ancient sediment); cite Salter 2014 for the reagent-contamination caveat
10. NEVER report read counts as abundance without mock-community calibration; reports as presence/absence, relative abundance with caveat, OR calibrated quantity
11. Fit site-occupancy models via `occumb` when replicates are available, citing Ficetola 2015
12. Report the "unassigned" fraction explicitly (typically 50-85% in non-vertebrate metabarcoding)

## Tips

- DADA2's `filterAndTrim` requires primer-free reads; running `filterAndTrim` on primer-bearing FASTQ produces a corrupted error model and ASV artifacts
- DADA2 `method = 'consensus'` chimera detection is conservative and standard; `method = 'pooled'` is aggressive and can over-merge real diversity
- `minBoot = 80` for genus-level taxonomic confidence; `minBoot = 50` for family-level; never accept species-level assignment from `minBoot < 80`
- Naive Bayes "confidence" reported by q2-feature-classifier is a scikit-learn calibration, NOT a true probability (Bokulich 2018); use phylogenetic placement (EPA-ng + gappa) for borderline assignments
- The Leray COI primer (mlCOIintF/jgHCO2198) has known coverage gaps for Platyhelminthes and some Nematoda (Wangensteen 2018); flag this in methods if those phyla are biologically relevant
- MiFish-U/E (Miya 2015) is the dominant fish eDNA marker globally; specifically targets a hypervariable 12S rRNA region (163-185 bp)
- For NovaSeq libraries (patterned flow cells), tag-jumping is ~10x higher than MiSeq; pre-2020 MiSeq-tuned filtering thresholds give false rare-species detections on NovaSeq data
- Decontam default threshold = 0.1 is over-aggressive in low-biomass eDNA (open ocean, ancient sediment); reduce to 0.05 and inspect each flagged ASV
- The "unassigned" fraction in metabarcoding is typically 50-85%; this reflects incomplete reference databases (Wangensteen 2018 demonstrated this for marine metabarcoding); report this honestly
- For COI metazoan metabarcoding, MIDORI2 LONGEST_NUC_GB259_CO1 (DADA2-formatted) is the standard reference; for fish 12S, MitoFish (Miya lab); for fungal ITS, UNITE 9.0+
- eDNA half-life depends strongly on temperature: ~4-15 hours at 20 deg C surface water; ~5-21 days at 4 deg C cold water (Strickler 2015); sample collection-to-storage time matters
- Read counts are NOT biomass: Elbrecht 2015 and Lamb 2019 documented weak, non-linear, primer-specific correlation; mock-community calibration is essential for quantitative inference
- For detection-probability correction with replicates, use `occumb` (multi-species occupancy) or `eDNAoccupancy` (single-species); cite Ficetola 2015
- Always include extraction blanks AND PCR-negative controls; without them, decontam frequency method is the only option (prevalence method requires controls)
- OBITools3 v3 vs v1 command rename: `obistat` -> `obi stats` (plural, with space); `obigrep` -> `obi grep`; v3 expects `.tar.gz` taxonomy archive, not a directory
- Salter 2014 documented that low-biomass samples can have a fictitious "microbiome" composed entirely of reagent-derived contaminant DNA; always run negative controls in proportion (~1 per 8-10 samples)

## Related Skills

- ecological-genomics/biodiversity-metrics - Diversity analysis from species occurrence tables (Hill numbers, Chao1 with singleton-bias caveat)
- ecological-genomics/community-ecology - Environmental gradient analysis of community composition (PERMANOVA + PERMDISP, JSDM)
- microbiome/amplicon-processing - 16S clinical microbiome alternative pipeline
- read-qc/quality-reports - Upstream read-quality assessment before primer trimming
- database-access/entrez-fetch - Retrieve reference sequences for custom taxonomy databases
