# eDNA Metabarcoding Pipeline - Usage Guide

## Overview
Complete environmental DNA metabarcoding workflow from raw amplicon sequences to biodiversity assessment and community ecology. Supports two processing paths: OBITools3 v3 (CLI-based, with the `obi stats` plural command convention and DMS-based ingestion) and DADA2 (R-based, single-nucleotide ASV resolution per Callahan 2017). Includes decontam contamination screening (Davis 2018; treated as classifier-not-truth), Hill-number diversity reported as effective species counts (Jost 2006) with coverage-based rarefaction (Chao & Jost 2012), beta-diversity decomposition with the mandatory PERMANOVA + PERMDISP pair (Anderson & Walsh 2013), and constrained ordination with vegan. Honest reporting of the species-level taxonomic-assignment gap (50-85% unassigned is typical) and the read-counts-are-not-abundance critique (Lamb 2019 meta-analysis).

## Prerequisites
```bash
# OBITools3 v3 (CLI path)
pip install obitools3
# or conda install -c bioconda obitools3

# General CLI tools
pip install cutadapt
# or conda install -c bioconda cutadapt fastqc multiqc

# R packages (DADA2 path + downstream analysis)
install.packages('BiocManager')
BiocManager::install(c('dada2', 'decontam', 'phyloseq'))
install.packages(c('iNEXT', 'vegan', 'indicspecies', 'metabaR'))

# Reference databases (download once, marker-specific):
# COI: BOLD reference library or MIDORI2
# 12S: MitoFish or 12S-seqdb
# ITS: UNITE (https://unite.ut.ee/)
# 16S/18S: SILVA (https://www.arb-silva.de/) or PR2 for 18S
```

**Input data:**
- Demultiplexed paired-end FASTQ files (one pair per sample)
- Sample metadata with sample type (field sample vs. negative control) AND DNA concentration (qPCR or Qubit) for decontam combined method
- Environmental variables for community analysis (e.g., temperature, depth, site)
- Primer sequences for the target marker

## Quick Start
Tell your AI agent what you want to do:
- "Process my COI amplicon data from raw FASTQ through species lists, using DADA2 ASVs and BOLD/MIDORI2 taxonomy"
- "Run the OBITools3 v3 pipeline (note: v3 uses `obi stats` plural and DMS ingestion)"
- "Apply decontam with combined method (concentration + negative controls) at threshold 0.1; treat output as screening, not classification"
- "Quantify tag-jumping rate; for NovaSeq libraries apply ~10x heavier filtering than MiSeq (patterned flow cells; the tag-jump phenomenon is Schnell 2015)"
- "Compare communities with PERMANOVA AND betadisper (mandatory pair per Anderson & Walsh 2013)"
- "Report Hill-number effective species counts (NOT raw Shannon) with coverage-based rarefaction at C=0.95"

## Example Prompts

### Full Pipeline
> "I have paired-end COI amplicon data from 50 water samples plus 5 negative controls. Process everything from raw FASTQ to Hill-number effective species counts at C=0.95 coverage, with combined decontam screening."

> "Run the complete eDNA pipeline: cutadapt linked-adapter primer removal, DADA2 ASVs with method='consensus' chimeras, decontam combined method at threshold 0.1, taxonomy against MIDORI2 with minBoot=80, Hill numbers q=0,1,2 with iNEXT coverage-based rarefaction."

### Specific Steps
> "I already have an ASV table from DADA2. Run contamination filtering with decontam combined method (negative controls + DNA concentration) at threshold 0.1; flag output as screening; manually verify each candidate ASV for biological plausibility."

> "Compare fish communities across upstream vs downstream sites with PERMANOVA + PERMDISP (BOTH must be reported; if betadisper is also significant, the location conclusion is not supported)."

### Marker-Specific
> "Process my fungal ITS2 eDNA data using DADA2 with the UNITE 9.0+ database."

> "Process MiFish 12S marine fish eDNA with the full pipeline; report Hill-number diversity at standardized coverage."

## What the Agent Will Do
1. Run FastQC and MultiQC on raw reads to assess quality
2. Remove primer sequences with Cutadapt (`--discard-untrimmed`); MANDATORY before DADA2 filterAndTrim (primer-bearing reads corrupt the DADA2 error model)
3. Validate: check reads per sample >1000, negative controls <100 reads
4. Merge paired ends and denoise (OBITools3 v3 with `obi stats`/`obi clean`/`obi ecotag`, OR DADA2 with `method='consensus'` chimeras)
5. Validate: chimera rate <20% (>30% indicates library-prep issues)
6. Apply decontam combined method (DNA concentration + negative controls; threshold 0.1 default, 0.05 for low-biomass); treat output as SCREENING and require biological-plausibility check per flagged ASV
7. Quantify residual tag-jumping rate from per-ASV cross-sample appearance; apply `metabaR::tagjumpslayer` with platform-appropriate threshold (~0.001-0.005 MiSeq; ~0.005-0.01 NovaSeq patterned flow cells)
8. Validate: assignment rate meets marker expectations (50-85% unassigned at species level is typical; report this gap honestly)
9. Compute Hill numbers q=0,1,2 as effective species counts (NOT raw Shannon/Simpson) with iNEXT coverage-based rarefaction at C=0.95; bound extrapolation at 2x reference sample size (Chao et al. 2014 doubling rule)
10. Validate: sample completeness >80%
11. Compare communities with PERMANOVA (`adonis2`) AND PERMDISP (`betadisper`); MANDATORY pair (Anderson & Walsh 2013). If betadisper is significant, location conclusion is not supported
12. Run constrained ordination (RDA on Hellinger-transformed data if DCA gradient <=3 SD; CCA if >3 SD)
13. Indicator species analysis with `multipatt(func='IndVal.g')` (group-equalized, NOT basic IndVal)
14. Export species table, Hill-number effective species counts, ordination plots; report PERMANOVA + PERMDISP together; document decontam manual-review decisions

## Tips
- Always include negative controls (extraction blanks AND PCR blanks) for decontam; without controls, only the frequency method is available (requires DNA concentration)
- Primer removal is MANDATORY before DADA2 filterAndTrim; primer-bearing reads corrupt the DADA2 error model
- OBITools3 v3 uses `obi stats` (plural; was `obistat` in v1) and `.tar.gz` taxonomy ingestion (NOT a directory)
- ASV (DADA2, UNOISE3) vs OTU (swarm v2): Callahan 2017 recommends ASVs for COI/12S/18S/fungi; Schloss 2021 critique applies for bacterial 16S where multi-copy rRNA may split single genomes into multiple ASVs
- Tag-jumping rate is ~10x higher on NovaSeq patterned flow cells than MiSeq (the tag-jump phenomenon is Schnell 2015); MiSeq-tuned thresholds underfilter NovaSeq data
- decontam output is SCREENING, not classification; the default threshold 0.1 is over-aggressive in low-biomass samples (reduce to 0.05 and inspect each flagged ASV)
- Read counts are NOT biomass; the correlation is weak, non-linear, and taxon-specific (Lamb 2019 meta-analysis). Report presence/absence or relative abundance with explicit caveat; quantitative claims require mock-community calibration
- Hill numbers are EFFECTIVE species counts (Jost 2006); report q=0 (richness), q=1 (exp(Shannon)), q=2 (1/Simpson) NOT raw Shannon/Simpson
- Coverage-based rarefaction (iNEXT type 3) is the postdoc-grade default; sample-size rarefaction is biased when assemblages differ in true diversity (Chao & Jost 2012)
- iNEXT extrapolation reliable only up to 2x the reference sample size; the default endpoint enforces this
- PERMANOVA + PERMDISP is mandatory (Anderson & Walsh 2013); never report PERMANOVA alone; if betadisper is significant, the location conclusion is not supported by dispersion-confounded PERMANOVA
- `adonis2()` is the modern PERMANOVA; `adonis()` was deprecated in vegan 2.6+
- For SES_MPD/MNTD, always specify the null-model choice explicitly (`taxa.labels`, `independentswap`, `richness`)
- Use DCA gradient length to choose between RDA (linear, <3 SD) and CCA (unimodal, >3 SD); apply Hellinger transformation before RDA
- For degraded eDNA (e.g., sediment cores), relax DADA2 maxEE to c(5,5) and reduce truncLen
- Increase sequencing depth for samples with low completeness rather than lowering quality thresholds
- For taxonomic ambiguity, phylogenetic placement (EPA-ng + gappa) is more accurate than closest-match assignment though 10-100x slower

## Related Skills
- ecological-genomics/edna-metabarcoding - Detailed eDNA processing with ASV/OTU decision and primer-bias caveats
- ecological-genomics/biodiversity-metrics - Hill-number diversity and Baselga/Podani beta partitions
- ecological-genomics/community-ecology - PERMANOVA + PERMDISP, JSDMs, RLQ + fourth-corner
- read-qc/quality-reports - Raw read quality assessment
- reporting/automated-qc-reports - Aggregate FastQC across samples with MultiQC
- microbiome/amplicon-processing - 16S clinical alternative
