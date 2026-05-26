# Variant Surveillance Usage Guide

## Overview

This skill assigns pathogen lineages and tracks variant frequencies over time for ongoing surveillance. It covers SARS-CoV-2 Pango lineage assignment via Pangolin UShER mode (pangoLEARN deprecated mid-2023), Nextclade clade + QC + mutation reporting, Nextstrain Augur surveillance pipelines, wastewater lineage deconvolution (Freyja, COJAC, alcov, lineagespot), lineage fitness modelling (Wenseleers / Bedford-Figgins multinomial logistic), recombinant detection (3SEQ, RDP4, Bolotie), and primer-scheme awareness for ARTIC V3 / V4.1 / V5.3.2 / Midnight 1200 amplicon surveillance. The decision-grade content is in the version discipline: pin pangolin-data, pin Nextclade dataset, pin Freyja barcode, pin ARTIC primer scheme per isolate, re-run the whole archive when versions update -- cross-time-point comparisons across dataset versions are silently invalid otherwise.

Distinct from `pathogen-typing` (which owns the typing call) and from `phylodynamics` (which owns R_e estimation from tree priors).

## Prerequisites

```bash
conda install -c bioconda pangolin nextclade augur freyja cojac usher matutils samtools lofreq ivar snakemake
pip install pandas biopython

R -e "BiocManager::install('lineagespot')"
```

## Quick Start

Tell the AI agent what surveillance question to answer:

- "Assign Pango lineages and Nextclade clades to these 1200 SARS-CoV-2 consensus sequences using UShER mode"
- "Deconvolve these wastewater samples into SARS-CoV-2 lineage frequencies with Freyja and tell me what JN.* sub-lineages are growing"
- "Detect emerging variants in wastewater early with COJAC co-occurrence on ARTIC V5.3.2 amplicons"
- "Build a Nextstrain Augur build for our regional SARS-CoV-2 surveillance from the last 6 months"
- "Compute lineage growth advantages using multinomial logistic regression and report covariance among lineages"
- "Pin pangolin-data and Nextclade dataset versions for this archive so calls are reproducible"
- "Check for recombinant candidates in this batch using Bolotie and reconcile with Pango designations"

## Example Prompts

### Pango lineage assignment with full version pinning

> "Twelve hundred SARS-CoV-2 consensus FASTAs from January through July 2025. Run Pangolin with `--analysis-mode usher` (UShER default since v4) and Nextclade run with the current sars-cov-2 dataset. Record `pangolin --all-versions` and the Nextclade dataset tag from `pathogen.json` for every batch. Reconcile any UShER vs pangoLEARN disagreement (UShER takes precedence). For any sample where Pango and Nextclade disagree on recombinant assignment (X-prefix vs parent lineage), flag for manual review. Pin pangolin-data version before producing the final report; if the version updated mid-analysis, re-run the entire batch."

### Wastewater Freyja deconvolution with barcode discipline

> "Twelve weekly wastewater samples from a regional treatment plant; targeted SARS-CoV-2 surveillance. For each sample: verify the Freyja barcode date postdates the sample collection date (if not, `freyja barcode-build` from the current UShER tree); `freyja variants` + `freyja demix` to produce per-lineage abundance. Report the `resid` column per sample (residual mass not assigned to known lineages); high resid indicates a novel lineage being missed. Cross-check ARTIC primer dropouts (V4.1 amplicons 64 / 76 / 88-90 known dropouts); mask failed amplicons before lineage assignment. Compare weekly trajectories of dominant lineages with a 14-day rolling smoother."

### COJAC for early variant detection in wastewater

> "Twenty wastewater samples; we suspect an emerging KP.3 sub-lineage. Run COJAC `cooc-mutbamscan` with KP.3.x defining-mutation YAML and the ARTIC V5.3.2 primer BED; report which samples show co-occurrence of at least 2 KP.3.x signature mutations on the same read pair. COJAC detected Alpha 13 days before clinical samples in Swiss data (Jahn 2022 *Nat Microbiol* 7:1151), so use for early detection; cross-confirm with Freyja for quantitative tracking."

### Nextstrain Augur regional build

> "Build a regional Nextstrain Augur build for SARS-CoV-2 surveillance covering one country and the last 6 months. Pull the latest pathogen build template from github.com/nextstrain/ncov; subsample to ~3000 genomes (document subsampling weights -- Hodcroft 2021 *Nature* 591:30 noted subsampling drives the result more than the data); align via `augur align`; tree via `augur tree`; refine with `augur refine --root oldest --timetree --coalescent skyline`; ancestral / translate / traits; export to Auspice JSON. Run a sensitivity analysis with alternative subsampling and report MRCA-date uncertainty."

### Lineage growth advantage with multinomial logistic regression

> "Eight months of weekly SARS-CoV-2 Pango lineage frequencies from a national surveillance dataset. Fit multinomial logistic regression in the Wenseleers / Bedford-Figgins framework (Abousamra 2024 *PLoS Comput Biol* 20:e1012443); report growth advantage per lineage with full covariance matrix (NOT marginal 95% CIs that hide covariance). Caveat that early growth estimates are systematically inflated and shrink as more time passes."

### Influenza vaccine-strain surveillance

> "H3N2 surveillance for the upcoming Northern Hemisphere flu season. Pull Nextstrain seasonal-flu build for H3N2; assign clades via Nextclade with the H3N2 dataset; report dominant antigenic-cluster shifts; flag egg-adapted mutations (L194P, T160K, G186V) in vaccine candidates. Compare to WHO biannual strain selection (February for NH; September for SH) and current circulating diversity."

## What the Agent Will Do

1. Identify the surveillance target (SARS-CoV-2, influenza, Mpox, RSV, H5N1, measles) and source type (consensus FASTA, paired-end FASTQ, wastewater BAM).
2. For consensus-based assignment, run Pangolin with `--analysis-mode usher` AND Nextclade with current datasets; record `pangolin --all-versions` and Nextclade dataset tag.
3. For wastewater pooled samples, verify Freyja barcode date postdates sample collection; `freyja barcode-build` if needed; `freyja variants` + `freyja demix`; report `resid`.
4. For ARTIC V4.1 / V5.3.2 / Midnight amplicon data, inspect per-amplicon coverage; mask low-coverage amplicons before lineage assignment.
5. Cross-check Pangolin UShER vs Nextclade per the reconciliation table; report both with provenance.
6. For published or regulatory output, pin EVERY tool version AND dataset version; re-run on update.
7. For longitudinal frequency tracking, build Nextstrain Augur pipeline or custom multinomial logistic; report covariance among lineages, not just marginal CIs.
8. For recombinant detection, use Bolotie / 3SEQ / RDP5; trust Pango-designation X-prefix assignments; submit novel candidates to cov-lineages issue tracker.
9. For Nextstrain Augur, document subsampling explicitly; run sensitivity analysis (Hodcroft 2021 *Nature* 591:30).

## Tips

- pangolin-data is updated weekly. Pin the version for any cross-time-point comparison. `pangolin --all-versions` outputs pangolin + pangolin-data + scorpio + constellations; record all four.
- pangoLEARN was deprecated mid-2023 (Pongmoragot 2024 *Virus Evol* 10:vead085). Pangolin `--analysis-mode usher` is the default since v4.
- Nextclade datasets are versioned independently from the executable; different dataset versions assign different mutations as "lineage-defining" because internal-node placement shifts as the tree grows.
- Freyja barcode date MUST postdate the sample collection date. Lineages designated after the barcode are silently invisible; appear as elevated parent-lineage abundance.
- `freyja barcode-build` is the current spelling (hyphen). Some legacy docs show `barcode_build` (underscore); that's wrong in current builds.
- For ARTIC SARS-CoV-2 sequencing, V4.1 amplicons 64 / 76 / 88-90 are known chronic dropouts; V5.3.2 fixed many but introduced others. Document the scheme version per isolate.
- Primer dropout produces N's that LOOK LIKE deletions. Mask failed amplicons; do not treat N-stretches as deletions in downstream lineage analysis.
- Recombinants get X-prefix Pango designations (XBB, XBC, XEC, XFG). Pangolin assigns recombinants to a parent lineage when no X designation exists yet; designation lag is days-to-months.
- Karthikeyan 2022 *Nature* 609:101 reported 11-day Omicron wastewater lead vs clinical, but subsequent analyses with different configurations produced detection lags from -5 to +3 days. The 11-day number is the most-favourable configuration.
- Wastewater RNA-to-cases conversion is NOT constant across variants. The wastewater literature documents variant-specific faecal shedding-rate variability; Omicron BA.1 shed less per case than Delta.
- For Nextstrain Augur, subsampling drives the result more than the data (Hodcroft 2021). Document subsampling configuration; run sensitivity.
- Wenseleers / Bedford-Figgins multinomial logistic-regression CI for one lineage hides covariance with all others. Early growth estimates are systematically too large; report rank-ordered with simultaneous CIs.
- GISAID weekly submission rate fell from ~500,000/week peak to ~5,000-20,000/week through 2024-2025. New-lineage detection lag has increased; emerging lineages often only identified at >5% prevalence in 2-3 countries.
- For H3N2 vaccine surveillance, antigenic drift between strain selection (February for NH) and circulation (October-March) is the recurring failure mode; egg-adapted mutations (L194P, T160K, G186V) alter HA antigenicity during manufacturing.
- MPXV Clade Ib (2024 DRC outbreak) required rapid re-tooling of pipelines built for SARS-CoV-2; pathogen-specific tooling does not transfer cleanly.

## Related Skills

- pathogen-typing - SARS-CoV-2 Pango / Nextclade assignment overlaps; this skill owns longitudinal frequency tracking
- phylodynamics - Lineage-stratified BDSKY / BICEPS R_e runs downstream
- transmission-inference - SARS-CoV-2 cluster definition combines lineage + 0-2 SNPs + epi link
- amr-surveillance - Antiviral drug-resistance mutation tracking
- phylogenetics/modern-tree-inference - Non-UShER topology
- phylogenetics/tree-io - Tree parsing for Augur output
- comparative-genomics/whole-genome-alignment - Reference-based alignment
- variant-calling/vcf-basics - VCF for lineage-defining mutations
- variant-calling/variant-calling - Variant calling for wastewater (lofreq, ivar)
- variant-calling/filtering-best-practices - Per-amplicon coverage filtering
- read-alignment/bwa-alignment - Read mapping upstream
- read-alignment/minimap2-alignment - Long-read for ARTIC-Midnight 1200
- read-qc/quality-reports - Sequencing QC upstream
- database-access/sra-data - SRA / INSDC retrieval
- data-visualization/multipanel-figures - Lineage frequency / wastewater plotting
- workflows/somatic-variant-pipeline - End-to-end orchestration patterns
