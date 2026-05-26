# Transmission Inference Usage Guide

## Overview

This skill infers who-infected-whom in outbreaks using pathogen genomics. It covers SNP-cluster definition (snp-dists + Gubbins for bacteria; HIV-TRACE for HIV) with pathogen-specific thresholds and explicit population caveats, Bayesian who-infected-whom inference (outbreaker2 with contact data; TransPhylo from a dated tree; phybreak; BadTrIP / SCOTTI for mixed-strain), source attribution (Mather 2013 islandR framework), within-host diversity and transmission bottleneck estimation (Sobel Leonard 2017 beta-binomial), and the conceptual distinctions that surveillance teams routinely conflate: cluster vs WIWS vs source-attribution; generation interval vs serial interval; pairwise SNP distance vs epidemiological linkage. The decision-grade content is in the caveats: there is NO universal SNP cutoff for transmission, direction of transmission cannot be inferred from pairwise SNP alone, unsampled intermediates inflate apparent R_e, and Walker 2013 thresholds derived from UK low-transmission settings inflate apparent recent-transmission rates 2-5x in high-burden settings.

Distinct from `phylodynamics` (which owns R_e estimation from tree priors) and from `pathogen-typing` (which owns cluster definition without WIWS inference).

## Prerequisites

```bash
conda install -c bioconda snippy snp-dists gubbins iqtree treetime hiv-trace beast2

R -e "install.packages(c('TransPhylo', 'outbreaker2', 'phybreak', 'BactDating', 'ape', 'igraph', 'transcluster'), repos='https://cran.r-project.org')"

pip install pandas biopython
```

## Quick Start

Tell the AI agent what transmission question to answer:

- "Define outbreak clusters from this MRSA cohort using the Coll 2017 <=15 SNP threshold after Gubbins"
- "Infer who infected whom in this densely sampled TB cluster using outbreaker2 with contact data"
- "Run TransPhylo on this dated Klebsiella tree to identify likely transmission pairs and unsampled intermediates"
- "Estimate the transmission bottleneck size from these donor-recipient deep-sequenced influenza pairs"
- "Cluster these HIV pol sequences with HIV-TRACE at the appropriate threshold for subtype C (NOT the US-CDC 1.5% default)"
- "Run source attribution for these Salmonella DT104 isolates against poultry / cattle / human reference panels"
- "Distinguish recent TB transmission from reactivation using TransPhylo within-host coalescent"

## Example Prompts

### Hospital MRSA outbreak with rich epi data

> "Forty methicillin-resistant *Staphylococcus aureus* isolates from a 6-month ICU surveillance project with full contact-tracing data (ward histories, overlap dates). Build core-genome alignment via snippy-core; mask recombination via Gubbins on `core.full.aln`; compute pairwise SNP distances via snp-dists; define clusters at Coll 2017 <=15 SNP threshold for hospital outbreaks. Then run outbreaker2 with the contact-tracing matrix `ctd`, MRSA generation-time prior, and explicit sampling proportion. Report posterior WIWS with credibility, unsampled-case count, and joint generation-interval posterior. Compare cluster definitions from raw SNP threshold vs outbreaker2 -- where they disagree, trust outbreaker2."

### TB outbreak with possible reactivation

> "Eighty *Mycobacterium tuberculosis* isolates from a multi-year clinical cohort with suspected reactivation cases. Snippy + snippy-core + Gubbins; build a dated tree with BactDating (Poisson model, not mixedgamma); pass to TransPhylo with TB-tuned generation-time prior (w.scale = months); within-host effective population size prior `startNeg` from literature; sampling-time prior `ws.*`. Run transcluster as a secondary check with TB-tuned SNP+time prior. Identify samples with implausibly early inferred transmission dates relative to clinical onset -- these are likely reactivations of years-old strains. Report transmission cluster membership AND reactivation flag per case."

### Source attribution for foodborne outbreak

> "One hundred *Salmonella* Typhimurium DT104 isolates from a national outbreak; reference collections from poultry (n=300), cattle (n=250), pigs (n=150), human (n=400). Run Mather 2013 *Science* 341:1514 Bayesian source attribution framework via islandR. Report per-source posterior probability per query isolate. Caveat: source-attribution circularity -- the model reproduces the host distribution of training data unless explicitly corrected. Run a rarified-reference sensitivity analysis and report multiple scenarios."

### HIV transmission cluster definition with subtype awareness

> "Two hundred HIV-1 pol sequences from a regional surveillance project; mixed subtype A and C from East / Southern Africa. Do NOT apply the US-CDC HIV-TRACE 1.5% threshold globally (subtype-B-tuned for MSM cohorts). Calibrate threshold per subtype against locally documented transmission pairs; if no local calibration exists, use the published subtype-specific threshold from UKHSA / ECDC and document the choice. Compare cluster definitions across subtype and discuss subtype-C diversity inflating cluster size at the US default."

### Influenza transmission bottleneck from deep sequencing

> "Ten household influenza A donor-recipient pairs with deep paired-end Illumina sequencing (>1000x median coverage). Call within-host minor variants with lofreq at MAF >=1%. Apply the Sobel Leonard 2017 *J Virol* 91:e00171-17 beta-binomial estimator to each pair to estimate transmission bottleneck size Nb. Report Nb posterior per pair; flag any pair with Nb >100 as inconsistent with the narrow-bottleneck literature (McCrone 2018 *eLife* 7:e35962 1-2 virions). Note that narrow bottleneck means consensus genomes are near-identical across transmission -- consensus-only transmission inference would lose this information."

### SARS-CoV-2 outbreak with ARTIC primer awareness

> "Thirty SARS-CoV-2 consensus sequences from a hospital outbreak; samples sequenced with mixed ARTIC primer schemes (some V3, some V4.1, some V5.3.2). Per-sample mask the dropout amplicons (V4.1 amplicons 64, 76, 88-90 in particular); compute pairwise SNP distance only on positions called in ALL samples; warn that primer-scheme heterogeneity makes naive cross-sample comparison invalid. Note that SARS-CoV-2 cluster definition is NOT defined by SNP alone; combine 0-2 SNPs + epi link + sampling window (Lythgoe 2021)."

## What the Agent Will Do

1. Disambiguate the transmission question: cluster definition vs who-infected-whom vs source attribution. Methods and data requirements differ.
2. For cluster definition, run snippy + snippy-core + Gubbins for bacteria, then snp-dists; apply pathogen-specific threshold (Walker 2013 TB; Coll 2017 MRSA; Eyre 2013 C. diff; Snitkin 2012 Klebs; EnteroBase / EFSA cgMLST for Salm / Listeria / E. coli; HIV-TRACE subtype-calibrated for HIV).
3. For dense outbreak with contact data, run outbreaker2 with `outbreaker_data(... ctd=contact_matrix)` and `create_config(n_iter=N)`.
4. For sparser outbreak or larger N, run TransPhylo on a BactDating- or BEAST-dated tree; supply pathogen-tuned generation-time and within-host Ne priors.
5. For mixed-strain infections (TB, HIV, prolonged carriage), use BadTrIP / SCOTTI to handle within-host diversity.
6. For longitudinal within-host samples, use BEASTLIER for identifiable within-host partition.
7. For source attribution, use Mather 2013 framework via islandR; rarify or inverse-weight reference panels to avoid sampling-intensity circularity.
8. For bottleneck size, use Sobel Leonard 2017 beta-binomial on deep-sequenced donor-recipient pairs.
9. Always document: SNP threshold AND its source population; whether `w_dens` encodes generation-interval or serial-interval; sampling proportion; primer scheme version per isolate (for amplicon surveillance).
10. Use "transmission consistent with genomics" not "transmission demonstrated" for direction claims based on pairwise SNP alone.

## Tips

- There is NO universal SNP cutoff. TB <=12 / MRSA <=15 / C. diff <=2 / Klebsiella <=21 / Salm cgMLST <=5 / HIV-TRACE 1.5%. Cite the source population for every threshold.
- Walker 2013 UK TB thresholds inflate apparent recent-transmission rates 2-5x in high-prevalence settings (Cape Town, Mumbai).
- Direction of transmission cannot be inferred from pairwise SNP distance alone (Worby, Lipsitch & Hanage 2014). Requires within-host samples, contact tracing, or both.
- Unsampled intermediates inflate inferred R_e and shorten inferred generation interval. outbreaker2 `pi` parameter or TransPhylo / SCOTTI model this explicitly.
- Generation interval (infection-to-infection) and serial interval (symptom-to-symptom) are NOT the same when incubation periods are highly variable or pre-symptomatic transmission is substantial. outbreaker2 and EpiNow2 want generation interval; using SI biases R_e.
- HIV-TRACE 1.5% is a US-CDC subtype B default. Under-clusters subtype C in southern Africa; over-clusters partial-genome longitudinal samples. Use locally validated thresholds.
- For TB / HIV / chronic infections, naive SNP cutoffs fail because of reactivation and within-host coalescence. TransPhylo / outbreaker2 with within-host-aware priors are mandatory.
- McCrone 2018 *eLife* 7:e35962 (influenza 1-2 virions) and Lythgoe 2021 *Science* 372:eabg0821 (SARS-CoV-2 <10 virions) showed narrow transmission bottlenecks. Consequence: donor and recipient consensus genomes are near-identical by default -- consensus-only transmission inference is WORSE than naive coalescent intuition predicts.
- Source-attribution circularity: the model reproduces training-data host distribution. Inverse-weight or rarify reference panels.
- For ARTIC SARS-CoV-2 cross-version transmission inference, primer dropouts produce N's that look like deletions. Mask failed amplicons; exclude positions in any sample's dropout regions.
- outbreaker2's iterations are set in the config object via `create_config(n_iter=N)`, NOT as `iters` to `outbreaker()`.
- TransPhylo `w.*` is generation-time Gamma prior; `ws.*` is sampling-time Gamma prior. Both must match the pathogen's biology.
- For ANY high-stakes claim (transmission direction, source attribution), supplement with forward simulation (SLiM, FAVITES, SEEDY) -- routinely under-done in published surveillance.

## Related Skills

- pathogen-typing - SNP-cluster / cgMLST cluster definition feeds WIWS inference
- phylodynamics - Time-scaled tree from BactDating / BEAST / TreeTime feeds TransPhylo
- amr-surveillance - Resistant-clone outbreak inference combines AMR + transmission
- variant-surveillance - Lineage assignment cross-checks transmission cluster boundaries
- phylogenetics/divergence-dating - General-purpose dating beyond outbreaks
- phylogenetics/bayesian-inference - BEAST mechanics beyond outbreak phylodynamics
- comparative-genomics/whole-genome-alignment - Core-genome alignment for SNP-typing
- variant-calling/vcf-basics - Per-isolate variant calls for SNP-typing
- variant-calling/variant-calling - SNP calling that feeds snp-dists
- read-alignment/bwa-alignment - Read mapping upstream
- data-visualization/network-visualization - Transmission tree visualisation
- workflows/somatic-variant-pipeline - End-to-end orchestration patterns
