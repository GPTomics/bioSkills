---
name: bio-epidemiological-genomics-transmission-inference
description: Infers person-to-person transmission from pathogen genomes using outbreaker2, TransPhylo, phybreak, BadTrIP, SCOTTI, BEASTLIER, and SNP-distance / cluster-picker approaches (HIV-TRACE for HIV; transcluster). Defines outbreak clusters using pathogen-specific SNP thresholds (NOT a universal cutoff -- TB <=12 SNPs; MRSA <=15; C. difficile <=2; Klebsiella <=21), models within-host diversity and transmission bottlenecks, integrates contact-tracing data, distinguishes generation from serial interval, and attributes source via Bayesian source attribution (islandR). Use when investigating outbreaks for who-infected-whom, defining SNP-cluster outbreak definitions, accounting for unsampled intermediates, choosing between outbreaker2 (rich epi data) and TransPhylo (genomic-only after a dated phylogeny), running source attribution between host populations, calling HIV-TRACE thresholds appropriate to the local subtype, or distinguishing recent transmission from reactivation in TB or chronic HIV.
tool_type: mixed
primary_tool: TransPhylo
---

## Version Compatibility

Reference examples tested with: TransPhylo 1.4+ (R), outbreaker2 1.2+ (R), phybreak 0.5+ (R), BadTrIP via BEAST 2.7+ package manager, BEASTLIER via BEAST 1.10+, transcluster 1.0+ (R), HIV-TRACE 1.5+, snp-dists 0.8+, ape 5.8+ (R), igraph 1.6+ (R), TreeTime 0.11+, BactDating 1.1+ (R), BEAST 2.7.6+, lofreq 2.1+, deepSNV via Bioconductor 3.18+, pandas 2.2+, BioPython 1.84+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('TransPhylo')`; `?inferTTree` to confirm arg names
- R: `packageVersion('outbreaker2')`; `?create_config` -- iteration count is set via `n_iter` in the `config` object, NOT as `iters` to `outbreaker()`
- Python: `pip show lofreq`; check whether deep variant calling supports the target MAF
- CLI: `snp-dists --help`; `hiv-trace --help`

If R rejects an argument, the function signature changed between minor releases; `?function_name` is authoritative.

# Transmission Inference

**"Who infected whom in this outbreak, and is this even an outbreak?"** -> Pick the question first (cluster definition vs WIWS who-infected-whom vs source attribution), then the method that fits the data (rich epi + dense sampling -> outbreaker2; sparse sampling + good dated tree -> TransPhylo; longitudinal within-host samples -> BEASTLIER / BadTrIP; rapid surveillance triage -> SNP-distance with pathogen-tuned threshold). Genomic distance is necessary but not sufficient for direction: two isolates 3 SNPs apart could be A->B, B->A, A->Unknown->B, or two-from-one common source. Direction inference requires temporal data, within-host diversity, contact-tracing data, or all three.

- R: `outbreaker2::outbreaker(data=outbreaker_data(dates=..., dna=..., w_dens=..., f_dens=..., ctd=...), config=create_config(n_iter=1e6))` -- dense outbreak with contact data
- R: `TransPhylo::inferTTree(ptree, mcmcIterations=1e5, w.shape=1.3, w.scale=10)` -- sparse outbreak from a dated tree
- CLI: `snp-dists -c gubbins.filtered_polymorphic_sites.fasta > pairwise.csv` -- pairwise SNP for cluster triage
- CLI: `hiv-trace --threshold 0.015` -- HIV cluster definition at the US-CDC default (subtype B); reconsider for non-B subtypes

## The Single Most Important Modern Insight -- There is no universal SNP cutoff for transmission

The pathogen-specific SNP threshold varies by 10x across taxa (TB <=12 SNPs, *C. difficile* <=2, MRSA <=15, *Salmonella* cgMLST <=5, *Klebsiella* <=21, SARS-CoV-2 not defined by SNP alone). Substitution rate, recombination, generation time, within-host diversity, and (for Mpox) APOBEC3 editing all vary by 100x. Walker 2013 *Lancet Infect Dis* 13:137 derived the TB <=12 SNP cutoff from UK Oxfordshire (low-transmission, contact-traced, household settings); applying the same threshold in Cape Town or Mumbai inflates apparent recent-transmission rates 2-5x because clonal isolates linked through long-past common ancestors get pooled with truly recent transmissions. Worby, Lipsitch & Hanage 2014 *PLoS Comput Biol* 10:e1003549 formally showed that within-host bacterial diversity puts an irreducible upper bound on the resolution of SNP-distance transmission-network reconstruction even with repeated sampling. Always cite the pathogen-specific source AND its derivation population; never apply a threshold outside its validated context without an explicit caveat. For TB / HIV / chronic infections, naive SNP cutoffs fail because of reactivation and within-host coalescence -- use TransPhylo or outbreaker2 with within-host-aware priors.

## Algorithmic Taxonomy

| Tool | Mechanism | Inputs | Output | Strength | Fails when |
|------|-----------|--------|--------|----------|------------|
| Pairwise SNP threshold (snp-dists; cluster picker) | Count SNPs between pairs; threshold + linkage | Core-SNP alignment | Adjacency at threshold | Fast; intuitive; standard for surveillance triage | Pathogen-specific cutoff; convergent evolution and recombination violate distance assumptions |
| HIV-TRACE (Kosakovsky Pond 2018 *Mol Biol Evol* 35:1812) | TN93 pairwise distance + threshold (default 1.5%) | HIV-1 pol or other gene | Cluster membership | CDC standard for US HIV surveillance | 1.5% threshold is US-CDC subtype B specific; under-clusters subtype C in southern Africa |
| outbreaker2 (Campbell 2018 *BMC Bioinformatics* 19:363) | MCMC; sequence + generation-interval + sampling-time + contact-tracing | Dated genomes + epi data | Posterior WIWS + unsampled intermediates + R_e | Integrates epi data explicitly; modular likelihood | ~100-200 cases practical limit; assumes one infection event per case (no within-host populations) |
| TransPhylo (Didelot 2017 *Mol Biol Evol* 34:997) | Coalescent within-host + birth-death between-host; colours a dated tree | Time-scaled tree + sampling dates | Posterior transmission tree + R_t + unsampled cases | Works from a tree, not raw genomes; scales to ~1000 tips; explicit within-host coalescence | Sensitive to within-host effective population size prior; requires good dated phylogeny |
| phybreak (Klinkenberg 2017 *PLoS Comput Biol* 13:e1005495) | Joint phylogeny + transmission inference via MCMC | Dated genomes | Posterior transmission tree | Proper within-host handling; fast for small outbreaks | <=100 cases; less benchmarked than outbreaker2/TransPhylo |
| BadTrIP (De Maio 2018 *PLoS Comput Biol* 14:e1006117) | Bayesian; explicit handling of multi-strain infections | Dated genomes | Posterior transmission tree with strain-level resolution | Handles within-host diversity / mixed infections (TB, HIV) | Slow; specialist tool |
| SCOTTI (De Maio 2016 *PLoS Comput Biol* 12:e1005130) | Structured-coalescent transmission inference (BEAST 2 package) | Dated genomes | Posterior transmission tree under structured coalescent | Sampling-aware; correctly models unsampled intermediates | Computationally heavy; specialist setup |
| BEASTLIER (Hall 2015 *PLoS Comput Biol* 11:e1004613) | Joint phylogeny + transmission partitioning | Dated genomes; ideally with multiple isolates per host | Posterior transmission tree with within-host partition | Postdoc-grade identifiability with within-host samples | Single-isolate-per-host data is under-identified |
| transcluster (Stimson 2019 *Mol Biol Evol* 36:587) | Per-pair posterior probability under SNP + time prior | Dated genomes | Per-pair cluster membership probability | Probabilistic; pathogen-tuned priors | Pair-level only; no full transmission tree |
| Sobel Leonard 2017 *J Virol* 91:e00171-17 beta-binomial bottleneck | Estimate transmission bottleneck size from donor-recipient deep sequencing | Donor + recipient deep-sequence allele frequencies | Bottleneck Nb posterior | Estimates an otherwise unobservable quantity | Requires deep-sequenced donor-recipient pairs |
| islandR / Bayesian source attribution (Mather 2013 *Science* 341:1514) | Bayesian per-population allele-frequency model | Reference collections per host source + query genome | Per-source posterior probability | Standard in Salmonella / Campylobacter food-safety surveillance | Source-attribution circularity: trained-on-distribution reproduces that distribution |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why wrong choices fail |
|----------|----------------------|------------------------|
| "Is this even an outbreak?" routine surveillance triage | `snp-dists` after Gubbins; pathogen-tuned threshold (Walker 2013 for TB, Eyre 2013 for C. diff, EFSA cgMLST <=5 for Salmonella); cross-check cgMLST distance | Universal SNP threshold across pathogens (10x variation) |
| Densely sampled outbreak with contact-tracing data | outbreaker2 with `ctd` contact matrix + generation-time prior + sampling-time prior | TransPhylo without epi data (loses information from contacts); naive SNP threshold (ignores within-host diversity) |
| Sparsely sampled, longer-time-scale outbreak | TransPhylo on a BactDating-derived dated tree | outbreaker2 (sampling-completeness assumption broken); SNP threshold inflates clusters with unsampled intermediates |
| TB outbreak with possible reactivation | TransPhylo + transcluster with TB-tuned priors; long within-host coalescent matters | SNP cutoff insufficient -- reactivation can have 0 SNPs from years-old strains |
| Hospital outbreak with possible mixed infection | BadTrIP / SCOTTI | Consensus-only methods (SNP distance, outbreaker2) ambiguous on mixed-strain |
| Multi-site outbreak with import suspected | TransPhylo + MASCOT-derived migration; source-attribution as separate analysis | Source attribution needs phylogeographic component beyond TransPhylo alone |
| Food-vehicle / environmental source attribution | islandR / Bayesian source attribution (Mather 2013 framework); manual cluster + phylogeographic plot | Naive phylogenetic placement loses the per-source priors |
| Sub-sampled outbreak (<50% cases sequenced) | outbreaker2 (handles unsampled cases explicitly with `pi` sampling parameter) | Raw SNP cutoff -- unsampled intermediates break SNP-distance reasoning |
| Recombining pathogen (S. pneumo, E. coli STEC, K. pneumoniae) | Gubbins / ClonalFrameML mask FIRST; then any of the above | Recombination inflates apparent SNP distance and creates false convergent transmission inference |
| HIV cluster definition | HIV-TRACE 1.5% for subtype B (US-CDC standard); reconsider for non-B subtypes | Applying 1.5% threshold globally without subtype caveat |
| Estimate transmission bottleneck | Sobel Leonard 2017 beta-binomial on deep-sequenced donor-recipient pairs | Consensus-only sequences cannot quantify bottleneck size |

Methodology evolves; before any high-stakes who-infected-whom claim, web-search "outbreak transmission inference benchmark <pathogen> 2025" for current best practice.

## outbreaker2 With Contact Data

**Goal:** Infer who-infected-whom posterior for a densely sampled outbreak with epi metadata, jointly estimating generation interval and unsampled-case proportion.

**Approach:** Build `outbreaker_data` with sampling dates, DNA alignment, generation-time density `w_dens`, sampling-time density `f_dens`, and contact-tracing matrix `ctd`; configure MCMC via `create_config(n_iter=N)`; run; summarise posterior over WIWS.

```r
library(outbreaker2)
library(ape)

dna <- read.dna('alignment.fasta', format='fasta')
dates <- read.csv('sampling_dates.csv')
ctd_matrix <- as.matrix(read.csv('contact_matrix.csv', row.names=1))

w_dens <- dgamma(1:30, shape=2.5, scale=2)  # generation time prior
f_dens <- dgamma(1:30, shape=2, scale=3)    # sampling-time prior

data <- outbreaker_data(dates=dates$collection_date, dna=dna,
                        w_dens=w_dens, f_dens=f_dens, ctd=ctd_matrix)

cfg <- create_config(n_iter=1e6, sample_every=200, find_import=TRUE)

res <- outbreaker(data=data, config=cfg)
summary(res)
```

`w_dens` is the generation-time distribution (time from infection of A to infection of B) -- NOT the serial interval (time between symptom onsets); using one in place of the other biases inference. Britton & Scalia Tomba *J R Soc Interface* 16:20180670 (2019) formalised the bias for emerging epidemics; for SARS-CoV-2 with substantial pre-symptomatic transmission (Ali 2020 *Science* 369:1106), the serial interval shortened from 7.8 to 2.2 days under NPI, and naive SI-based inference was biased.

## TransPhylo From a Dated Tree

**Goal:** Infer transmission tree posterior from a time-scaled phylogeny when raw genomes are not directly usable or when the outbreak is too large for outbreaker2 (>200 cases).

**Approach:** Time-scale the tree first (BactDating after Gubbins for bacteria; BEAST or TreeTime for viruses); convert to TransPhylo `ptree` with `ptreeFromPhylo`; run `inferTTree` with generation-time prior and within-host effective population size prior; summarise via `medTTree` (medoid transmission tree) and posterior probabilities per WIWS pair.

```r
library(TransPhylo)
library(ape)

tree <- read.nexus('dated_tree.nexus')
date_last_sample <- 2024.95

ptree <- ptreeFromPhylo(tree, dateLastSample=date_last_sample)

w.shape <- 1.3
w.scale <- 10
ws.shape <- 1.1
ws.scale <- 7
neg <- 0.5

res <- inferTTree(ptree, mcmcIterations=1e5,
                  w.shape=w.shape, w.scale=w.scale,
                  ws.shape=ws.shape, ws.scale=ws.scale,
                  startNeg=neg, dateT=date_last_sample + 0.1)

med_tree <- medTTree(res)
pairs <- extractTTree(med_tree)$ttree
```

`w.*` is the generation-time Gamma prior; `ws.*` is the sampling-time Gamma prior. Both must reflect the pathogen's biology (e.g., TB w.scale = months; SARS-CoV-2 w.scale = days). Wrong priors silently bias the transmission-tree posterior.

## SNP-Cluster Definition With Pathogen-Specific Thresholds

**Goal:** Define outbreak clusters from a recombination-masked core-SNP alignment using the published pathogen-specific threshold, with the threshold's source population caveated.

**Approach:** Snippy -> snippy-core -> Gubbins on `core.full.aln` for bacteria -> snp-dists -> single-linkage clustering at the pathogen-specific threshold; cite Walker 2013 (TB), Eyre 2013 (C. diff), Coll 2017 (MRSA), Snitkin 2012 (Klebsiella) per organism; flag any extrapolation outside the threshold's validation population.

```bash
snippy-core --ref reference.fa --prefix core snippy_out/*
run_gubbins.py --prefix gubbins core.full.aln
snp-dists -c gubbins.filtered_polymorphic_sites.fasta > pairwise.csv
```

```python
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

dist = pd.read_csv('pairwise.csv', index_col=0)
condensed = dist.values[np.triu_indices(len(dist), k=1)]

THRESHOLD_TB = 12   # Walker 2013 Lancet Infect Dis 13:137 -- UK low-transmission
THRESHOLD_MRSA = 15  # Coll 2017 Clin Infect Dis 65:1781
THRESHOLD_CDIFF = 2  # Eyre 2013 NEJM 369:1195
THRESHOLD_KPNEUMO = 21  # Snitkin 2012 Sci Transl Med 4:148ra116

linkage_matrix = linkage(condensed, method='single')
clusters = fcluster(linkage_matrix, t=THRESHOLD_TB, criterion='distance')
```

## Per-Method Failure Modes

### Pairwise SNP threshold applied outside its validation population

**Trigger:** Walker 2013 UK 5/12-SNP TB threshold applied to Cape Town or Mumbai high-transmission settings.

**Mechanism:** Walker 2013 *Lancet Infect Dis* 13:137 calibrated the 5/12 SNP threshold on Oxfordshire community / household contact-traced data (low-transmission). In high-prevalence settings, clonal isolates linked through long-past common ancestors fall within the threshold without recent direct transmission.

**Symptom:** Country-level Mtb genomic-epi report shows 60-80% of cases in "transmission clusters", far exceeding clinical contact-tracing rates.

**Fix:** Cite the threshold's source population; for high-prevalence settings, derive a local threshold from epidemiologically-anchored case pairs in the local cohort rather than importing a UK-low-transmission cutoff. For transmission-direction claims, supplement with TransPhylo / outbreaker2.

### Direction of transmission asserted from pairwise SNP distance alone

**Trigger:** Outbreak report concluding "A -> B" because A has earlier sampling date and 3 SNPs from B.

**Mechanism:** A 3-SNP pairwise difference is consistent with A->B, B->A, Unknown->both, or A->Unknown->B. Worby, Lipsitch & Hanage 2014 *PLoS Comput Biol* 10:e1003549 formalised the irreducible uncertainty. Earlier sampling date does not establish earlier infection date because of within-host evolution and asymptomatic carriage.

**Symptom:** Outbreak conclusions claim directionality without within-host data or contact tracing; reviewers from the Didelot / Worby groups push back.

**Fix:** Use "transmission consistent with genomics" not "transmission demonstrated". For direction claims, require within-host samples (BEASTLIER), contact-tracing data (outbreaker2 with `ctd`), or both. Cite Worby 2014 as the upper bound on what SNP distance can establish.

### Unsampled intermediates collapsed into A->B direct links

**Trigger:** Outbreak with <50% sequencing coverage; transmission inference assumes all cases sampled.

**Mechanism:** When sampling is incomplete, inferred A->B "direct" transmissions are routinely A->Unknown->B chains. This systematically inflates inferred R_e (longer chains compressed), underestimates generation interval, and biases topology toward bushy trees.

**Symptom:** Inferred R_e is implausibly high (each "tip" appears to spawn extra children once unsampled intermediates collapse into apparent direct links); generation interval estimate is implausibly short; topology appears bushier than expected.

**Fix:** Use outbreaker2 with explicit `pi` (sampling proportion) parameter, or TransPhylo / SCOTTI which model unsampled intermediates explicitly. Cite the unsampled-intermediates caveat in every transmission-inference report.

### Narrow transmission bottleneck makes consensus-only inference WORSE than coalescent intuition predicts

**Trigger:** Consensus-genome transmission-pair inference for a pathogen with documented narrow bottleneck (influenza 1-2 virions per McCrone 2018 *eLife* 7:e35962; SARS-CoV-2 <10 virions per Lythgoe 2021 *Science* 372:eabg0821).

**Mechanism:** When the transmission bottleneck is narrow, donor and recipient consensus genomes are near-identical *by default* -- the bottleneck strips most within-host diversity. Near-identity therefore does NOT discriminate direct transmission from infection by an unsampled intermediate or from a shared common source. Naive coalescent intuition predicts that "more transmissions = more divergence"; the opposite is true under a narrow bottleneck.

**Symptom:** Most pairs in a dense outbreak appear identical or 1 SNP apart; SNP-distance-based cluster definitions become uninformative; transmission-direction claims based on consensus difference are unfalsifiable.

**Fix:** For narrow-bottleneck pathogens, supplement consensus-based methods with deep within-host variant calling (lofreq / deepSNV / VarScan2 at MAF >= 1%) on donor-recipient pairs; estimate bottleneck size explicitly via Sobel Leonard 2017 *J Virol* 91:e00171-17 beta-binomial estimator; report transmission claims as "consistent with" rather than "demonstrated by" consensus identity. Pair-level resolution requires within-host data; without it, claim only cluster membership, not direction.

### Generation interval and serial interval used interchangeably

**Trigger:** outbreaker2 / EpiNow2 / similar tools fed the serial-interval distribution (`w_dens` set from symptom-to-symptom data) when the model wants generation-interval (infection-to-infection).

**Mechanism:** Generation interval = time from infection of A to infection of B; serial interval = time from symptom onset of A to symptom onset of B. They differ when incubation periods vary or pre-symptomatic transmission is substantial. Britton & Scalia Tomba 2019 *J R Soc Interface* 16:20180670 formalised the bias for emerging epidemics; Ali 2020 *Science* 369:1106 showed for SARS-CoV-2 the SI shortened from 7.8 to 2.2 days under NPI.

**Symptom:** Inferred R_e is biased; comparison to case-based R_t (also often SI-based) shows compounding bias.

**Fix:** Document which distribution `w_dens` actually encodes. For SARS-CoV-2 with substantial pre-symptomatic transmission, generation interval is ~5 days in the ancestral-strain literature; serial interval was ~4-5 days early but shortened to 2-3 under NPI. Cite Britton 2019.

### HIV-TRACE 1.5% threshold applied to non-subtype-B HIV

**Trigger:** HIV-TRACE run on subtype C sequences from southern Africa with the default 1.5% TN93 threshold.

**Mechanism:** Kosakovsky Pond et al 2018 *Mol Biol Evol* 35:1812 documented HIV-TRACE methodology; the 1.5% threshold is the US-CDC default tuned for subtype B in MSM cohorts. Subtype C in southern Africa has higher diversity per unit time and more recent epidemics; the 1.5% threshold under-clusters there.

**Symptom:** Cluster definitions in southern African subtype C HIV surveillance under-detect transmission; comparison to US surveillance literature shows incompatible cluster sizes.

**Fix:** Tune threshold for the local subtype and population; cite the local validation. UKHSA / ECDC use different thresholds; document which.

### Source attribution circularity

**Trigger:** Bayesian source attribution model (Mather 2013 *Science* 341:1514 framework) trained on a reference collection that over-represents one host population.

**Mechanism:** Source-attribution models reproduce the host-distribution of their training data unless explicitly corrected. If 80% of training isolates are from cattle, the model will tend to attribute new isolates to cattle even when the true source is poultry.

**Symptom:** Source attribution reproduces the sampling intensity of the reference collection; conclusions are circular.

**Fix:** Weight by inverse sampling intensity per source category; use rarefied reference collections; report attribution alongside the reference-collection composition as a caveat.

### Primer-scheme dropout misread as real divergence

**Trigger:** SARS-CoV-2 outbreak comparison across samples sequenced with different ARTIC primer schemes (V3 / V4 / V4.1 / V5.3.2); "differences" concentrated in one amplicon are interpreted as real SNPs.

**Mechanism:** ARTIC primer dropouts produce N's or reference-derived consensus calls in failed amplicons (Itokawa 2020 *PLoS ONE* 15:e0239403); these LOOK LIKE deletions or reference matches in downstream analysis but are missing data. Cross-scheme comparison without masking failed amplicons produces spurious transmission differences.

**Symptom:** Cluster definitions differ implausibly between ARTIC-V3 and ARTIC-V4.1 samples; "differences" cluster in known dropout amplicons (V4.1 amplicons 64, 76, 88-90).

**Fix:** Mask failed amplicons per sample (samtools depth + per-amplicon coverage); document primer scheme version per isolate; for transmission inference, exclude positions in any sample's dropout regions.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| outbreaker2 and TransPhylo disagree on WIWS | Different sampling-completeness assumptions; outbreaker2 expects ~dense sampling, TransPhylo handles sparse | Pick the method whose assumption matches the data; cite the choice |
| SNP threshold cluster and outbreaker2 cluster differ | SNP threshold ignores temporal data and contacts | Trust outbreaker2 (integrates more evidence); SNP cluster is triage only |
| Two consecutive Pangolin versions give different lineage for a "transmission pair" | Lineage definitions revised | Re-run both samples against a single Pango / pangolin-data version |
| TB cluster definition flips between 5 and 12 SNP threshold | Walker 2013 ambiguous range | Run TransPhylo for transmission-direction posterior; report SNP-distance with cluster picker certainty |
| HIV cluster differs between HIV-TRACE 1.5% and 2.0% | Threshold sensitivity at boundary | Subtype-specific calibration; cite the chosen threshold's validation |
| Source attribution differs between islandR runs with different reference panels | Sampling-intensity bias | Re-run with rarefied or inverse-weighted reference; report multiple scenarios |

## Quantitative Thresholds

| Pathogen | "Outbreak cluster" threshold | Source / rationale |
|----------|------------------------------|--------------------|
| *Mycobacterium tuberculosis* (whole-genome core SNP) | <=12 SNPs (likely transmission); <=5 SNPs (recent transmission) | Walker 2013 *Lancet Infect Dis* 13:137 (UK low-transmission setting) |
| *Staphylococcus aureus* (core genome) | <=15 SNPs (within hospital outbreak); <=40 SNPs (broader temporal cluster) | Coll 2017 *Clin Infect Dis* 65:1781 |
| *Klebsiella pneumoniae* (KPC outbreak) | <=21 SNPs | Snitkin 2012 *Sci Transl Med* 4:148ra116 |
| *Salmonella enterica* (cgMLST EnteroBase) | <=5 allelic differences (cluster); <=7 (extended cluster) | EnteroBase / EFSA harmonised |
| *Listeria monocytogenes* (PulseNet cgMLST) | <=4 allelic differences | PulseNet protocol convention |
| *E. coli* (cgMLST, EnteroBase) | <=10 allelic differences (STEC outbreak) | EnteroBase convention |
| *Neisseria gonorrhoeae* | <=25 core SNPs (transmission) | UKHSA STI framework |
| *Clostridioides difficile* (core SNP, recombination-masked) | <=2 SNPs (likely direct); <=10 (plausible within 6 months) | Eyre 2013 *NEJM* 369:1195 |
| SARS-CoV-2 (whole-genome) | No fixed cutoff; 0-2 SNPs + epi link + sampling window | Lythgoe 2021 *Science* 372:eabg0821 |
| HIV-1 subtype B (TN93 distance) | 1.5% genetic distance (HIV-TRACE default; US-CDC standard) | Kosakovsky Pond 2018 *Mol Biol Evol* 35:1812 |
| Mpox clade IIb | <=2 SNPs cluster threshold; APOBEC3 editing inflates apparent distance | Mpox 2022 outbreak APOBEC3-editing literature |
| Transmission bottleneck -- influenza | ~1-2 virions (narrow) | McCrone 2018 *eLife* 7:e35962 |
| Transmission bottleneck -- SARS-CoV-2 | <10 virions (tight) | Lythgoe 2021 *Science* 372:eabg0821 |
| Generation interval -- SARS-CoV-2 ancestral | ~5 days | SARS-CoV-2 ancestral-strain literature |

CRITICAL: a number from one pathogen does NOT transfer to another. Always cite the source population.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| outbreaker2 rejects `iters` arg | Iterations set via `n_iter` in the `config` object | `create_config(n_iter=N)` |
| TransPhylo MCMC fails to converge | Within-host Ne prior misspecified; bad input tree | Tune `startNeg`; verify tree dating quality |
| Cluster definition flips between linkage methods | Single-linkage vs complete-linkage on borderline pairs | Document; sensitivity analysis |
| outbreaker2 estimates implausible R_e | Sampling proportion mis-specified | Set `pi` based on epi knowledge or estimate within outbreaker2 |
| Transmission inferred between two distant lineages | Recombination unmasked | Run Gubbins on `core.full.aln` first |
| HIV-TRACE clusters incompatible across labs | Different subtype calibration | Document subtype; use locally validated threshold |
| Source attribution always pointing at one host | Reference-collection bias | Re-weight or rarify reference panel |
| `snp-dists -t` rejected | `-t` flag doesn't exist; default IS tab; `-c` for CSV | Use `-c` for CSV; default for TSV |
| Snippy outputs disagree across samples | Different reference; reference mismatch silently shifts SNP coordinates | Always document reference; use same reference cross-lab |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "What SNP threshold and on what population?" | Cite Walker 2013 / Eyre 2013 / Coll 2017 per pathogen; caveat the population if extrapolating |
| "Were unsampled intermediates handled?" | outbreaker2 `pi` parameter or TransPhylo / SCOTTI explicit modelling; never a raw SNP-distance method on sub-sampled data |
| "Direction of transmission inference?" | Within-host samples + contact tracing required for direction; otherwise "consistent with" phrasing |
| "Generation interval vs serial interval?" | Documented `w_dens` source; cite Britton 2019 if SI used as approximation for GI |
| "Why TransPhylo / outbreaker2 / phybreak?" | Decision tree based on sampling completeness, dataset size, contact-tracing availability |
| "Was within-host diversity considered?" | TransPhylo's within-host coalescent OR BadTrIP for mixed-strain; bottleneck size from Sobel Leonard 2017 if relevant |
| "HIV-TRACE 1.5% threshold outside subtype B?" | Acknowledged US-CDC subtype B origin; either use locally validated threshold or document caveat |
| "Source attribution sampling-intensity bias?" | Re-weighted reference collection or rarified; cite Mather 2013 limitation |
| "Was forward simulation run as a sanity check?" | SLiM / FAVITES / SEEDY if claims are high-stakes; routinely under-done in published transmission inference |

## References

- Worby CJ, Lipsitch M, Hanage WP (2014) Within-host bacterial diversity hinders accurate reconstruction of transmission networks from genomic distance data. *PLoS Comput Biol* 10(3):e1003549. doi:10.1371/journal.pcbi.1003549
- Campbell F, Didelot X, Fitzjohn R, Ferguson N, Cori A, Jombart T (2018) outbreaker2: a modular platform for outbreak reconstruction. *BMC Bioinformatics* 19(Suppl 11):363. doi:10.1186/s12859-018-2330-z
- Didelot X, Fraser C, Gardy J, Colijn C (2017) Genomic infectious disease epidemiology in partially sampled and ongoing outbreaks. *Mol Biol Evol* 34(4):997-1007. doi:10.1093/molbev/msw275
- Klinkenberg D, Backer JA, Didelot X, Colijn C, Wallinga J (2017) Simultaneous inference of phylogenetic and transmission trees in infectious disease outbreaks. *PLoS Comput Biol* 13(5):e1005495. doi:10.1371/journal.pcbi.1005495
- De Maio N, Worby CJ, Wilson DJ, Stoesser N (2018) Bayesian reconstruction of transmission within outbreaks using genomic variants. *PLoS Comput Biol* 14(4):e1006117. doi:10.1371/journal.pcbi.1006117
- De Maio N, Wu CH, Wilson DJ (2016) SCOTTI: efficient reconstruction of transmission within outbreaks with the structured coalescent. *PLoS Comput Biol* 12(9):e1005130. doi:10.1371/journal.pcbi.1005130
- Hall M, Woolhouse M, Rambaut A (2015) Epidemic reconstruction in a phylogenetics framework: transmission trees as partitions of the node set. *PLoS Comput Biol* 11(12):e1004613. doi:10.1371/journal.pcbi.1004613
- Stimson J, Gardy J, Mathema B et al (2019) Beyond the SNP threshold: identifying outbreak clusters using inferred transmissions. *Mol Biol Evol* 36(3):587-603. doi:10.1093/molbev/msy242
- Walker TM, Ip CLC, Harrell RH et al (2013) Whole-genome sequencing to delineate Mycobacterium tuberculosis outbreaks: a retrospective observational study. *Lancet Infect Dis* 13(2):137-146. doi:10.1016/S1473-3099(12)70277-3
- Coll F, Harrison EM, Toleman MS et al (2017) Longitudinal genomic surveillance of MRSA in the UK reveals transmission patterns in hospitals and the community. *Clin Infect Dis* 65(11):1781-1789. doi:10.1093/cid/cix645
- Eyre DW, Cule ML, Wilson DJ et al (2013) Diverse sources of C. difficile infection identified on whole-genome sequencing. *N Engl J Med* 369(13):1195-1205. doi:10.1056/NEJMoa1216064
- Snitkin ES, Zelazny AM, Thomas PJ et al (2012) Tracking a hospital outbreak of carbapenem-resistant Klebsiella pneumoniae with whole-genome sequencing. *Sci Transl Med* 4(148):148ra116. doi:10.1126/scitranslmed.3004129
- Lythgoe KA, Hall M, Ferretti L et al (2021) SARS-CoV-2 within-host diversity and transmission. *Science* 372(6539):eabg0821. doi:10.1126/science.abg0821
- McCrone JT, Woods RJ, Martin ET et al (2018) Stochastic processes constrain the within and between host evolution of influenza virus. *eLife* 7:e35962. doi:10.7554/eLife.35962
- Sobel Leonard A, Weissman DB, Greenbaum B, Ghedin E, Koelle K (2017) Transmission bottleneck size estimation from pathogen deep-sequencing data, with an application to human influenza A virus. *J Virol* 91(14):e00171-17. doi:10.1128/JVI.00171-17
- Britton T, Scalia Tomba G (2019) Estimation in emerging epidemics: biases and remedies. *J R Soc Interface* 16(150):20180670. doi:10.1098/rsif.2018.0670
- Ali ST, Wang L, Lau EHY et al (2020) Serial interval of SARS-CoV-2 was shortened over time by nonpharmaceutical interventions. *Science* 369(6507):1106-1109. doi:10.1126/science.abc9004
- Kosakovsky Pond SL, Weaver S, Leigh Brown AJ, Wertheim JO (2018) HIV-TRACE (TRAnsmission Cluster Engine): A tool for large-scale molecular epidemiology of HIV-1 and other rapidly evolving pathogens. *Mol Biol Evol* 35(7):1812-1819. doi:10.1093/molbev/msy016
- Mather AE, Reid SWJ, Maskell DJ et al (2013) Distinguishable epidemics of multidrug-resistant Salmonella Typhimurium DT104 in different hosts. *Science* 341(6153):1514-1517. doi:10.1126/science.1240578
- Itokawa K, Sekizuka T, Hashino M, Tanaka R, Kuroda M (2020) Disentangling primer interactions improves SARS-CoV-2 genome sequencing by multiplex tiling PCR. *PLoS ONE* 15(9):e0239403. doi:10.1371/journal.pone.0239403

## Related Skills

- pathogen-typing - SNP-cluster / cgMLST cluster definition feeds transmission inference
- phylodynamics - Time-scaled tree from BactDating / BEAST / TreeTime feeds TransPhylo
- amr-surveillance - Resistant-clone outbreak inference combines AMR + transmission
- variant-surveillance - Lineage assignment cross-checks transmission cluster boundaries
- phylogenetics/divergence-dating - Calibrated trees for non-pathogen contexts
- phylogenetics/bayesian-inference - BEAST mechanics beyond outbreak phylodynamics
- comparative-genomics/whole-genome-alignment - Core-genome alignment for SNP-typing
- variant-calling/vcf-basics - Per-isolate variant calls for SNP-typing
- variant-calling/variant-calling - SNP calling that feeds snp-dists
- read-alignment/bwa-alignment - Read mapping upstream
- data-visualization/network-visualization - Transmission tree visualisation
- workflows/somatic-variant-pipeline - End-to-end orchestration patterns
