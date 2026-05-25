---
name: bio-ecological-genomics-species-delimitation
description: Delimits putative species boundaries from molecular data within the de Queiroz 2007 unified-lineage framework using ASAP (Puillandre 2021 successor to ABGD), mPTP C++ (Kapli 2017 successor to bPTP; bPTP is Python NOT R), GMYC single/multi-threshold (Pons 2006; Fujisawa 2013), multilocus BPP v4 with prior calibration from data (NOT defaults; Yang 2015), SNAPP + BFD* for SNP delimitation, DELINEATE (Sukumaran 2021) speciation-process modeling to address Sukumaran & Knowles 2017 PNAS critique that MSC delimits structure not species, integrative-taxonomy congruence (Padial 2010; Carstens 2013), Dsuite for introgression testing before sister claims (Malinsky 2021), and Meyer & Paulay 2005 barcoding-gap-absence caveat. Use when delineating species from DNA barcoding data, resolving cryptic complexes, choosing among ASAP/mPTP/BPP/DELINEATE, calibrating BPP priors, distinguishing introgression from ILS, or applying the Sukumaran-Knowles oversplitting correction.
tool_type: mixed
primary_tool: ASAP
---

## Version Compatibility

Reference examples tested with: ASAP (current CLI), mPTP (current C++), bPTP (Python from github.com/zhangjiajie/PTP), splits 1.0+ for GMYC, BPP 4.7+, SNAPP/BEAST 2.7+, DELINEATE (current Python), Dsuite 0.5+, ape 5.7+, fossil 0.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Species Delimitation

**"Delineate species boundaries from my DNA barcoding or genomic data"** -> Apply distance-based (ASAP) and tree-based (mPTP) methods for primary delimitation, multilocus coalescent (BPP) for genomic confirmation, DELINEATE for speciation-process modeling to address MSC oversplitting, integrative-taxonomy validation across multiple lines of evidence, and Dsuite to test introgression before claiming sister relationships.
- CLI: `asap` (downloadable C binary) OR web (https://bioinfo.mnhn.fr/abi/public/asap/) for single-locus primary delimitation
- CLI: `mptp` (https://github.com/Pas-Kapli/mptp) for tree-based multi-rate PTP (faster than bPTP)
- Python: `bPTP` (https://github.com/zhangjiajie/PTP) for legacy tree-based delimitation — NOTE: Python not R
- R: `splits::gmyc()` for GMYC on ultrametric trees
- CLI: `bpp` (https://bpp.github.io) for multilocus MSC delimitation; control-file driven
- CLI: `Dsuite` for D-statistics / f4-ratio / f-branch introgression testing

## The Single Most Important Modern Insight -- MSC Methods Delimit Structure, NOT Species

Sukumaran & Knowles 2017 *PNAS* 114(7):1607-1612 established that BPP, BFD*, and other multispecies-coalescent methods **delimit genetic structure, not species**. The mathematical reason: the MSC model treats every panmictic population as a "species" in the parameterization. Applied to data with substructure (isolated demes within a species), MSC methods partition that structure into "species" — leading to systematic oversplitting in published literature, especially in lizards, frogs, geckos, and insects.

**Modern best practice (post-2021):**
1. Run MSC methods (BPP, BFD*) as PRELIMINARY enumeration of candidate population lineages
2. Apply DELINEATE (Sukumaran 2021 *PLoS Comput Biol* 17:e1008924) to test which lineages have completed speciation vs are in the speciation process
3. Validate via integrative taxonomy (Padial 2010; Carstens 2013) — morphological, ecological, geographic congruence

A second cornerstone: **the barcoding gap is OFTEN absent in real data** (Meyer & Paulay 2005 *PLoS Biol* 3:e422). ABGD and ASAP are theoretically grounded in the gap; when the gap is absent, both fail gracefully but their output is unreliable. Always inspect the pairwise distance histogram.

A third: **Dsuite-checked introgression must precede sister-species claims.** D-statistics, f4-ratio, and f-branch (Malinsky 2021) distinguish admixture from incomplete lineage sorting; without this check, ILS-driven discordance is misinterpreted as gene flow.

## Algorithmic Taxonomy

| Method | Input | Approach | Strength | Fails when |
|--------|-------|----------|----------|------------|
| ASAP (Puillandre 2021) | Aligned sequences | Hierarchical-clustering distance partitioning with new scoring | Fast; modern; web + CLI | Single locus; works only if barcoding gap exists |
| ABGD (Puillandre 2012) | Aligned sequences | Distance gap-detection across priors | Legacy; superseded by ASAP | Less robust when gap is weak |
| GMYC (Pons 2006) | Ultrametric tree | Yule-coalescent transition threshold | Theoretical foundation | Requires time-calibrated tree |
| Multi-threshold GMYC (Fujisawa & Barraclough 2013) | Ultrametric tree | Per-lineage transition thresholds | Heterogeneous rates | Same; not always more powerful |
| bPTP (Zhang 2013) | Rooted phylogeny | Bayesian PTP MCMC | Posterior support per partition | Single intraspecific rate assumption; slower than mPTP |
| mPTP (Kapli 2017) | Rooted phylogeny | Multi-rate PTP, C++ | Faster (5+ orders); per-species intraspecific rates | Heterogeneity-friendly; for shallow phylogenies bPTP may be more powerful |
| BPP A10/A11 (Yang & Rannala) | Multi-locus alignments | Bayesian multispecies coalescent | Rigorous, multilocus | Computationally heavy; oversplits per Sukumaran-Knowles |
| SNAPP (Bryant 2012) | Biallelic SNP data | Coalescent species-tree bypassing gene trees | High-power genomic | Cubic complexity; not for very large datasets |
| BFD* (Leache 2014) | SNP data via SNAPP | Bayes-factor delimitation | Genomic-scale species hypothesis testing | Same oversplitting risk |
| DELINEATE (Sukumaran 2021) | BPP output + guide tree | Speciation-process modeling | Distinguishes structure from species | Requires preliminary delimitation as input |
| Dsuite (Malinsky 2021) | VCF + populations | D-statistic, f4-ratio, f-branch | Tests introgression vs ILS | Tree must accurately reflect history |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Primary species hypothesis from single locus | ASAP first; cross-check with mPTP | ASAP is the modern successor to ABGD; faster and better-scoring |
| Tree-based delimitation, heterogeneous intraspecific rates | mPTP | Per-species rate model; 5+ orders faster than bPTP |
| Tree-based delimitation, shallow phylogeny | bPTP (single-rate) or mPTP | When rates are similar, single-rate model may be more powerful |
| Time-calibrated tree from BEAST / chronos | GMYC (single or multi-threshold) | Requires ultrametric; both threshold variants worth comparing |
| Multilocus genomic delimitation | BPP A10/A11 + DELINEATE | BPP alone oversplits per Sukumaran-Knowles 2017 |
| SNP-based species hypothesis test | SNAPP + BFD* | Coalescent species-tree from biallelic SNPs |
| Confirm BPP result is not oversplitting | DELINEATE (Sukumaran 2021) | Speciation-process modeling distinguishes structure from species |
| Test introgression before sister-species claim | Dsuite (D, f4-ratio, f-branch) | Distinguishes admixture from ILS |
| Cryptic species complex | Integrative taxonomy: genetic + morphological + ecological congruence | Single-method conclusion unreliable |
| ABGD/ASAP returns no clear partition | Inspect pairwise distance histogram; do not force partition | Barcoding gap may be absent (Meyer & Paulay 2005) |
| Recently-diverged populations (Ne*t << 1) | Caution; MSC methods very prone to oversplit | Incomplete lineage sorting indistinguishable from shallow structure |
| Conservation management unit definition | Moritz 1994 ESU/MU framework, NOT MSC delimitation | Sukumaran-Knowles caveat applies |

## ASAP — The Modern Successor to ABGD

**Goal:** Primary species hypothesis from a single-locus alignment via hierarchical clustering of pairwise distances, scored by a new asap-score (Puillandre 2021).

**Approach:** Run ASAP via the web interface (https://bioinfo.mnhn.fr/abi/public/asap/) OR the downloadable C binary. ASAP ranks candidate partitions by asap-score (lower = better); always inspect the top 5-10 partitions, not just the best, and check for large score gaps that signal robust partitioning. Output: ranked species partitions with asap-score, p-value, and threshold distance per partition.

```bash
# ASAP CLI (downloadable C binary)
# Substitution model: K2P standard for COI; p-distance for very closely related; JC for general
./asap -d K2P -o asap_results/ aligned_sequences.fasta

# Inspect output
# - asap_output.csv: ranked partitions by asap-score (lower = better)
# - Always look at top 5-10 partitions; large score gaps = robust partitioning
# - The "best" partition is the top-ranked; report secondary partitions too
```

```python
# Python ASAP-style analysis for inspection (NOT a substitute for ASAP)
# Useful for understanding the distance landscape
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
import numpy as np
import matplotlib.pyplot as plt

alignment = AlignIO.read('aligned_sequences.fasta', 'fasta')
calc = DistanceCalculator('identity')  # use K2P-friendly model in practice
dm = calc.get_distance(alignment)

names = [r.id for r in alignment]
n = len(names)
dist_array = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        dist_array[i][j] = dm[names[i], names[j]]

# Inspect distance histogram for barcoding gap (often ABSENT per Meyer & Paulay 2005)
upper = dist_array[np.triu_indices(n, k=1)]
plt.hist(upper, bins=50)
plt.xlabel('Pairwise distance')
plt.ylabel('Frequency')
plt.title('Distance histogram — look for bimodality (gap)')
plt.savefig('barcoding_gap_histogram.pdf')

# Hierarchical clustering at threshold scan
condensed = squareform(dist_array)
Z = linkage(condensed, method='average')
thresholds = np.arange(0.01, 0.10, 0.005)
for t in thresholds:
    clusters = fcluster(Z, t=t, criterion='distance')
    print(f'Threshold {t:.3f}: {len(set(clusters))} groups')
```

## mPTP — Multi-Rate PTP for Heterogeneous Intraspecific Rates

**Goal:** Delimit species from a phylogeny by detecting the transition from speciation to coalescent branching rates, allowing per-species intraspecific rates (Kapli 2017).

**Approach:** Run mPTP on a rooted ML tree (e.g., from RAxML or IQ-TREE). mPTP's multi-rate model is more flexible than bPTP's single-rate; for shallow phylogenies with similar rates, bPTP may be competitive. mPTP is C++ (5+ orders of magnitude faster than bPTP), available from https://github.com/Pas-Kapli/mptp.

```bash
# mPTP installation: build from C source
# git clone https://github.com/Pas-Kapli/mptp
# cd mptp && ./autogen.sh && ./configure && make

# Run mPTP ML inference
mptp --ml --tree_file rooted_tree.nwk --output_file mptp_ml

# MCMC for posterior support
mptp --mcmc 1000000 --mcmc_sample 1000 --mcmc_burnin 100000 \
     --tree_file rooted_tree.nwk --output_file mptp_mcmc

# Output: species partition with branch annotations
# > 0.95 posterior support: strong; 0.80-0.95: moderate; < 0.80: uncertain
```

## bPTP — Legacy Tree-Based Delimitation (Python, NOT R)

**Goal:** Bayesian Poisson Tree Processes delimitation with single intraspecific rate (Zhang 2013).

**Approach:** Run bPTP via the Python package from https://github.com/zhangjiajie/PTP (NOTE: this is a Python tool, NOT an R package; `install.packages('PTP')` does not exist). For heterogeneous-rate data, prefer mPTP.

```bash
# bPTP Python install (NOT R)
pip install git+https://github.com/iTaxoTools/PTP-pyqt5

# Run bPTP MCMC
python -m PTP.PTP -t rooted_tree.nwk -o bptp_results \
    --type bayesian --ngen 100000 --burnin 0.1 --seed 42

# Output: posterior support per species partition (> 0.95 strong)
# bPTP over-splits when populations have strong geographic substructure
# Cross-check with ASAP and/or mPTP for consensus
```

## GMYC with Mandatory Ultrametric Tree

**Goal:** Detect the threshold on an ultrametric phylogeny where branching transitions from interspecific (Yule) to intraspecific (coalescent).

**Approach:** Convert ML tree to ultrametric with `ape::chronos` (penalized likelihood) or BEAST (Bayesian time-calibration), then run GMYC. Single-threshold vs multi-threshold; multi-threshold handles heterogeneous rates across the tree.

```r
library(splits)
library(ape)

tree <- read.tree('rooted_tree.nwk')

# Convert to ultrametric (REQUIRED for GMYC)
# chronos: penalized likelihood; lambda=1 mid-clock
# For rigorous time-calibration, use BEAST2 instead
ultrametric_tree <- chronos(tree, lambda = 1)
class(ultrametric_tree) <- 'phylo'
stopifnot(is.ultrametric(ultrametric_tree))

# Single-threshold GMYC
gmyc_single <- gmyc(ultrametric_tree, method = 'single')
summary(gmyc_single)
cat('GMYC species (single-threshold):', gmyc_single$entity[1], '\n')
cat('LR test p-value:', gmyc_single$p.value[1], '\n')

# Multi-threshold for heterogeneous rates
gmyc_multi <- gmyc(ultrametric_tree, method = 'multiple')
summary(gmyc_multi)

# Compare single vs multi-threshold; report both
species_single <- spec.list(gmyc_single)
species_multi <- spec.list(gmyc_multi)
```

## BPP A10/A11 — Multilocus MSC Delimitation with Calibrated Priors

**Goal:** Bayesian multilocus species delimitation under the multispecies coalescent (Yang & Rannala 2010; current Flouri 2018 BPP v4).

**Approach:** Configure BPP A10 (joint delimitation + species tree estimation) with priors on theta and tau ESTIMATED FROM DATA — NOT defaults. Yang 2015 *Curr Zool* 61:854-865 gives explicit guidance: `thetaprior = G(2, 2000)` for small/moderate data; tighter for large genomic data. Wrong priors dominate the analysis.

```text
# BPP control file for A10 analysis
# CRITICAL: priors below are CALIBRATED to data, not defaults
# Defaults are tuned for primate-mammal data; insects/fish/plants need adjustment

seed = 12345

seqfile = alignment.phy
Imapfile = imap.txt                   # map individuals -> putative species
outfile = bpp_results.txt
mcmcfile = bpp_mcmc.txt

speciesdelimitation = 1 0 2 0.5       # 1 = delimitation mode; rjMCMC parameters
speciestree = 1                       # 1 = estimate species tree
speciesmodelprior = 1                 # 1 = uniform prior on rooted trees

species&tree = 4  speciesA speciesB speciesC speciesD
                  10 8 12 6
                  ((speciesA, speciesB), (speciesC, speciesD));

# CALIBRATED theta prior (per-locus nucleotide diversity)
# G(2, 2000) -> mean = 0.001 (calibrate from observed per-locus pi)
# Yang 2015 Curr Zool 61:854 has taxon-specific recommendations
thetaprior = 2 2000

# CALIBRATED tau prior (divergence times)
# Calibrate from observed maximum genetic distance / mutation rate
tauprior = 2 2000

# MCMC
nsample = 100000
sampfreq = 2
burnin = 50000

# Run multiple independent chains; check convergence
```

```bash
# Run BPP with multiple independent chains
for i in 1 2 3 4; do
    bpp --cfile bpp_control.bpp --seed $((1000 * i)) --out bpp_run_$i.txt &
done
wait

# Compare chains; posterior probabilities must agree across chains
# PP > 0.95 for a delimitation: strong support
# PP 0.50-0.95: moderate; integrate other evidence
```

## DELINEATE — Addressing the Sukumaran-Knowles Oversplitting

**Goal:** Test whether candidate lineages from BPP represent fully-formed species OR incomplete-speciation structure (Sukumaran 2021).

**Approach:** Provide DELINEATE with a guide tree of population lineages and a hypothesis about which are species. The output: probability that each population lineage is a true species vs an incipient lineage. THE post-2021 best practice complement to BPP.

```bash
# DELINEATE (https://github.com/jsukumaran/delineate)
# Python tool
pip install delineate

# Configure delineate input
# - guide tree: rooted Newick with population lineages as tips
# - constraint file: which tips are confirmed species, which to test
delineate-estimate partitions \
    --tree-file population_lineage_tree.nwk \
    --config-file delineate_config.json \
    --output-prefix delineate_run

# Output: probability that each lineage is a species
# Use ALONGSIDE BPP; do not claim species from BPP alone post-2021
```

## Dsuite — Introgression vs Incomplete Lineage Sorting

**Goal:** Test for introgression before claiming sister-species relationships, distinguishing admixture from ILS.

**Approach:** Compute D-statistic (ABBA-BABA), f4-ratio (admixture fraction), and f-branch (Malinsky 2021 *Mol Ecol Resour* 21:584-595) on a VCF with population assignments. Non-zero D indicates introgression between H3 and one of H1/H2.

```bash
# Dsuite (https://github.com/millanek/Dsuite)
# Sets file: 3 columns (sample_id, population/species, optional outgroup tag)
# Tree file: Newick with population names matching the sets file

# All-trio D-statistic
Dsuite Dtrios -t species_tree.nwk sets.txt input.vcf.gz

# F-branch decomposition (more interpretable than raw D)
Dsuite Fbranch species_tree.nwk Dtrios_output_tree.txt > fbranch.txt

# Plot f-branch matrix
dtools.py fbranch.txt species_tree.nwk

# D > 0 with |Z| > 3: significant ABBA-BABA imbalance -> introgression
# f4-ratio gives the admixture fraction estimate
```

## Per-Method Failure Modes

### MSC delimitation oversplits population structure as species

**Trigger:** Running BPP, BFD*, or other multispecies-coalescent methods on data with within-species substructure (isolated demes), then claiming the delimited lineages are species.

**Mechanism:** The MSC model treats every panmictic population as a "species" in its parameterization (Sukumaran & Knowles 2017). When applied to data with population substructure, MSC partitions that structure into "species."

**Symptom:** Many "species" delimited; each with low sample size; lineages defined by geography rather than morphology; BPP posterior support is high but biological reality is unclear.

**Fix:** Run DELINEATE (Sukumaran 2021) on the BPP output to test which lineages have completed speciation. Validate via integrative taxonomy (Padial 2010): require congruence across genetic + morphological + ecological evidence before publishing species.

### ABGD/ASAP forced to partition data without a barcoding gap

**Trigger:** Reporting "best partition" from ABGD or ASAP when the pairwise distance histogram is unimodal (no gap).

**Mechanism:** ABGD and ASAP detect a gap in the pairwise-distance distribution that separates intra- from inter-specific divergences. Meyer & Paulay 2005 *PLoS Biol* 3:e422 demonstrated the gap is frequently absent.

**Symptom:** Distance histogram is unimodal (no clear bimodality); ASAP partitions have similar scores across many K values; ABGD priors disagree.

**Fix:** Inspect the distance histogram first; if no clear gap, do NOT report ASAP/ABGD results as definitive; integrate other lines of evidence.

### bPTP oversplit on populations with strong substructure

**Trigger:** Running bPTP on a tree where the same biological species has multiple isolated demes contributing genetic structure.

**Mechanism:** bPTP detects branching-rate shifts as species boundaries; isolated demes within a species show coalescent-rate branching that bPTP partitions as separate species.

**Symptom:** More bPTP species than morphologically recognized; species correspond to geographic regions.

**Fix:** Cross-check with ASAP (more conservative) and with mPTP at multi-rate option; treat bPTP/mPTP as preliminary; require congruence with other lines of evidence.

### BPP with default theta and tau priors

**Trigger:** Running BPP without calibrating priors to the data; using priors copied from a primate/mammal example file.

**Mechanism:** BPP results depend critically on theta (per-locus nucleotide diversity) and tau (divergence time) priors. Default priors are typically tuned for mammal data; for insects, fish, plants with different mutation-rate and Ne contexts, defaults can dominate the posterior.

**Symptom:** Posterior support clusters at 1.0 OR 0.0; results inconsistent with morphological / ecological data; sensitivity to alternative priors is very high.

**Fix:** Calibrate priors from data per Yang 2015 *Curr Zool* 61:854: estimate per-locus heterozygosity for theta mean; estimate divergence times from observed maximum genetic distances and mutation rate.

### Reporting D-statistic without f-branch

**Trigger:** Reporting raw D-statistics for many taxa without decomposing introgression patterns.

**Mechanism:** D-statistic measures asymmetric site-pattern counts; with taxon-rich datasets, D > 0 can be confounded by ancient introgression in the outgroup or by ghost lineages. f-branch (Malinsky 2018) decomposes admixture across the tree more interpretably.

**Symptom:** Many trios with D > 0; difficult to identify which specific introgression event drives the pattern.

**Fix:** Report f-branch decomposition alongside D; check for consistent patterns across phylogenetic levels.

## Quantitative Thresholds

| Threshold | Value | Source / rationale |
|-----------|-------|-------------------|
| COI K2P "species" cutoff | typically 2-3% pairwise distance (variable) | Animal barcoding convention; verify per taxon |
| ASAP score gap signaling robust partition | Large drop between consecutive scores | Puillandre 2021; inspect top 5-10 partitions |
| bPTP/mPTP posterior support | > 0.95 strong; 0.80-0.95 moderate; < 0.80 uncertain | Standard Bayesian inference |
| GMYC LR test p-value | < 0.05 detects significant transition | Pons 2006 |
| BPP posterior probability | > 0.95 strong; 0.50-0.95 moderate; < 0.50 uncertain | Yang 2015 |
| BPP minimum independent chains | 4 with different seeds | Convergence check |
| Dsuite D significance | |Z| > 3 | Standard test statistic |
| Sukumaran-Knowles validation rule | DELINEATE OR integrative-taxonomy congruence required post-2021 | Modern best practice |
| Integrative taxonomy congruence | Genetic + morphological + ecological all support hypothesis | Padial 2010; Carstens 2013 |

## Common errors

| Error | Cause | Solution |
|-------|-------|----------|
| ABGD/ASAP returns no clear partition | Barcoding gap absent in data | Inspect distance histogram; do not force partition |
| bPTP install error in R | bPTP is Python, NOT R | `pip install git+https://github.com/iTaxoTools/PTP-pyqt5` |
| GMYC error "not ultrametric" | ML tree directly fed without chronos | Run `ape::chronos()` first or use BEAST output |
| BPP all PP = 1.0 (oversplit) | Default priors too informative; population substructure | Calibrate priors from data; run DELINEATE alongside |
| BPP no convergence across chains | Insufficient burnin or sampling | Increase nsample and burnin; multiple chains |
| Dsuite format error | Sets file column ordering | First column = sample ID, second = species/population |
| mPTP not found | Build from source github.com/Pas-Kapli/mptp | `pip install mptp` does not exist |
| ASAP web upload "too many sequences" | > 10^4 sequences | Use downloadable C binary instead of web |

## References

- de Queiroz K (2007) Species concepts and species delimitation. *Syst Biol* 56(6):879-886. doi:10.1080/10635150701701083
- Pons J, Barraclough TG, Gomez-Zurita J et al. (2006) GMYC original. *Syst Biol* 55(4):595-609. doi:10.1080/10635150600852011
- Fujisawa T, Barraclough TG (2013) Revised GMYC. *Syst Biol* 62(5):707-724. doi:10.1093/sysbio/syt033
- Puillandre N, Lambert A, Brouillet S, Achaz G (2012) ABGD original. *Mol Ecol* 21(8):1864-1877. doi:10.1111/j.1365-294X.2011.05239.x
- Puillandre N, Brouillet S, Achaz G (2021) ASAP. *Mol Ecol Resour* 21(2):609-620. doi:10.1111/1755-0998.13281
- Zhang J, Kapli P, Pavlidis P, Stamatakis A (2013) PTP / bPTP. *Bioinformatics* 29(22):2869-2876. doi:10.1093/bioinformatics/btt499
- Kapli P, Lutteropp S, Zhang J et al. (2017) mPTP. *Bioinformatics* 33(11):1630-1638. doi:10.1093/bioinformatics/btx025
- Yang Z, Rannala B (2010) BPP original. *Proc Natl Acad Sci USA* 107(20):9264-9269. doi:10.1073/pnas.0913022107
- Yang Z (2015) BPP review with prior guidance. *Curr Zool* 61(5):854-865. doi:10.1093/czoolo/61.5.854
- Flouri T, Jiao X, Rannala B, Yang Z (2018) BPP v4 genomic-scale. *Mol Biol Evol* 35(10):2585-2593. doi:10.1093/molbev/msy147
- Bryant D, Bouckaert R, Felsenstein J et al. (2012) SNAPP. *Mol Biol Evol* 29(8):1917-1932. doi:10.1093/molbev/mss086
- Leache AD, Fujita MK, Minin VN, Bouckaert RR (2014) BFD*. *Syst Biol* 63(4):534-542. doi:10.1093/sysbio/syu018
- Sukumaran J, Knowles LL (2017) MSC delimits structure, not species. *Proc Natl Acad Sci USA* 114(7):1607-1612. doi:10.1073/pnas.1607921114
- Sukumaran J, Holder MT, Knowles LL (2021) DELINEATE (speciation-process modeling). *PLoS Comput Biol* 17(5):e1008924. doi:10.1371/journal.pcbi.1008924
- Padial JM, Miralles A, De la Riva I, Vences M (2010) Integrative future of taxonomy. *Front Zool* 7:16. doi:10.1186/1742-9994-7-16
- Carstens BC, Pelletier TA, Reid NM, Satler JD (2013) How to fail at species delimitation. *Mol Ecol* 22(17):4369-4383. doi:10.1111/mec.12413
- Malinsky M, Matschiner M, Svardal H (2021) Dsuite. *Mol Ecol Resour* 21(2):584-595. doi:10.1111/1755-0998.13265
- Malinsky M, Svardal H, Tyers AM et al. (2018) Whole-genome sequences of Malawi cichlids reveal multiple radiations interconnected by gene flow (f-branch original). *Nat Ecol Evol* 2(12):1940-1955. doi:10.1038/s41559-018-0717-x
- Meyer CP, Paulay G (2005) Barcoding-gap-absence critique. *PLoS Biol* 3(12):e422. doi:10.1371/journal.pbio.0030422

## Related Skills

- ecological-genomics/edna-metabarcoding - Generate barcode sequences from environmental samples for primary delimitation
- ecological-genomics/conservation-genetics - Population-level genetic assessment (use ESU/MU framework, NOT species delimitation, for management units)
- phylogenetics/tree-io - Tree input/output for tree-based delimitation methods
- phylogenetics/modern-tree-inference - ML and Bayesian tree construction for delimitation input
- database-access/entrez-fetch - Retrieve barcode sequences from GenBank for reference comparison
