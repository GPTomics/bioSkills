# Species Delimitation Usage Guide

## Overview

Delimits putative species boundaries from molecular data within the de Queiroz 2007 unified-lineage framework (species are separately evolving metapopulation lineages; methods test operational criteria, not competing definitions) using complementary approaches: distance-based partitioning via ASAP (Puillandre 2021 successor to ABGD), tree-based mPTP (Kapli 2017 multi-rate, C++; replaces bPTP for most uses), GMYC single- and multi-threshold on ultrametric trees (Pons 2006; Fujisawa 2013), Bayesian multispecies-coalescent BPP v4 with prior calibration from data (NOT default priors), SNAPP + BFD* for SNP-based delimitation, DELINEATE (Sukumaran 2021) for speciation-process modeling that addresses the Sukumaran & Knowles 2017 oversplitting critique, integrative-taxonomy congruence following Padial 2010 and Carstens 2013 ("How to fail at species delimitation"), and Dsuite D-statistics + f-branch for testing introgression vs incomplete lineage sorting before claiming sister species (Malinsky 2021). Includes the Meyer & Paulay 2005 caveat that the barcoding gap is frequently absent in real data.

## Prerequisites

- ASAP web interface (https://bioinfo.mnhn.fr/abi/public/asap/) OR downloadable C binary
- mPTP: build from C source at https://github.com/Pas-Kapli/mptp (faster successor to bPTP)
- bPTP (Python, NOT R): `pip install git+https://github.com/iTaxoTools/PTP-pyqt5`
- R with splits for GMYC (`install.packages('splits', repos = 'http://R-Forge.R-project.org')`)
- ape for tree manipulation (`install.packages('ape')`)
- phytools for `force.ultrametric` and tree visualization (`install.packages('phytools')`)
- BPP v4 (https://github.com/bpp/bpp; control-file driven)
- SNAPP via BEAST 2 (https://www.beast2.org/snapp/)
- DELINEATE (https://github.com/jsukumaran/delineate; Python, current speciation-process modeling)
- Dsuite (https://github.com/millanek/Dsuite) for D-statistic and f-branch
- fossil for partition comparison (`install.packages('fossil')`)

## Quick Start

Tell your AI agent what you want to do:

- "Run primary species delimitation with ASAP (modern successor to ABGD) on my COI alignment"
- "Apply mPTP for tree-based delimitation with heterogeneous intraspecific rates"
- "Run BPP A10 with priors CALIBRATED to data (not defaults); follow with DELINEATE to test for oversplitting"
- "Test introgression with Dsuite (D, f4-ratio, f-branch) BEFORE claiming sister-species relationships"
- "Compare ASAP / mPTP / GMYC partitions for multi-method consensus; report congruence per Carstens 2013"
- "Inspect pairwise distance histogram for barcoding-gap presence (often ABSENT per Meyer & Paulay 2005)"

## Example Prompts

### Primary Delimitation

> "I have an aligned FASTA of 200 COI barcode sequences. Run ASAP with K2P distances. Report the top 5-10 partitions ranked by asap-score; look for large score gaps signaling robust partitioning. Also inspect the pairwise distance histogram for barcoding-gap presence."

### Tree-Based Delimitation

> "Run mPTP on my rooted ML tree (faster and supports heterogeneous rates compared to bPTP). For shallow phylogenies with similar rates, also cross-check with bPTP (Python, not R)."

### MSC Oversplitting Caveat

> "Run BPP A10 with priors CALIBRATED from data: estimate theta prior mean from observed per-locus heterozygosity; estimate tau prior from observed maximum divergence and mutation rate. Run multiple independent chains. THEN apply DELINEATE on the BPP output to test which lineages represent fully-formed species vs incomplete-speciation structure (Sukumaran 2021)."

### Introgression vs ILS

> "Before claiming sister-species relationships among my candidate species, run Dsuite Dtrios for all-trio D-statistics and Fbranch for tree-aware admixture decomposition (Malinsky 2021). |Z| > 3 with D > 0 indicates introgression beyond ILS."

### Integrative Taxonomy

> "Compare ASAP + mPTP + GMYC partitions; compute adjusted Rand index; report individuals assigned to the same species by >= 2 methods as 'congruent'. Recommend morphological / ecological / geographic validation for any candidate species per Padial 2010 and Carstens 2013."

### Conservation Management Context

> "Do NOT use BPP / BFD* / multispecies-coalescent methods for management unit (MU) or evolutionarily significant unit (ESU) definition. Use the Moritz 1994 framework (reciprocal monophyly + nuclear divergence). MSC methods oversplit per Sukumaran-Knowles 2017."

## What the Agent Will Do

1. Inspect the pairwise distance histogram FIRST; check for a clear barcoding gap (often absent per Meyer & Paulay 2005); if no gap, flag the data as unsuitable for distance-based partition methods
2. Run ASAP as the primary single-locus delimitation method (modern successor to ABGD); report top 5-10 partitions ranked by asap-score
3. Run mPTP for tree-based delimitation as the modern default (faster, supports heterogeneous rates); only fall back to bPTP for very shallow phylogenies
4. Convert ML trees to ultrametric (chronos or BEAST) BEFORE running GMYC; verify with `is.ultrametric()`
5. For multilocus / genomic data: run BPP A10 with priors CALIBRATED from data (NOT defaults); run multiple independent chains for convergence check
6. Apply DELINEATE alongside BPP to address the Sukumaran-Knowles 2017 oversplitting critique
7. For sister-species claims: run Dsuite Dtrios for D-statistics and Fbranch for admixture decomposition; |Z| > 3 with D > 0 signals introgression
8. Compute adjusted Rand index across methods; report congruent vs discordant assignments
9. Recommend integrative-taxonomy validation (Padial 2010; Carstens 2013) for any species hypothesis: genetic + morphological + ecological congruence
10. For conservation management questions: use ESU / MU frameworks (Moritz 1994), NOT species-delimitation methods

## Tips

- The Sukumaran & Knowles 2017 critique is the most important modern insight: MSC methods (BPP, BFD*) delimit genetic STRUCTURE, not species; pure-MSC species claims post-2017 are considered methodologically incomplete by serious systematics reviewers
- DELINEATE (Sukumaran 2021) is the modern complement to BPP for genomic-scale delimitation; addresses the oversplitting problem by modeling speciation as an extended process
- ASAP supersedes ABGD as the modern distance-based primary delimitation method; faster, less prior-dependent, ranks partitions by asap-score
- mPTP supersedes bPTP for most uses (5+ orders of magnitude faster, multi-rate model); for shallow phylogenies with similar intraspecific rates, bPTP may still be useful
- bPTP is a PYTHON tool, not R; `install.packages('PTP')` does not exist; `pip install git+https://github.com/iTaxoTools/PTP-pyqt5`
- ASAP is a STANDALONE C binary OR web service (https://bioinfo.mnhn.fr/abi/public/asap/); `install.packages('ASAP')` does not exist for the species delimitation tool
- The COI barcoding gap is typically at 2-3% pairwise distance for animals BUT often absent (Meyer & Paulay 2005); always inspect the distance histogram first
- BPP results depend critically on theta and tau priors; defaults are tuned for mammal data and may be wrong for insects, fish, or plants
- Yang 2015 *Curr Zool* gives taxon-specific prior recommendations; calibrate theta prior mean from per-locus heterozygosity estimates
- GMYC requires strictly ultrametric trees; use `ape::chronos()` for quick conversion or BEAST2 for rigorous time-calibration; verify with `is.ultrametric()` (which has small floating-point tolerance)
- For sister-species claims with possible hybridization: run Dsuite D-statistic and f-branch BEFORE concluding sister relationship; ILS-driven discordance can be mistaken for introgression and vice versa
- Integrative taxonomy (Padial 2010) requires congruence across multiple lines of evidence: genetic, morphological, ecological, geographic; single-method species claims are unreliable
- For conservation management, the ESU / MU framework (Moritz 1994) is more appropriate than species delimitation; cite the Sukumaran-Knowles caveat as the reason for not using BPP/BFD* for management units
- Recently-diverged populations (Ne*t << 1) are particularly prone to MSC oversplitting; for shallow divergences, weight DELINEATE and integrative-taxonomy evidence higher than raw BPP posterior

## Related Skills

- ecological-genomics/edna-metabarcoding - Generate barcode sequences from environmental samples
- ecological-genomics/conservation-genetics - Population-level genetic assessment (use ESU/MU framework not species delimitation)
- phylogenetics/tree-io - Tree input/output for tree-based delimitation
- phylogenetics/modern-tree-inference - ML and Bayesian tree construction
- database-access/entrez-fetch - Retrieve barcode sequences from GenBank for reference
