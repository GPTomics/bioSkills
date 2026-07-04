# Introgression Detection - Usage Guide

## Overview

Tests for inter-population admixture combine site-pattern statistics (ABBA-BABA, f4), tree-topology methods (Dsuite f-branch, QuIBL, Twisst), explicit networks (PhyloNet, qpGraph), and haplotype-based methods (sprime, Relate). The fundamental confounder is **incomplete lineage sorting (ILS)**: a significant Patterson D-statistic is consistent with introgression OR with asymmetric ILS OR with ancestral structure OR with ghost lineages (Eriksson & Manica 2012 PNAS 109:13956).

For publication-grade introgression claims, combine D-statistic with **Fbranch** mapping (Malinsky 2018), **Twisst / QuIBL** topology weighting at the locus level, and at least one network method.

## Prerequisites

```bash
# Modern Dsuite
git clone https://github.com/millanek/Dsuite && cd Dsuite && make
# Or: conda install -c bioconda dsuite

# Population-genetics ABBA-BABA + qpAdm
git clone https://github.com/DReichLab/AdmixTools && cd AdmixTools && make
# Modern AdmixTools v2 (R wrapper)
Rscript -e "remotes::install_github('uqrmaie1/admixtools')"

# TreeMix
conda install -c bioconda treemix

# Locus-level topology
git clone https://github.com/miriamtnzr/QuIBL  # QuIBL
git clone https://github.com/simonhmartin/twisst # Twisst

# Network methods
git clone https://github.com/NakhlehLab/PhyloNet && cd PhyloNet && mvn package

# Haplotype-based
conda install -c bioconda sprime

# Hybridization-specific
git clone https://github.com/pblischak/HyDe && cd HyDe && pip install -e .

# Optimal migration count
Rscript -e "install.packages('OptM')"
```

## Quick Start

Tell the AI agent what to do:
- "Compute Patterson's D-statistic between these populations using Dsuite Dtrios"
- "Apply Fbranch to assign admixture signal to specific tree branches"
- "Run TreeMix with migration edges and select optimal number with OptM"
- "Use sprime to identify archaic introgressed haplotype tracts in modern human populations"

## Example Prompts

### Inter-Species Admixture Test

> "Test for introgression between Heliconius melpomene and H. cydno using Dsuite Dtrios with H. erato as outgroup. Apply Fbranch mapping. Cross-validate with Twisst topology weighting in 50 kb windows across the genome to identify specific introgressed regions."

### Population Admixture Modeling

> "Test if the modern Argentinean population can be modeled as admixed from Spanish + Native American sources using qpAdm rotation. Cross-validate with TreeMix migration edges. Use OptM (Evanno method) to select optimal number of migration edges."

### Archaic Introgression

> "Identify Neanderthal introgressed haplotype tracts in modern non-African human samples using sprime. Verify that tract lengths cluster around expected size (>50 kb typical for archaic). Compare introgression rates across populations."

## What the Agent Will Do

1. **Validate inputs**: VCF format; sample-population mapping; outgroup choice; population sample sizes
2. **Compute D-statistic** with Dsuite Dtrios across all population trios
3. **Apply Fbranch** to assign admixture to specific tree branches
4. **Cross-validate locus level** with Twisst (topology weighting per window) or QuIBL (per gene tree)
5. **Build migration model** with TreeMix; select migration count via OptM
6. **qpAdm rotation tests** for proposed source populations
7. **For ancient admixture**: qpGraph manual topology design
8. **For archaic introgression**: sprime per-individual analysis
9. **Report**: D + Z-score per trio, Fbranch attributions, Twisst topology weights, qpAdm rotation pass/fail
10. **Caveats**: ILS vs introgression, ghost lineages, ancestral structure, outgroup distance

## Tips

- Always interpret D > 0 as "consistent with introgression OR ILS OR ancestral structure" (Eriksson & Manica 2012)
- Use Fbranch (Malinsky 2018) to assign signal to specific branches; sharper than raw D
- Cross-validate with Twisst / QuIBL at locus level to distinguish ILS from introgression
- For multi-population analysis, qpAdm + TreeMix combination is gold standard
- TreeMix migration count selection: use OptM Evanno method, not raw likelihood
- For divergent species, use ABBAclustering (Koppetsch-Malinsky-Matschiner 2024 Syst Biol)
- Outgroup distance: < 100 Myr ideal; > 200 Myr introduces convergent-evolution artifacts
- One individual per population is standard (Dsuite default); for population-genetic, 5+ samples
- Always report sample sizes; small N inflates D variance
- Multiple-testing correction across many trios: FDR or Bonferroni
- For sprime, restrict to populations where archaic introgression source is known
- For ancient admixture, qpGraph with explicit topology > automatic methods
- ILS confound dominates short internodes; report internode length context
- Cytonuclear discordance: mtDNA introgression may differ from autosomal

## Related Skills

comparative-genomics/hgt-detection - Distinguish HGT from hybridization in microbes
comparative-genomics/gene-tree-species-tree-reconciliation - Reconciliation alternative for introgression
comparative-genomics/synteny-analysis - Per-window topology analysis
population-genetics/population-structure - PCA / ADMIXTURE precedes introgression testing
population-genetics/selection-statistics - Selection on introgressed regions
phylogenetics/species-trees - Coalescent species tree
variant-calling/joint-calling - Multi-sample VCF input
read-alignment/bwa-alignment - Alignment underlies variant calling
