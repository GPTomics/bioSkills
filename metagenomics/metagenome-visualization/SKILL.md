---
name: bio-metagenomics-visualization
description: Turns a shotgun profiler table (MetaPhlAn relative abundance, Bracken counts, HUMAnN function tables) into honest figures and defensible community statistics with phyloseq, vegan, microViz, and Python. Covers why an ordination/bar/diversity number is a modeling choice that can manufacture a result, the MetaPhlAn-percent-vs-Bracken-counts fork that decides everything, CLR/Aitchison vs Bray-Curtis, Hill numbers and why shotgun richness is a database readout, pairing PERMANOVA with betadisper, and the multi-tool differential-abundance consensus. Use when plotting taxonomic/functional profiles, computing alpha/beta diversity, running ordination/PERMANOVA, or testing differential abundance. For amplicon/QIIME2 stats see the microbiome category; for compositional theory see abundance-estimation.
tool_type: mixed
primary_tool: phyloseq
---

## Version Compatibility

Reference examples tested with: phyloseq 1.46+, vegan 2.6+, microViz 0.12+, ALDEx2 1.34+, pandas 2.2+, scikit-bio 0.6+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `ktImportTaxonomy` (no args) to confirm Krona column flags on the installed build

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The input profiler decides the toolchain: MetaPhlAn gives relative abundance (cannot rarefy; ships an SGB tree so UniFrac is available); Bracken gives counts (rarefaction and count models valid; no tree); HUMAnN gives gene-family/pathway features. Record the profiler, its filtering (Kraken `--confidence`, Bracken threshold, MetaPhlAn `--stat_q`), and every modeling choice below - they are not recoverable from the figure.

# Metagenome Visualization

**"Show me how my communities differ."** -> Choose a transform, a distance, and a test - each a modeling choice that can create or erase the difference - then declare them and show the conclusion survives them.
- R: phyloseq + vegan + microViz on a parsed profiler table
- Python: pandas + scikit-bio + matplotlib (plotting and wrangling; R is primary for community stats)

Scope: visualizing and testing a shotgun profiler table. Input generation -> kraken-classification, metaphlan-profiling, abundance-estimation, functional-profiling. Compositional theory and absolute load -> abundance-estimation. Amplicon/QIIME2 stats -> the microbiome category. Generic plotting primitives -> data-visualization.

## The Single Most Important Modern Insight -- A Figure Is a Modeling Choice, Not an Observation

The distance metric, the transform, the rarefaction depth, the confidence filter, and the top-N cutoff are knobs turned before the answer is seen; turning them differently gives a different paper. A shotgun table is a compositional, depth-confounded, false-positive-laden estimate, and the figure inherits all of it. The job is to declare the modeling choices and show they did not fabricate the conclusion. Memory hooks:

- A stacked bar of relative abundance hides absolute load and visually inflates whatever is already dominant - it shows the relative race, not whether the community bloomed or collapsed.
- Bray-Curtis on relative abundance is the field default and is compositionally incoherent; Aitchison (Euclidean on CLR) is correct and almost nobody runs it.
- Shotgun richness is a readout of the database and the confidence threshold, not biology.
- A significant PERMANOVA may be a difference in variability, not composition.
- The differentially abundant taxa reported depend more on which DA tool was run than on biology.

The fork that decides everything downstream: **MetaPhlAn = percent, Bracken = counts.** Rarefaction is only meaningful for counts; CLR needs a pseudocount on percentages; richness estimators (Chao1/Observed) require integer counts and are meaningless on MetaPhlAn relative abundances.

## Honest Composition Plots

Krona gives the full drillable hierarchy - often the honest answer to "what is in it," versus a bar that pre-collapses to the top 10:

```bash
kreport2krona.py -r kraken.kreport -o krona.txt && ktImportText krona.txt -o krona.html
# Confirm ktImportTaxonomy column flags with no-arg usage; they drift across versions.
```

For a stacked bar: collapse rare taxa to "Other" but label how many taxa and what percent that hides; state n per group (never stack one representative sample); use a colorblind-safe palette of <=12 colors; show absolute load alongside if total-biomass data exist. microViz `comp_barplot()` handles Other and palettes.

## Alpha Diversity: Hill Numbers and Richness Honesty

**Goal:** Report within-sample diversity in interpretable units without letting shotgun richness masquerade as biology.

**Approach:** Use Hill numbers (effective species) at q=0/1/2; prefer evenness-weighted q=1/q=2 for shotgun because richness (q=0) is dominated by database size and false positives; only compute richness estimators on integer counts.

```r
library(phyloseq); library(vegan)
# estimate_richness Observed/Chao1 assume INTEGER COUNTS - valid for Bracken, meaningless on MetaPhlAn %.
alpha <- estimate_richness(ps_counts, measures = c('Shannon', 'InvSimpson'))
hill_q1 <- exp(alpha$Shannon)        # effective species (q=1), interpretable units
hill_q2 <- alpha$InvSimpson          # q=2, dominance-weighted; robust to rare-taxon noise
```

Richness moves drastically with the classifier's confidence/threshold and 25-70% of shotgun species can be false positives - report the filtering, show a rarefaction or coverage curve, and gate richness behind unique-k-mer evidence (KrakenUniq). The rarefaction debate is unresolved: McMurdie & Holmes 2014 (*PLoS Comput Biol* 10:e1003531) call rarefying inadmissible for differential abundance; Schloss 2024 (*mSphere* 9:e00355-23) defends it for diversity. Decide per analysis - rarefy/coverage-standardize for diversity, model-based for DA - and show the curve either way.

## Beta Diversity and Ordination

| Distance | Compositionally coherent | Needs tree | When |
|----------|--------------------------|------------|------|
| Bray-Curtis | no | no | field default, intuitive; label it incoherent |
| Aitchison (Euclidean on CLR) | yes | no | the compositional-correct choice; needs zero handling |
| Robust Aitchison / RPCA (DEICODE) | yes | no | handles sparsity without a pseudocount |
| Weighted/unweighted UniFrac | partial | yes | only if a tree exists (MetaPhlAn SGB tree; not Bracken/HUMAnN) |

PCoA decomposes any distance; PCA on CLR is the compositional-coherent ordination and gives taxon loadings (which taxa drive an axis) that PCoA cannot; NMDS reports stress (not variance); UMAP does not preserve global distances and is for spotting clusters, not measuring dissimilarity. Pick the geometry for a stated reason and show the conclusion survives at least one alternative.

```r
library(vegan)
dist_bc <- vegdist(otu_matrix, method = 'bray')   # samples as rows
pm <- adonis2(dist_bc ~ group, permutations = 999, by = 'terms')
# ALWAYS pair PERMANOVA with a dispersion test - a significant adonis2 can be a spread difference,
# not a location shift (Anderson & Walsh 2013), especially for unbalanced designs.
bd <- betadisper(dist_bc, group); permutest(bd, permutations = 999)
```

If betadisper is significant, the PERMANOVA is ambiguous - report both.

## Differential Abundance: Consensus, Not a Single Tool

DA methods disagree wildly across datasets (Nearing 2022 *Nat Commun* 13:342); the taxa called significant depend more on the tool than on biology. Run at least two compositionally aware methods, report the intersect as high-confidence and the union as exploratory, and name every tool. ALDEx2 and ANCOM-II were the most conservative/consistent in Nearing; LinDA and ANCOM-BC were well FDR-controlled in Yang & Chen 2022. Prevalence-filter first and BH-correct across taxa.

| Tool | Model | Citation |
|------|-------|----------|
| ALDEx2 | Dirichlet Monte-Carlo -> CLR | Fernandes 2014 *Microbiome* 2:15 |
| ANCOM-BC | log-linear + sampling-fraction bias correction | Lin & Peddada 2020 *Nat Commun* 11:3514 |
| MaAsLin2 | general linear model + covariates | Mallick 2021 *PLoS Comput Biol* 17:e1009442 |
| LinDA | linear model on CLR + bias correction | Zhou 2022 *Genome Biol* 23:95 |

The anti-pattern to forbid: an uncorrected Wilcoxon or t-test on raw relative abundances - wrong on compositionality and on multiple testing in one line. Most tools want counts (Bracken); for MetaPhlAn percentages use methods that accept proportions (MaAsLin2/LinDA) or convert to pseudo-counts.

## Compositional Ordination in Python

**Goal:** Produce a compositionally coherent ordination in Python instead of the incoherent StandardScaler-then-PCA on raw relative abundance.

**Approach:** Replace zeros, CLR-transform, then PCA on the CLR coordinates (this is Aitchison-PCA); the loadings are interpretable as taxa.

```python
import pandas as pd
from skbio.stats.composition import clr, multi_replace   # renamed from multiplicative_replacement in skbio 0.6
from sklearn.decomposition import PCA

ab = pd.read_csv('merged_abundance.txt', sep='\t', index_col=0)
ab = ab[ab.index.str.contains(r'\|s__') & ~ab.index.str.contains(r'\|t__')]  # species rows only
proportions = (ab.T.values / ab.T.values.sum(axis=1, keepdims=True))
clr_mat = clr(multi_replace(proportions))                # zeros replaced, then CLR (NOT StandardScaler on relab)
pca = PCA(n_components=2).fit(clr_mat)
coords = pca.transform(clr_mat)
```

## Per-Method Failure Modes

### The transform/metric/rarefaction triad manufactures the result
**Trigger:** defaulting to Bray-Curtis PCoA and presenting it as "the" answer. **Mechanism:** metric, transform, and rarefaction each change the geometry. **Symptom:** separation that vanishes under Aitchison/RPCA, or appears only at one rarefaction depth. **Fix:** state the choice; show the conclusion survives a second reasonable choice; report % variance / stress.

### PERMANOVA dispersion confound
**Trigger:** "communities differed (adonis2 p<0.001)" with no dispersion check. **Mechanism:** pseudo-F responds to within-group spread, not only centroid location (unbalanced designs especially). **Symptom:** a "difference" that is really higher variability in one group. **Fix:** pair every adonis2 with betadisper + permutest; report both.

### Single-tool differential abundance
**Trigger:** reporting the DA tool that gives the prettiest story. **Mechanism:** tools disagree (Nearing 2022). **Symptom:** findings that do not replicate. **Fix:** consensus of >=2 compositional tools; intersect = confident; name them; BH-correct.

### Richness as biology
**Trigger:** plotting Observed/Chao1 from a k-mer classifier as diversity. **Mechanism:** richness tracks database size and false positives; estimators assume integer counts. **Symptom:** "richness" differences driven by depth/filtering; meaningless Chao1 on MetaPhlAn percentages. **Fix:** prefer Hill q=1/q=2; gate richness behind confidence filtering and unique-k-mer evidence.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Prevalence filter (e.g. present in >=10% of samples) | DA best practice | reduces multiple-testing and unstable zero-dominated taxa |
| BH FDR across taxa | standard | many taxa = many tests; report q-values |
| NMDS stress < 0.2 usable, < 0.1 good | Clarke 1993 *Aust J Ecol* 18:117 | above 0.2 the configuration is suspect |
| Hill q=0,1,2 reported together | Jost 2007; Chao 2014 | characterize the richness-evenness spectrum, not richness alone |
| Pair adonis2 with betadisper | Anderson & Walsh 2013 *Ecol Monogr* 83:557 | dispersion can masquerade as a location difference |
| Consensus of >=2 DA tools | Nearing 2022 *Nat Commun* 13:342 | single-tool hits do not replicate |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Chao1/Observed nonsensical | run on MetaPhlAn relative abundances | compute richness on Bracken counts only |
| PCA shows no compositional structure | StandardScaler on raw relab | CLR-transform first, then PCA (Aitchison-PCA) |
| "Significant" PERMANOVA challenged in review | no dispersion test | add betadisper + permutest |
| DA hits do not replicate | single tool, uncorrected test | >=2 compositional tools, BH correction |
| UniFrac errors on Bracken data | no phylogeny for NCBI taxonomy | UniFrac needs a tree (MetaPhlAn SGB tree only) |
| Krona flags rejected | version-fragile column flags | run `ktImportTaxonomy` with no args to confirm |

## References

- McMurdie PJ, Holmes S. 2013. phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS One* 8:e61217.
- Ondov BD, Bergman NH, Phillippy AM. 2011. Interactive metagenomic visualization in a web browser. *BMC Bioinformatics* 12:385.
- Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. 2017. Microbiome datasets are compositional: and this is not optional. *Front Microbiol* 8:2224.
- Nearing JT, Douglas GM, Hayes MG, et al. 2022. Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun* 13:342.
- Anderson MJ, Walsh DCI. 2013. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions. *Ecol Monogr* 83:557-574.
- McMurdie PJ, Holmes S. 2014. Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Comput Biol* 10:e1003531.
- Schloss PD. 2024. Waste not, want not: revisiting the analysis that called into question the practice of rarefaction. *mSphere* 9:e00355-23.
- Fernandes AD, Reid JN, Macklaim JM, et al. 2014. Unifying the analysis of high-throughput sequencing datasets. *Microbiome* 2:15.
- Lin H, Peddada SD. 2020. Analysis of compositions of microbiomes with bias correction. *Nat Commun* 11:3514.
- Mallick H, Rahnavard A, McIver LJ, et al. 2021. Multivariable association discovery in population-scale meta-omics studies. *PLoS Comput Biol* 17:e1009442.
- Barnett DJM, Arts ICW, Penders J. 2021. microViz: an R package for microbiome data visualization and statistics. *J Open Source Softw* 6:3201.

## Related Skills

- metaphlan-profiling - Generates the relative-abundance table (with the SGB tree)
- kraken-classification - Generates Kraken/Bracken count input
- abundance-estimation - Compositional theory, normalization, and absolute load
- functional-profiling - HUMAnN function tables tested with the same DA logic
- microbiome/diversity-analysis - Amplicon/QIIME2 diversity and differential abundance for ASV input
- data-visualization/ggplot2-fundamentals - Generic plotting primitives
- workflows/metagenomics-pipeline - End-to-end shotgun analysis
