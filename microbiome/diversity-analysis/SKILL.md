---
name: bio-microbiome-diversity-analysis
description: Alpha and beta diversity of an amplicon (16S/ITS) ASV/OTU community table - observed features, Shannon, Pielou evenness, Faith PD, Bray-Curtis, Jaccard, weighted/unweighted/generalized UniFrac, Aitchison/RPCA - via QIIME2 core-metrics-phylogenetic, phyloseq/vegan, and scikit-bio. Covers the three knobs that set the answer before it is seen (rarefaction sampling depth, the tree, the metric), why core-metrics silently deletes samples below the sampling depth, why de novo trees lose to SEPP fragment-insertion and Greengenes2, why unweighted and weighted UniFrac can flip the story, why observed features is an ASV count not a species count, the QIIME2-log2 vs R-ln Shannon mismatch, and pairing PERMANOVA (adonis2) with betadisper. Use when summarizing whole-community richness/evenness or testing group differences in community structure. Per-taxon testing -> differential-abundance. Shotgun tables -> metagenomics/metagenome-visualization. Shared CoDA/rarefaction theory -> metagenomics/abundance-estimation.
tool_type: mixed
primary_tool: phyloseq
---

## Version Compatibility

Reference examples tested with: phyloseq 1.46+, vegan 2.6+, picante 1.8+, GUniFrac 1.8+, scikit-bio 0.6+, QIIME2 2024.2+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `qiime <plugin> <action> --help` to confirm flags
- Python: `pip show scikit-bio` then `help(skbio.diversity.beta_diversity)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

scikit-bio 0.6.0 renamed OTU to taxon across the API and drifted metric kwargs (`otu_ids=` vs newer forms) - discover names with `skbio.diversity.get_beta_diversity_metrics()` before hard-coding. UniFrac/Faith PD results inherit the tree (de novo vs SEPP vs Greengenes2 reference build) AND the chosen sampling depth - record both alongside the QIIME2 release that produced the `.qza` artifacts.

# Diversity Analysis

**"Compare microbial diversity across my samples"** -> Summarize within-sample richness/evenness (alpha) and between-sample dissimilarity (beta) - but only after declaring the rarefaction depth, the tree, and the metric, because each is a knob that sets the answer before it is seen.
- CLI: `qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth N --m-metadata-file md.tsv --output-dir cm/`
- R: `estimate_richness(ps_rare)` for alpha; `UniFrac(ps_rare, weighted=)` / `vegdist()` then `adonis2()` + `betadisper()` for beta

Scope: whole-community summary (a number or ordination per sample) of an amplicon ASV/OTU table plus a tree. Per-taxon between-group testing -> differential-abundance. Shotgun profiler tables (MetaPhlAn/Bracken) -> metagenomics/metagenome-visualization. The shared CoDA and rarefaction-debate theory lives in metagenomics/abundance-estimation; the Hill-number and PERMANOVA-dispersion theory in metagenomics/metagenome-visualization - cross-referenced here, not re-derived. Tree handling -> phylogenetics/tree-io.

## The Single Most Important Modern Insight -- A Diversity Number Is the Output of Three Knobs Turned Before the Answer Appears

An alpha or beta diversity value is not a measurement of the community; it is the output of three choices made before the number appears - the rarefaction DEPTH, the TREE, and the METRIC. Turn them differently and the conclusion can change. The job is to declare all three and show the result survives a second reasonable choice, not to run `core-metrics-phylogenetic` and read the p-value. The quietest and most dangerous knob is the depth:

1. **`--p-sampling-depth` is a sample-deletion knob in a normalization costume.** core-metrics rarefies every sample to the depth and, per the QIIME2 docs, silently drops every sample whose total count is below it - no warning, just fewer points in the PCoA. The dropped samples are the lowest-yield ones (the lowest-biomass swab, the sickest patient, the failed extraction), so the loss is almost never random. Pick the depth from the feature-table summary plus the alpha-rarefaction plateau, report the depth AND the dropped samples, and confirm the conclusion at a nearby depth. Rarefying to `min(sample_sums)` is the worst of both worlds - one tiny library drags everyone to noise.
2. **UniFrac and Faith PD are only as real as the tree, and the tree is a model not a property of the data.** A de novo MAFFT+FastTree tree from ~250 bp reads is poorly resolved and arbitrarily midpoint-rooted (Janssen 2018); SEPP fragment-insertion into a full-length reference, or Greengenes2 placement, gives stable topology and correct associations - and SEPP/GG2 align 16S with shotgun (McDonald 2024). SEPP also drops fragments that fail to insert, a second silent table-shrink.
3. **Rarefy for diversity, never for differential abundance.** Rarefaction-to-even-depth is defensible for alpha/beta (Schloss 2024); for DA it discards count information a compositional model needs (McMurdie 2014). Keep the raw counts; rarefy only into the diversity branch; route DA to differential-abundance on the unrarefied table.

## Tool / Metric Taxonomy

| Metric / tool | Citation | What it measures / does | When |
|---------------|----------|-------------------------|------|
| Observed features | - | ASV richness (Hill q=0); most depth-sensitive; an ASV count, not species | richness, but report denoising params; prefer Hill q1/q2 |
| Shannon | - | entropy = richness+evenness (Hill q=1 = exp(H')); QIIME2 log2/bits, R ln/nats | balanced diversity; report exp(H') to dodge the base |
| Pielou evenness | Pielou 1966 *J Theor Biol* 13:131 | H'/ln(S); 0-1; isolates evenness from richness | when evenness is the question |
| Faith PD | Faith 1992 *Biol Conserv* 61:1 | sum of branch lengths spanning observed taxa; phylogenetic q=0 | amplicon-native richness; needs a tree |
| Jaccard | - | presence/absence dissimilarity; no tree | membership turnover; depth/rare-ASV sensitive |
| Bray-Curtis | - | abundance dissimilarity; no tree; compositionally incoherent | abundance default; intuitive, label the caveat |
| Unweighted UniFrac | Lozupone 2005 *Appl Environ Microbiol* 71:8228 | branch length unique to one community (presence/absence) | rare/divergent lineages + topology; needs a tree |
| Weighted UniFrac | Lozupone 2007 *Appl Environ Microbiol* 73:1576 | branch length weighted by abundance difference | abundant-lineage shifts; needs a tree |
| Generalized UniFrac | Chen 2012 *Bioinformatics* 28:2106 | alpha in [0,1] interpolating unweighted-weighted | alpha=0.5 compromise; powerful for moderately abundant lineages |
| Aitchison / RPCA | Martino 2019 *mSystems* 4:e00016-19 | CLR + matrix completion; ordination with feature loadings | compositionally coherent; sparse data; no pseudocount |
| SEPP insertion | Janssen 2018 *mSystems* 3:e00021-18 | places ASVs into a full-length reference tree | the preferred tree for short reads |
| Greengenes2 | McDonald 2024 *Nat Biotechnol* 42:715 | unified genome+16S reference tree | makes 16S UniFrac comparable to shotgun |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Need a phylogenetic metric (UniFrac, Faith PD) | SEPP-into-reference or Greengenes2 tree | de novo from short reads is unstable (Janssen 2018) |
| De novo tree is the only option | treat unweighted UniFrac with suspicion | topology noise on ~250 bp reads dominates it |
| Change is in rare/low-abundance lineages | unweighted UniFrac, observed features | presence/absence + topology see rare taxa |
| Change is a bloom of dominant taxa | weighted UniFrac, Bray-Curtis | abundance-weighted metrics see dominant shifts |
| Do not want to metric-shop | generalized UniFrac alpha=0.5 + report both un/weighted | Chen 2012 compromise; single-metric hit is tentative |
| Richness vs evenness question | observed/Faith (q0) AND Shannon-exp (q1) / InvSimpson (q2) | span the richness-evenness spectrum |
| Compositional, want axis-driving taxa | RPCA (DEICODE/gemelli) | CLR ordination with interpretable loadings |
| Picking a rarefaction depth | feature-table summarize + alpha-rarefaction plateau | depth must retain samples AND saturate richness |
| Per-taxon "which bug changed" | -> differential-abundance | diversity is whole-community; DA is per-feature |
| Shotgun profiler table, not amplicon | -> metagenomics/metagenome-visualization | no per-feature tree; different idiom |

## Choosing the Sampling Depth (the biggest lever)

**Goal:** Pick a rarefaction depth that saturates richness while retaining an acceptable fraction of samples, and know exactly which samples were dropped.

**Approach:** Read the per-sample frequency distribution from the feature-table summary, find where the alpha-rarefaction curve plateaus, set the depth there, then declare the depth and the dropped-sample list.

```bash
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv   # per-sample frequencies; the depth lives here

qiime diversity alpha-rarefaction \
    --i-table table.qza --i-phylogeny rooted-tree.qza \
    --p-max-depth 20000 \   # set near the median sample depth; the curve panel shows survivors per depth
    --m-metadata-file metadata.tsv --o-visualization alpha-rarefaction.qzv

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza --i-table table.qza \
    --p-sampling-depth 10000 \   # on the observed-features plateau; SILENTLY DROPS samples below this
    --m-metadata-file metadata.tsv --output-dir core-metrics-results
```

core-metrics-phylogenetic rarefies the table, computes the four alpha vectors (`faith_pd_vector`, `observed_features_vector`, `shannon_vector`, `evenness_vector`) and four beta matrices (`unweighted_unifrac_`, `weighted_unifrac_`, `jaccard_`, `bray_curtis_distance_matrix`), and produces a PCoA + Emperor plot for each beta metric. The non-phylogenetic twin `qiime diversity core-metrics` drops Faith PD and both UniFracs and needs no tree.

## Building the Tree (a modeling choice, not a fixed step)

**Goal:** Obtain a phylogeny over the ASVs that does not inject topology noise into UniFrac/Faith PD.

**Approach:** Prefer SEPP fragment-insertion into a full-length reference (or Greengenes2 placement) over a de novo build from short reads; for de novo, mask the alignment and accept that unweighted UniFrac will be shaky.

```bash
qiime fragment-insertion sepp \
    --i-representative-sequences rep-seqs.qza \
    --i-reference-database sepp-refs-gg-13-8.qza \
    --p-threads 4 \
    --o-tree insertion-tree.qza --o-placements insertion-placements.qza

qiime fragment-insertion filter-features \
    --i-table table.qza --i-tree insertion-tree.qza \
    --o-filtered-table table-sepp.qza --o-removed-table removed-table.qza   # fragments that failed to insert are DROPPED
```

De novo is `qiime phylogeny align-to-tree-mafft-fasttree` (MAFFT align -> mask -> FastTree2 -> midpoint root) - acceptable only when no reference package fits the marker/region, and unweighted UniFrac on it must be treated as suspect.

## Alpha Diversity in R (counts on the rarefied table)

**Goal:** Compute richness and evenness per sample and test for a group difference without confounding by sequencing depth.

**Approach:** Rarefy to a chosen depth, estimate Hill-spanning metrics, test with a non-parametric test (escalate to a linear/mixed model for covariates), and report effective species exp(H').

```r
library(phyloseq); library(vegan)

ps_rare <- rarefy_even_depth(ps, sample.size = chosen_depth, rngseed = 42, replace = FALSE)
alpha <- estimate_richness(ps_rare, measures = c('Observed', 'Shannon', 'InvSimpson'))   # q0, exp gives q1, q2
alpha$Group <- sample_data(ps_rare)$Group
alpha$Shannon_eff <- exp(alpha$Shannon)   # effective species; base-invariant in interpretation (Hill q=1)

kruskal.test(Shannon ~ Group, data = alpha)   # non-parametric; escalate to lme4/nlme for covariates or repeated measures
```

Faith PD in R uses `picante::pd(otu_matrix, tree, include.root = TRUE)`. The Shannon from `estimate_richness` is in natural log (nats); QIIME2 reports log2 (bits) - report `exp(Shannon)` to compare across the two.

## Beta Diversity in R (report weighted AND unweighted)

**Goal:** Quantify between-sample dissimilarity with phylogenetic and abundance-weighted views, then test the group effect while ruling out a dispersion artifact.

**Approach:** Compute both UniFrac variants (and generalized UniFrac alpha=0.5), ordinate by PCoA, run adonis2 for location, and ALWAYS pair it with betadisper for spread.

```r
wu  <- UniFrac(ps_rare, weighted = TRUE)    # abundant-lineage view
uwu <- UniFrac(ps_rare, weighted = FALSE)   # rare-lineage + topology view
# generalized UniFrac alpha=0.5 (Chen 2012 compromise):
gu  <- as.dist(GUniFrac::GUniFrac(t(as(otu_table(ps_rare), 'matrix')), phy_tree(ps_rare), alpha = 0.5)$unifracs[, , 'd_0.5'])

meta <- data.frame(sample_data(ps_rare))
adonis2(wu ~ Group, data = meta, permutations = 999)   # >=999 permutations; significance = LOCATION
permutest(betadisper(wu, meta$Group))                  # MANDATORY: is it dispersion, not location?
```

If betadisper is significant the adonis2 result is ambiguous (location vs spread) - state it. The PERMANOVA-dispersion theory is shared; see metagenomics/metagenome-visualization. For a compositionally coherent ordination with feature loadings use RPCA (DEICODE `qiime deicode rpca` / gemelli). The Python engine is scikit-bio (`skbio.diversity.beta_diversity`, `skbio.stats.ordination.pcoa`, `skbio.stats.distance.permanova`).

## Per-Method Failure Modes

### Sampling-depth sample-massacre
**Trigger:** a `--p-sampling-depth` higher than some samples' totals. **Mechanism:** core-metrics drops every sample below the depth with no warning. **Symptom:** fewer points in the PCoA than samples in the metadata; the lost ones skew low-biomass. **Fix:** pick the depth from the rarefaction plateau, report the dropped-sample list, confirm at a nearby depth.

### De novo tree noise
**Trigger:** UniFrac/Faith PD on a MAFFT+FastTree tree from short reads. **Mechanism:** ~250 bp reads give an unstable topology and arbitrary midpoint root. **Symptom:** unweighted-UniFrac separation that vanishes under SEPP insertion or weighted UniFrac. **Fix:** use SEPP-into-reference or Greengenes2; treat de novo unweighted UniFrac as suspect.

### Unweighted-vs-weighted flip
**Trigger:** reporting only the UniFrac variant that gives p<0.05. **Mechanism:** unweighted listens to rare/short branches, weighted to abundant lineages. **Symptom:** the two disagree and the chosen one is the significant one. **Fix:** report both plus generalized alpha=0.5; state which lineage axis each implicates.

### Rarefy-then-reuse-for-DA
**Trigger:** feeding the rarefied table to a differential-abundance tool. **Mechanism:** rarefaction discards count information the DA model needs. **Symptom:** underpowered or distorted DA. **Fix:** keep raw counts; rarefy only into the diversity branch; route DA to differential-abundance.

### Observed-features-as-species
**Trigger:** comparing raw ASV counts across runs/studies as "richness". **Mechanism:** ASV count tracks DADA2 truncation/maxEE/pooling and intragenomic 16S copy variants, not just biology. **Symptom:** richness shifts with denoising settings. **Fix:** prefer Hill q1/q2; report observed features with the denoising parameters stated.

### Shannon base mismatch
**Trigger:** comparing a QIIME2 Shannon to an R Shannon. **Mechanism:** QIIME2 uses log2 (bits), R `diversity`/`estimate_richness` natural log (nats). **Symptom:** numbers differ by a constant factor and look like a real effect. **Fix:** state the base, convert, or report `exp(H')`.

### PERMANOVA dispersion
**Trigger:** a significant adonis2 read as a composition shift. **Mechanism:** pseudo-F responds to within-group spread, not only centroid location (shared theory; metagenomics/metagenome-visualization). **Symptom:** significant adonis2 with significant betadisper. **Fix:** always run betadisper/permutest alongside; report both.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Sampling depth on the observed-features plateau | Janssen 2018; QIIME2 docs | depth must saturate richness while retaining samples; report dropped list |
| Do NOT use `min(sample_sums)` as the depth | McMurdie 2014 | one tiny library drags every sample to under-saturated noise |
| Generalized UniFrac alpha = 0.5 | Chen 2012 *Bioinformatics* 28:2106 | most powerful for moderately abundant lineages; beats running un/weighted jointly |
| Report Hill q = 0, 1, 2 together | (shared; metagenomics/metagenome-visualization) | spans richness (q0) -> evenness-weighted (q2) |
| PERMANOVA permutations >= 999 | vegan docs | resolution floor for p ~ 0.001; use 9999 for publication |
| Pair adonis2 with betadisper | Anderson & Walsh 2013 (shared) | distinguishes a location shift from a dispersion difference |
| Rarefy for diversity, not for DA | McMurdie 2014; Schloss 2024 | per-analysis decision, not a global switch |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| PCoA has fewer points than samples | `--p-sampling-depth` dropped low-count samples | lower the depth or report the loss; never assume zero drops |
| `UniFrac` errors / Faith PD missing | no `phy_tree` slot in the phyloseq object | attach a SEPP/GG2 (preferred) or de novo tree |
| Unweighted UniFrac significant, weighted not | change is in rare lineages, or de novo tree noise | report both; verify the tree; treat single-metric hit as tentative |
| R and QIIME2 Shannon disagree | log base differs (nats vs bits) | report `exp(H')`; convert by `log2(e)` |
| adonis2 p<0.001 but groups visually overlap | dispersion difference, not location | run betadisper; report it |
| scikit-bio `otu_ids=` deprecation warning | 0.6 renamed OTU to taxon; `otu_ids=` kept as a deprecated alias | `get_beta_diversity_metrics()` and `help()` to find current kwargs |
| Diversity tracks host/plant content | host mitochondria/chloroplast 16S not removed | filter Mitochondria/Chloroplast features (see taxonomy-assignment) before computing diversity |
| "Community" in a near-sterile/low-biomass sample | reagent kitome not removed | sequence controls + run decontam upstream (amplicon-processing; metagenomics/contamination-controls) |

## References

- Faith DP. 1992. Conservation evaluation and phylogenetic diversity. *Biol Conserv* 61:1-10.
- Pielou EC. 1966. The measurement of diversity in different types of biological collections. *J Theor Biol* 13:131-144.
- Lozupone C, Knight R. 2005. UniFrac: a new phylogenetic method for comparing microbial communities. *Appl Environ Microbiol* 71:8228-8235.
- Lozupone CA, Hamady M, Kelley ST, Knight R. 2007. Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. *Appl Environ Microbiol* 73:1576-1585.
- Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD, Collman RG, Bushman FD, Li H. 2012. Associating microbiome composition with environmental covariates using generalized UniFrac distances. *Bioinformatics* 28:2106-2113.
- Janssen S, McDonald D, Gonzalez A, et al. 2018. Phylogenetic placement of exact amplicon sequences improves associations with clinical information. *mSystems* 3:e00021-18.
- Mirarab S, Nguyen N, Warnow T. 2012. SEPP: SATe-enabled phylogenetic placement. *Pac Symp Biocomput* 2012:247-258.
- McDonald D, Jiang Y, Balaban M, et al. 2024. Greengenes2 unifies microbial data in a single reference tree. *Nat Biotechnol* 42:715-718.
- Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A, Knight R, Zengler K. 2019. A novel sparse compositional technique reveals microbial perturbations. *mSystems* 4:e00016-19.
- McDonald D, Vazquez-Baeza Y, Koslicki D, et al. 2018. Striped UniFrac: enabling microbiome analysis at unprecedented scale. *Nat Methods* 15:847-848.
- McMurdie PJ, Holmes S. 2014. Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Comput Biol* 10:e1003531.
- Schloss PD. 2024. Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere* 9:e00354-23.
- McMurdie PJ, Holmes S. 2013. phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS One* 8:e61217.

## Related Skills

- amplicon-processing - Generate the ASV table and representative sequences upstream
- taxonomy-assignment - Label the ASVs summarized here
- differential-abundance - Per-taxon between-group testing on the unrarefied counts
- qiime2-workflow - The QIIME2 CLI home for core-metrics and tree building
- phylogenetics/tree-io - Read, write, and root the UniFrac/Faith PD tree
- metagenomics/abundance-estimation - Shared CoDA and rarefaction-debate theory
- metagenomics/metagenome-visualization - Shared Hill-number and PERMANOVA-dispersion theory; diversity/ordination on shotgun profiler tables
- data-visualization/ggplot2-fundamentals - Custom ordination and diversity plots
