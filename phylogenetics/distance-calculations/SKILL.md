---
name: bio-phylo-distance-calculations
description: Build model-corrected evolutionary distance matrices and distance trees (NJ, BIONJ, FastME, UPGMA) with Biopython Bio.Phylo plus R ape/phangorn/FastME. Covers why a distance is a model-corrected estimate of substitutions per site that undercounts raw because of multiple/back/parallel hits (saturation); why the matrix discards the per-site information ML keeps; the LogDet/paralinear fix for compositional heterogeneity; the UPGMA molecular-clock trap; and the Bio.Phylo landmine that DistanceCalculator offers only identity/matrix distances, not JC/K80/TN93. Use when computing a distance matrix, building a fast NJ/FastME tree, seeding an ML search, barcoding, or testing substitution saturation before a deep tree. Routes ML and starting-tree work to modern-tree-inference, alignment quality to alignment/alignment-io, and tree I/O to tree-io.
tool_type: mixed
primary_tool: Bio.Phylo.TreeConstruction
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+. Model-corrected alternatives: ape 5.8+, phangorn, FastME, scikit-bio.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(module.function)` to check signatures
- R: `packageVersion('ape')` then `?dist.dna` to verify model strings
- CLI: `fastme --version` then `fastme --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Bio.Phylo `DistanceCalculator` offers only identity / substitution-matrix distances and has NO model-based JC/K80/TN93; for model-corrected DNA distances use ape `dist.dna(model=...)` (introspect its model list and `gamma=` support if the version differs).

# Distance Calculations -- Distance Erases the Information ML Keeps; the Correction Is the Model

**"Build a tree from my alignment"** -> Estimate model-corrected pairwise distances, then cluster them into a tree, knowing the correction is where the biology lives and the matrix discards what ML uses.
- Python: `DistanceCalculator(...).get_distance(aln)` then `DistanceTreeConstructor().nj(dm)` (Bio.Phylo, identity-only)
- R: `dist.dna(x, model='TN93')` then `fastme.bal(d)` / `nj(d)` (ape, model-corrected)

Scope: computing evolutionary distance matrices and building/bootstrapping distance trees (NJ, BIONJ, FastME, UPGMA), and testing saturation before trusting them. ML inference, model selection, and where NJ seeds an ML search -> modern-tree-inference. Alignment quality and trimming, which gate every distance -> alignment/alignment-io. Reading/writing/converting the trees these methods emit -> tree-io. Rooting and pruning -> tree-manipulation.

## The Single Most Important Modern Insight

A distance-based phylogeny is built in two lossy steps and both losses are silent. First, a distance is a model-corrected estimate of expected substitutions per site, not the observed proportion of differing sites; second, the matrix collapses each sequence pair to one scalar and discards the site-pattern information that maximum likelihood evaluates. Distance methods are therefore fast but doubly limited -- by how good the correction is, and by what a matrix of scalars can encode -- and they are structurally incapable of using the per-site signal that makes ML the gold standard for hard problems. Three load-bearing facts:

1. **The correction is the whole biological content.** The raw p-distance (proportion of differing sites) systematically UNDERCOUNTS true divergence because every site that mutated twice (back-substitution), in parallel on both lineages, or to a third state (multiple hit) is missed, and the gap grows without bound toward saturation. The algorithm downstream (NJ, FastME) is just arithmetic on the numbers handed to it; feed it p-distances on divergent data and it returns a confidently wrong tree, fast.
2. **The matrix throws away what ML keeps.** Two completely different alignments can produce the same distance matrix and therefore the same tree. ML never collapses the data, so it cannot be fooled this way; once the matrix exists, which sites drove the divergence, among-site rate heterogeneity, and competing site-pattern signal are gone.
3. **Consistency is conditional, not a guarantee.** NJ is statistically consistent -- it returns the true tree as sequence length grows -- but only if the input distances are additive and correctly estimated. A consistent algorithm fed biased (saturated, misspecified) distances converges confidently on the WRONG tree. "I used NJ so it is consistent" is not a defense; the correction is exactly what most users get wrong.

## Distance Corrections

Every nucleotide correction inverts the probability that a site differs under a model, so as p approaches the equilibrium difference (0.75 for equal base frequencies under JC) the correction diverges to infinity and its variance explodes -- the mathematical face of saturation.

| Correction | Corrects for | When to use | ape model |
|---|---|---|---|
| p-distance (raw / Hamming) | nothing | very shallow / barcoding; saturation plots | `'raw'`, `'N'` |
| Jukes-Cantor (JC69) | multiple hits, one rate | minimal correction, quick sanity tree | `'JC69'` |
| Kimura 2-param (K80) | + transition/transversion bias | default DNA quick distance | `'K80'` (default) |
| Tamura-Nei (TN93) | + two transition rates + unequal base freqs | mtDNA, richest closed-form | `'TN93'` |
| LogDet / paralinear | COMPOSITIONAL heterogeneity (non-stationary, general Markov) | GC drift across taxa; compositional attraction suspected | `'logdet'`, `'paralin'` |
| Gamma-corrected | among-site rate variation (ASRV) | alpha small (< ~1); long branches | `dist.dna(..., gamma=alpha)` |
| Protein matrix (LG/WAG/JTT/PAM) | amino-acid replacement, multiple hits | protein data (use phangorn `dist.ml`) | n/a (LG/WAG modern; PAM legacy) |

LogDet / paralinear (Lockhart et al. 1994; Lake 1994) is the non-obvious one: JC/K80/TN93 all assume stationarity (constant base composition across the tree), and when composition drifts (thermophiles, AT-rich insect mtDNA, GC-rich chloroplasts) every stationary correction groups taxa by base composition rather than ancestry -- compositional attraction, a sibling of long-branch attraction that fools NJ and stationary-model ML alike. LogDet is computed from the determinant of the 4x4 pairwise divergence matrix and is consistent under the general Markov model, at the cost of needing long sequences and not easily taking gamma. If a topology flips between JC and LogDet, suspect compositional heterogeneity.

## Algorithm Decision

| Algorithm | Criterion | Clock? | Cost | Use / pitfall |
|---|---|---|---|---|
| UPGMA | average-linkage clustering | ASSUMES strict clock (ultrametric) | O(n^2) | avoid for molecular phylo; wrong TOPOLOGY under rate variation; serotyping/dendrograms only |
| NJ (Saitou & Nei 1987) | greedy balanced minimum evolution (Q-matrix) | no | O(n^3) | fast default; consistent ONLY if distances correct/additive; barcoding, large n |
| BIONJ (Gascuel 1997) | NJ + variance-weighted reduction | no | O(n^3) | better than NJ on large/noisy distances; a common ML starting-tree generator (IQ-TREE2/RAxML-NG default to parsimony starts) |
| FastME (Lefort et al. 2015) | balanced minimum evolution + NNI/SPR search | no | ~NJ-speed | the modern best distance tree; actively searches, not a single greedy pass |

UPGMA is the key warning: it forces every tip equidistant from the root, which molecular data essentially never satisfies because lineages evolve at different rates. A fast-evolving lineage gets pulled toward the tips and grouped by total divergence rather than true ancestry -- a clustering-by-rate artifact analogous to long-branch attraction. Flat rule: do not use UPGMA for molecular phylogeny unless a clock is independently established. NJ does not assume a clock because its Q-matrix corrects each pairwise distance for that taxon's average divergence to everyone else. NJ/BIONJ/FastME can seed ML searches: some ML programs build a BIONJ or ME starting tree in O(n^3) before likelihood NNI/SPR (though IQ-TREE2 and RAxML-NG default to parsimony starting trees), so this machinery still has a place in an all-ML pipeline.

## Tool Taxonomy

| Tool (lang) | Distance step | Tree step | Note |
|---|---|---|---|
| Bio.Phylo.TreeConstruction (Py) | `DistanceCalculator` identity / BLOSUM / PAM matrices ONLY | `DistanceTreeConstructor().nj()` / `.upgma()` | NO JC/K80/TN93; `'identity'` is p-distance, named matrices are score distances |
| ape (R) | `dist.dna(model=, gamma=)` -- full menu incl. LogDet, gamma | `nj()`, `bionj()`, `fastme.bal()`, `fastme.ols()` | the reference engine for real corrections |
| phangorn (R) | `dist.ml(model=)` -- DNA JC69/F81 only, but all protein matrices | `upgma()`, `NJ()` | the natural choice for protein distances |
| FastME (CLI) | JC/K2P/F84/TN93/LogDet/protein | NJ/BIONJ/BME/OLS-ME + NNI/SPR + bootstrap | standalone for very large n; `-m B` = balanced ME |
| scikit-bio (Py) | none (requires a pre-corrected matrix) | `skbio.tree.nj(dm)` | pure-Python NJ; `neg_as_zero=` since 0.6.3 |

Loud flag, the second flat rule: Bio.Phylo and scikit-bio do NOT correct for multiple hits -- their built-in distance is p-distance or score-based. Presenting Bio.Phylo's `'identity'` distance as a "Jukes-Cantor tree" is a common and wrong shortcut. For a real correction, compute the matrix in ape/FastME (or compute it yourself) and pass it in.

## Build a Model-Corrected Distance Matrix and NJ Tree

**Goal:** Turn an alignment into a tree, contrasting Bio.Phylo's identity-only path against ape's model-corrected path.

**Approach:** In Python use Bio.Phylo for the NJ/UPGMA algorithm and pure-Python I/O, but treat its distance as an uncorrected p-distance; for any divergent data move the distance step to ape `dist.dna`, then cluster with FastME (best) or NJ.

```python
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

aln = AlignIO.read('alignment.fasta', 'fasta')
calc = DistanceCalculator('identity')        # identity-only: this is a p-distance, NOT a JC/K80 correction
dm = calc.get_distance(aln)                   # multiple/back/parallel hits are NOT corrected here
tree = DistanceTreeConstructor().nj(dm)       # NJ algorithm is correct; the DISTANCES are the limitation
# For a model-corrected distance, build the matrix in ape (below) and pass it to skbio/Bio.Phylo, or stay in R.
```

```r
library(ape)
aln <- read.dna('alignment.fasta', format = 'fasta')
d <- dist.dna(aln, model = 'TN93', gamma = 0.5)   # model-corrected; gamma applies ASRV (alpha < 1 = strong)
tree <- fastme.bal(d, nni = TRUE, spr = TRUE)     # balanced minimum evolution: the modern best distance tree
# nj(d) / bionj(d) are the faster single-pass alternatives; bionj seeds ML searches
```

## Pre-flight: Test Substitution Saturation

**Goal:** Decide whether the data retain phylogenetic signal before trusting any deep distance tree.

**Approach:** Compute Xia's index of substitution saturation Iss and compare it to the simulation-derived critical value Iss.c; a saturation plot of transitions against a corrected distance is the visual companion.

```r
library(ape)                    # the formal entropy-based Iss vs Iss.c test lives in DAMBE; this is the ape/base-R saturation plot
d_jc  <- dist.dna(aln, model = 'JC69')
ts_tv <- dist.dna(aln, model = 'TS')   # transitions; plot against d_jc -- a PLATEAU means saturated, signal erased
plot(d_jc, ts_tv)                       # unsaturated = roughly linear; bent-over transition curve = drop those sites
# Interpretation gate (Xia 2003): Iss < Iss.c => signal retained (usable); Iss >= Iss.c => substantially saturated, do not use.
```

## Bootstrap a Distance Tree

**Goal:** Quantify how reproducible each clade is under resampling -- precision, not accuracy.

**Approach:** Resample alignment columns with replacement, recompute the matrix with the SAME correction, rebuild with the SAME algorithm, and summarize clade frequencies on a majority-rule consensus.

```r
library(ape)
boot <- boot.phylo(tree, aln, function(x) fastme.bal(dist.dna(x, model = 'TN93')), B = 500)
# 100-1000 replicates standard; support measures SAMPLING stability only -- it cannot detect a bias in the distances
```

In Python, Bio.Phylo's `bootstrap_consensus(aln, 100, DistanceTreeConstructor(calc, 'nj'), majority_consensus)` does the same on the identity distance (same correction caveat).

## When Distance Is Legitimate vs a Trap

Legitimate or preferred: a quick exploratory / sanity tree before a long ML run; very large n (thousands+ tips) where ML is infeasible (large-scale barcoding, OTU/pangenome trees); barcoding and population-level shallow data, where saturation is negligible and the per-site information ML keeps adds little (NJ on K2P distances is the literal DNA-barcoding standard); and as the starting tree for ML/Bayesian search.

A trap, do not: publication-grade deep phylogeny or formal hypothesis testing (dating, selection, contested deep nodes) -- use ML or Bayesian with model selection; any dataset that fails a saturation test (Xia Iss >= Iss.c), which no algorithm rescues; data with strong compositional heterogeneity unless using LogDet/paralinear; and mistaking consistency for a guarantee -- NJ is consistent only with correct distances, and saturation/misspecification break it exactly as they break ML.

## Per-Method Failure Modes

### Saturation Makes Distances Plateau and Mislead at Depth
**Trigger:** Deep divergences with many pairwise p-distances above ~0.5; transitions exhausted while transversions still climb.
**Mechanism:** So many superimposed substitutions accumulate that observed differences approach the random expectation; the correction inflates violently near its singularity and carries no remaining signal about deep splits.
**Symptom:** Corrected distances explode and become unstable; the saturation plot's transition curve flattens; deep nodes are unstable across models.
**Fix:** Run the Xia Iss test; exclude transitions, third codon positions, or saturated partitions; do not build a deep distance (or ML) tree on saturated data.

### UPGMA Returns the Wrong Topology Under Rate Variation
**Trigger:** Lineages evolving at different rates analyzed with UPGMA.
**Mechanism:** UPGMA forces an ultrametric tree (clock), so a fast-evolving lineage is pushed toward the tips and grouped by total divergence, not ancestry.
**Symptom:** Topology differs from an NJ tree on the same matrix; fast lineages cluster together.
**Fix:** Use NJ or FastME; never UPGMA for molecular phylogeny without an independently established clock.

### Identity / Uncorrected Distance Model on Divergent Data
**Trigger:** Bio.Phylo `'identity'` or ape `'raw'` on anything beyond shallow divergence; calling it a "JC tree".
**Mechanism:** No multiple-hit correction, so divergence is undercounted nonlinearly and worst on the longest branches, amplifying long-branch artifacts.
**Symptom:** Branch lengths and topology shift when a real correction is applied; long-branch taxa attract.
**Fix:** Use ape `dist.dna(model='TN93')` or a gamma correction; for compositional skew use LogDet; for protein use phangorn `dist.ml(model='LG')`.

### Bootstrap on a Method That Discarded the Per-Site Information
**Trigger:** High distance-NJ bootstrap on saturated or compositionally biased data.
**Mechanism:** Resampling removes sampling noise, not a systematic bias baked into the distances, so the wrong split reproduces every replicate.
**Symptom:** Confident, reproducible support on a clade that moves when the correction or saturated sites change.
**Fix:** Bootstrap measures precision, not accuracy; fix the distances (model + saturation) first.

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|---|---|---|
| p-distance saturation onset | unstable as nucleotide p -> ~0.5-0.6; singularity at 0.75 (equal base freqs) | many pairwise p > ~0.5 = red flag, run a saturation test |
| Xia substitution-saturation test | Iss < Iss.c => usable; Iss >= Iss.c (esp. asymmetric Iss.c, the stricter bar) => do not use | Xia et al. 2003 (DAMBE) |
| JC/correction singularity | log argument -> 0 as p -> 0.75 | d and its variance diverge to infinity |
| Gamma shape | alpha < ~1 = strong ASRV, materially changes distances/topology | ignoring ASRV inflates long-branch artifacts |
| ts/tv ratio | commonly 2-15 (~15 primate mtDNA control region) | large ratios make K80/TN93 over JC matter (Tamura & Nei 1993) |
| Bootstrap replicates | 100-1000; clades < ~70% conventionally unsupported | precision not truth (Hillis & Bull 1993) |

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| NJ tree on `model='identity'` called a Jukes-Cantor tree | DistanceCalculator does NOT do JC/K80/TN93 | compute the corrected matrix in ape/FastME, pass it as a matrix |
| Corrected distances explode / `NaN` | p near or above 0.75 singularity (saturation) | test saturation; drop saturated sites; do not trust deep distances |
| UPGMA tree disagrees with NJ | UPGMA clock assumption violated by rate variation | use NJ/FastME for molecular data |
| Unrelated GC-rich taxa group together; topology flips JC vs LogDet | non-stationary composition under a stationary correction | use LogDet/paralinear or a non-stationary ML model |
| LogDet returns `NaN` / undefined | short sequences drive det(F) <= 0 | use longer alignments or a stationary correction |
| High bootstrap on a wrong clade | resampling masks a systematic bias in the distances | fix the model/saturation; support is precision, not accuracy |

## References

Saitou N, Nei M. 1987. The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution* 4(4):406-425.
Gascuel O. 1997. BIONJ: an improved version of the NJ algorithm based on a simple model of sequence data. *Molecular Biology and Evolution* 14(7):685-695.
Lefort V, Desper R, Gascuel O. 2015. FastME 2.0: a comprehensive, accurate, and fast distance-based phylogeny inference program. *Molecular Biology and Evolution* 32(10):2798-2800.
Xia X, Xie Z, Salemi M, Chen L, Wang Y. 2003. An index of substitution saturation and its application. *Molecular Phylogenetics and Evolution* 26(1):1-7.
Lockhart PJ, Steel MA, Hendy MD, Penny D. 1994. Recovering evolutionary trees under a more realistic model of sequence evolution. *Molecular Biology and Evolution* 11(4):605-612.
Lake JA. 1994. Reconstructing evolutionary trees from DNA and protein sequences: paralinear distances. *PNAS* 91(4):1455-1459.
Tamura K, Nei M. 1993. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees. *Molecular Biology and Evolution* 10(3):512-526.
Schliep KP. 2011. phangorn: phylogenetic analysis in R. *Bioinformatics* 27(4):592-593.
Hillis DM, Bull JJ. 1993. An empirical test of bootstrapping as a method for assessing confidence in phylogenetic analysis. *Systematic Biology* 42(2):182-192.

## Related Skills

- modern-tree-inference - ML inference, model selection, and where a BIONJ/ME distance tree seeds the ML search
- tree-manipulation - rooting and pruning the unrooted trees NJ/FastME emit
- tree-io - reading, writing, and converting the trees these methods produce
- alignment/alignment-io - the alignment whose quality gates every distance; filter ambiguous blocks first
