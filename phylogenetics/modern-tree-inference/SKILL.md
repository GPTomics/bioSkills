---
name: bio-phylo-modern-tree-inference
description: Infers maximum-likelihood phylogenetic trees with IQ-TREE2 and RAxML-NG -- model selection (ModelFinder), branch support (UFBoot2, SH-aLRT), concordance factors (gCF/sCF), partitioning, topology tests, and long-branch-attraction control. Covers why an ML tree inherits every flaw of the assumed model and the fixed alignment, why reported support measures repeatability under resampling and not correctness, why UFBoot uses a >=95 cutoff and not the bootstrap-70 rule, and why a node with UFBoot 100 but gCF ~35 is essentially unresolved ILS rather than a clade. Use when inferring an ML tree, selecting a substitution or partition model, choosing or interpreting support measures, testing an a-priori topology, or diagnosing LBA. Routes model-free distance trees to distance-calculations, posteriors to bayesian-inference, and species trees under ILS to species-trees.
tool_type: cli
primary_tool: IQ-TREE2
---

## Version Compatibility

Reference examples tested with: IQ-TREE 2.2+ / 2.3+, RAxML-NG 1.2+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `iqtree2 --version` then `iqtree2 --help` to confirm flags
- CLI: `raxml-ng --version` then `raxml-ng --help` to confirm flags

If code throws an unrecognized-argument or model-parse error, introspect the installed tool and adapt the example to match the actual API rather than retrying.

IQ-TREE2 uses single-dash documented forms (`-alrt`, `-bnni`, `-B`); `-B`/`-T` are v2.x (v1.x used `-bb`/`-nt`). Do NOT write `--alrt`. The likelihood site-concordance flag `--scfl` requires IQ-TREE 2.2.2+ (older builds have only the parsimony `--scf`).

# Modern ML Tree Inference -- ML Support Measures Repeatability, Not Correctness

**"Build a maximum-likelihood tree with support from my alignment"** -> Select a substitution model, search topology and branch lengths that maximize the likelihood, then attach support that quantifies repeatability and concordance that quantifies genealogical agreement.
- CLI: `iqtree2 -s aln.fasta -m MFP -B 1000 -bnni -alrt 1000` (model selection + dual support, all built in)
- CLI: `raxml-ng --all --msa aln.fasta --model GTR+G --bs-metric fbp,tbe` (very large trees, transfer bootstrap, precise branch lengths)

Scope: ML estimation of topology, branch lengths, model selection, branch support, concordance factors, partitioning, topology tests, and LBA control. Model-free distance/NJ trees and distance correction -> distance-calculations. Posterior distributions, MCMC, and CAT-GTR -> bayesian-inference. Per-locus gene trees summarized into a species tree under ILS -> species-trees. Time-scaled trees -> divergence-dating.

## The Single Most Important Modern Insight

An ML tree is the topology and branch lengths that maximize the likelihood under an ASSUMED substitution model, conditioned on a FIXED multiple-sequence alignment. It inherits every flaw of both: a misaligned column is a fabricated character the model dutifully fits, and a misspecified model biases the point estimate toward a wrong topology that more data only sharpens. The reported support measures REPEATABILITY under resampling, not correctness. Three load-bearing facts:

1. **High support is consistent with being wrong.** Bootstrap, UFBoot, SH-aLRT, and aBayes all ask whether a branch reappears when the data or tree is perturbed. Under model misspecification every replicate reproduces the same bias, so support climbs toward 100% precisely as the inference becomes more wrong. At genome scale, sampling error vanishes and UFBoot ~100 on most branches is the default, not a signal.
2. **More data fixes variance, not bias.** Adding sites shrinks sampling error and makes the estimate more confident, but the systematic error from saturation, compositional heterogeneity, and ILS is bias that CONCENTRATES with scale. The cure for confident-wrong trees is better models, better alignments, and better diagnostics -- never more bootstrap replicates.
3. **Concordance factors are the honest measure at genome scale.** gCF/sCF ask what FRACTION of genes or sites actually contains a branch. A node with UFBoot 100 and gCF 35 means the concatenated likelihood is certain but only ~35% of loci endorse that branch -- biologically unresolved (ILS or introgression), and the bootstrap answered a question nobody should have asked. Report CFs on every phylogenomic tree; treat bootstrap as necessary-not-sufficient.

## Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| IQ-TREE2 | Minh 2020 | ML search + ModelFinder + UFBoot2 + SH-aLRT + gCF/sCF + AU test + C60/PMSF, all built in | the default for almost all work |
| RAxML-NG | Kozlov 2019 | ML search, transfer bootstrap, terrace-aware, MPI/checkpointing | very large trees, TBE, precise branch lengths, long HPC runs |
| PhyML | Guindon 2010 | ML search, origin of SH-aLRT | legacy/teaching; SH-aLRT is now in IQ-TREE2 |
| FastTree | Price 2010 | approximate ML, single-pass NNI/SPR | a fast first-pass tree on thousands of sequences; not for final support |

IQ-TREE2 vs RAxML-NG (both hill-climb the same likelihood; they differ in built-in features and scaling, not correctness):

| Need | Use |
|------|-----|
| Model selection, UFBoot2, SH-aLRT, concordance factors, AU test, mixture/PMSF | IQ-TREE2 (built in; nothing else bundles all of this) |
| Very large trees (thousands of taxa), low memory, MPI | RAxML-NG |
| Transfer bootstrap (TBE) for rogue-taxon mega-trees | RAxML-NG (`--bs-metric tbe`) |
| Most precise branch lengths for downstream dating | RAxML-NG |
| Robust checkpoint/restart on long runs | RAxML-NG |

Common production pattern: model-select, compute CFs, and topology-test in IQ-TREE2; do heavy bootstrap on a mega-tree in RAxML-NG with `--bs-metric fbp,tbe`.

## Model Selection

ModelFinder (Kalyaanamoorthy 2017) scores the substitution matrix and the rate-heterogeneity model jointly, including FreeRate categories the old jModelTest/ProtTest generation never tested, and ranks by BIC (the default; k*ln(n) penalty favors simpler models that generalize on large n).

- `-m MFP` -- ModelFinder Plus: test all models by BIC, then search with the winner. The standard. (`-m MF` selects only; avoid `-m TEST`/`-m TESTONLY`, the legacy jModelTest-style limited set.)
- **+G vs +R.** `+G` (discrete Gamma, one alpha) is the workhorse for single short genes. `+R` (FreeRate) freely estimates each category's rate AND weight, capturing the non-Gamma, often multimodal rate distributions of concatenated phylogenomic data; expect `+R3` to `+R6` to win on BIC at scale.
- **The +I+G trap.** The invariant-sites proportion (+I) and the Gamma shape (+G) describe overlapping rate distributions, so the likelihood surface has a flat ridge: the two estimates are individually near-meaningless and start-dependent. Prefer `+R` (its slowest category absorbs near-invariant sites); use `+G` alone unless BIC genuinely demands +I+G.
- **Partition model selection.** `-m MFP+MERGE` fits per-partition models AND greedily merges partitions that fit the same model (BIC-chosen scheme, the successor to PartitionFinder greedy). Pair with `-rcluster 10` (relaxed clustering: only test the top 10% most-similar pairs) for many partitions.
- **Site-heterogeneous mixtures for deep data.** Empirical matrices (LG, WAG) assume one residue-frequency vector for the whole alignment; real proteins do not, and that across-site compositional heterogeneity is the chief driver of deep LBA. The ML answer is the C10..C60 profile-mixture series (`LG+C60+F+G`), made tractable by PMSF (Wang 2018): a guide-tree pass computes one posterior-mean profile per site, and the real search uses those frozen profiles. PMSF is the standard recommendation for deep / LBA-prone protein phylogenomics.

## Branch Support -- the Heart

Five measures, three different questions; the cardinal sin is cross-comparing their cutoffs.

| Metric | What it perturbs / measures | Strong cutoff | Tool / flag | Failure mode |
|--------|------------------------------|---------------|-------------|--------------|
| Standard bootstrap (FBP) | resample sites; clade frequency (binary) | >=70 (folklore) | RAxML-NG `--bs-metric fbp`; IQ-TREE `-b` | slow; crushed by rogue taxa in big trees |
| Ultrafast bootstrap 2 (UFBoot) | RELL-resampled log-Ls; ~unbiased clade prob | **>=95 (NOT 70)** | IQ-TREE2 `-B 1000` (+`-bnni`) | inflates under model violation -> use `-bnni` |
| SH-aLRT | local NNI likelihood ratio (no resampling) | >=80 | IQ-TREE2 `-alrt 1000` | conservative on very short branches |
| aBayes | posterior from 3 NNI Ls, flat prior | >=0.95 | IQ-TREE2 `-abayes` | anti-conservative; never the sole criterion |
| Transfer bootstrap (TBE) | gradual transfer distance under resampling | no fixed cutoff; > FBP | RAxML-NG `--bs-metric tbe` | permissive; "fuzzy" branch identity |

UFBoot2 (Hoang 2018) uses the RELL trick (resample site log-likelihoods, reuse a candidate tree set) to run hundreds of times faster than the Felsenstein bootstrap (1985), and its values are CLOSER to unbiased clade probabilities -- which is exactly why the strong-support cutoff is 95, not the conservative-bootstrap 70. Treating UFBoot 70 as "good" is a category error. `-bnni` re-optimizes each replicate tree by NNI to rein in the inflation that model violation causes; use `-B 1000 -bnni` routinely. UFBoot is on a DIFFERENT scale from the standard bootstrap -- never read it with the BP-70 rule.

SH-aLRT (Guindon 2010) does not resample data; for each branch it tests whether the ML likelihood beats its two best NNI rearrangements. UFBoot (data perturbation) and SH-aLRT (tree perturbation) have different failure modes, so the community-standard joint criterion requires both:

> A branch is strongly supported iff SH-aLRT >= 80% AND UFBoot >= 95%.

For large rogue-taxon-prone trees, the binary Felsenstein bootstrap lets a single wandering tip crush an otherwise-recovered deep branch; transfer bootstrap (TBE, Lemoine 2018, RAxML-NG `--bs-metric tbe`) replaces the in/out indicator with a gradual transfer distance and rescues those branches, at the cost of being more permissive.

## Concordance Factors

Bootstrap quantifies statistical confidence given the concatenated data; concordance factors (Minh 2020) quantify how much of the actual data carries a branch.

- **gCF (gene concordance factor):** the percentage of decisive single-locus gene trees that contain the exact branch. Needs per-locus gene trees.
- **sCF (site concordance factor):** the percentage of decisive sites supporting the branch, from sampled quartets; works on a single concatenated alignment with no gene trees.
- **sCFL (likelihood sCF, Mo 2023):** uses ancestral-state likelihoods rather than parsimony quartet counting, substantially reducing (not abolishing) homoplasy and taxon-sampling bias. Prefer `--scfl` over the old `--scf` on IQ-TREE 2.2.2+.

```bash
# one gene tree per locus from a directory of locus alignments (-S = separate, no concatenation)
iqtree2 -S loci_dir -m MFP -B 1000 -T AUTO --prefix loci   # -B 1000 = UFBoot per gene tree, needed to contract weak branches before ASTRAL

# gene + likelihood site concordance against a fixed concatenated tree (-te fixes the tree)
iqtree2 -te concat.treefile -s concat.fasta --gcf loci.treefile --scfl 100 -T 4 --prefix concord
#  --gcf loci.treefile   per-locus gene trees for gCF
#  --scfl 100            100 sampled quartets per branch for likelihood sCF (higher = more stable)
```

Outputs `concord.cf.tree` (Newick with gCF/sCF labels) and `concord.cf.stat` (per-branch gCF, gDF1, gDF2, gDFP, sCF). A node with UFBoot 100 but gCF ~35 (gDF1 ~33, gDF2 ~30) is genes split three ways: the concatenated point estimate barely edges the alternatives and the node is biologically unresolved -- the signature of ILS or introgression, not a clade. Report CFs alongside support on every phylogenomic tree.

## Topology Tests

For testing an a-priori hypothesis ("can I reject that X and Y are monophyletic?") against the ML tree, by comparing a set of fixed trees. Build the constrained tree with `-g constraint.tree`, then evaluate both trees:

```bash
# trees.nex holds the unconstrained ML tree + the constrained/alternative trees
iqtree2 -s aln.fasta -m <model> -z trees.nex -n 0 -zb 10000 -au --prefix autest
#  -z trees.nex   trees to compare        -n 0   no fresh search, just evaluate
#  -zb 10000      RELL replicates (>=1000) -au    add the AU test (must accompany -zb)
```

The AU test (Shimodaira 2002) uses multiscale bootstrap resampling to correct both the selection bias of the KH test (invalid on the data-selected ML tree) and the over-conservatism of the SH test (which rejects less as the candidate set is padded). Use AU by default. **p-AU < 0.05 means that tree is REJECTED**; p-AU >= 0.05 means it is in the 95% confidence set (failure to reject is not acceptance -- weak data fails to reject many trees). Report SH/KH only for completeness.

## Partitioning

When splitting an alignment into partitions, the branch-length linkage choice is what people get wrong:

| Mode | Flag | Branch lengths | Use when |
|------|------|----------------|----------|
| Edge-equal | `-q part.nex` | identical across partitions | partitions share rate (rare, restrictive) |
| Edge-linked proportional | `-p part.nex` | shared topology, per-partition rate multiplier | DEFAULT -- genes evolve at different speeds, share history |
| Edge-unlinked | `-Q part.nex` | fully independent per partition | genuine heterotachy; parameter-hungry, overfits |

`-p` (edge-linked proportional) is the standard: one rate multiplier per partition over a shared topology. Over-partitioning spends degrees of freedom without bias reduction and inflates variance; the antidote is to start fine (gene x codon position) and let `-m MFP+MERGE -rcluster 10` find the coarsest BIC-justified scheme. Prefer a merged scheme over a hand-picked maximal one; prefer `-p` over `-Q`.

## Per-Method Failure Modes

### Long-Branch Attraction
**Trigger:** Two or more independently fast-evolving lineages on long branches separated by a short internode.
**Mechanism:** Convergent/homoplastic substitutions on the long branches look like shared ancestry; a site-homogeneous model cannot separate convergence from homology and groups them -- bias that GROWS with more sites.
**Symptom:** Fast taxa group with the outgroup or each other at 100% bootstrap; the grouping collapses under a better model or when a long-branch taxon is removed.
**Fix:** Site-heterogeneous model (C60/PMSF) first; remove the fastest sites and watch the node; drop the long-branch taxon or use a closer outgroup; cross-check with SR4/Dayhoff recoding. Believe a deep node only when it survives all of these, not when it merely has UFBoot 100 under LG+G.

### Model Underspecification
**Trigger:** A single inadequate model on heterogeneous data; the best site-homogeneous model by BIC still inadequate at depth.
**Mechanism:** Wrong matrix, missing rate heterogeneity, or no partitioning biases the topology while support stays high.
**Symptom:** Biologically implausible nodes with full support that move under a richer model.
**Fix:** `-m MFP` (+MERGE for multi-locus), FreeRate `+R`, and at amino-acid depth a C60/PMSF mixture. Best-by-BIC among site-homogeneous models is not sufficient deep in the tree.

### Over-Partitioning
**Trigger:** Hundreds of hand-defined partitions, especially under `-Q`.
**Mechanism:** Each partition's model is estimated from too little data; parameter and branch-length estimates get noisy without reducing bias.
**Symptom:** Slow runs, noisy estimates, degraded support.
**Fix:** `-m MFP+MERGE -rcluster 10` to the coarsest BIC-justified scheme; use `-p`, not `-Q`.

### The Support-Accuracy Gap
**Trigger:** UFBoot/bootstrap ~100 everywhere, including implausible or conflicting nodes.
**Mechanism:** Support measures repeatability of a possibly-biased estimate; concatenation pools conflicting gene signals so the likelihood is certain while the loci disagree.
**Symptom:** Full support but gCF ~33 / sCF near its ~33% floor on contested nodes.
**Fix:** Compute gCF/sCFL; treat UFBoot 100 + gCF ~33 as UNRESOLVED; require SH-aLRT >=80 AND UFBoot >=95; use `-bnni`. If most genes reject the ML resolution, route to a coalescent species tree -> species-trees.

## Quantitative Thresholds

| Quantity | Threshold | Source |
|----------|-----------|--------|
| UFBoot strong support | >=95 (NOT the bootstrap 70) | Hoang 2018 |
| SH-aLRT strong support | >=80 | Guindon 2010 |
| Joint rule | SH-aLRT >=80 AND UFBoot >=95 | community standard |
| Standard bootstrap "good" | >=70 (different metric, folklore) | Felsenstein 1985 / Hillis 1993 |
| aBayes strong | >=0.95 (anti-conservative; never alone) | Anisimova 2011 |
| Bootstrap replicates | UFBoot `-B` >=1000; SH-aLRT `-alrt` >=1000; AU `-zb` 10000 | IQ-TREE docs |
| gCF reading | >~75 agree; ~50 conflicted; ~33 effective polytomy; <33 with higher alternative = possible wrong resolution | Minh 2020 |
| sCF floor | ~33% (three quartet resolutions); ~33 = no site signal | Minh 2020 |
| AU test | p-AU < 0.05 => tree REJECTED | Shimodaira 2002 |
| Model selection | rank by BIC (`-m MFP` default), k*ln(n) penalty | Kalyaanamoorthy 2017 |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Unknown argument --alrt` | wrote the GNU double-dash form | IQ-TREE2 uses single-dash `-alrt`, `-bnni`, `-B` |
| `-bb` / `-nt` not recognized | v1.x flags on a v2.x binary | use `-B` (bootstrap) and `-T` (threads) in 2.x |
| Reading UFBoot 80 as "supported" | applied the bootstrap-70 rule to a different scale | use UFBoot >=95 AND SH-aLRT >=80 |
| Fully-supported deep node distrusted by reviewer | no concordance factors reported | compute gCF/sCFL; treat high-support/low-CF as unresolved |
| `--scfl` unrecognized | IQ-TREE older than 2.2.2 | upgrade, or fall back to parsimony `--scf` |
| Concatenated tree confidently wrong on a rapid radiation | ILS; concatenation is inconsistent in the anomaly zone | infer per-locus gene trees and a coalescent species tree -> species-trees |
| AU test "fails to reject" the alternative | weak data, or candidate set padded | do not pad the set; failure to reject is not acceptance |

## References

Felsenstein J. 1978. Cases in which parsimony or compatibility methods will be positively misleading. *Systematic Zoology* 27(4):401-410.
Felsenstein J. 1985. Confidence limits on phylogenies: an approach using the bootstrap. *Evolution* 39(4):783-791.
Guindon S, Dufayard J-F, Lefort V, Anisimova M, Hordijk W, Gascuel O. 2010. New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. *Systematic Biology* 59(3):307-321.
Price MN, Dehal PS, Arkin AP. 2010. FastTree 2: approximately maximum-likelihood trees for large alignments. *PLoS ONE* 5(3):e9490.
Anisimova M, Gil M, Dufayard J-F, Dessimoz C, Gascuel O. 2011. Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. *Systematic Biology* 60(5):685-699.
Shimodaira H. 2002. An approximately unbiased test of phylogenetic tree selection. *Systematic Biology* 51(3):492-508.
Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14(6):587-589.
Wang H-C, Minh BQ, Susko E, Roger AJ. 2018. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. *Systematic Biology* 67(2):216-235.
Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. 2018. UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35(2):518-522.
Lemoine F, Domelevo Entfellner J-B, Wilkinson E, Correia D, Davila Felipe M, De Oliveira T, Gascuel O. 2018. Renewing Felsenstein's phylogenetic bootstrap in the era of big data. *Nature* 556(7702):452-456.
Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A. 2019. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics* 35(21):4453-4455.
Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution* 37(5):1530-1534.
Minh BQ, Hahn MW, Lanfear R. 2020. New methods to calculate concordance factors for phylogenomic datasets. *Molecular Biology and Evolution* 37(9):2727-2733.
Mo YK, Lanfear R, Hahn MW, Minh BQ. 2023. Updated site concordance factors minimize effects of homoplasy and taxon sampling. *Bioinformatics* 39(1):btac741.

## Related Skills

- distance-calculations - model-corrected distances and fast NJ trees as a model-free alternative
- bayesian-inference - posteriors, MCMC convergence, and CAT-GTR site-heterogeneous models
- species-trees - coalescent species-tree estimation when concordance factors reveal ILS
- divergence-dating - time-scaled trees from the ML topology
- tree-manipulation - rooting, pruning, and collapsing low-support nodes
- alignment/alignment-io - the alignment whose homology assumption the ML tree trusts as fixed
