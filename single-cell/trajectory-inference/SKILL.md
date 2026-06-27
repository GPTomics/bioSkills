---
name: bio-single-cell-trajectory-inference
description: Infers developmental trajectories, pseudotime, RNA velocity, and directed fate probabilities from single-cell data using PAGA, Slingshot, Monocle3, DPT, Palantir, scVelo, and CellRank 2. Use when ordering cells along a differentiation continuum, choosing a trajectory method by topology, rooting pseudotime, estimating RNA velocity direction, computing fate probabilities near a bifurcation, or judging whether an inferred trajectory is real.
tool_type: mixed
primary_tool: Monocle3
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, scVelo 0.3+, CellRank 2.0+, Monocle3 1.3+, Slingshot 2.x

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Trajectory Inference

**"Find the developmental trajectory in my data"** -> Order cells along a continuous manifold and assign a pseudotime, fate probabilities, or velocity direction.
- Python: `sc.tl.paga`/`sc.tl.dpt` (scanpy), `palantir`, `scvelo`, `cellrank` (kernels + GPCCA)
- R: `slingshot`, `monocle3`, `tradeSeq`

## Governing Principle

A snapshot is not a movie. Every method substitutes transcriptomic similarity for temporal adjacency: it takes a single static sample and imposes an ordering or a Markov process on the kNN graph. The output is meaningful only when the population is a genuine continuum of asynchronously-progressing cells, sampled densely enough to bridge intermediate states. Five rules follow and they drive every downstream decision.

1. Pseudotime is geometry, not a clock. The axis is monotone in progression at best; equal pseudotime intervals are NOT equal real-time intervals, and rates differ across lineages. Reading an interval as a duration is a category error.
2. The existence of a continuum is a BIOLOGICAL judgment made BEFORE ordering. Any method returns numbers on discrete cell types or on noise; no algorithm tests whether a continuum exists. Decide topology first with PAGA connectivity, then order.
3. Root-cell choice flips every gene trend. Pseudotime is defined only up to an origin; the sign of every trend and which cells are "early" invert with the root. Anchor the root with orthogonal evidence (known marker, real sampling time, velocity, or a stemness score), never by eye on a UMAP, which distorts global geometry. When no progenitor marker and no real timepoints are available and velocity is invalid (mature or non-cycling tissue), fall back to a model-free stemness score (CytoTRACE/CytoTRACE2) to nominate the root, and treat the resulting direction as a hypothesis rather than an established origin.
4. Near bifurcations use fate PROBABILITIES, not hard branch labels. A multipotent progenitor's fate is genuinely undetermined, so a hard assignment to one lineage is biologically false. Represent each cell as a distribution over terminal fates (Palantir, CellRank).
5. Transcriptomic state does not fully predict fate (Weinreb 2020). Sister cells in an indistinguishable state systematically diverge in fate, so a state-based branch call can be systematically wrong at exactly the decision point it most wants to resolve. Frame fate calls as predictions, not measurements, and validate with orthogonal evidence.

Snapshot dynamics are formally non-identifiable (Weinreb 2018, gauge freedom): a single snapshot constrains a family of dynamics, not one, and a unique answer requires extra assumptions (a potential field, known birth-death rates, or real timepoints). Treat any single-method pseudotime as a hypothesis to cross-validate, never a measurement.

## Method Decision Table

Saelens 2019 benchmarked 45 methods across 110 real + 229 synthetic datasets: no single method wins across all topologies, so selection is topology-first.

| Method | Model / assumption | Use when | Fails when |
|--------|--------------------|----------|------------|
| PAGA (Wolf 2019) | cluster-graph connectivity = observed vs expected inter-cluster edges | deciding IF a continuum exists; unknown, disconnected, or cyclic topology; the mandatory first step | gives connectivity not pseudotime/direction; threshold is manual; resolution-dependent |
| Slingshot (Street 2018) | MST on cluster centroids + simultaneous principal curves | known tree/bifurcating topology; smooth per-lineage curves; tradeSeq DE | depends entirely on input clustering; no cycles/disconnection; scales poorly |
| Monocle3 (Cao 2019) | principal graph (reversed graph embedding) in UMAP space | tree, cyclic, or disconnected topology without pre-clustered lineages; Moran's-I trajectory DE | graph learned in UMAP inherits global distortion; resolution/seed-sensitive |
| DPT (Haghverdi 2016) | reversible diffusion random walk; diffusion distance from root | linear or simple-branching; fast native-scanpy scalar after PAGA fixes topology | undirected (reversible model for an irreversible process); root-sensitive; weak at branching |
| Palantir (Setty 2019) | directed Markov chain on diffusion graph oriented by an early cell | branching fate; needs fate probabilities + differentiation-potential entropy | root-sensitive; entropy is a model-internal proxy; auto terminal states can be spurious |
| CellRank 2 (Weiler 2024) | non-reversible Markov chain; direction from pluggable kernels; GPCCA macrostates | directed fate mapping; multiview evidence; millions of cells; uncertainty-aware fate probabilities | kernel-quality dependent (garbage direction in, confident states out); n_states selection; metastability assumption |

Choosing by topology: linear -> Slingshot/DPT; bifurcating/multifurcating -> Slingshot/Palantir/CellRank 2; tree (>2 branches) -> Slingshot/Monocle3/PAGA-Tree; cyclic -> PAGA (Slingshot/Monocle2 cannot); disconnected/unknown -> PAGA first, then per-component pseudotime. Methodology evolves; verify current best practice against installed package docs before committing to one method, and report multi-method concordance.

### Decide Topology First With PAGA

**Goal:** Test whether putative branches are truly connected before committing to a continuous model.
**Approach:** Partition the kNN graph, compute PAGA connectivity, prune weak edges by threshold, then seed a global-faithful UMAP from PAGA.

```python
import scanpy as sc, numpy as np
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
sc.tl.leiden(adata, resolution=1.0, flavor='igraph', n_iterations=2, directed=False)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, threshold=0.03, color='leiden')   # prune low-connectivity (likely spurious) edges
sc.tl.umap(adata, init_pos='paga')                  # global topology preserved, local detail kept
```

The `threshold` in `sc.pl.paga` is the key judgment call: isolated clusters with no surviving edges are discrete cell types, not trajectory branches, and must not be forced into one ordering.

### Diffusion Pseudotime From an Anchored Root

**Goal:** Assign a scalar pseudotime once topology is fixed.
**Approach:** Run a diffusion map, set the root as a positional index on `adata.uns['iroot']`, then run DPT.

```python
sc.tl.diffmap(adata, n_comps=15)
adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type'] == 'HSC')[0]   # root anchored by a known marker, not by eye
sc.tl.dpt(adata, n_dcs=10, n_branchings=0)          # n_branchings=0 -> pure pseudotime; branch mode is fragile
```

`iroot` is a positional integer into `adata.obs_names`, set on `adata.uns` BEFORE `dpt`. The entire ordering and the sign of every gene trend flip with this choice.

### Fate Probabilities With Palantir

**Goal:** Represent each cell as a distribution over terminal fates near bifurcations.
**Approach:** Build a diffusion-map multiscale space, run Palantir from an early cell, and read pseudotime, entropy, and branch probabilities.

```python
import palantir
dm_res = palantir.utils.run_diffusion_maps(adata, n_components=5)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
pr_res = palantir.core.run_palantir(ms_data, early_cell='HSC_cell_id', terminal_states=None, num_waypoints=1200)
# pr_res.pseudotime, pr_res.entropy (differentiation potential), pr_res.branch_probs
```

Entropy of the fate-probability vector is the differentiation-potential proxy: high near multipotent cells, falling toward 0 as cells commit. Auto terminal-state detection can miss real fates or invent spurious ones, so verify terminals against markers.

### Directed Fate Mapping With CellRank 2

**Goal:** Infer initial states, terminal states, and uncertainty-aware fate probabilities from any directional evidence source.
**Approach:** Build a directed transition matrix from one or more kernels, combine with a connectivity kernel for smoothing, then coarse-grain into macrostates with GPCCA.

```python
import cellrank as cr
pk = cr.kernels.PseudotimeKernel(adata, time_key='dpt_pseudotime').compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
combined = 0.8 * pk + 0.2 * ck                      # weights are a researcher choice; sweep them

g = cr.estimators.GPCCA(combined)
g.compute_macrostates(n_states=10, cluster_key='leiden')   # n_states from the Schur/eigenvalue spectral gap
g.predict_terminal_states(method='stability')
g.predict_initial_states(n_states=1)
g.compute_fate_probabilities()
g.compute_lineage_drivers()
```

Kernels decouple WHERE direction comes from (RealTime when timepoints exist, Pseudotime/CytoTRACE otherwise, Velocity only when trustworthy, Connectivity for smoothing) from WHAT is computed (GPCCA macrostates + fate probabilities). Prefer the RealTimeKernel for time courses. Fate probabilities are a deterministic function of the transition matrix, so a wrong kernel yields confidently wrong, well-formed probabilities with no internal warning; check that conclusions survive dropping the velocity kernel.

### Slingshot and Monocle3 (R)

**Goal:** Fit smooth lineage curves (Slingshot) or a principal graph (Monocle3) and order cells.
**Approach:** Slingshot needs user-supplied dimred + cluster labels + a start cluster; Monocle3 learns its own graph and roots by node.

```r
library(slingshot)
sce <- slingshot(sce, clusterLabels='seurat_clusters', reducedDim='UMAP', start.clus='HSC')
pt  <- slingPseudotime(sce)     # cells x lineages; NA off-lineage; a trunk cell scores in every descendant lineage
```

```r
library(monocle3)
cds <- cluster_cells(cds)                           # produces clusters AND partitions
cds <- learn_graph(cds, use_partition = TRUE)       # TRUE allows disconnected trajectories
cds <- order_cells(cds, root_pr_nodes = root_node)  # root via graph node name, anchored by biology
graph_test_res <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)   # Moran's I trajectory DE
```

`start.clus` is mandatory in practice for Slingshot; downstream DE goes through tradeSeq (`fitGAM` then `associationTest` for any-variation-along-pseudotime or `startVsEndTest` for endpoint contrasts), not Slingshot itself. Monocle3's own trajectory DE is `graph_test` above. Monocle3's principal graph is learned in UMAP space, so loops and branches can be embedding artifacts.

## RNA Velocity

RNA velocity infers the time derivative of the spliced-mRNA state from the lag between unspliced (nascent) and spliced mRNA: velocity ds/dt = beta*u - gamma*s. It is a model-based extrapolation on a timescale of hours, and every downstream claim inherits the model's assumptions.

| Mode (`mode=`) | Model | Use when | Fails when |
|----------------|-------|----------|------------|
| `'deterministic'` | La Manno steady-state regression on extreme quantiles | quick first pass; well-separated induction/repression | assumes common splicing rate and that data spans both steady states; transient populations mis-fit |
| `'stochastic'` (default) | adds 2nd-moment treatment; GLS on both moments | a more robust gamma without the dynamical EM cost | still steady-state; same constant-rate assumption |
| `'dynamical'` | full likelihood EM; per-gene alpha/beta/gamma + latent time | transient states; needs gene-shared latent time | `recover_dynamics` dominates runtime; can still mis-fit multi-kinetics genes |

**Goal:** Estimate velocity direction and a latent-time ordering.
**Approach:** Compute moments, recover dynamics (dynamical only), compute velocity, build the velocity graph, then sanity-check confidence and phase portraits before any embedding plot.

```python
import scvelo as scv
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata)                      # dynamical only
scv.tl.velocity(adata, mode='dynamical')            # DEFAULT is 'stochastic'; pass 'dynamical' explicitly
scv.tl.velocity_graph(adata)
scv.tl.velocity_confidence(adata)                   # inspect BEFORE trusting the stream plot
scv.pl.velocity(adata, var_names=['GATA1'])         # per-gene phase portrait, not just the embedding
```

Bergen 2021 failure modes are the DEFAULT expectation, not edge cases. Velocity is unreliable or invalid in mature/terminal/non-dividing systems (adult neurons, steady-state tissue), where little net du/dt means noise dominates and arrows can point backward; under heterogeneous kinetics, one global gamma per gene mis-fits multi-branch systems; and a clean 2D stream plot can manufacture coherence the high-dimensional field lacks. Deeper still (Gorin 2022), the velocity ODE is a deterministic reduction of a stochastic process, intronic reads are a biased proxy for nascent RNA (internal priming, intron retention, 3' and length bias all corrupt gamma), and confidence/coherence metrics reward the kNN smoothing of the moments step rather than correspondence to truth (Zheng 2023). Do not consume raw arrows: feed velocity into CellRank 2 as ONE kernel, validate against known markers or metabolic labeling, and gate interpretation with uncertainty (veloVI `get_directional_uncertainty`).

Quantifier disagreement is first-order, not a detail (Soneson 2021): velocyto vs kb-python (`nac`) vs alevin-fry (USA mode) vs STARsolo (`--soloFeatures Gene Velocyto`) use different intron models and ambiguous-read rules, which shift the unspliced/spliced ratio, change gamma, and can flip velocity sign on borderline genes. Single-nucleus data is intron-rich; use the nascent/mature (`nac`/spliceu) framing. A direction that is not stable across at least two quantifiers is a pipeline artifact, not a finding.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Smooth pseudotime axis through what are actually discrete cell types | no real continuum; kNN bridges islands with spurious edges | run PAGA first; if clusters have no surviving connectivity edges, do not order them |
| Every gene trend reverses between runs | root chosen by eye / on a UMAP; ordering flips with origin | anchor `iroot`/`root_pr_nodes` with a known marker, real time, velocity, or stemness |
| Branch assignment unstable across parameters | hard-assigning progenitors whose fate is genuinely undetermined | report fate PROBABILITIES (Palantir branch_probs, CellRank), do not hard-assign near bifurcations |
| Trajectory passes through a near-empty region | rare/fast-traversed intermediate state is unsampled; graph interpolates a void | check cell density along the path; treat the gap as missing data, not a real intermediate |
| Velocity stream looks clean but points backward | mature/terminal/non-cycling system; little net du/dt, noise dominates | velocity is invalid here; do not interpret arrows; validate with markers/lineage or drop velocity |
| Velocity direction flips when the quantifier changes | intron model / ambiguous-read handling differs across tools | re-run with a second quantifier; only trust direction stable across both (Soneson 2021) |
| high `velocity_confidence` but biologically wrong arrows | metric rewards kNN smoothing, not truth (Zheng 2023) | sweep `n_neighbors`; require orthogonal validation, not the confidence score alone |
| CellRank invents discrete macrostates from a smooth flow | metastability assumption violated; GPCCA forced to partition a continuum | show the Schur/eigenvalue spectrum; justify n_states by a real gap or treat states as coarse-graining artifacts |
| Pseudotime intervals reported as durations | pseudotime is monotone in progression, not time | only RealTimeKernel/WOT exploit actual time; do not read intervals as elapsed hours |

## Related Skills

- single-cell/clustering - Leiden clusters and the kNN graph that PAGA, DPT, and the velocity moments step all depend on
- single-cell/preprocessing - normalization, HVG selection, and PCA whose choices the inferred axis inherits
- single-cell/lineage-tracing - orthogonal lineage ground truth that tests whether a state-based trajectory predicts fate
- single-cell/cell-communication - downstream signaling analysis along the inferred trajectory
- differential-expression/deseq2-basics - pseudobulk DE between trajectory endpoints or branches

## References

Haghverdi L, Buttner M, Wolf FA, Buettner F, Theis FJ (2016). Diffusion pseudotime robustly reconstructs lineage branching. Nat Methods 13(10):845-848.
Street K, Risso D, Fletcher RB, et al. (2018). Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics 19:477.
Cao J, Spielmann M, Qiu X, et al. (2019). The single-cell transcriptional landscape of mammalian organogenesis (Monocle3). Nature 566(7745):496-502.
Wolf FA, Hamey FK, Plass M, et al. (2019). PAGA: graph abstraction reconciles clustering with trajectory inference. Genome Biology 20:59.
Setty M, Kiseliovas V, Levine J, et al. (2019). Characterization of cell fate probabilities in single-cell data with Palantir. Nat Biotechnol 37:451-460.
Saelens W, Cannoodt R, Todorov H, Saeys Y (2019). A comparison of single-cell trajectory inference methods. Nat Biotechnol 37(5):547-554.
Lange M, Bergen V, Klein M, et al. (2022). CellRank for directed single-cell fate mapping. Nat Methods 19(2):159-170.
Weiler P, Lange M, Klein M, Pe'er D, Theis FJ (2024). CellRank 2: unified fate mapping in multiview single-cell data. Nat Methods 21(7):1196-1205.
La Manno G, Soldatov R, Zeisel A, et al. (2018). RNA velocity of single cells. Nature 560:494-498.
Bergen V, Lange M, Peidli S, Wolf FA, Theis FJ (2020). Generalizing RNA velocity to transient cell states through dynamical modeling (scVelo). Nat Biotechnol 38(12):1408-1414.
Bergen V, Soldatov RA, Kharchenko PV, Theis FJ (2021). RNA velocity - current challenges and future perspectives. Mol Syst Biol 17(8):e10282.
Gayoso A, Weiler P, Lotfollahi M, et al. (2024). Deep generative modeling of transcriptional dynamics for RNA velocity analysis (veloVI). Nat Methods 21:50-59.
Weinreb C, Wolock S, Tusi BK, Socolovsky M, Klein AM (2018). Fundamental limits on dynamic inference from single-cell snapshots. PNAS 115(10):E2467-E2476.
Weinreb C, Rodriguez-Fraticelli A, Camargo FD, Klein AM (2020). Lineage tracing on transcriptional landscapes links state to fate (LARRY). Science 367(6479):eaaw3381.
Gorin G, Fang M, Chari T, Pachter L (2022). RNA velocity unraveled. PLoS Comput Biol 18(9):e1010492.
Zheng SC, Stein-O'Brien G, Boukas L, Goff LA, Hansen KD (2023). Pumping the brakes on RNA velocity by understanding and interpreting RNA velocity estimates. Genome Biology 24(1):246.
Soneson C, Srivastava A, Patro R, Stadler MB (2021). Preprocessing choices affect RNA velocity results for droplet scRNA-seq data. PLoS Comput Biol 17(1):e1008585.
