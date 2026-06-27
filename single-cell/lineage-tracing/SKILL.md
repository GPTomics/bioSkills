---
name: bio-single-cell-lineage-tracing
description: Reconstructs single-cell lineage trees and clonal relationships from CRISPR/Cas9 scars, static expressed barcodes (LARRY/CellTag), or somatic mtDNA mutations using Cassiopeia, Startle, and CoSpar. Use when building a phylogeny from barcode scars, choosing a tree-reconstruction solver, handling homoplasy and dropout, grouping clones from mtDNA, integrating clone with transcriptomic state, or judging whether a state-based fate call is trustworthy.
tool_type: python
primary_tool: Cassiopeia
---

## Version Compatibility

Reference examples tested with: Cassiopeia 2.0+, CoSpar 0.3+, scanpy 1.10+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Lineage Tracing

**"Reconstruct cell lineage from barcodes"** -> Read heritable marks across single cells and ask which cells share which marks to recover an ontogenetic phylogeny or clonal grouping.
- Python: `cassiopeia` (scar-tree reconstruction), `cospar` (clone + state integration), mtDNA variant callers (mgatk/MAESTER pipelines)

## Governing Principle

Transcriptomic state does NOT fully predict fate. Weinreb 2020 (LARRY) showed sister cells in an indistinguishable transcriptomic state systematically diverge in fate, so the information that decides a bifurcation is heritable but invisible to the measured transcriptome. Three consequences drive every decision here.

1. Lineage is orthogonal to expression, not redundant with it. A purely state-based trajectory method is systematically wrong about commitment for such populations, which is the empirical justification for every "barcode + transcriptome in the same cell" assay and for integrative tools like CoSpar.
2. Reconstruction is phylogenetics on error-prone scars. A scar/barcode tree inherits every pathology of molecular phylogenetics, sharpened by CRISPR peculiarities: homoplasy (independent cells acquire the identical scar), dropout (missing vs unedited confusion), saturation of a finite editable array, and non-clock editing rates.
3. The state->fate map can be one-to-many. A state-based branch call can be confident precisely because it is blind to the heritable variable that actually decides fate. Prospective tracing (engineered barcodes installed BEFORE the process) is the only design that can measure fate independently of state and thus test state->fate; retrospective tracing (mtDNA read out after) recovers ancestry but cannot by itself establish what state preceded a fate.

Two hard errors a reviewer presses on. Missing-as-unedited: an uncaptured edit recorded as the "0" state is indistinguishable from a site that genuinely never edited, and heritable excision dropout removes a whole character across an entire clade, biasing topology, not merely adding noise. Homoplasy: Cas9 indel outcomes are highly non-uniform, so a handful of indels dominate and parsimony falsely fuses unrelated lineages that share a frequent scar. Topology error also compounds toward the root, where the deepest, most consequential splits rest on the fewest characters.

## Assay Decision Table

Choose the recording technology by the question, not by availability.

| Assay | Mark / model | Use when | Fails when |
|-------|--------------|----------|------------|
| CRISPR scar array (GESTALT/scGESTALT/ScarTrace/LINNAEUS) | cumulative irreversible Cas9 indels = phylogenetic characters | deep tree topology in an engineered organism; scar + scRNA-seq in the same cell | saturation caps depth; homoplasy fuses lineages; heritable dropout deletes clades; not usable in native human tissue |
| Static expressed barcode (LARRY, Weinreb 2020) | one unique inherited lentiviral barcode per founder = a flat clone | clean state->fate maps via split-and-profile; proving state underdetermines fate | gives clonal membership, NOT division-order topology; library << founders causes barcode collisions |
| Combinatorial/sequential barcode (CellTag, Biddy 2018) | combination of expressed tags + nested timepoints | extra clonal resolution and coarse multi-level nesting | shallow trees; collisions; not a resolved phylogeny |
| Somatic mtDNA (Ludwig 2019; mtscATAC-seq; MAESTER) | drifting heteroplasmy of somatic mtDNA variants | retrospective tracing in primary human tissue with no engineering | low mutation rate -> few informative variants; hotspot homoplasy; coverage/dropout; heteroplasmy drift + selection; gives clonal grouping not deep trees |

## Reconstruction Solver Decision Table

For scar data, the solver is a separate choice from the assay (Cassiopeia, Jones 2020; Startle, Sashittal 2023).

| Solver | Model | Use when | Fails when |
|--------|-------|----------|------------|
| VanillaGreedySolver | top-down parsimony, split on most-frequent mutation | fast first pass; 10^4-10^5 cells | greedy split errors propagate; sensitive to homoplasy |
| ILPSolver | integer-LP Steiner tree (Gurobi); near-optimal | small clades needing accuracy | expensive; does not scale to large trees |
| HybridSolver | greedy top + ILP on small subclades | the practical default for large data | inherits greedy errors at the top split |
| NeighborJoiningSolver | distance-based (weighted Hamming) | quick comparison baseline; non-character distances | less accurate than parsimony on scar characters |
| Startle (Startle-ILP / Startle-NNI) | star-homoplasy: a character mutates at most once per root-to-leaf path | severe homoplasy/dropout breaks parsimony | ILP cost; NNI is heuristic at scale |

Run a panel of solvers, not one: it is rare for a single solver to be optimal over all parts of a tree, and agreement across solvers is the practical certainty signal. Always weight indels by their formation probability (down-weight frequent low-information scars) and report robustness to homoplasy and dropout. Methodology evolves; verify the current solver API and recommended defaults against the installed Cassiopeia docs.

### Build a Character Matrix and Reconstruct a Tree

**Goal:** Turn scar calls into a maximum-parsimony lineage tree with missing data modeled explicitly.
**Approach:** Load a cells x sites character matrix (0 = unedited, 1+ = distinct scars, -1 = missing), assess missingness and informativeness, then solve.

```python
import cassiopeia as cas
import numpy as np

tree = cas.data.CassiopeiaTree(character_matrix=char_matrix, cell_meta=cell_meta)
print(f'cells {tree.n_cell}  characters {tree.n_character}  missing {(char_matrix == -1).mean():.2%}')

solver = cas.solver.VanillaGreedySolver()
solver.solve(tree, collapse_mutationless_edges=True)   # collapse edges with no supporting mutation
newick = tree.get_newick()
```

The `-1` missing state must stay distinct from the `0` unedited state: collapsing missing into unedited is the single most consequential preprocessing error, since heritable dropout is tree-correlated and silently erases real structure.

### Compare Solvers and Score Tree Robustness

**Goal:** Quantify how much the topology depends on solver choice and on homoplasy/dropout.
**Approach:** Solve with several solvers, then compare the resulting trees with Robinson-Foulds and the depth-stratified triplets-correct metric.

```python
hybrid = cas.solver.HybridSolver(top_solver=cas.solver.VanillaGreedySolver(), bottom_solver=cas.solver.ILPSolver(), cell_cutoff=200)
nj = cas.solver.NeighborJoiningSolver(dissimilarity_function=cas.solver.dissimilarity_functions.weighted_hamming_distance)
for s in (hybrid, nj):
    s.solve(tree)                                       # solve independent copies in practice
rf, rf_max = cas.critique.robinson_foulds(tree_a, tree_b)
triplet_acc = cas.critique.triplets_correct(tree_a, tree_b)
```

Triplets-correct is depth-stratified, so it exposes the field's hard truth: deep (near-root) splits are the least certain and the most consequential, while well-supported leaf structure is often the least interesting biologically.

### Build a Character Matrix From Raw Reads

**Goal:** Go from aligned barcode reads to an allele table and character matrix.
**Approach:** Resolve UMIs, align to the reference, call alleles, group cells into clonal populations, then convert the allele table.

```python
umi_table = cas.pp.resolve_umi_sequence(molecule_table, output_directory='.', min_umi_per_cell=10)
aligned = cas.pp.align_sequences(umi_table, ref_filepath='barcode_reference.fa')
alleles = cas.pp.call_alleles(aligned, ref_filepath='barcode_reference.fa')
alleles = cas.pp.call_lineage_groups(alleles, output_directory='.')
char_matrix, priors, state_map = cas.pp.convert_alleletable_to_character_matrix(alleles)
```

`convert_alleletable_to_character_matrix` returns indel priors alongside the matrix; pass those priors to the solver so frequent low-information scars are down-weighted against homoplasy.

### Integrate Clones With State Using CoSpar

**Goal:** Recover early fate bias from sparse clonal barcodes rather than assuming the manifold encodes fate.
**Approach:** Fit a transition map jointly from clonal observations and transcriptomic similarity, then read fate bias and fate maps.

```python
import cospar as cs
adata = cs.hf.read('lineage_traced.h5ad')
adata = cs.pp.initialize_adata_object(adata, X_clone=adata.obsm['X_clone'], time_info=adata.obs['time_info'])
adata = cs.tmap.infer_Tmap_from_multitime_clones(adata, smooth_array=[15, 10, 5], sparsity_threshold=0.1)
cs.tl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'])
cs.pl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'])
```

CoSpar operationalizes Weinreb 2020: it propagates fate probabilities onto cells lacking clonal labels and is robust to severe downsampling of lineage data, but it needs paired clone + state and does NOT build a phylogenetic tree (clones are flat). CoSpar needs MULTIPLE independent clones to be lineage-informed; with effectively one clone the constraint is vacuous and the transition map degenerates to transcriptomic similarity, the state-only answer CoSpar exists to correct. For tree topology from scars, use Cassiopeia or Startle.

## Threshold and Parameter Rationale

| Parameter | Typical value | Rationale |
|-----------|---------------|-----------|
| min_umi_per_cell | ~10 | below this, allele calls are dominated by sequencing noise |
| missing fraction per cell | drop > ~0.5 | cells missing most characters carry little phylogenetic signal and inflate ambiguity |
| informative character | states in > 1 cell | a scar seen in one cell cannot group lineages; uninformative for topology |
| indel prior weighting | from empirical indel frequencies | frequent microhomology-driven indels are high-homoplasy, low-information; down-weight them |
| barcode library complexity | >> number of founders | small libraries cause collisions (two founders share a barcode -> phantom merged clone) |
| HybridSolver cell_cutoff | ~200 | subclades below the cutoff are solved exactly by ILP; above it, greedily |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Distinct lineages collapse into one clade | missing data coded as the unedited 0 state | keep -1 missing distinct from 0; model dropout, never treat it as unedited |
| Parsimony fuses unrelated cells | homoplasy: independent cells share a frequent indel | weight indels by formation probability; use Startle's star-homoplasy model under heavy convergence |
| Late divisions are unresolved near the leaves | editable array saturated; recording stopped early | use inducible/paced recorders; report the per-site edit fraction distribution |
| Two founders appear as one giant clone | barcode library too small relative to founders -> collision | use library complexity >> cell number; estimate collisions empirically |
| Impossible chimeric clones / character vectors | doublets carry two barcode/scar sets | run doublet detection and barcode-consistency filtering before reconstruction |
| Deep splits flip between solvers | early splits rest on the fewest, most-overwritten characters | report branch support; trust leaf structure more than the root; run a solver panel |
| mtDNA "tree" is actually clonal blobs | low somatic mutation rate; hotspot homoplasy; heteroplasmy drift and selection | claim clonal grouping not deep ordered trees; blacklist NUMTs/RNA-edit/hotspot sites |
| State-based branch call confidently wrong | state underdetermines fate (Weinreb 2020); map is one-to-many | frame fate as a prediction; validate with prospective lineage data, integrate with CoSpar |

## Related Skills

- single-cell/trajectory-inference - state-based pseudotime/velocity that lineage data tests and corrects
- single-cell/preprocessing - QC, doublet handling, and normalization upstream of barcode and clone calls
- single-cell/clustering - cell-type labels annotated onto tree leaves and clones
- phylogenetics/modern-tree-inference - general phylogenetic inference, parsimony vs ML, and branch support

## References

Weinreb C, Rodriguez-Fraticelli A, Camargo FD, Klein AM (2020). Lineage tracing on transcriptional landscapes links state to fate during differentiation (LARRY). Science 367(6479):eaaw3381.
McKenna A, Findlay GM, Gagnon JA, Horwitz MS, Schier AF, Shendure J (2016). Whole-organism lineage tracing by combinatorial and cumulative genome editing (GESTALT). Science 353(6298):aaf7907.
Raj B, Wagner DE, McKenna A, et al. (2018). Simultaneous single-cell profiling of lineages and cell types in the vertebrate brain (scGESTALT). Nat Biotechnol 36(5):442-450.
Alemany A, Florescu M, Baron CS, Peterson-Maduro J, van Oudenaarden A (2018). Whole-organism clone tracing using single-cell sequencing (ScarTrace). Nature 556(7699):108-112.
Spanjaard B, Hu B, Mitic N, et al. (2018). Simultaneous lineage tracing and cell-type identification using CRISPR-Cas9-induced genetic scars (LINNAEUS). Nat Biotechnol 36:469-473.
Biddy BA, Kong W, Kamimoto K, et al. (2018). Single-cell mapping of lineage and identity in direct reprogramming (CellTag). Nature 564:219-224.
Wang SW, Herriges MJ, Hurley K, Kotton DN, Klein AM (2022). CoSpar identifies early cell fate biases from single-cell transcriptomic and lineage information. Nat Biotechnol 40:1066-1074.
Ludwig LS, Lareau CA, Ulirsch JC, et al. (2019). Lineage tracing in humans enabled by mitochondrial mutations and single-cell genomics. Cell 176(6):1325-1339.
Lareau CA, Ludwig LS, Muus C, et al. (2021). Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling (mtscATAC-seq). Nat Biotechnol 39:451-461.
Miller TE, Lareau CA, Verga JA, et al. (2022). Mitochondrial variant enrichment from high-throughput single-cell RNA sequencing resolves clonal populations (MAESTER). Nat Biotechnol 40:1030-1034.
Jones MG, Khodaverdian A, Quinn JJ, et al. (2020). Inference of single-cell phylogenies from lineage tracing data using Cassiopeia. Genome Biology 21:92.
Sashittal P, Schmidt H, Chan M, Raphael BJ (2023). Startle: a star homoplasy approach for CRISPR-Cas9 lineage tracing. Cell Systems 14(12):1113-1121.
