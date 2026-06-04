---
name: bio-imaging-mass-cytometry-phenotyping
description: Assign cell types from marker expression in IMC/MIBI data using clustering (PhenoGraph/FlowSOM/Leiden/Pixie), marker-based probabilistic classifiers (Astir), or image-context CNNs (CellSighter), covering the double-positive segmentation artifact, lineage-vs-state markers, the two spillover types, and why a "cell type" in imaging is conditioned on a segmentation guess. Use when phenotyping segmented IMC cells, choosing clustering vs classification, diagnosing implausible double-positive populations, separating lineage from functional markers, or transferring labels across a cohort.
tool_type: python
primary_tool: scanpy
---

## Version Compatibility

Reference examples tested with: scanpy 1.10+, anndata 0.10+, astir 0.1.4+, numpy 1.26+, scikit-learn 1.4+, FlowSOM 2.10+ (R)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: arcsinh cofactor for IMC single-cell means is ~1, not the suspension-CyTOF 5 -- do not hard-code 5. Astir assigns a per-cell probability and routes below-threshold cells (default 0.7) to "Unknown" rather than forcing a call. CellSighter consumes raw multi-channel image crops + masks (not a mean matrix). FlowSOM consensus metaclustering can override `set.seed()` via ConsensusClusterPlus.

# Cell Phenotyping for IMC

**"Assign cell types to my segmented IMC cells"** -> Map each cell's marker profile to an identity, while distinguishing real co-expression from segmentation/spillover artifacts.
- Python: `scanpy.tl.leiden` (cluster then annotate), `astir` (marker-dictionary classifier)
- R: `FlowSOM` (self-organizing-map clustering)

## The Single Most Important Modern Insight -- a "cell type" in imaging is an inference conditioned on a segmentation guess

In suspension CyTOF each event is one physically isolated cell; in imaging, every cell-by-marker row is the integral of pixels inside a polygon a segmentation algorithm drew, and that polygon is wrong at a non-trivial fraction of cells -- so the most dangerous phenotypes are not biology but boundary artifacts. The canonical case is the CD3+CD20+ ("T/B") double-positive, also CD3+CD68+ and panCK+CD45+. It arises by two distinct mechanisms that are indistinguishable in the mean matrix: segmentation merging (one polygon spans a T cell and a B cell) and lateral spillover (a neighbor's membrane bleeds across the boundary even with perfect masks). The diagnostic tell that separates artifact from biology is spatial: artifactual double-positives localize to cell BORDERS and to high-density regions, so a suspect population must be mapped back onto the image before it is believed (CellSighter authors state the matrix cannot separate the two). The asymmetry that drives method choice: clustering CREATES the artifact as a named population, while a marker-dictionary classifier (Astir) REFUSES it -- a true double-positive vector matches no defined type and is quarantined as "Unknown" rather than crowned a new lineage. This is why imaging-aware groups increasingly prefer (semi-)supervised phenotyping for the lineage layer, and why mean-expression clustering imported wholesale from CyTOF inherits none of the spatial information that would let it notice the polygon was wrong.

## Phenotyping Approach Taxonomy

| Approach | Tools | Input | Robust to bad segmentation? | Failure signature |
|----------|-------|-------|-----------------------------|-------------------|
| Unsupervised clustering | PhenoGraph, FlowSOM, Leiden | cell x marker mean matrix | No -- averages spilled signal into a fake type | phantom double-positive clusters; resolution-dependent type count |
| Pixel-then-cell clustering | Pixie (ark-analysis) | pixel x marker, then cell | More -- avoids committing to a segmentation mean early | parameter-sensitive; still unsupervised |
| Marker-based probabilistic | Astir | mean matrix + marker->type YAML | Partially -- ambiguous cells -> "Unknown" | high Unknown rate if dictionary/markers wrong |
| Image-context CNN | CellSighter, MAPS | raw image crops + masks + labels | Yes -- sees where the signal sits | needs representative labels not harvested from clustering |
| Segmentation-aware mixture | STARLING | mean matrix + doublet prior | Yes -- models a cell as a mixture of two | newer; verify priors |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Can write marker->celltype rules, no training labels | Astir (lineage layer) | deterministic, fast, "Unknown" for ambiguous, separates type from state |
| Have expert-labeled cells, segmentation/spillover is a known problem | CellSighter | image context rejects border/spillover double-positives |
| Severe segmentation doubt | STARLING | explicitly models doublet/contamination mixtures |
| Annotated reference cohort, want label transfer | STELLAR | graph model using neighborhood + expression |
| Exploratory, no priors, accept manual annotation | Pixie (most robust) or Leiden/FlowSOM (least) | always run the double-positive image-diagnostic first |
| Any across-condition comparison of the resulting types | hand off to differential-analysis | phenotyping and statistical-unit choice are orthogonal |

## Load and Transform

**Goal:** Build the single-cell matrix on the correct count scale.

**Approach:** Arcsinh with cofactor ~1 for IMC means (not 5), and keep raw counts available. Treat zeros as genuine low ion counts plus Poisson noise, not technical dropout -- scRNA-style imputation hallucinates expression.

```python
import scanpy as sc
import anndata as ad
import numpy as np

adata = ad.read_h5ad('imc_segmented.h5ad')
adata.layers['counts'] = adata.X.copy()
adata.X = np.arcsinh(adata.X / 1.0)   # cofactor ~1 for IMC single-cell means, not 5
```

## Marker-Based Classification with Astir

**Goal:** Assign lineage with a principled abstention instead of a forced call.

**Approach:** Encode marker->celltype rules in a YAML with separate `cell_type` and `cell_state` blocks; Astir returns a per-cell probability and labels below-threshold cells "Unknown". The Unknown rate is itself QC -- 40% Unknown means the dictionary or panel is mis-specified, not that the cells are exotic.

```python
from astir.data import from_anndata_yaml

# inputs are PATHS: an .h5ad and a marker YAML with a cell_type block (CD3->T, CD20->B,
# CD68->Macrophage; no type is both) and an optional cell_state block (Ki67, PD-1)
ast = from_anndata_yaml('imc_segmented.h5ad', 'markers.yaml')
ast.fit_type()
celltypes = ast.get_celltypes(threshold=0.7)   # per-cell labels; < 0.7 -> 'Unknown' (information, not failure)
```

## Cluster on Lineage Markers Only

**Goal:** Discover structure without splitting one type into activation states.

**Approach:** Cluster on lineage markers only; mixing continuous state markers (Ki67, PD-1) fragments one type into proliferating/resting pseudo-types. Validating clusters with the same markers used to cluster is circular -- confirm with held-out evidence (spatial context, independent markers).

```python
lineage = ['CD45', 'CD3', 'CD8', 'CD4', 'CD20', 'CD68', 'E-cadherin']   # lineage only, no Ki67/PD-1
sub = adata[:, lineage]
sc.pp.pca(sub, n_comps=min(15, len(lineage)))
sc.pp.neighbors(sub, n_neighbors=15)
sc.tl.leiden(sub, resolution=0.5)
adata.obs['leiden'] = sub.obs['leiden']
# report cluster stability across resolutions/seeds rather than one hand-picked setting
```

## Diagnose Double-Positive Populations

**Goal:** Decide whether an implausible co-expressing population is biology or artifact.

**Approach:** A real co-expressing cell has the second marker over its own membrane/cytoplasm; an artifact has it concentrated on the border adjacent to a donor neighbor. Quantify how often the suspect cells sit next to a cell of the donor type -- border + donor-adjacency means spillover/merge, not a lineage.

```python
import squidpy as sq

sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
suspect = adata.obs['cell_type'] == 'CD3+CD20+?'
# if suspect cells are overwhelmingly adjacent to true B cells (the CD20 donor), the CD20
# is spillover/merge, not endogenous -- treat the population as a QC failure, not a discovery
```

## Per-Method Failure Modes

### Clustering -- arbitrary resolution invents types
**Trigger:** tuning Leiden resolution / FlowSOM metacluster count until clusters match expectation. **Mechanism:** the resolution directly sets the type count; it is an identifiability hole, not a tuning knob. **Symptom:** unreproducible type counts; clusters drift across samples. **Fix:** fix the type set with a dictionary/classifier, or report stability across resolutions and seeds.

### FlowSOM -- seed override
**Trigger:** `set.seed()` then consensus metaclustering, expecting reproducibility. **Mechanism:** ConsensusClusterPlus resets the seed internally. **Symptom:** cluster identities differ between runs. **Fix:** set the seed inside the consensus call; assess label stability across runs.

### "I compensated, so no double-positives"
**Trigger:** running CATALYST channel compensation and assuming spatial spillover is handled. **Mechanism:** channel/isotope spillover and lateral/optical spillover are different physical problems. **Symptom:** double-positives persist after channel compensation. **Fix:** channel compensation early (pixel level), REDSEA boundary compensation after segmentation; neither fixes a merged segment -- improve segmentation first.

### Imputing IMC zeros
**Trigger:** scRNA-style dropout imputation on the count matrix. **Mechanism:** IMC zeros are largely genuine low counts, not a capture-dropout mechanism. **Symptom:** hallucinated expression, inflated positivity. **Fix:** model low counts as low counts; do not impute.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| arcsinh cofactor ~1 (IMC means) | Hunter 2024 *Cytometry A* 105:36 | preserves positive/negative separation; 5 over-compresses |
| Astir assignment threshold 0.7 (package default) | Geuenich 2021 *Cell Syst* 12:1173 | principled abstention; the Unknown rate is a QC metric |
| ~40 markers, no redundancy | panel design | one channel can decide a fate -- verify the load-bearing channel per type |
| CellSighter labels NOT from clustering | Amitay 2023 *Nat Commun* 14:4302 | clustering-derived labels re-import the double-positive artifact |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Tidy CD3+CD20+ cluster reported as a lineage | clustering legitimized a segmentation/spillover artifact | diagnose border-localization on the image; treat as QC failure |
| One T-cell type split into two clusters | state markers (Ki67) mixed into lineage clustering | cluster lineage on lineage markers; profile state within type |
| ~40% of cells "Unknown" in Astir | mis-specified dictionary or missing type | inspect Unknown cells; iterate the YAML; tune 0.7 consciously |
| Cluster identities drift between analyses | stochastic clustering / FlowSOM seed | pin seeds, assess stability; do not assume "cluster 7" is stable |
| "Disease has more Tregs" with p~0 | cell-level testing (pseudoreplication) | aggregate to per-patient proportions; see differential-analysis |

## References

- Levine JH, Simonds EF, Bendall SC, et al. 2015. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. *Cell* 162(1):184-197. — PhenoGraph.
- Van Gassen S, Callebaut B, Van Helden MJ, et al. 2015. FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data. *Cytometry A* 87(7):636-645. — FlowSOM.
- Traag VA, Waltman L, van Eck NJ. 2019. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep* 9(1):5233. — Leiden.
- Geuenich MJ, Hou J, Lee S, et al. 2021. Automated assignment of cell identity from single-cell multiplexed imaging and proteomic data. *Cell Syst* 12(12):1173-1186.e5. — Astir.
- Amitay Y, Bussi Y, Feinstein B, Bagon S, Milo I, Keren L. 2023. CellSighter: a neural network to classify cells in highly multiplexed images. *Nat Commun* 14:4302. — image-context classification.
- Liu CC, Greenwald NF, et al. 2023. Robust phenotyping of highly multiplexed tissue imaging data using pixel-level clustering. *Nat Commun* 14:4618. — Pixie.
- Brbic M, Cao K, Hickey JW, et al. 2022. Annotation of spatially resolved single-cell data with STELLAR. *Nat Methods* 19(11):1411-1418. — label transfer.
- Campbell KR, et al. 2025. Segmentation aware probabilistic phenotyping of single-cell spatial protein expression data. *Nat Commun* 16:389. — STARLING.
- Bai Y, Zhu B, Rovira-Clave X, et al. 2021. Adjacent Cell Marker Lateral Spillover Compensation and Reinforcement for Multiplexed Images. *Front Immunol* 12:652631. — REDSEA.
- Chevrier S, Crowell HL, Zanotelli VRT, et al. 2018. Compensation of Signal Spillover in Suspension and Imaging Mass Cytometry. *Cell Syst* 6(5):612-620.e5. — channel spillover.
- Hunter B, Nicorescu I, Foster E, et al. 2024. OPTIMAL: An OPTimized Imaging Mass cytometry AnaLysis framework for benchmarking segmentation and data exploration. *Cytometry A* 105(1):36-53. — arcsinh cofactor 1 for IMC single-cell means.

## Related Skills

- cell-segmentation - double-positives diagnose the segmentation/spillover that phenotyping inherits
- data-preprocessing - arcsinh cofactor and channel spillover compensation
- spatial-analysis - phenotype labels feed neighborhood and niche analysis
- differential-analysis - comparing cell-type proportions across conditions at the patient level
- interactive-annotation - mapping clusters back onto tissue to confirm they are real
- flow-cytometry/clustering-phenotyping - FlowSOM/PhenoGraph background for suspension data
