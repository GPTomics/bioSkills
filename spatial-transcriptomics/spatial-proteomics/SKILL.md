---
name: bio-spatial-transcriptomics-spatial-proteomics
description: Analyzes multiplexed antibody-imaging data (CODEX/PhenoCycler, MIBI-TOF, IMC, CyCIF, Opal/Vectra mIF) as continuous protein intensity rather than transcript counts, using scimap and squidpy. Use when choosing an intensity transform/normalization (arcsinh cofactor vs z-score vs percentile -- NOT log1p-of-counts) and correcting channel spillover and antibody-batch effects; deciding whether to phenotype by gating or by clustering on intensities; recognizing that a bounded antibody panel makes marker absence uninformative; treating whole-cell segmentation (Mesmer) as the dominant error source; and knowing which platform applies and when to defer to the imaging-mass-cytometry skills for the IMC pipeline.
tool_type: python
primary_tool: scimap
---

## Version Compatibility

Reference examples tested with: scimap 2.0+, scanpy 1.10+, anndata 0.10+, squidpy 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Proteomics Analysis

**"Analyze my CODEX/MIBI/IMC multiplexed-imaging data"** -> Turn a cell-by-marker protein-intensity matrix into phenotyped cells and spatial neighborhoods, while treating intensity as a continuous, confounded signal.
- Python: per-marker `arcsinh`/z-score transform -> spillover/batch correction -> `scimap.tl.phenotype_cells()` (gating) or `scimap.tl.cluster()` (clustering) -> `squidpy.gr.nhood_enrichment()`

## Governing Principle

Protein intensity is continuous with antibody, batch, and staining confounds -- it is not a molecule count, and treating it like one is a category error that propagates through every downstream result.

A multiplexed-imaging measurement is the reporter signal (fluorescence photons or secondary-ion/metal counts) for an antibody bound to its epitope. That signal scales with antibody affinity, conjugation efficiency, staining-day conditions, fixation, and detector response -- none of which are the abundance of the protein, and all of which differ between markers and between samples. Consequences that separate this from RNA-based spatial omics: there is no Poisson/NB count model, so the variance-stabilizing transform is arcsinh (or z-score/percentile), NEVER log1p applied as if the values were UMIs; metal/fluor channels leak into each other (spillover) and must be compensated; and antibody-batch and staining variation must be normalized before any cross-sample comparison. Hickey 2021 (*Front Immunol* 12:727626) made the cost concrete: crossing 5 normalizations x 4 clustering methods on ONE CODEX dataset produced 20 different cell-type annotations -- the normalization-and-clustering choice, not the biology, dominated the phenotype calls.

The antibody panel is targeted, so absence is uninformative. A 20-100 marker panel is chosen a priori; the phenotype space is bounded by it exactly as an imaging RNA panel (Xenium/MERFISH/CosMx) bounds detectable transcripts. A cell type whose defining markers are off-panel is invisible or silently mis-assigned to the nearest panel-defined type -- there is no de-novo discovery. "Marker X is absent" usually means "X was not stained," not "the protein is not there."

Segmentation is the dominant downstream error source -- the per-cell intensity vector is only as good as the cell mask. There is no native cell in an image; a cell-by-marker matrix exists only after an algorithm draws boundaries. Lateral spillover of membrane/cytoplasmic signal into neighboring masks fabricates phantom double-positive cells (a CD3+CD20+ "cell" is usually a T cell touching a B cell), and every neighborhood, niche, and proximity result inherits that error. Whole-cell segmentation on a membrane/boundary stain (Mesmer/DeepCell, trained on the ~1M-cell TissueNet, is the multiplexed-imaging standard) is the highest-leverage decision; nucleus-only loses cytoplasmic signal and nucleus-expansion assumes round equal cells. The IMC/MIBI pipeline mechanics live in the imaging-mass-cytometry category -- this skill owns the platform breadth and the intensity reframe; defer the deep pipeline there.

## The Platform-Breadth Decision

This skill owns BREADTH across antibody-based platforms; IMC pipeline DEPTH lives in imaging-mass-cytometry. The first question is which platform produced the data, because chemistry sets the confounds.

| Platform | Chemistry | Markers | Strengths | Dominant confounds |
|----------|-----------|---------|-----------|--------------------|
| CODEX / PhenoCycler (Goltsev 2018) | DNA-barcoded antibodies, iterative fluorescent reporter cycles | ~50-60+ | High-plex on a standard fluorescence microscope; sub-um optical | Cycle-to-cycle registration drift, photobleaching/tissue degradation over many cycles; continuous intensity |
| MIBI-TOF (Angelo 2014; Keren 2019) | Lanthanide-metal antibodies, ion beam + TOF mass spec | ~40 metal channels | High resolution (~260-500 nm); low autofluorescence | Slow, small FOV; semi-quantitative (secondary-ion yield, detector); isotopic/channel crosstalk (spillover) |
| IMC (Giesen 2014) | Metal-isotope antibodies, UV laser ablation + CyTOF | ~40 metal channels | Metal multiplexing, no autofluorescence | ~1 um, slow ablation; metal-channel spillover (Chevrier 2018); conjugation-efficiency bias. Pipeline -> imaging-mass-cytometry |
| CyCIF / t-CyCIF (Lin 2018) | Cyclic IF: stain ~4 dyes, image, bleach, restain | up to ~60 | Conventional optical microscope, accessible | Bleach/restain degrades antigenicity; registration drift; autofluorescence |
| Opal / Vectra mIF (Parra 2017) | Tyramide-amplified multispectral IF | ~6-8 | Clinical-grade, FFPE-validated | Low plex; spectral unmixing artifacts; amplification nonlinearity |

CODEX/CyCIF/Opal yield continuous FLUORESCENCE intensity; MIBI/IMC yield semi-quantitative METAL counts (still not transcript counts -- they carry detector and spillover effects, not Poisson sampling). All five share the targeted-panel and segmentation traps above. When the data is specifically IMC or MIBI and the question is the end-to-end processing workflow (spillover compensation, segmentation execution, FlowSOM phenotyping on metal channels), defer to imaging-mass-cytometry rather than reimplementing it here.

## Transform and Normalize Intensities

**Goal:** Put marker intensities on a comparable, variance-stabilized scale and remove antibody-batch and staining confounds before phenotyping -- without imposing a count model.

**Approach:** Apply arcsinh with a per-dataset-tuned cofactor (or z-score/percentile), then correct channel spillover and batch; choose the transform deliberately, because this choice dominates the cell-type calls.

| Transform | Form | Best when | Fails / caveat |
|-----------|------|-----------|----------------|
| arcsinh (cofactor) | `arcsinh(x / cofactor)` | Mass-cytometry-like intensities (CyTOF/IMC/MIBI); compresses high values, near-linear near zero | Cofactor ~5 is a CyTOF CONVENTION (Bendall 2011), NOT auto-optimal for imaging -- too small over-expands near-zero noise into spurious populations; tune and sanity-check per dataset |
| z-score (per marker) | `(x - mean) / sd` | Cross-marker comparability for clustering | Sensitive to outliers; assumes roughly symmetric post-transform spread |
| percentile / min-max (per marker) | clip to e.g. 1st-99th pct, scale 0-1 | Gating-style cutoffs; robust to extreme bright pixels | Throws away absolute scale; per-image rescaling can erase real cross-sample differences |
| log1p-of-counts | `log(1 + x)` with NB/Poisson tooling | RNA UMI counts | WRONG for intensity -- there is no count process; imposes a model the data does not follow |

scimap's `pp.rescale` fits a per-marker two/three-component Gaussian mixture to set the 0-1 gating scale (an intensity-aware step, distinct from log1p-as-counts); for clustering, an explicit arcsinh or z-score on `adata.X` is the transparent choice.

```python
import numpy as np
import scimap as sm

# Tune the cofactor: start at 5 (CyTOF convention) but verify the near-zero
# population is not split into a phantom 'positive' cluster for each marker.
cofactor = 5
adata.layers['intensity'] = adata.X.copy()              # stash raw intensities
adata.X = np.arcsinh(adata.X / cofactor)                # variance-stabilize; NOT log1p-of-counts

# Antibody/staining-batch correction across images or staining days.
# Intensity differences between batches masquerade as biology -- correct before merging.
sm.pp.combat(adata, batch_key='imageid')                # batch_key names the confound column in .obs
```

Channel spillover (isotopic impurity and oxide/abundance-sensitivity crosstalk for metals; spectral bleed for fluorophores) creates false double-positive cells. Estimate a spillover matrix from single-stain bead controls and correct by non-negative least squares (Chevrier 2018, implemented in CATALYST/spillR). For IMC/MIBI specifically, run compensation through the imaging-mass-cytometry/data-preprocessing skill rather than reimplementing the matrix here.

## Phenotype Cells: Gating vs Clustering

**Goal:** Assign each cell a cell-type label from its marker-intensity vector.

**Approach:** Choose GATING (flow-cytometry-style positive/negative thresholds encoded as a marker workflow) when the panel has canonical lineage markers and the types are known a priori, or CLUSTERING (Leiden/PhenoGraph/FlowSOM on transformed intensities) for unsupervised discovery within the bounded panel; gating and clustering can give materially different calls, and cluster boundaries shift with the transform, cofactor, k, and segmentation spillover.

scimap gating expects a phenotype-workflow DataFrame, not a dict: first column = group, second = cell-type name, remaining columns = marker names holding `pos`/`neg`/`allpos`/`allneg`/`anypos`/`anyneg`.

```python
import pandas as pd
import scimap as sm

# Build the gating workflow (or load a CSV). 'allpos' = all listed markers must clear the gate.
workflow = pd.DataFrame([
    ['lineage', 'T_cell',     'allpos', 'allpos', 'neg',    'neg'],
    ['lineage', 'B_cell',     'allpos', 'neg',    'allpos', 'neg'],
    ['lineage', 'Macrophage', 'allpos', 'neg',    'neg',    'allpos'],
    ['lineage', 'Tumor',      'neg',    'neg',    'neg',    'neg'],
], columns=['group', 'phenotype', 'CD45', 'CD3', 'CD20', 'CD68'])

sm.pp.rescale(adata, gate=None, method='by_image')       # per-marker GMM sets the 0-1 scale; 'by_image' rescales each image separately
sm.tl.phenotype_cells(adata, phenotype=workflow, gate=0.5, label='phenotype')   # gate=0.5 after rescale
```

```python
# Unsupervised alternative: cluster the transformed intensities, then annotate clusters by marker means.
sm.tl.cluster(adata, method='leiden', resolution=1.0, label='leiden')
# A 'protein absent' cluster may simply lack the marker on the panel -- annotate against the panel, not the transcriptome.
```

Phenotyping on a targeted panel cannot discover a type whose markers are off-panel; an unexpected "negative-for-everything" cluster is often an unstained type, not a novel state. Audit phantom double-positives (e.g. CD3+CD20+) as likely segmentation spillover before treating them as biology.

## Spatial Neighborhood and Interaction Analysis

**Goal:** Quantify which phenotypes are spatial neighbors more or less than chance, and summarize recurrent cellular neighborhoods.

**Approach:** Build a spatial graph on cell centroids, then run a permutation-based neighborhood-enrichment test; for niches, summarize each cell's k-nearest-neighbor window by composition and cluster the windows (Schurch 2020 cellular-neighborhoods logic). Co-occurrence is not communication, and a "neighborhood" inherits every segmentation/normalization error upstream.

```python
import squidpy as sq

sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=10)   # imaging cells are a point cloud, not a grid -> 'generic'
sq.gr.nhood_enrichment(adata, cluster_key='phenotype')              # label-permutation null; z-scores in adata.uns
sq.pl.nhood_enrichment(adata, cluster_key='phenotype')
```

```python
# Recurrent cellular neighborhoods (niches): per-cell composition of the local window, then cluster.
sm.tl.spatial_count(adata, phenotype='phenotype', method='knn', knn=10, label='neighborhood_counts')
sm.tl.cluster(adata, method='kmeans', k=8, use_raw=False, label='neighborhood')   # k is a biological choice; report k +/- 1 sensitivity
```

The neighbor count k and the upstream phenotype calls both define the result; report the window size and show sensitivity. On a 20-100 marker panel the relevant ligand AND receptor AND cofactors are rarely all present, so multiplexed-imaging "communication" is almost always cell-type PROXIMITY (niche co-occurrence), not measured ligand-receptor co-localization -- a weaker inference than transcriptomic LR, because the LR pair was never measured.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Spurious "positive" populations near zero; clusters that split noise | Treated intensity as counts with `log1p`, or used an untuned tiny arcsinh cofactor | Use arcsinh with a per-dataset-tuned cofactor (start ~5, verify near-zero is not over-expanded), or z-score/percentile -- never log1p-of-counts |
| Phantom CD3+CD20+ (or any lineage-incompatible) double-positive cells | Channel spillover and/or segmentation lateral spillover between adjacent cells | Compensate spillover (NNLS, Chevrier 2018) and audit segmentation; treat double-positives as artifacts until proven |
| Cell types differ wildly between two runs of the same data | Normalization x clustering choice dominates calls (Hickey 2021: 20 annotations from one CODEX dataset) | Fix and report the transform, cofactor, normalization, and clustering; do not present one pipeline's calls as ground truth |
| Cross-sample comparison shows a "batch" cell type | Antibody-lot/staining/fixation intensity differences not normalized | Correct batch (combat or per-image rescale) before merging or comparing samples |
| Concluded a cell type or marker is "absent" | Read panel absence as biological absence on a bounded antibody panel | State that absence on a targeted panel is uninformative; the marker was likely not stained |
| Neighborhood/interaction result looks strong but is not reproducible | Built on bad segmentation masks; enrichment inherits the mask error | Validate segmentation (membrane-stain whole-cell, Mesmer) before trusting any spatial result |
| `sm.tl.spatial_cluster` returns one cluster or nonsense | Ran it before building the neighborhood matrix it reads | Compute `sm.tl.spatial_count` (or `sm.tl.spatial_lda`) first, then point `spatial_cluster(df_name=...)` at that result |
| `phenotype_cells` errors or mislabels everything | Passed a dict instead of the workflow DataFrame, or did not `rescale` first | Pass a group/phenotype/marker DataFrame with pos/neg/allpos codes; run `sm.pp.rescale` before phenotyping |

## Related Skills

- image-analysis - whole-cell segmentation upstream of every per-cell intensity vector (the dominant error source)
- imaging-mass-cytometry/cell-segmentation - Mesmer/Cellpose segmentation execution and error propagation for IMC/MIBI
- imaging-mass-cytometry/phenotyping - FlowSOM/Phenograph phenotyping on metal channels and the double-positive artifact
- spatial-transcriptomics/spatial-multiomics - integrating spatial proteomics with matched spatial transcriptomics (ADT/CytAssist)
- spatial-transcriptomics/spatial-statistics - permutation nulls, neighborhood enrichment, and co-occurrence shared with squidpy

## References

- Goltsev Y, Samusik N, Kennedy-Darling J, et al. (2018) Deep profiling of mouse splenic architecture with CODEX multiplexed imaging. Cell 174(4):968-981. DOI 10.1016/j.cell.2018.07.010
- Angelo M, Bendall SC, Finck R, et al. (2014) Multiplexed ion beam imaging of human breast tumors. Nature Medicine 20(4):436-442. DOI 10.1038/nm.3488
- Keren L, Bosse M, Thompson S, et al. (2019) MIBI-TOF: a multiplexed imaging platform relates cellular phenotypes and tissue structure. Science Advances 5(10):eaax5851. DOI 10.1126/sciadv.aax5851
- Giesen C, Wang HAO, Schapiro D, et al. (2014) Highly multiplexed imaging of tumor tissues with subcellular resolution by mass cytometry. Nature Methods 11(4):417-422. DOI 10.1038/nmeth.2869
- Lin J-R, Izar B, Wang S, et al. (2018) Highly multiplexed immunofluorescence imaging of human tissues and tumors using t-CyCIF and conventional optical microscopes. eLife 7:e31657. DOI 10.7554/eLife.31657
- Parra ER, Uraoka N, Jiang M, et al. (2017) Validation of multiplex immunofluorescence panels using multispectral microscopy for immune-profiling of FFPE human tumor tissues. Scientific Reports 7(1):13380. DOI 10.1038/s41598-017-13942-8
- Greenwald NF, Miller G, Moen E, et al. (2022) Whole-cell segmentation of tissue images with human-level performance using large-scale data annotation and deep learning (Mesmer/DeepCell). Nature Biotechnology 40(4):555-565. DOI 10.1038/s41587-021-01094-0
- Stringer C, Wang T, Michaelos M, Pachitariu M (2021) Cellpose: a generalist algorithm for cellular segmentation. Nature Methods 18(1):100-106. DOI 10.1038/s41592-020-01018-x
- Chevrier S, Crowell HL, Zanotelli VRT, et al. (2018) Compensation of signal spillover in suspension and imaging mass cytometry. Cell Systems 6(5):612-620. DOI 10.1016/j.cels.2018.02.010
- Bendall SC, Simonds EF, Qiu P, et al. (2011) Single-cell mass cytometry of differential immune and drug responses across a human hematopoietic continuum. Science 332(6030):687-696. DOI 10.1126/science.1198704
- Hickey JW, Tan Y, Nolan GP, Goltsev Y (2021) Strategies for accurate cell type identification in CODEX multiplexed imaging data. Frontiers in Immunology 12:727626. DOI 10.3389/fimmu.2021.727626
- Schurch CM, Bhate SS, Barlow GL, et al. (2020) Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front. Cell 182(5):1341-1359. DOI 10.1016/j.cell.2020.07.005
