---
name: bio-imaging-mass-cytometry-interactive-annotation
description: Interactive cell annotation and image QC for IMC/MIBI using napari, napari-imc, Mantis Viewer, and cytomapper, covering the pixels-to-cell-table bridge, overlaying masks to catch segmentation/spillover artifacts, inter-annotator variability as the accuracy ceiling, contrast-as-threshold, and building class-balanced ground-truth label sets. Use when manually labeling cells, generating training data for a classifier, QC-ing segmentation on the image, confirming clusters are spatially real, or choosing an annotation viewer.
tool_type: python
primary_tool: napari
---

## Version Compatibility

Reference examples tested with: napari 0.4.18+, napari-imc 0.7+, numpy 1.26+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: core napari cannot open `.mcd` -- the `napari-imc` plugin reads raw Fluidigm/Standard BioTools files in machine coordinates. napari `add_labels` edits integer masks (segmentation QC); `add_points` with a `features` table holds per-cell categorical labels. Mantis Viewer is a standalone Electron app from CANDELbio/Parker Institute (NOT the Bodenmiller group). napari has no journal paper -- cite the Zenodo DOI.

# Interactive Annotation

**"Manually annotate cell types in my IMC data"** -> Look at the image with masks overlaid to label cells, generate ground truth, and catch the artifacts a cell table hides.
- Python: `napari.Viewer` with `napari-imc` (raw `.mcd`), `add_labels`/`add_points`
- R: `cytomapper::cytomapperShiny` (gate then see the gate on tissue)

## The Single Most Important Modern Insight -- annotation is the pixels-to-cell-table bridge, and the image vetoes the table

The IMC pipeline collapses a multi-GB pixel stack into a cell-by-marker table via a segmentation mask and mean-intensity extraction, and that table has thrown away three things only the image still contains: whether the mask boundary matches a real cell, where a marker's signal physically sits (its spillover provenance), and morphology/context. None of these surface as an outlier in the table -- a merged CD3+CD20+ "cell" looks like a plausible rare double-positive, not an error -- so the table is seductive precisely because it is tidy and statistical. Annotation is the only QC step that lets pixels veto the table. The expert habit is distrust of any cell-table claim (a cluster, a double-positive population, a rare type) until it has been seen on the image with its mask overlaid. Three concrete bridge operations follow: overlay the mask on the DNA + summed-membrane channels and walk the image (the irreplaceable segmentation QC); paint clusters back onto tissue, where a real cluster forms coherent structures (a tumor nest, a T-cell zone) and an artifact cluster scatters as salt-and-pepper haze along the boundary between two real clusters (spatial incoherence is the tell); and gate in image space, not only expression space, because a biaxial gate drawn on the table is a thresholding decision made blind to the pixels. Two hard limits frame all of it: manual labels are not ground truth (expert-vs-expert concordance is only ~86%, and a substantial share of disagreements reflect genuinely ambiguous cells, so chasing >90% classifier accuracy is chasing noise above human agreement), and the display contrast limit IS a positivity threshold, so auto-scaling per image silently moves what counts as "positive" field to field.

## Viewer Landscape

| Tool | Stack | What it is for |
|------|-------|----------------|
| napari + napari-imc | Python | open raw `.mcd`, overlay channels + masks + annotation layers, scriptable |
| napari-steinpose | napari plugin | human-in-the-loop Cellpose segmentation inside napari |
| Mantis Viewer | Electron (CANDELbio/Parker) | dedicated ground-truth label + region/population curation on huge images |
| cytomapper / cytomapperShiny | R/Bioconductor | gate on up to ~24 markers and see the gated cells painted on tissue |
| cytoviewer | R/Bioconductor | interactive image + mask overlay colored by cell metadata |
| TissUUmaps 3 | browser/WebGL | whole-slide QC of marker/point overlays at 10^7+ points |
| QuPath | Java | whole-slide pathology; map clusters/cells back to tissue for validation |

## Decision Tree by Scenario

| Task | Tool | Why |
|------|------|-----|
| Open and inspect a raw `.mcd` | napari + napari-imc | only it reads IMC in machine coordinates |
| Generate ground-truth labels for a classifier | Mantis Viewer or napari points | built for population curation across large images |
| Gate a population and confirm it on tissue | cytomapperShiny | the gate-then-see-on-tissue loop |
| Whole-slide point-overlay QC | TissUUmaps | scales to 10^7 points |
| Segmentation mask QC/edit | napari `add_labels` / napari-steinpose | paint/fill integer masks |
| Confirm a cluster is real | paint cluster back on tissue (any viewer) | spatial incoherence = artifact |

## Open Raw IMC and Overlay the Mask

**Goal:** Put the raw channels, the segmentation mask, and an annotation layer in one canvas.

**Approach:** Use napari-imc for the `.mcd`, overlay the mask outlines on DNA + a summed membrane channel, and add a labels or points layer for annotation. Fix the contrast limits explicitly so "positive" means the same thing across fields.

```python
import napari

viewer = napari.Viewer()
viewer.open('slide.mcd', plugin='napari-imc')          # reads acquisitions + panoramas
viewer.add_labels(cell_masks, name='masks')            # outlines to audit boundaries
# fix contrast (a display limit IS a positivity threshold) -- record it in the protocol
viewer.add_image(dna, name='DNA', contrast_limits=[0, 20], colormap='gray', blending='additive')
napari.run()
```

## Paint Clusters Back onto Tissue

**Goal:** Decide whether a data-driven cluster is real biology or an artifact.

**Approach:** Color the mask by cluster and look: a real cluster forms coherent spatial structures; an artifact cluster scatters along the boundary between two real clusters.

```python
import numpy as np

def cluster_label_image(masks, cell_ids, cluster_of_cell):
    # paint each cell's cluster id back onto its mask; view in napari and demand spatial
    # coherence (nests/zones/sheets). Salt-and-pepper haze along a boundary = artifact.
    out = np.zeros_like(masks)
    lut = dict(zip(cell_ids, cluster_of_cell))
    for cid, cl in lut.items():
        out[masks == cid] = cl + 1
    return out
```

## Build a Class-Balanced Ground-Truth Set

**Goal:** Produce a training set that learns rare types and generalizes across batches.

**Approach:** Tissue is wildly imbalanced (hundreds vs tens-of-thousands of cells per class), so deliberately over-sample rare types when choosing what to label, and spread annotation across multiple patients/compartments -- label breadth beats label depth. Aim for hundreds-to-low-thousands of confident cells per class.

```python
import numpy as np

def sample_cells_to_annotate(cell_ids, cell_types, per_class=300, rng=None):
    # over-sample rare classes toward a target count; annotating random fields lets rare
    # types stay unlearnably sparse
    rng = rng or np.random.default_rng(0)
    picks = []
    for ct in np.unique(cell_types):
        pool = cell_ids[cell_types == ct]
        picks.append(rng.choice(pool, size=min(per_class, len(pool)), replace=False))
    return np.concatenate(picks)
```

## Per-Trap Failure Modes

### Trusting a table double-positive
**Trigger:** a CD3+CD20+ population from the cell table. **Mechanism:** under-segmentation or lateral spillover makes a chimeric vector that looks like a plausible rare type. **Symptom:** a "novel doublet lineage". **Fix:** overlay mask + both channels on those exact cells; two abutting nuclei means a segmentation artifact.

### Per-image auto-contrast while annotating
**Trigger:** letting the viewer auto-scale each field. **Mechanism:** the contrast limit is a positivity threshold; auto-scaling moves it per field. **Symptom:** "positive" drifts image to image; inconsistent labels. **Fix:** fix and record contrast limits; apply the same transform to every annotator and image.

### Chasing >90% classifier accuracy
**Trigger:** treating manual labels as ground truth. **Mechanism:** expert-vs-expert concordance is ~86%, and many disagreements reflect genuinely ambiguous cells. **Symptom:** overfitting to one annotator's noise. **Fix:** use multi-annotator consensus for the evaluation set; report inter-annotator agreement as the ceiling.

### Annotating one ROI deeply
**Trigger:** labeling many cells in a single field. **Mechanism:** batch/staining variation between images is a top failure mode. **Symptom:** the classifier generalizes only to that ROI. **Fix:** spread annotation across patients/images and compartments; breadth over depth.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Inter-annotator concordance ~86% | Amitay 2023 *Nat Commun* 14:4302 | the realistic accuracy ceiling; many disagreements are genuinely ambiguous |
| Hundreds-to-low-thousands confident cells/class | Amitay 2023; Shaban 2024 | enough to train; rare classes and inter-image variation bind, not total count |
| Over-sample rare classes + Poisson-resample augmentation | Amitay 2023 *Nat Commun* 14:4302 | tissue is imbalanced; signal is ion counts |
| Labels NOT harvested from clustering | annotation hygiene | clustering-derived labels re-import the double-positive artifact |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `.mcd` will not open in napari | core napari has no IMC reader | install and use the `napari-imc` plugin |
| Annotation layer behaves unexpectedly | wrong layer type | `add_labels` for masks, `add_points` (with `features`) for per-cell labels |
| Cluster looks tight but is biologically odd | spillover/segmentation artifact | paint it on tissue; demand spatial coherence |
| Rare cell type unlearnable | random-field annotation | deliberately over-sample rare types; augment |
| Classifier plateaus below expectation | exceeding inter-annotator agreement | accept the ~86% ceiling; consensus-label the evaluation set |

## References

- napari contributors. 2019. napari: a multi-dimensional image viewer for Python. Zenodo. doi:10.5281/zenodo.3555620. — no journal paper exists; cite the Zenodo DOI.
- Amitay Y, Bussi Y, Feinstein B, Bagon S, Milo I, Keren L. 2023. CellSighter: a neural network to classify cells in highly multiplexed images. *Nat Commun* 14:4302. — inter-annotator concordance; ground-truth and augmentation.
- Shaban M, et al. 2024. MAPS: pathologist-level cell type annotation from tissue images through machine learning. *Nat Commun* 15:28. — annotation scale and class imbalance.
- Geuenich MJ, Hou J, Lee S, et al. 2021. Automated assignment of cell identity from single-cell multiplexed imaging and proteomic data. *Cell Syst* 12(12):1173-1186.e5. — marker-prior labels as expert annotation.
- Bankhead P, Loughrey MB, Fernandez JA, et al. 2017. QuPath: Open source software for digital pathology image analysis. *Sci Rep* 7:16878. — whole-slide validation.
- Pielawski N, Andersson A, Avenel C, et al. 2023. TissUUmaps 3: Improvements in interactive visualization, exploration, and quality assessment of large-scale spatial omics data. *Heliyon* 9(5):e15306. — whole-slide point QC.
- Chevrier S, Crowell HL, Zanotelli VRT, et al. 2018. Compensation of Signal Spillover in Suspension and Imaging Mass Cytometry. *Cell Syst* 6(5):612-620.e5. — channel spillover (distinct from lateral spillover caught visually).

## Related Skills

- cell-segmentation - mask overlay is the irreplaceable segmentation QC
- phenotyping - annotation supplies labels/priors and confirms clusters are real
- quality-metrics - the image catches artifacts that table statistics cannot
- data-preprocessing - contrast/transform choices mirror preprocessing thresholds
- spatial-analysis - spatial coherence of a painted cluster validates it
