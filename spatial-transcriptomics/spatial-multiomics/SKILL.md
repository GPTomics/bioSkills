---
name: bio-spatial-transcriptomics-spatial-multiomics
description: Integrates spatial RNA with a second modality (protein, ATAC, or histone marks) on spatial CITE-seq, DBiT-seq, spatial-ATAC, or Visium CytAssist data. Use when deciding vertical (same-pixel co-profiling -> WNN/MOFA joint factors) versus diagonal (serial adjacent sections -> registration via PASTE/STalign) integration; recognizing that modalities from serial sections are DIFFERENT cells so joint same-cell methods do not apply; handling a bounded antibody/feature panel where absence is uninformative; or treating a pixel/spot as a multi-cell mixture rather than a single cell.
tool_type: python
primary_tool: muon
---

## Version Compatibility

Reference examples tested with: muon 0.1+, mudata 0.2+, scanpy 1.10+, anndata 0.10+, mofapy2 0.7+, paste-bio 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Spatial Multi-omics Integration

**"Integrate my spatial RNA with protein / ATAC / a histone mark"** -> Decide whether the modalities are co-measured on the SAME pixels or come from DIFFERENT sections, then either build one joint representation or register two distinct cell populations.
- Same pixel (vertical): `muon`/`mudata` container -> `muon.tl.mofa` joint factors, or WNN-style per-modality weighting.
- Different sections (diagonal): `paste`/`paste2` (optimal-transport alignment), `STalign` (diffeomorphic registration), `GPSA` (common coordinate).

## Governing Principle

Serial-section modalities are DIFFERENT cells -- integration across them is registration, not coupling. A z-step of even 5-10 um lands on a different cell population, so RNA on section N and ATAC on section N+1 share NO cell. Aligning them (PASTE/STalign/GPSA) produces an approximate coordinate correspondence between distinct populations, never a cell-to-cell correspondence. Any cross-modality "coupling" read off a registered pair is a statistical imputation across non-identical cells, not a measurement. The single most-abused move in this subfield is presenting diagonal, registration-based integration with the confident language of same-cell co-measurement (Vandereyken 2023 *Nat Rev Genet*).

Vertical (same-pixel co-profiling) versus diagonal (different cells, must register or anchor) is therefore THE decision, and it is categorical, not a tuning choice. Same-pixel platforms (the Fan-lab DBiT family, spatial-CITE-seq, Visium CytAssist Gene+Protein) justify joint same-cell methods -- WNN, MOFA, totalVI-style models -- because pairing is real. Applying WNN/MOFA diagonally, across serial sections, silently violates the matching assumption: the math runs and returns factors and weights that look meaningful but encode registration artifacts as if they were joint biology.

Two further traps ride on top of the fork. First, every spatial multi-omics platform measures PIXELS or spots (10-55 um) or sub-micron DNBs that must be binned -- "single-cell multi-omics in space" is almost always a downstream binning/segmentation CLAIM, not a property of the measurement, and a multi-cell pixel mixing a sender and a receiver can manufacture apparent within-cell cross-modal coupling. Second, the protein (or ATAC, or histone-mark) side is a TARGETED panel: the antibody set is chosen a priori, so a protein "absent" was likely never in the panel, exactly as a targeted RNA panel bounds what RNA can be detected -- absence is uninformative on either side.

## The Integration Regime Fork

| Regime | What is shared | Example assay | Correct approach | Fails when |
|---|---|---|---|---|
| Vertical / same-pixel | Same pixels (real pairing) | spatial-CITE-seq, DBiT-seq, Visium CytAssist, spatial epigenome-transcriptome | Joint factor / weighted-graph: MOFA, WNN, totalVI-style | Treating a multi-cell pixel as one cell; ignoring panel bound |
| Diagonal / serial-section | Nothing (distinct cells) | RNA on section N + ATAC on section N+1; two technologies on adjacent slices | Spatial REGISTRATION: PASTE/PASTE2, STalign, GPSA | Feeding registered pairs into WNN/MOFA as if same-cell |

When the regime is ambiguous (a vendor markets "co-profiling" but the modalities were run on adjacent sections), default to diagonal: assume different cells until same-pixel co-capture is documented.

## Spatial Multi-omics Platform Table

| Platform | Modalities co-measured | Same-pixel vs serial | Resolution | Integration approach | Best when / fails when |
|---|---|---|---|---|---|
| spatial-CITE-seq (Liu/Fan 2023) | ~100s proteins (ADTs) + whole transcriptome | Same pixel | 20-25 um pixels | Vertical: MOFA/WNN on MuData | Whole-transcriptome RNA + real protein pairing; fails as single cells (pixel = several cells), protein bounded by ADT panel |
| DBiT-seq (Liu/Fan 2020) | mRNA + protein (antibody DNA tags) | Same pixel | 10/25/50 um pixels | Vertical: MOFA/WNN | True co-capture; even 10 um pixel is 1-several cells, no segmentation |
| spatial-ATAC-seq (Deng/Fan 2022) | Chromatin accessibility (single modality) | Same pixel grid | 20/50 um | Integrate with RNA section by REGISTRATION (diagonal) | ATAC-only per run; sparse per-pixel; pair to RNA only across sections |
| Spatial epigenome-transcriptome (Zhang/Fan 2023) | (ATAC or one histone mark) + RNA, same pixel | Same pixel | 20 um (near-single-cell) | Vertical: joint factors | Genuinely paired per pixel; one mark per run; pixels still not segmented cells |
| Visium CytAssist Gene+Protein (10x, commercial) | Whole transcriptome + ~31-35 protein panel | Same spot | 55 um spots | Vertical BUT deconvolve: spot = many cells | Same-spot pairing; protein limited to validated panel; no peer-reviewed primary paper |
| Stereo-CITE-seq (BGI, preprint) | mRNA + protein on Stereo-seq array | Same array | sub-micron DNBs (must bin) | Vertical after binning | Flag as not peer-reviewed; protein sensitivity under-characterized |

Methods and platforms move fast here; verify the current co-capture claim and the recommended joint method against the vendor and tool docs before committing, especially whether a "multiomics" product co-captures on one pixel or runs adjacent sections.

## Vertical: Build a Same-Pixel MuData (RNA + Protein)

**Goal:** Assemble one MuData whose RNA and protein modalities are indexed on the SAME pixels, the precondition for any same-cell joint method.

**Approach:** Wrap each modality as an AnnData, share one pixel index and the spatial coordinates, then intersect observations so every pixel carries both modalities before joint modeling.

```python
import muon as mu
import mudata as md
import scanpy as sc

rna = sc.read_h5ad('spatial_cite_rna.h5ad')        # pixels x genes, whole transcriptome
prot = sc.read_h5ad('spatial_cite_adt.h5ad')       # SAME pixels x bounded ADT panel

# both modalities MUST be indexed on identical pixel barcodes -- this is what makes pairing real
mdata = md.MuData({'rna': rna, 'prot': prot})
mdata.obsm['spatial'] = rna.obsm['spatial']        # one shared coordinate frame for both modalities

sc.pp.normalize_total(mdata['rna']); sc.pp.log1p(mdata['rna'])
sc.pp.highly_variable_genes(mdata['rna'])
# ADT is a targeted panel: CLR across pixels, not log1p-of-counts; absence of a marker is uninformative
mu.prot.pp.clr(mdata['prot'])

mu.pp.intersect_obs(mdata)                          # keep only pixels present in BOTH modalities
```

## Vertical: Joint MOFA Factors Across Modalities

**Goal:** Learn interpretable latent factors shared across RNA and protein, with per-modality variance explained, on co-measured pixels.

**Approach:** Run MOFA on the MuData, then read the joint embedding and inspect which factors are RNA-driven versus protein-driven before clustering or spatial mapping.

```python
# MOFA assumes the modalities are matched on the same pixels -- valid ONLY because this is same-pixel data
mu.tl.mofa(mdata, n_factors=10, use_var='highly_variable', outfile=None)   # writes mdata.obsm['X_mofa']

sc.pp.neighbors(mdata, use_rep='X_mofa')
sc.tl.leiden(mdata)

# mu.tl.mofa trains in place (returns None); write the model with outfile= and
# read per-modality variance explained back with mofax. A factor that is ~100%
# one modality is not "multi-omic" structure -- it is that modality's own signal.
# import mofax; m = mofax.mofa_model('mofa_model.hdf5'); m.get_variance_explained()
```

WNN is the alternative joint method when interpretable factors are not needed: build a per-modality reduction, then learn per-pixel modality weights. A single modality dominating the weights is the same red flag as in single-cell CITE-seq (see single-cell/multimodal-integration). The protein panel is bounded, so a pixel scoring "negative" for a marker may simply lack that antibody, not the protein.

## Diagonal: Register Serial Sections (Different Cells)

**Goal:** Place two adjacent-section datasets in a common coordinate frame so that NEARBY does not require SAME-CELL.

**Approach:** Compute an optimal-transport alignment between the two slices over expression plus physical distance, then stack them on shared coordinates -- the output is a coordinate map between DISTINCT cell populations, never a cell correspondence.

```python
import paste as pst

# sliceA and sliceB are AnnData from ADJACENT sections -- DIFFERENT cells, not the same tissue plane
pi = pst.pairwise_align(sliceA, sliceB)            # transport plan between spots, NOT a cell-to-cell map

# stack onto a shared frame for joint visualization / neighborhood comparison only
new_slices = pst.stack_slices_pairwise([sliceA, sliceB], [pi])
# any cross-modality "coupling" inferred here is imputed across non-identical cells -- label it as such
```

PASTE2 handles partial overlap between sections; STalign performs diffeomorphic (LDDMM) registration; GPSA learns a Gaussian-process common coordinate. All produce a coordinate correspondence, not a cell correspondence -- do not feed the registered pair into WNN/MOFA as though the pixels were paired.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| WNN/MOFA "joint" factors that look biological but irreproducible across replicate sections | Ran a same-cell joint method across SERIAL sections (different cells) | Use registration (PASTE/STalign/GPSA); reserve WNN/MOFA for same-pixel data only |
| Cross-modal "coupling" reported as same-cell co-regulation from adjacent slices | Treated a registered coordinate map as a cell-to-cell correspondence | State the regime is diagonal; report coupling as imputed across non-identical cells, not measured |
| Marker called "absent" / cell type "missing" in the protein modality | Antibody was never in the bounded ADT panel | Treat panel absence as uninformative; check the panel manifest before any negative claim |
| A factor or weight is ~100% one modality but called "multi-omic" | Per-modality variance not inspected; one denser modality dominates | Report per-modality variance explained / per-pixel weights; down-weight or denoise the dominant modality |
| ADT/protein values explode after log1p | Modeled a targeted protein panel as RNA counts | Use CLR (or arcsinh) across pixels for protein, not log1p-of-UMIs |
| "Single-cell multi-omics" conclusions from pixel data | Treated a 20-55 um pixel/spot as one cell | Bin/segment explicitly and disclose it; a multi-cell pixel can manufacture within-cell cross-modal coupling |
| MOFA / WNN errors on pixel mismatch | RNA and protein indexed on non-identical pixel barcodes | `mu.pp.intersect_obs(mdata)` so every pixel carries both modalities |

## Related Skills

- spatial-transcriptomics/spatial-proteomics - protein-intensity (not count) handling, arcsinh, segmentation-dominated panels for the protein side
- spatial-transcriptomics/spatial-data-io - load each modality with the correct reader before assembling a MuData
- single-cell/multimodal-integration - WNN/totalVI/MOFA+ mechanics, ADT denoising, and the anchor-structure fork in non-spatial data
- multi-omics-integration/mofa-integration - MOFA factor interpretation and likelihood choice across modalities

## References

Liu Y, DiStasio M, Su G, et al. High-plex protein and whole transcriptome co-mapping at cellular resolution with spatial CITE-seq. Nat Biotechnol 41(10):1405-1409 (2023).
Liu Y, Yang M, Deng Y, et al. High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue (DBiT-seq). Cell 183(6):1665-1681 (2020).
Deng Y, Bartosovic M, Kukanja P, et al. Spatial profiling of chromatin accessibility in mouse and human tissues. Nature 609(7926):375-383 (2022).
Zhang D, Deng Y, Kukanja P, et al. Spatial epigenome-transcriptome co-profiling of mammalian tissues. Nature 616(7955):113-122 (2023).
Zeira R, Land M, Strzalkowski A, Raphael BJ. Alignment and integration of spatial transcriptomics data (PASTE). Nat Methods 19(5):567-575 (2022).
Argelaguet R, Arnol D, Bredikhin D, et al. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome Biol 21:111 (2020).
Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data (WNN). Cell 184(13):3573-3587 (2021).
Vandereyken K, Sifrim A, Thienpont B, Voet T. Methods and applications for single-cell and spatial multi-omics. Nat Rev Genet 24(8):494-515 (2023).
Marconato L, Palla G, Yamauchi KA, et al. SpatialData: an open and universal data framework for spatial omics. Nat Methods 22(1):58-62 (2025).
