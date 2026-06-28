---
name: bio-spatial-transcriptomics-spatial-communication
description: Maps cell-cell communication and ligand-receptor co-expression in spatial transcriptomics (Visium, Xenium, MERFISH, CosMx, Slide-seq) with Squidpy ligrec, COMMOT, stLearn, CellChat-spatial, and NicheNet. Use when choosing a method by whether spatial distance is actually modeled (squidpy ligrec is space-blind cluster-permutation vs COMMOT optimal-transport is distance-aware vs stLearn neighborhood vs CellChat-spatial filter) and by secreted-vs-contact-dependent range; choosing the ligand-receptor database knowingly because it drives the result as much as the algorithm; guarding against segmentation-spillover circularity that fabricates short-range hits; treating every ligand-receptor score as a co-expression hypothesis on a confidence ladder, not validated signaling; correcting for thousands of pair-by-cell-type-pair permutation tests; and recognizing that a targeted imaging panel rarely contains the relevant ligands and receptors so a "no communication" call is uninformative.
tool_type: python
primary_tool: squidpy
---

## Version Compatibility

Reference examples tested with: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, commot 0.0.3+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

COMMOT is the distance-aware alternative shown below; LIANA+ wraps and benchmarks several methods/databases; CellChat v2 (R) and NicheNet (R) are noted in the decision table but not coded here.

# Spatial Cell-Cell Communication

**"Map cell-cell communication in my spatial data"** -> Score ligand mRNA in a sender population against receptor mRNA in a receiver population, optionally weighted by spatial proximity, against a permutation null.
- Python: `squidpy.gr.ligrec` (CellPhoneDB engine, space-blind) or `commot.tl.spatial_communication` (optimal transport, distance-aware)
- R: CellChat v2 spatial mode, or NicheNet for downstream receiver-response linkage

## Governing Principle

Co-expression is not communication, and segmentation spillover fabricates exactly the short-range signal these tools reward.

Every standard tool (squidpy ligrec, COMMOT, CellChat, CellPhoneDB, stLearn, SpaTalk, NicheNet) ultimately computes some function of ligand mRNA in a sender cell type and receptor mRNA in a receiver cell type, against a permutation null. None measures protein, binding, secretion, diffusion, or a downstream response. A "significant" ligand-receptor pair means "ligand mRNA in A and receptor mRNA in B are spatially co-expressed above a permutation null, under database D and radius r" -- a co-expression HYPOTHESIS, full stop. The field routinely reports these correlative pairs in the language of validated signaling ("cell type A signals to B via pathway X"); a careful analyst refuses that wording and reads every call as a testable hypothesis (Armingol 2021 *Nat Rev Genet* 22:71-88).

Spatial proximity is a weak filter, not evidence. A distance constraint removes the absurd long-range calls a non-spatial method would make, but two cells being adjacent and co-expressing a pair is nowhere near sufficient for signaling. Adding a radius converts an implausible call into a plausible-looking hypothesis -- that is all it does.

The spillover circularity is the spatial-specific trap. In imaging platforms, distance-dependent transcript mis-assignment between adjacent cells (segmentation spillover) bleeds a sender's ligand transcripts into a touching receiver and vice versa, manufacturing precisely the short-range ligand-receptor co-occurrence these methods detect. A "hit" between two touching cell types can therefore be pure segmentation artifact, and spillover is strongest between adjacent heterotypic cells -- the exact pairs a communication analysis is built to find (Mitchel 2026 *Nat Genet* 58:434). Validate every short-range call against segmentation quality before believing it (see spatial-transcriptomics/image-analysis).

The database is result-determining, as much as the algorithm. CellPhoneDB, CellChatDB, and other resources differ in content, complex/subunit handling, and curation; swapping the resource changes the inferred network as much as or more than swapping the method (Dimitrov 2022 *Nat Commun* 13:3224). Two tools agreeing is often two tools sharing a database, not independent confirmation. Report the database and version as a primary methods parameter.

On capture/spot platforms there is a second spatial trap distinct from imaging spillover: a Visium spot is a 1-10-cell MIXTURE, so the "cell types" fed to a communication tool are deconvolution estimates, and a ligand and its receptor can co-reside within the SAME multi-cell spot with no inter-cellular signaling implied at all. Running ligrec on spot clusters treats regions/niches as cell types and compounds deconvolution error into the L-R call. Prefer single-cell-resolution data, or deconvolve first and restrict the analysis to spots where the sender and receiver types are estimated to be present, and read spot-level calls as the weakest rung of the confidence ladder.

## Method Decision Table

The first question is not "which tool" -- it is "does this method actually model spatial distance, and does it distinguish secreted from contact-dependent range?" Most do not.

| Method | Spatial mechanism | L-R database | Best when | Fails when |
|--------|-------------------|--------------|-----------|------------|
| squidpy `ligrec` (CellPhoneDB engine) | NONE by default -- permutes cluster labels; any cluster can "talk" to any cluster | CellPhoneDB (via omnipath) | Fast scanpy-native CellPhoneDB baseline on large data | Misused as "spatial" -- it is space-blind unless cells are pre-restricted; sensitive to clustering granularity |
| COMMOT (Cang & Nie) | Collective optimal transport; distance COST with a per-pathway cutoff; isotropic, diffusion-like | CellPhoneDB/CellChat-derived, built in | True spatial data; want competition among L-R species + sender/receiver direction maps | One characteristic length per pathway -- cannot separate secreted vs contact; sensitive to the cutoff; transport is not flux |
| stLearn cci | L-R co-expression within local neighborhoods; two-level (label + position) permutation | User-supplied (CellPhoneDB-style) | Spot/imaging hotspot maps of where a pair co-occurs | Tests spatial co-expression enrichment, not signaling; depends on neighborhood radius |
| CellChat v2 (spatial mode, R) | Distance constraint applied as a FILTER on an expression-driven mass-action score | CellChatDB (cofactor-aware; differs from CellPhoneDB) | Pathway-level aggregation, cofactor/antagonist modeling, hierarchy summaries | Mass-action over group-averaged expression is not kinetics; spatial mode still filters a non-spatial score |
| NicheNet (R) | NONE -- not spatial; uses analyst-defined sender/receiver sets | Curated integrated prior network | Linking a ligand to DOWNSTREAM target-gene response in the receiver | Not a detector; the prior is fixed/generic; a "top ligand" is a prior-weighted hypothesis |
| MISTy (R) | Multi-view random forests over juxta/para radii; reports view importances | NONE -- marker-to-marker, no L-R DB | Highly-multiplexed imaging; dissecting spatial co-variation without an L-R DB | Models correlative spatial structure, not signaling flux |

Secreted vs contact-dependent ligands need DIFFERENT ranges, and most tools apply ONE cutoff to all pairs. A juxtacrine pair (Notch-Delta, contact-only) and a diffusible chemokine have categorically different interaction lengths; a single global radius over-calls one class and misses the other. Paracrine spread is a reaction-diffusion process, not a hard radius -- interrogate any fixed cutoff (the ~500 um conventions in the literature are conveniences, not biology). Because methods and databases genuinely compete here, verify current best practice against the latest LIANA+/benchmark docs before committing.

## The Confidence Ladder

A communication claim earns confidence by climbing, not by a low p-value:

co-expression (bare L-R) < proximity-conditioned co-expression < downstream receiver-response support (NicheNet target DE up in neighboring receivers) < orthogonal protein co-localization (the ligand AND receptor protein imaged together) < perturbation (block the ligand/receptor, measure the receiver).

Almost nothing in the spatial literature reaches the perturbation tier. The single most defensible computational move is to require BOTH spatial proximity AND a coherent downstream transcriptional response in the receiver: co-expression proposes, receiver-response disposes.

## Space-Blind Baseline (squidpy ligrec)

**Goal:** Rank ligand-receptor pairs that are co-expressed across annotated cell-type pairs as a fast CellPhoneDB-style baseline.

**Approach:** Run the permutation engine over cluster labels; recognize this is space-blind -- it tests "is this pair unusual for these two cell types," not "is this pair unusually co-located." It needs cell-type annotations (see spatial-transcriptomics/spatial-domains) and fetches the database from omnipath (internet required).

```python
import squidpy as sq

adata = sq.datasets.seqfish()                            # built-in single-cell-resolution fixture with celltype labels
res = sq.gr.ligrec(
    adata,
    cluster_key='celltype_mapped_refined',
    n_perms=1000,                                        # permutation null; more = stabler p-values, slower
    threshold=0.01,                                      # min FRACTION of cells in a cluster expressing the gene -- NOT a p-value
    use_raw=False,                                       # seqfish has no .raw; default True errors here
    copy=True,
)
pvalues = res['pvalues']                                 # MultiIndex columns = (cluster_1, cluster_2), index = (ligand, receptor)
means = res['means']
```

The `threshold` argument is the expression-fraction floor inside the engine, not a significance cutoff -- mislabeling it as a p-value is a common error. The default fetches all omnipath interactions; pass `interactions=<DataFrame with 'source'/'target'>` to pin a known database/version.

## Honest Significance and Multiple Testing

**Goal:** Extract co-expression hypotheses without manufacturing a network from nominal p-values.

**Approach:** The test space is thousands of L-R pairs times every ordered cell-type pair; correct over the full space and treat survivors as ranked hypotheses, not findings.

```python
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

flat = pvalues.stack([0, 1], future_stack=True).rename('pval').reset_index()
flat = flat.dropna(subset=['pval'])
flat['padj'] = multipletests(flat['pval'].values, method='fdr_bh')[1]   # correct over the WHOLE pair x celltype-pair space
hits = flat[flat['padj'] < 0.05].sort_values('padj')
print(f'{len(hits)} co-expression hypotheses survive BH-FDR out of {len(flat)} tests')
```

Reporting top-ranked pairs at nominal p without an honest corrected null is how interaction networks get manufactured. Label permutation and position permutation answer different questions; neither asks "is there signaling."

## Distance-Aware Inference (COMMOT)

**Goal:** Score communication with spatial distance actually in the model, and respect a finite signaling range rather than letting any cluster talk to any cluster.

**Approach:** Optimal transport moves ligand "mass" to receptor "mass" across real coordinates under a per-pathway distance cost; set `dis_thr` to the signaling length scale and handle heteromeric complexes explicitly. Use micron coordinates, not pixels.

```python
import commot as ct

# database= identifiers drift across commot releases -- verify against the installed version
df_ligrec = ct.pp.ligand_receptor_database(database='CellPhoneDB_v4.0', species='human')
ct.tl.spatial_communication(
    adata,
    database_name='cellphonedb',
    df_ligrec=df_ligrec,
    dis_thr=200,                                         # signaling range in COORDINATE UNITS (um) -- one length per pathway, isotropic
    heteromeric=True,                                    # respect multi-subunit complexes (e.g. TGFBR1_TGFBR2)
)
# sender/receiver signaling stored in adata.obsm['commot-cellphonedb-sum-sender'] / '-receiver'
```

`dis_thr` is a single characteristic length applied isotropically -- it cannot distinguish a contact-only pair from a diffusing cytokine, so set it per pathway when secreted and juxtacrine pairs are both in play, and report it. Optimal transport gives a directional, competition-aware map, but "transport" is a model device, not measured flux.

## Visualize and Audit

**Goal:** Inspect top pairs while keeping the spillover and panel caveats visible.

**Approach:** Plot the ligrec dotplot for chosen sender/receiver groups, then overlay the actual ligand and receptor expression in space to eyeball whether a "hit" sits exactly at a cell-type boundary (the spillover signature).

```python
sq.pl.ligrec(res, source_groups='Endothelium', alpha=0.05, swap_axes=True)
```

If a short-range hit localizes to the seam between two touching cell types, suspect transcript spillover before signaling: re-check the segmentation, or test whether the pair survives on a re-segmented (Baysor/proseg) matrix.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "Cell type A signals to B via pathway X" written as a finding | Treating an L-R score as validated signaling | Report it as a co-expression hypothesis; climb the confidence ladder (receiver-response, protein co-localization, perturbation) |
| Short-range hit between two touching cell types that vanishes after re-segmentation | Distance-dependent transcript spillover manufactured the co-occurrence | Validate against segmentation quality; re-run on a Baysor/proseg matrix (spatial-transcriptomics/image-analysis) |
| Two tools "confirm" the same interaction | They share the same L-R database, not independent evidence | Report database + version as a primary parameter; vary the resource (Dimitrov 2022) |
| Juxtacrine pair over-called or cytokine missed | One global distance cutoff applied to secreted and contact-dependent pairs alike | Set range per signaling class; interrogate any fixed radius (~500 um is a convention) |
| `sq.gr.ligrec` results look space-aware but are not | ligrec permutes cluster labels -- it is space-blind by default | Use COMMOT/stLearn for distance-modeled inference, or pre-restrict cells to a neighborhood |
| L-R call between two cell types inside one Visium spot | A spot is a 1-10-cell mixture; "cell types" are deconvolution estimates and both genes can live in the same spot | Prefer single-cell-resolution data, or deconvolve then restrict to spots where both types are present; treat spot-level calls as the weakest evidence |
| Hundreds of "significant" pairs at nominal p | No correction over thousands of pair x cell-type-pair tests | Apply BH-FDR over the FULL test space; treat survivors as ranked hypotheses |
| Almost no genes match the database; "no communication found" | Targeted imaging panel (Xenium/MERFISH/CosMx) lacks the relevant ligands/receptors/cofactors | Absence on a panel is uninformative; check panel coverage before concluding |
| `threshold` filters nothing / errors as a p-value | `threshold` is the expression-fraction floor, not a significance cutoff | Set it as a fraction (e.g. 0.01-0.1); filter significance on `pvalues` afterward |
| `ValueError` about `.raw` in ligrec | `use_raw=True` default with no `.raw` present | Pass `use_raw=False` |

## Related Skills

- spatial-transcriptomics/image-analysis - the segmentation-spillover circularity that fabricates short-range L-R hits; validate here first
- spatial-transcriptomics/spatial-neighbors - build the spatial graph that distance-aware methods inherit
- spatial-transcriptomics/spatial-statistics - neighborhood enrichment and permutation-null co-occurrence for which-types-co-occur questions
- spatial-transcriptomics/spatial-domains - annotate the cell types that define senders and receivers
- single-cell/cell-communication - the non-spatial CellPhoneDB/CellChat/NicheNet baseline and database choice
- pathway-analysis/go-enrichment - enrich downstream receiver-response programs (note: enrichment on predicted L-R lists is circular)

## References

- Cang Z, Zhao Y, Almet AA, et al. (2023) Screening cell-cell communication in spatial transcriptomics via collective optimal transport (COMMOT). Nature Methods 20(2):218-228. DOI 10.1038/s41592-022-01728-4
- Dimitrov D, Turei D, Garrido-Rodriguez M, et al. (2022) Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nature Communications 13:3224. DOI 10.1038/s41467-022-30755-0
- Palla G, Spitzer H, Klein M, et al. (2022) Squidpy: a scalable framework for spatial omics analysis. Nature Methods 19(2):171-178. DOI 10.1038/s41592-021-01358-2
- Jin S, Guerrero-Juarez CF, Zhang L, et al. (2021) Inference and analysis of cell-cell communication using CellChat. Nature Communications 12:1088. DOI 10.1038/s41467-021-21246-9
- Browaeys R, Saelens W, Saeys Y (2020) NicheNet: modeling intercellular communication by linking ligands to target genes. Nature Methods 17(2):159-162. DOI 10.1038/s41592-019-0667-5
- Pham D, Tan X, Balderson B, et al. (2023) Robust mapping of spatiotemporal trajectories and cell-cell interactions in healthy and diseased tissues (stLearn). Nature Communications 14:7739. DOI 10.1038/s41467-023-43120-6
- Tanevski J, Ramirez Flores RO, Gabor A, et al. (2022) Explainable multiview framework for dissecting spatial relationships from highly multiplexed data (MISTy). Genome Biology 23:97. DOI 10.1186/s13059-022-02663-5
- Mitchel J, Gao T, Petukhov V, et al. (2026) Impact and correction of segmentation errors in spatial transcriptomics. Nature Genetics 58:434-444. DOI 10.1038/s41588-025-02497-4
- Armingol E, Officer A, Harismendy O, Lewis NE (2021) Deciphering cell-cell interactions and communication from gene expression. Nature Reviews Genetics 22(2):71-88. DOI 10.1038/s41576-020-00292-x
