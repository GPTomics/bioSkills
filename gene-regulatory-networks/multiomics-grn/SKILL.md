---
name: bio-gene-regulatory-networks-multiomics-grn
description: Build enhancer-driven gene regulatory networks (eGRNs) by integrating single-cell RNA-seq and ATAC-seq using SCENIC+, CellOracle base GRNs, Pando, FigR, DIRECT-NET, TRIPOD, and scMEGA. Covers the accessibility-defines-enhancers principle, peak-to-gene linking and its cell-composition confound, the paired-vs-unpaired decision, and TF-region-gene eRegulon triplets. Use when analyzing 10x multiome or paired/unpaired scRNA+scATAC to infer cis-regulatory GRNs. For RNA-only regulons see scenic-regulons; for in silico TF perturbation see perturbation-simulation.
tool_type: python
primary_tool: SCENIC+
---

## Version Compatibility

Reference examples tested with: SCENIC+ (current Snakemake workflow), pycisTopic 2.0+, pycistarget 1.0+, scanpy 1.10+, MACS3 3.0+; FigR/Signac/ArchR (R).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

SCENIC+ has undergone major API churn: the manual-object API (cisTopicObject, separate pycistarget, SCENICPLUS objects) is superseded by a Snakemake workflow (`scenicplus init_snakemake`). Pre-2024 tutorials are stale; verify against scenicplus.readthedocs.io before coding.

# Multiomics GRN Inference

**"Build an enhancer-driven gene regulatory network from my multiome data"** -> Integrate scRNA-seq and scATAC-seq to identify eRegulons: transcription-factor -> enhancer-region -> target-gene triplets that link TF motif occupancy in accessible chromatin to gene expression.
- Python: SCENIC+ Snakemake pipeline (pycisTopic -> pycistarget -> eRegulons)
- Python: CellOracle base GRN (motif scan in Cicero co-accessible regions); R: Pando / FigR

## The Single Most Important Modern Insight -- Accessibility Defines the Candidate Enhancers, but Every Inference Arrow Leaks

The advance of multiomic GRN inference over expression-only methods is **where the candidate regulatory regions come from**: expression-only tools can only search promoter-proximal motifs, while scATAC nominates the actual distal enhancers active in these cells -- and most cell-type-specific regulatory information is distal. So an eGRN edge is a **triplet (TF -> region -> gene)**, with the region as the mechanistic anchor available for validation (ChIP/CUT&Tag, CRISPRi, reporter). But the reasoning chain leaks at every step: motif-present does not mean the TF binds (motifs are short, degenerate, and shared across a family -- the model cannot tell GATA1 from GATA2); accessible does not mean this TF holds the peak open; and peak-gene correlation does not mean the peak controls the gene. Treat an eGRN as a **prioritized hypothesis list**, not a wiring diagram.

The hardest and least-appreciated step is **peak-to-gene linking, which is confounded by cell-type composition**: across a heterogeneous dataset, any peak open in a cell type correlates with any gene expressed in that same type, whether or not the peak regulates it -- cell identity is a massive shared latent factor. Genome-wide peak-gene correlation therefore mostly recovers co-marker pairs. Mitigations (distance window, within-cell-type or GC/accessibility-matched null, requiring a motif, requiring the TF to be co-expressed) reduce but never eliminate it. Metacell aggregation, needed to beat scATAC sparsity, then introduces pseudo-replication: metacell-derived p-values are not calibrated significance and should be treated as ranking scores.

## eGRN Method Taxonomy

| Method | Citation | Approach | Data regime | Note |
|--------|----------|----------|-------------|------|
| SCENIC+ | Bravo Gonzalez-Blas 2023 *Nat Methods* | topics -> motif enrichment -> region-to-gene & TF-to-gene GBM -> eRegulons | paired or separate | reference eGRN method; heavy (Snakemake, cluster job) |
| CellOracle base GRN | Kamimoto 2023 *Nature* | motif scan in Cicero co-accessible regions; prebuilt base GRNs | scRNA alone OK | base GRN is a prior; feeds perturbation-simulation |
| Pando | Fleck 2023 *Nature* | regression with TF x peak-accessibility interaction term | paired (Seurat) | regions = peaks intersect conserved/annotated CREs |
| FigR | Kartha 2022 *Cell Genomics* | DORCs (genes with many correlated peaks) -> TF-DORC scores | SHARE-seq/paired | `pairCells` for unpaired (pairing caps quality) |
| DIRECT-NET | Zhang 2022 *Sci Adv* | XGBoost CRE-gene importance -> TF via motif | paired or scATAC alone | |
| TRIPOD | Jiang 2022 *Cell Syst* | nonparametric TF-peak-gene trio test, matched/conditional | paired | strong false-positive control |
| scMEGA | Li 2023 *Bioinform Adv* | integrate -> trajectory -> TF-gene network | paired, trajectory | lighter-weight |
| GLUE | Cao & Gao 2022 *Nat Biotechnol* | graph-linked latent embedding | **unpaired/diagonal** | run first to integrate, then a paired method |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Paired 10x Multiome / SHARE-seq, want full eGRN | SCENIC+ | reference method; eRegulon triplets + activity |
| Paired data, lighter/faster | Pando or DIRECT-NET | single-regression / XGBoost, less infrastructure |
| Rigorous false-positive control on trios | TRIPOD | conditional/matched testing removes the composition confound |
| scRNA-seq only (no ATAC) | -> CellOracle prebuilt base GRN | base GRN is a prior; no paired data needed |
| Unpaired scRNA + scATAC (separate experiments) | GLUE to integrate first, then FigR/SCENIC+ | computed pairing caps all downstream link confidence |
| RNA-only regulons, no enhancers needed | -> scenic-regulons | promoter-proximal motif pruning suffices |
| Goal is in silico TF perturbation | -> perturbation-simulation | CellOracle/Dynamo simulate; build the base GRN here |

## SCENIC+ Pipeline (current Snakemake workflow)

**Goal:** Assemble eRegulons (TF -> enhancer -> gene) from paired or separate scRNA + scATAC.

**Approach:** Topic-model the ATAC with pycisTopic, call consensus peaks from per-cell-type pseudobulk (so cell-type labels are needed before peak calling), run pycistarget motif enrichment, then region-to-gene and TF-to-gene GBM regression to build triplets; orchestrate with Snakemake.

```python
# Initialize and run the Snakemake workflow (the supported modern entry point).
# init_snakemake scaffolds a config.yaml pointing at the scATAC fragments, the scRNA
# AnnData, the motif/cisTarget databases, and the cell-type annotation used for pseudobulk.
# scenicplus init_snakemake --out_dir scenicplus_run
# (edit scenicplus_run/Snakemake/config/config.yaml, then run from inside that dir:)
# cd scenicplus_run/Snakemake && snakemake --cores 16

# Region-to-gene search space defaults to min 1kb / max 150kb from the gene, capped at the
# nearest neighboring gene's promoter -- narrower than the +/-500kb used by ArchR/Signac.
```

```python
# Inspect the resulting eRegulons (TF -> region -> gene triplets). The output filename and
# directory are set in config.yaml (output_data); the direct (high-confidence) and extended
# (motif-similarity-inferred) tables are written there. The exact spelling has varied across
# versions, so resolve it by glob rather than hard-coding.
import glob, pandas as pd
ereg_file = glob.glob('scenicplus_run/**/eRegulon*direct*.tsv', recursive=True)[0]
eregulons = pd.read_csv(ereg_file, sep='\t')
summary = (eregulons.groupby('TF')
           .agg(n_regions=('Region', 'nunique'), n_genes=('Gene', 'nunique'))
           .sort_values('n_genes', ascending=False))
# Direct vs extended is a motif-to-TF annotation CONFIDENCE distinction, not topology:
# direct = curated/orthology; _extended adds motif-similarity-inferred (larger, noisier).
```

## CellOracle Base GRN (works without paired multiome)

**Goal:** Build a base GRN -- the candidate TF -> gene scaffold -- from accessibility, as a prior for context-specific modeling.

**Approach:** Define active regions by Cicero co-accessibility (in R), scan them for TF motifs, and format the result as the base GRN; or load a prebuilt base GRN and skip ATAC entirely.

```python
import celloracle as co
import pandas as pd

# Custom base GRN: peaks already filtered to Cicero co-accessible, promoter-linked regions.
peaks = pd.read_parquet('processed_peak_file.parquet')   # columns: peak_id, gene_short_name
tfi = co.motif_analysis.TFinfo(peak_data_frame=peaks, ref_genome='hg38')
tfi.scan(fpr=0.02)                                       # motif FPR
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary()
base_grn = tfi.to_dataframe()

# Or skip ATAC: prebuilt base GRN as a prior (CellOracle ships ~10 species + mouse atlas).
# base_grn = co.data.load_mouse_scATAC_atlas_base_GRN()
```

The base GRN constrains which TF -> gene edges are allowed; the context-specific weights are then learned per cluster in perturbation-simulation.

## FigR (paired, DORC-based, R)

**Goal:** Identify domains of regulatory chromatin (DORCs) and the TFs that regulate them.

**Approach:** Correlate peaks to genes, call DORCs (genes with an unusually large number of significant peaks), then score TF-DORC associations from motif enrichment plus TF expression correlation.

```r
library(FigR)
# Step 1: peak-gene correlations (smoothed over KNN metacells for sparsity)
cisCor <- runGenePeakcorr(ATAC.se = atac_se, RNAmat = rna_mat,
                          genome = 'hg38', nCores = 8, p.cut = 0.05)
# Step 2: DORC scores then TF-DORC regulation scores (signed: activator/repressor)
dorcMat <- getDORCScores(atac_se, cisCor, geneList = unique(cisCor$Gene), nCores = 8)
figR <- runFigRGRN(ATAC.se = atac_se, dorcTab = cisCor, dorcMat = dorcMat,
                   rnaMat = rna_mat, genome = 'hg38', nCores = 8)
```

## Per-Method Failure Modes

### Cell-composition confound in peak-gene linking
**Trigger:** genome-wide peak-gene correlation with no within-type control. **Mechanism:** any peak open in a cell type correlates with any gene expressed in it. **Symptom:** "links" that are just cell-type co-markers. **Fix:** restrict to the distance window, use a matched null (Signac) or within-type correlation, and require a motif.

### Metacell p-value laundering
**Trigger:** reporting astronomically small p-values from KNN-smoothed metacells. **Mechanism:** metacells are not independent (overlapping cells). **Symptom:** millions of "significant" links. **Fix:** treat metacell p-values as ranking scores; apply FDR honest about pseudo-replication.

### Motif-family overclaiming
**Trigger:** naming a specific TF (GATA1) from a family motif. **Mechanism:** paralogs share near-identical motifs. **Symptom:** confident single-TF claims where only a family is detectable. **Fix:** report the motif/family and require orthogonal evidence to single out a member.

### Unstated/wrong pairing in unpaired data
**Trigger:** computing "joint" correlations on separately-measured scRNA and scATAC. **Mechanism:** cells were computationally matched, not co-measured. **Symptom:** no pairing method or accuracy reported. **Fix:** integrate with GLUE/anchors first; report pairing quality; treat links as upper-bounded by it.

### Treating eGRN edges as ground truth
**Trigger:** a published wiring diagram with no orthogonal validation. **Mechanism:** accessibility != binding != regulation, and motif/proximity validation is circular. **Symptom:** no ChIP/CRISPRi/perturbation check; validation uses the same motif DB used to build the net. **Fix:** validate against independent ChIP-seq or perturbation data.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Region-to-gene window: 1kb-150kb (SCENIC+) | Bravo Gonzalez-Blas 2023 | capped at nearest gene promoter; narrower than +/-500kb conventions |
| Peak-gene window +/-500kb (ArchR/Signac/Cicero) | tool defaults | distal enhancer reach; pair with a matched null |
| MACS3 `--keep-dup all`, BEDPE/shift-extsize | ATAC convention | fragment-based peak calling for the region universe |
| Motif FPR ~0.02 (CellOracle scan) | CellOracle default | motif-match false-positive rate |
| eRegulon: direct vs _extended | SCENIC+ | direct = curated motif2TF; extended adds inferred (noisier) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| stale SCENIC+ API errors | following pre-2024 manual-object tutorials | use the Snakemake workflow; check scenicplus.readthedocs.io |
| consensus peak step fails | no cell-type labels for pseudobulk | annotate cell types before peak calling |
| feather DB format error | region-based vs gene-based DB mismatch / old ctxcore | use region-based DBs for ATAC; align versions |
| huge spurious link set | composition confound + no matched null | window + within-type/matched-null + motif requirement |
| empty eRegulons | species/assembly/motif-collection mismatch | match genome, motif DB, and gene IDs |

## References

- Bravo Gonzalez-Blas C, et al. 2023. SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. *Nat Methods* 20(9):1355-1367.
- Kamimoto K, et al. 2023. Dissecting cell identity via network inference and in silico gene perturbation (CellOracle). *Nature* 614(7949):742-751.
- Fleck JS, et al. 2023. Inferring and perturbing cell fate regulomes in human brain organoids (Pando). *Nature* 621:365-372.
- Kartha VK, et al. 2022. Functional inference of gene regulation using single-cell multi-omics (FigR). *Cell Genomics* 2(9):100166.
- Zhang L, Zhang J, Nie Q. 2022. DIRECT-NET. *Sci Adv* 8(22):eabl7393.
- Jiang Y, et al. 2022. TRIPOD: nonparametric single-cell multiomic trio characterization. *Cell Syst* 13(9):737-751.e4.
- Li Z, et al. 2023. scMEGA. *Bioinform Adv* 3(1):vbad003.
- Cao ZJ, Gao G. 2022. Multi-omics integration and regulatory inference with graph-linked embedding (GLUE). *Nat Biotechnol* 40(9):1458-1466.
- Pliner HA, et al. 2018. Cicero: cis-regulatory DNA interactions from single-cell accessibility. *Mol Cell* 71(5):858-871.e8.

## Related Skills

- scenic-regulons - RNA-only regulon inference (promoter-proximal motif pruning)
- perturbation-simulation - in silico TF perturbation built on a CellOracle base GRN
- coexpression-networks - undirected co-expression baseline
- single-cell/scatac-analysis - scATAC preprocessing with Signac/ArchR
- atac-seq/atac-peak-calling - peak calling for the region universe
- chip-seq/motif-analysis - motif enrichment and TF binding for validation
