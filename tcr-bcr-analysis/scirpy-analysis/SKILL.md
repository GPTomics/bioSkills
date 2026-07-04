---
name: bio-tcr-bcr-analysis-scirpy-analysis
description: Integrates single-cell paired TCR/BCR (10x VDJ, AIRR, dandelion, BD Rhapsody) with gene expression in an AnnData/MuData object using scirpy - chain-pairing QC, clonotype definition, clonal expansion, diversity, repertoire overlap, V(D)J usage, and VDJdb specificity. Operates on the awkward-array AIRR model (adata.obsm['airr'], accessed via get.airr after pp.index_chains), not legacy per-chain obs columns. Use when deciding clonotype definition for TCR (exact CDR3-nt identity via define_clonotypes) versus BCR (nucleotide distance clustering via define_clonotype_clusters with normalized_hamming plus same_v_gene/same_j_gene, because somatic hypermutation shatters identity clonotypes); tuning receptor_arms (all vs any), dual_ir, and within_group; filtering chain_qc categories (multichain doublets, orphan dropout, extra-VJ dual-TCR) without biasing clonal-expansion and diversity estimates; and overlaying clonality onto the transcriptomic UMAP.
tool_type: python
primary_tool: scirpy
---

## Version Compatibility

Reference examples tested with: scirpy 0.24+, scanpy 1.10+, anndata 0.10+, mudata 0.3+, awkward 2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: since scirpy 0.13 the AIRR receptor data lives as an awkward array in `adata.obsm['airr']`, NOT in per-chain `adata.obs['IR_VJ_1_*']` columns. Fields are read with `scirpy.get.airr(...)` after `scirpy.pp.index_chains(...)`. Legacy (pre-0.13) objects must be migrated with `scirpy.io.upgrade_schema()`. Paired GEX+AIRR is held in a MuData with modalities `gex` and `airr`; tool/plot functions take the MuData and namespace their output obs columns as `airr:<column>`.

# scirpy Analysis

**"Analyze my single-cell paired TCR/BCR alongside gene expression"** -> Ingest AIRR receptor records, QC chain pairing, define clonotypes, and overlay clonality on the transcriptomic embedding, all in one AnnData/MuData object.
- Python: `scirpy.io.read_10x_vdj()` / `read_airr()` / `from_dandelion()`, `scirpy.pp.index_chains()`, `scirpy.tl.chain_qc()`, `scirpy.tl.define_clonotypes()` (TCR) or `scirpy.tl.define_clonotype_clusters()` (BCR)

## The governing principles

Two choices silently bias every downstream expansion, diversity, and overlap number. State both explicitly whenever reporting a result.

Principle 1 - the clonotype definition is a choice, and TCR and BCR need different ones. TCR does not somatically hypermutate, so all progeny of a founding T cell share the exact CDR3 nucleotide sequence: `tl.define_clonotypes` (implicit `metric='identity', sequence='nt'`) on CDR3-nt plus V and J is correct and defensible. B cells DO hypermutate during affinity maturation, so lineage members are NOT identical - `tl.define_clonotypes` shatters one true BCR lineage into dozens of fake singletons and destroys expansion and diversity estimates (Gupta 2015 *Bioinformatics* 31:3356). BCR requires `tl.define_clonotype_clusters` with a distance metric (`normalized_hamming`), `sequence='nt'` (SHM acts on nucleotides), and `same_v_gene=True, same_j_gene=True` to approximate clonal lineages. Even then scirpy returns clonal CLUSTERS, not germline-rooted phylogenies - hand off to Immcantation/dandelion/Dowser for true lineages, mutation calling, and selection.

Principle 2 - single-cell chain QC is the domain-specific hard part, and blanket filtering biases clonality upward. `tl.chain_qc` labels each cell (single pair, orphan VJ/VDJ, extra VJ/VDJ, two full chains, multichain, ambiguous). Multichain and TCR+BCR-ambiguous cells are likely doublets and are excluded from clonotype definition regardless. But dropping ALL orphan and extra-chain cells is not free: large clones capture both chains more often, so orphans are enriched for singletons, and deleting them preferentially removes small clones - inflating apparent clonal expansion and deflating diversity. Match the filter to the question, and report it.

## Clonotype definition: which function

**Goal:** Pick the clonotyping approach that matches the receptor's biology.

| Question | Function | metric / sequence | Best when | Fails when |
|----------|----------|-------------------|-----------|------------|
| TCR clonal identity | `define_clonotypes` | identity / nt (implicit) | TCR (no SHM); exact founder-lineage identity | Applied to BCR - SHM fragments lineages |
| BCR clonal lineage (approx.) | `define_clonotype_clusters` | normalized_hamming / nt + same_v_gene + same_j_gene | BCR; group SHM-diverged members of one lineage | Threshold not data-derived -> chains/merges clones |
| Convergent/functional TCR clusters | `define_clonotype_clusters` | tcrdist or alignment / aa | Antigen-convergent TCRs (different nt, same specificity) | Interpreted as recombination-event counts |
| Reconcile with bulk beta/heavy-only | `define_clonotype_clusters` | receptor_arms='VDJ' | Matching single-cell to bulk TRB/IGH repertoires | Paired-chain specificity is discarded |

`pp.ir_dist` computes and caches the VJ/VDJ distance matrices; the subsequent `define_*` call MUST use the SAME `metric` and `sequence` or the cached distances silently mismatch the grouping. For BCR the `cutoff` for `normalized_hamming` is a PERCENT distance, not a nucleotide count: `cutoff=15` means 15% mismatch (~85% identity). Set it from the bimodal distance-to-nearest-neighbor histogram (within-clone mode near 0 vs between-clone mode).

## Clustering parameters and their biology

| Parameter | Options | Meaning / when to change |
|-----------|---------|--------------------------|
| `receptor_arms` | all / any / VJ / VDJ | `all` (default): BOTH VJ (alpha/light) and VDJ (beta/heavy) must match - stringent, high specificity. `any` rescues single-arm dropout but can merge distinct clones sharing only a beta (beta convergence is real). `VDJ` mimics bulk beta/heavy-only clonotyping. |
| `dual_ir` | any / primary_only / all | Handles two chains of one arm. ~30% of T cells carry two productive TRA (allelic inclusion, Padovan 1993 *Science* 262:422) - so extra-VJ is real dual-TCR, not junk. `primary_only` uses the highest-UMI chain; `any` links cells sharing any chain; `all` requires both to correspond. |
| `same_v_gene` / `same_j_gene` | False / True | Require identical V (and J) gene, not just CDR3. Two cells can convergently share a CDR3 from different V genes; requiring same V/J enforces common ancestry. Turn ON for BCR lineage stringency. |
| `within_group` | 'receptor_type' (default) / obs col | Never merge clonotypes across this grouping. Default stops a B cell and a T cell joining one clonotype; set to sample/patient to forbid cross-sample clonotypes. |

## Load VDJ and build the joint object

**Goal:** Ingest receptor contigs and pair them with gene expression in one MuData.

**Approach:** `read_10x_vdj` (or `read_airr` / `from_dandelion`) returns an AIRR AnnData; wrap it with the GEX AnnData in a MuData keyed `gex`/`airr`, then index chains before any QC or clonotyping.

```python
import scirpy as ir
import scanpy as sc
import mudata as mu

adata_gex = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
adata_airr = ir.io.read_10x_vdj('filtered_contig_annotations.csv')  # returns an AnnData, does NOT modify in place
# ir.io.read_airr(['tra.tsv', 'trb.tsv'])  # AIRR TSV from dandelion/Immcantation/airrflow
# ir.io.from_dandelion(dandelion_obj)      # round-trip a dandelion Dandelion object
# ir.io.upgrade_schema(legacy_adata)       # migrate a pre-0.13 obs-column object first

mdata = mu.MuData({'gex': adata_gex, 'airr': adata_airr})
ir.pp.index_chains(mdata)  # REQUIRED before QC/clonotyping; builds obsm['chain_indices']
```

## Chain QC and question-aware filtering

**Goal:** Categorize chain pairing and remove doublets without silently biasing clonality.

**Approach:** Run `chain_qc`, always drop multichain and TCR+BCR-ambiguous doublets, and decide orphan/extra retention by the downstream question - keep orphans for pure GEX overlay, drop them only for paired-clonotype/specificity work.

```python
ir.tl.chain_qc(mdata)  # writes obs: airr:receptor_type, airr:receptor_subtype, airr:chain_pairing
print(mdata.obs['airr:chain_pairing'].value_counts())

# Always exclude likely doublets from clonotype definition.
drop = ['multichain']
keep_types = mdata.obs['airr:receptor_type'].isin(['TCR', 'BCR'])  # exclude 'ambiguous' (TCR+BCR doublet)
paired = mdata[keep_types & ~mdata.obs['airr:chain_pairing'].isin(drop)].copy()

# For paired-clonotype/specificity analysis also require a complete receptor (drop orphans),
# but note this preferentially deletes small clones -> inflates expansion, deflates diversity.
complete = paired[paired.obs['airr:chain_pairing'].isin(['single pair', 'extra VJ', 'extra VDJ'])].copy()
```

## Define clonotypes - TCR (identity)

**Goal:** Group T cells sharing an exact CDR3-nucleotide founder rearrangement.

**Approach:** Cache identity distances, then partition; identity on CDR3-nt plus matching arms is the correct TCR clonotype.

```python
ir.pp.ir_dist(complete, metric='identity', sequence='nt', cutoff=0)
ir.tl.define_clonotypes(complete, receptor_arms='all', dual_ir='primary_only')  # writes airr:clone_id
print('TCR clonotypes:', complete.obs['airr:clone_id'].nunique())
```

## Define clonotypes - BCR (distance clusters)

**Goal:** Group SHM-diverged B cells of one lineage that identity clonotyping would shatter.

**Approach:** Use normalized Hamming distance on nucleotides within same-V/same-J partitions - this approximates a clonal lineage; identity clonotyping is WRONG for BCR.

```python
# cutoff=15 is a PERCENT distance for normalized_hamming (15% mismatch ~= 85% identity), NOT 15 nt;
# confirm from the distance-to-nearest-neighbor histogram (bimodal trough).
ir.pp.ir_dist(complete, metric='normalized_hamming', sequence='nt', cutoff=15)
ir.tl.define_clonotype_clusters(
    complete,
    sequence='nt', metric='normalized_hamming',
    receptor_arms='all', dual_ir='any',
    same_v_gene=True, same_j_gene=True,   # enforce common ancestry for lineage-grade clones
)  # writes airr:cc_nt_normalized_hamming (a clonotype-cluster id column)
# For true germline-rooted lineages, SHM, and selection: hand off to immcantation-analysis / dandelion / Dowser.
```

## Clonal expansion and diversity

**Goal:** Quantify how expanded each clone is and how diverse each group's repertoire is.

**Approach:** Bin cells by clone size, then compute per-group diversity - but remember both numbers depend entirely on the QC filter and clonotype definition above, so report them alongside.

```python
# target_col is resolved WITHIN the airr modality, so pass the bare name 'clone_id', not 'airr:clone_id'.
ir.tl.clonal_expansion(mdata, target_col='clone_id')  # bins per cell: singleton / 2 / >= 3 (breakpoints=(1, 2))
ir.pl.clonal_expansion(mdata, target_col='clone_id', groupby='airr:receptor_subtype')

# Alpha diversity per group; groupby names a full mdata.obs column, so it keeps its modality prefix.
ir.tl.alpha_diversity(mdata, groupby='gex:sample', target_col='clone_id', metric='normalized_shannon_entropy')

# Pairwise repertoire sharing (public/expanded clones, trafficking) - depth-sensitive; compare at equal depth.
ir.tl.repertoire_overlap(mdata, groupby='gex:sample', target_col='clone_id')
ir.pl.repertoire_overlap(mdata, groupby='gex:sample')
```

## Overlay clonality on the transcriptome

**Goal:** See which cell states the expanded clones occupy.

**Approach:** Cluster on GEX independently (never on receptor sequence), then color the transcriptomic UMAP by a clonality column pushed into the GEX modality.

```python
# GEX pipeline lives on mdata['gex']: normalize -> HVG -> PCA -> neighbors -> leiden -> umap (see single-cell/clustering).
mdata['gex'].obs['clonal_expansion'] = mdata.obs['airr:clonal_expansion']
sc.pl.umap(mdata['gex'], color='clonal_expansion')

# clonotype_modularity tests whether a clone's cells are more transcriptionally connected than random
# (needs sc.pp.neighbors on the GEX modality first) - distinguishes a coherent functional clone from scatter.
ir.tl.clonotype_modularity(mdata, target_col='clone_id')
```

## Gene usage and specificity

**Goal:** Summarize V(D)J segment usage and annotate antigen specificity.

**Approach:** Plot usage/spectratype directly; for specificity, match receptors against a reference database by sequence distance (not ML prediction).

```python
ir.pl.vdj_usage(mdata, full_combination=False)      # V-D-J segment flow (Sankey/ribbon)
ir.pl.spectratype(mdata, chain='VDJ_1', color='airr:receptor_subtype')  # CDR3-length distribution

# Antigen specificity by sequence match to a reference DB (reuses the ir_dist machinery).
vdjdb = ir.datasets.vdjdb()
ir.tl.ir_query(mdata, vdjdb, metric='identity', sequence='aa')
ir.tl.ir_query_annotate(mdata, vdjdb, include_ref_cols=['antigen.species', 'antigen.epitope'])
# For deeper TCR specificity modelling leave scirpy for tcrdist3 / CoNGA (see specificity-annotation).
```

## Export AIRR

**Goal:** Hand the receptor table to a bulk/interchange tool.

**Approach:** Write the AIRR modality as a standard rearrangement TSV; do NOT reconstruct it from stale per-chain obs columns (they no longer exist).

```python
ir.io.write_airr(mdata['airr'], 'scirpy_airr.tsv')
# Pull specific fields for a custom table with the get accessor, not obs indexing:
junction_vj = ir.get.airr(mdata, 'junction_aa', 'VJ_1')   # pandas Series
with ir.get.airr_context(mdata, 'junction_aa', ['VJ_1', 'VDJ_1']):
    pass  # AIRR fields temporarily materialized into obs for grouping/plotting
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| BCR lineages appear as hundreds of singletons; no expansion | Identity `define_clonotypes` used on B cells; SHM makes members non-identical | Use `define_clonotype_clusters` with `metric='normalized_hamming', sequence='nt', same_v_gene=True, same_j_gene=True` |
| `KeyError: 'IR_VJ_1_junction_aa'` / obs receptor columns missing | Pre-0.13 schema assumed; AIRR now lives in `obsm['airr']` | Access via `ir.get.airr(...)` after `pp.index_chains`; migrate legacy objects with `io.upgrade_schema()` |
| Expansion looks high, diversity looks low vs a collaborator | Blanket-filtered all orphan/extra-chain cells, deleting small clones | Keep orphans for GEX overlay; only drop them for paired-clonotype work, and report the filter |
| `define_*` gives grouping that ignores the chosen metric | `pp.ir_dist` metric/sequence differ from the `define_*` call | Match `metric` and `sequence` between `ir_dist` and `define_clonotype_clusters` |
| Clonotype/QC functions error or return nothing | `pp.index_chains` not run before QC/clonotyping | Run `ir.pp.index_chains(mdata)` immediately after building the MuData |
| Real T/B cells look receptor-negative | GEX-only cells (contig dropout) treated as VDJ-negative | Keep GEX-only cells with `NaN` clonotype for cell-state analysis; only drop VDJ-only cells failing GEX QC |
| Spurious shared/secondary chains in a hyperexpanded sample | Ambient VDJ mRNA from a dominant clone mis-assigned to droplets | Start from CellRanger `filtered_contig_annotations` (is_cell/high_confidence/productive), then drop secondary chains with very low UMI support (`duplicate_count`/`consensus_count`, e.g. < 2-3) before clonotyping |
| BCR distance clusters look degraded even with correct settings | CellRanger BCR contigs are not IMGT-numbered and include partial/nonproductive contigs | Reannotate with IgBLAST (dandelion/airrflow) before `define_clonotype_clusters`, or hand off to Immcantation |

## Related Skills

- mixcr-analysis - Process raw single-cell VDJ FASTQ
- immcantation-analysis - Proper BCR clonal lineages and SHM downstream
- specificity-annotation - Antigen-specificity clustering on single-cell clonotypes
- single-cell/data-io - Load and manage the GEX AnnData/MuData
- single-cell/clustering - Cell-state clustering to overlay clonality
- single-cell/doublet-detection - Corroborate multichain doublet calls

## References

- Sturm G, Szabo T, Fotakis G, Haider M, Rieder D, Trajanoski Z, Finotello F. Scirpy: a Scanpy extension for analyzing single-cell T-cell receptor-sequencing data. *Bioinformatics* 2020;36(18):4817-4818. doi:10.1093/bioinformatics/btaa611.
- Suo C, Polanski K, Dann E, et al. Dandelion uses the single-cell adaptive immune receptor repertoire to explore lymphocyte developmental origins. *Nature Biotechnology* 2024;42:40-51. doi:10.1038/s41587-023-01734-7.
- Gupta NT, Vander Heiden JA, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH. Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. *Bioinformatics* 2015;31(20):3356-3358. doi:10.1093/bioinformatics/btv359.
- Padovan E, Casorati G, Dellabona P, Meyer S, Brockhaus M, Lanzavecchia A. Expression of two T cell receptor alpha chains: dual receptor T cells. *Science* 1993;262:422-424. doi:10.1126/science.8211163.
- Vander Heiden JA, Marquez S, Marthandan N, et al. AIRR Community standardized representations for annotated immune repertoires. *Frontiers in Immunology* 2018;9:2206. doi:10.3389/fimmu.2018.02206.
