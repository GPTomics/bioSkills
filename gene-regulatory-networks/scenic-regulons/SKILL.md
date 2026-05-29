---
name: bio-gene-regulatory-networks-scenic-regulons
description: Infer transcription factor regulons from single-cell RNA-seq with pySCENIC by combining GRNBoost2 co-expression, cisTarget motif-enrichment pruning, and AUCell per-cell activity scoring. Covers the motif-pruning-as-directionality principle, regulon specificity scoring, run-to-run stability, and database/species matching. Use when identifying TF regulons, scoring TF activity per cell, finding master regulators of cell identity, or comparing regulon activity across conditions. For enhancer-driven multiomic GRNs see multiomics-grn; for bulk inference and VIPER protein-activity see grn-inference.
tool_type: python
primary_tool: pySCENIC
---

## Version Compatibility

Reference examples tested with: pySCENIC 0.12+, ctxcore 0.2+, arboreto 0.1.6+, scanpy 1.10+, loompy 3.0+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The motif-DB machinery lives in `ctxcore`; a ctxcore/feather-format version mismatch is the most common silent failure. pySCENIC is most reliable on a dedicated Python 3.10 environment.

# SCENIC Regulons

**"Identify transcription factor regulons and score TF activity from my scRNA-seq data"** -> Run the pySCENIC three-step pipeline: infer TF-target co-expression with GRNBoost2, prune to direct targets by cis-regulatory motif enrichment with cisTarget, then score per-cell regulon activity with AUCell.
- CLI: `pyscenic grn` -> `pyscenic ctx` -> `pyscenic aucell`
- Python: `arboreto_with_multiprocessing.py` for the GRN step (avoids the dask breakage)

## The Single Most Important Modern Insight -- Motif Pruning Is What Converts Co-expression into Directed Regulation

Step 1 (GRNBoost2) produces **undirected co-expression only** -- it is no better than WGCNA and inherits all of co-expression's confounding (indirect edges, batch, cell-cycle). The entire conceptual payload of SCENIC is **Step 2 (cisTarget)**: for each module it asks whether the candidate TF's binding motif is significantly enriched (NES >= 3.0) in the cis-regulatory space of the module's targets, and keeps only the targets in the motif's leading edge. This (a) imposes a mechanistic prior -- the TF can physically bind near its retained targets, (b) breaks the symmetry of co-expression into a TF -> target direction, and (c) discards indirect targets. **A "regulon" is by definition only the post-cisTarget TF plus its direct targets.** Modules that were never pruned are co-expression modules, and calling them regulons misuses the word.

The second non-obvious consequence is **AUCell: regulon activity is not TF expression.** AUCell ranks genes within each cell and computes the area under the recovery curve for the regulon's gene set, so activity can be high even when the TF's own mRNA is dropout-zero (TF transcripts are sparse). Showing TF *expression* in place of regulon AUC -- or "validating" activity by its correlation with TF expression -- misses the method's point and is circular. SCENIC regulons remain motif-supported co-expression: a strong, directed hypothesis worth a knockdown, not proof of causal regulation.

## Pipeline Taxonomy

| Step | Tool | Produces | Key parameter | Watch out for |
|------|------|----------|---------------|---------------|
| 1. GRN | GRNBoost2 (or GENIE3) | TF-target co-expression adjacencies | `--seed`, `--num_workers` | stochastic; not reproducible without a fixed seed |
| 2. Prune | cisTarget (ctxcore) | regulons (direct targets) | `--nes_threshold 3.0`, `--rank_threshold 5000` | feather DB + motif2TF version must match |
| 3. Score | AUCell | per-cell regulon activity (AUC) | `--auc_threshold 0.05` | this is the top-fraction, NOT the binarization cut |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| scRNA-seq, want TF regulons + per-cell activity | pySCENIC grn/ctx/aucell | the canonical workflow |
| GRN step hangs / KilledWorker | `arboreto_with_multiprocessing.py` | arboreto's dask backend breaks on newer dask |
| Need reproducible regulons | run GRN 10-100x, keep links recurring >80% | GRNBoost2/GENIE3 are stochastic |
| Which regulons mark a cell type | Regulon Specificity Score (RSS) | JSD-based specificity, not just magnitude |
| Paired scRNA + scATAC available | -> multiomics-grn (SCENIC+) | accessibility defines enhancers; eRegulons add the region layer |
| Bulk RNA-seq / want protein activity | -> grn-inference (ARACNe + VIPER) | SCENIC is single-cell; VIPER reads TF activity from bulk |
| Compare activity across conditions | run SCENIC once on the integrated object | raw AUC is population-relative; batch survives into regulons |

## Required Databases

cisTarget needs three matched resources: ranking database(s), motif-to-TF annotations, and the TF list -- all the same species/assembly/symbol namespace. Download from `resources.aertslab.org/cistarget/`.

```bash
# Human hg38 gene-based rankings (~1.5 GB each). Run ctx with BOTH search-space DBs
# (500bp+100bp around TSS, and TSS +/-10kb) so the leading-edge logic pools them.
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
# The ranking-DB version (mc9nr / v10) and the motif2tf annotation version MUST match.
```

## Step 1: GRN Inference (use the multiprocessing wrapper)

**Goal:** Infer TF-target co-expression adjacencies as candidate regulatory modules.

**Approach:** Run GRNBoost2 via the bundled multiprocessing script (single-node, stable) rather than the dask backend, and fix the seed so the stochastic boosting is reproducible.

```bash
# arboreto's dask backend breaks on dask>=2.x (silent hangs, KilledWorker).
# The bundled multiprocessing wrapper is the supported workaround.
python arboreto_with_multiprocessing.py \
    filtered.loom allTFs_hg38.txt \
    --method grnboost2 --output adj.tsv \
    --num_workers 8 --seed 42
```

## Step 2: Prune to Regulons by Motif Enrichment

**Goal:** Keep only TF-target links whose target genes are enriched for the TF's binding motif -- the step that confers directness and direction.

**Approach:** Load the ranking databases and motif2TF annotations, build candidate modules from the adjacencies, and run cisTarget pruning; targets surviving motif enrichment (NES >= 3.0) form the regulon.

```python
import glob, pickle, pandas as pd
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from ctxcore.rnkdb import FeatherRankingDatabase

adjacencies = pd.read_csv('adj.tsv', sep='\t')
expr = pd.read_csv('expr.csv', index_col=0)            # cells x genes
modules = list(modules_from_adjacencies(adjacencies, expr))

dbs = [FeatherRankingDatabase(f, name=f) for f in glob.glob('*.genes_vs_motifs.rankings.feather')]
# rank_threshold=5000 matches the CLI default (the prune2df Python default is 1500).
df = prune2df(dbs, modules, 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl', rank_threshold=5000)
regulons = df2regulons(df)                              # TF + direct targets only

with open('regulons.pkl', 'wb') as fh:
    pickle.dump(regulons, fh)
```

CLI equivalent for steps 1-2 (`pyscenic grn`, then `pyscenic ctx adj.tsv DB.feather --annotations_fname motifs.tbl --expression_mtx_fname filtered.loom -o reg.csv`). ctx verified defaults: `--rank_threshold 5000`, `--auc_threshold 0.05`, `--nes_threshold 3.0`, `--min_genes 20`. `--mask_dropouts` now defaults to False (matching R SCENIC); it changes the TF-target correlation sign that splits activating `(+)` from repressing `(-)` regulons, so report the setting used.

## Step 3: AUCell Per-Cell Activity

**Goal:** Score each regulon's activity in every cell, robustly to dropout.

**Approach:** Rank genes within each cell, integrate the recovery curve over the top fraction (`auc_threshold`, default 0.05 = top 5%), and emit a cell-by-regulon AUC matrix.

```python
from pyscenic.aucell import aucell

# auc_threshold = top 5% of the ranking integrated for the AUC -- NOT a binarization cut.
auc_mtx = aucell(expr, regulons, auc_threshold=0.05, num_workers=8)
auc_mtx.to_csv('auc_matrix.csv')
```

## Interpretation: Specificity and Binarization

**Goal:** Surface the regulons that define each cell type and convert activity to on/off states for clustering.

**Approach:** Use the Regulon Specificity Score (Jensen-Shannon divergence vs an idealized cell-type-specific distribution) for identity regulators, and binarize the AUC distribution (bimodal -> density threshold) for state heatmaps.

```python
from pyscenic.rss import regulon_specificity_scores
from pyscenic.binarization import binarize

cell_types = pd.read_csv('cell_types.csv', index_col=0)['cell_type']
rss = regulon_specificity_scores(auc_mtx, cell_types)     # high RSS = identity regulator
binary_mtx, thresholds = binarize(auc_mtx)                # per-regulon on/off
```

RSS (rewards specificity) and a per-cluster AUC z-score (rewards magnitude) can disagree; prefer RSS for "which regulon marks this cluster."

## Per-Method Failure Modes

### Calling unpruned modules "regulons"
**Trigger:** skipping ctx, or dropping the NES threshold to admit everything. **Mechanism:** without motif enrichment the output is co-expression, not direct regulation. **Symptom:** no motif DB/version reported; implausibly large "regulons." **Fix:** always run cisTarget; report DB + motif2TF versions and the search-space windows.

### Dask hang in the GRN step
**Trigger:** native arboreto on dask>=2.x. **Mechanism:** scheduler incompatibility. **Symptom:** silent hang or KilledWorker. **Fix:** use `arboreto_with_multiprocessing.py` (single-node, stable).

### Species / assembly mismatch
**Trigger:** mouse genes against an hg38 ranking DB, or HGNC vs MGI symbol mismatch. **Mechanism:** gene IDs do not map into the database. **Symptom:** near-empty regulon set. **Fix:** match expression IDs, ranking DB, and motif2TF to one species/assembly/namespace.

### Cross-condition AUC comparison without batch control
**Trigger:** comparing raw AUC across separately-run SCENIC analyses or strong batches. **Mechanism:** AUC is relative to the population it was ranked within; batch-driven co-expression can pass motif enrichment by chance. **Symptom:** a "condition-specific regulator" that tracks the batch. **Fix:** run SCENIC once on the integrated object; sanity-check condition regulons against batch.

### Over-reading _extended or _- regulons
**Trigger:** using `_extended` regulons for direct-binding claims, or building a story on `(-)` repressor activity. **Mechanism:** `_extended` adds orthology/similarity-inferred (low-confidence) motif annotations; negative regulons are sparse and weakly enriched. **Symptom:** direct-regulation claims from low-confidence edges. **Fix:** default to high-confidence positive regulons; treat `_extended`/`(-)` as hypotheses.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| NES >= 3.0 (motif enrichment) | Aibar 2017 / iRegulon (Janky 2014) | recovery-curve enrichment cutoff defining a supported motif |
| auc_threshold = 0.05 (top 5%) | pySCENIC default | fraction of the ranking integrated for the AUC |
| GRN reruns: keep links recurring >80% of runs | Van de Sande 2020 | GRNBoost2/GENIE3 are stochastic; recurrence = high confidence |
| min_genes = 20 per regulon | pySCENIC default | smaller target sets give unstable AUC |
| >= a few hundred cells per cell type | practical | rare clusters and doublets inflate spurious regulons |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "not a cisTarget Feather database in v1 or v2 format" | ctxcore/DB version mismatch | download current DB; align `ctxcore` version |
| empty regulon set | species/assembly or symbol mismatch | match gene IDs to the DB namespace |
| different regulons each run | unset seed in GRN step | fix `--seed`; run multiple seeds and intersect |
| activity != TF expression confuses the reader | conflating regulon AUC with TF mRNA | report AUCell activity; that independence is the point |
| ctx returns nothing | missing/mismatched `--annotations_fname` | supply matching motif2TF; check DB is gene-based (not region-based) |

## References

- Aibar S, et al. 2017. SCENIC: single-cell regulatory network inference and clustering. *Nat Methods* 14(11):1083-1086.
- Van de Sande B, et al. 2020. A scalable SCENIC workflow for single-cell gene regulatory network analysis. *Nat Protoc* 15(7):2247-2276.
- Moerman T, et al. 2019. GRNBoost2 and Arboreto. *Bioinformatics* 35(12):2159-2161.
- Janky R, et al. 2014. iRegulon: cisTarget ranking-and-recovery framework. *PLoS Comput Biol* 10(7):e1003731.
- Suo S, et al. 2018. Revealing critical regulators of cell identity (Regulon Specificity Score). *Cell Rep* 25(6):1436-1445.e3.
- Huynh-Thu VA, et al. 2010. GENIE3. *PLoS ONE* 5(9):e12776.

## Related Skills

- multiomics-grn - enhancer-driven eRegulons from paired scRNA+scATAC (SCENIC+)
- grn-inference - bulk GRN inference and VIPER TF protein-activity (the Califano lineage)
- coexpression-networks - undirected co-expression modules (what step 1 produces alone)
- single-cell/clustering - cluster cells before regulon and RSS analysis
- single-cell/preprocessing - QC, doublet removal, and normalization of scRNA-seq inputs
- single-cell/doublet-detection - remove doublets that inflate spurious regulons
