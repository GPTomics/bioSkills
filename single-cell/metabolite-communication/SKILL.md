---
name: bio-single-cell-metabolite-communication
description: Infers metabolite-mediated cell-cell communication from scRNA-seq by scoring enzyme-to-sensor pairs (MEBOCOST), with metabolic flux (scFEA), FBA state (Compass), and neurotransmitter (NeuronChat) alternatives. Use when studying metabolic crosstalk between cell types, predicting metabolite secretion and sensing, or deciding which metabolic-communication method fits and how speculative the result is.
tool_type: python
primary_tool: MeboCost
---

## Version Compatibility

Reference examples tested with: mebocost 1.0+, scanpy 1.10+, anndata 0.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolite-Mediated Cell Communication

**"Find which cell types exchange metabolites"** -> Infer a metabolite's presence from the cells expressing its synthesizing enzymes, then score communication to cells expressing its sensor or transporter.
- Python: `mebocost.create_obj()` -> `infer_commu()` (metabolite-sensor scoring)

## Governing Principle

Metabolite-mediated communication is a DOUBLE inference and the most speculative layer of cell-cell communication. scRNA-seq never measures metabolites; their levels are inferred from the expression of synthesizing enzymes, then a chain of further assumptions is stacked: enzyme mRNA -> enzyme protein -> enzyme ACTIVITY -> metabolic FLUX -> intracellular metabolite POOL -> SECRETION/export -> extracellular concentration in space -> import/SENSING by the receiver. Every arrow is an assumption and none is measured. On top of this sit all the ligand-receptor caveats (proxy, no spatial geometry in dissociated data, abundance and depth confounds, ambient RNA), so the output is hypothesis-generation ONLY. A defensible claim is never "cell A produces metabolite X" but "cell A expresses the machinery consistent with producing X". Validation is non-optional and means metabolomics (LC-MS), mass-spectrometry imaging (MALDI/DESI), isotope tracing, or perturbation of the enzyme or sensor - not another expression-based method. State the enzyme->flux->level->sensing chain explicitly whenever reporting a result.

## Method Decision Table

| Method | Tests what / null | Use when | Fails when |
|--------|-------------------|----------|------------|
| MEBOCOST | Metabolite SENDER->RECEIVER communication; estimates an extracellular metabolite level from producing-enzyme expression, scores enzyme->sensor pairs, permutation FDR over shuffled labels | The question is metabolite crosstalk between cell types (enzyme-sensor), analogous to CellPhoneDB for L-R | Synthase mRNA present but substrate/cofactor absent; transporter "sensor" is bidirectional/promiscuous so sender/receiver direction is wrong; no spatial geometry |
| scFEA | Per-cell metabolic FLUX through modules; graph neural network solver enforcing flux balance (in approximately out) | The question is relative flux per cell to FEED metabolite reasoning, not direct communication | Read as absolute mol/s (fluxes are relative); needs the matching module/stoichiometry files for the species; not a CCC tool itself |
| Compass | Per-cell metabolic STATE via flux-balance analysis; reaction penalty inversely proportional to enzyme expression over Recon2, outputs a score per reaction per cell | Comparing metabolic state between conditions (e.g. pathogenic vs non-pathogenic cells) | Used to claim a SECRETED metabolite communicates (output is reaction favorability, not secretion); assumes steady state, questionable for differentiated non-proliferating cells; heavy compute, micropool first |
| NeuronChat | Neurotransmitter/neuromodulator communication; vesicular release machinery and synthesis enzymes vs target receptor abundance | Neural systems specifically (glutamate, GABA, dopamine, serotonin, neuropeptides) | Applied outside neural tissue; same expression-proxy and geometry limits |

Methods and their curated databases evolve; before committing, verify current best practice, the metabolite-sensor database version, and required config/species files against the installed package docs.

## Confounds That Mimic Communication

| Confound | How it manufactures a fake signal | Mitigation |
|----------|-----------------------------------|------------|
| Compounded inference | Enzyme mRNA is a poor proxy for metabolite concentration (post-transcriptional control, substrate availability, allostery, compartmentalization), so a "secreted" metabolite may never be made | State the enzyme->flux->level->sensing chain; require metabolomics/MSI/tracing before any production claim |
| Bidirectional transporters | A transporter labeled a "sensor" may export rather than import, and many move several metabolites, so sender/receiver direction can be inverted | Treat transporter-based calls as lower-confidence than dedicated-receptor calls; check transport directionality literature |
| Ambient RNA | Soup of highly expressed transcripts inflates enzyme/sensor "expression" in clusters that do not transcribe them | Decontaminate (SoupX/DecontX/CellBender) before inference |
| Cell-type abundance and depth | Larger clusters tighten the permutation null and deeper cells detect more genes, inflating significance independent of biology | Down-sample, run on integrated counts, do not compare raw counts across conditions |
| No spatial geometry | Metabolites diffuse and degrade, but dissociated data has no coordinates, so a "communication" may be between cells never co-located | Validate proximity with spatial metabolomics/MSI; do not claim neighbor exchange from dissociated data |

## Run MEBOCOST

**Goal:** Score metabolite sender->receiver communication between cell types with permutation significance.

**Approach:** Build a MEBOCOST object from a log-normalized AnnData with cell-type labels and a config file pointing at the metabolite-sensor database, then run permutation inference; results carry both the communication score and an FDR.

```python
from mebocost import mebocost
import scanpy as sc

adata = sc.read_h5ad('adata_annotated.h5ad')   # log-normalized, gene SYMBOLS not Ensembl IDs

# config_path points to mebocost.conf listing the metabolite-enzyme-sensor database paths
# cutoff_prop=0.15: a gene must be expressed in >=15% of a group to count (dropout floor)
# species MUST match the data: mouse data against the human enzyme/sensor DB returns almost nothing
mebo = mebocost.create_obj(adata=adata, group_col='cell_type', condition_col=None,
                           met_est='mebocost', config_path='./mebocost.conf', species='human',
                           cutoff_exp='auto', cutoff_met='auto', cutoff_prop=0.15,
                           sensor_type='All', thread=8)

# n_shuffle=1000: label-permutation null for FDR; min_cell_number=10 drops tiny groups
commu_res = mebo.infer_commu(n_shuffle=1000, seed=12345, Return=True,
                             min_cell_number=10, pval_method='permutation_test_fdr',
                             pval_cutoff=0.05, thread=None)
```

## Filter and Summarize Results

**Goal:** Extract the significant, defensible metabolite communications.

**Approach:** Filter on the permutation FDR (not the raw p-value), then summarize by metabolite and by sender->receiver pair; column names are capitalized in the result table.

```python
sig = commu_res[commu_res['permutation_test_fdr'] < 0.05].copy()

# Result columns: Sender, Receiver, Metabolite_Name, Sensor, Annotation (Transporter/Enzyme),
# Commu_Score, Norm_Commu_Score, met_in_sender, sensor_in_receiver, permutation_test_fdr
sig['pair'] = sig['Sender'] + ' -> ' + sig['Receiver']
top_metabolites = sig['Metabolite_Name'].value_counts().head(10)
top_pairs = sig['pair'].value_counts().head(10)

# Transporter-based sensors are lower-confidence (bidirectional); separate them
transporter_calls = sig[sig['Annotation'] == 'Transporter']
```

## Compare Conditions

**Goal:** Find metabolite communications that differ between conditions (e.g. tumor vs normal).

**Approach:** Either pass `condition_col` to a single object, or run MEBOCOST separately per condition subset and compare the significant sets; never compare raw interaction counts across conditions without controlling for cell number and depth.

```python
results = {}
for cond in adata.obs['condition'].unique():
    sub = adata[adata.obs['condition'] == cond].copy()
    obj = mebocost.create_obj(adata=sub, group_col='cell_type', met_est='mebocost',
                              config_path='./mebocost.conf', species='human',
                              cutoff_prop=0.15, thread=8)
    results[cond] = obj.infer_commu(n_shuffle=1000, seed=12345, Return=True,
                                    min_cell_number=10, pval_cutoff=0.05)
# A "differential" metabolite call (significant in one condition only) is a HYPOTHESIS for metabolomics
```

## Threshold and Permutation Rationale

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `cutoff_prop` | 0.15 | A chosen dropout floor (MEBOCOST tutorials use 0.15-0.25, not a fixed package default): a gene expressed in <15% of a group is mostly dropout, but real low-abundance signaling is also discarded, so tune per dataset |
| `cutoff_exp` / `cutoff_met` | 'auto' | MEBOCOST data-derived thresholds for calling a gene/metabolite present; set manually only with a documented reason |
| `n_shuffle` | 1000 | Stable label-permutation FDR; the FDR is about label shuffling, not actual metabolite flux |
| `min_cell_number` | 10 | Groups under ~10 cells give unstable mean expression and inflated scores |
| `permutation_test_fdr` cutoff | 0.05 | Filter on the FDR, not the raw permutation p-value; significance is statistical, not a measured concentration |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `KeyError: 'metabolite'` / `'pval'` | Result columns are capitalized (Sender, Receiver, Metabolite_Name, Commu_Score, permutation_test_fdr) | Use the exact column names from `commu_res.columns` |
| `AttributeError: module 'mebocost' has no attribute 'create_obj'` | Wrong import; `create_obj` lives in the submodule | `from mebocost import mebocost` then `mebocost.create_obj(...)` |
| Almost no metabolites detected | Genes are Ensembl IDs, data is not log-normalized, `config_path` database is missing, or `species` does not match the data (mouse data run against the human enzyme/sensor DB) | Convert to gene symbols, log-normalize, point `config_path` at a valid mebocost.conf, set `species` to match the organism |
| A cell type "secretes" a metabolite implausibly | Synthase mRNA present but substrate/cofactor absent, or ambient RNA inflated the enzyme | Decontaminate ambient RNA; treat as "machinery consistent with", validate with metabolomics |
| Sender/receiver direction looks reversed | Sensor is a bidirectional/promiscuous transporter | Check `Annotation == 'Transporter'` calls separately; confirm transport direction |
| More communications in condition B than A | Counts scale with cell number and depth | Compare score magnitudes or matched subsets, not raw counts |

## Related Skills

- single-cell/cell-communication - Ligand-receptor CCC; the single-inference counterpart this skill mirrors at one extra remove
- single-cell/cell-annotation - Cell-type labels define metabolite senders and receivers
- single-cell/preprocessing - Log-normalization, gene-symbol mapping, and ambient-RNA decontamination happen here, before inference
- metabolomics/pathway-mapping - Places inferred metabolites in pathway context and informs which to prioritize
- metabolomics/isotope-tracing - Orthogonal flux validation that a producing cell actually makes the metabolite
- systems-biology/flux-balance-analysis - Genome-scale FBA underlying Compass-style per-cell metabolic state

## References

- Zheng R, et al. MEBOCOST maps metabolite-mediated intercellular communications using single-cell RNA-seq. Nucleic Acids Res 53(12):gkaf569 (2025). PMID 40568942.
- Alghamdi N, et al. A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data [scFEA]. Genome Res 31(10):1867-1884 (2021).
- Wagner A, et al. Metabolic modeling of single Th17 cells reveals regulators of autoimmunity [Compass]. Cell 184(16):4168-4185 (2021).
- Zhao W, et al. Inferring neuron-neuron communications from single-cell transcriptomics through NeuronChat. Nat Commun 14(1):1128 (2023).
- Dimitrov D, et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13:3224 (2022). [discordance framing]
- Young MD, Behjati S. SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience 9(12):giaa151 (2020).
