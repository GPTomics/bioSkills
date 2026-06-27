# Cell-Cell Communication - Usage Guide

## Overview

Cell-cell communication (CCC) analysis ranks ligand-receptor interactions between cell types from scRNA-seq. Every output is a co-expression proxy, not proof of signaling, and competing methods disagree because they estimate different quantities (specificity vs magnitude vs probability). The defensible workflow is consensus-first (LIANA), with a resource-sensitivity check and orthogonal validation, reserving NicheNet for the distinct mechanistic question of which ligand drives a receiver's response.

## Prerequisites

```bash
# Consensus + specificity (Python)
pip install liana cellphonedb scanpy anndata
```

```r
# Pathway probability + downstream activity (R)
install.packages('Seurat')                         # CellChat and nichenetr build on Seurat objects
devtools::install_github('jinworks/CellChat')      # repo moved from sqjin
devtools::install_github('saeyslab/nichenetr')     # plus NicheNet model files from Zenodo
```

## Quick Start

Tell your AI agent what you want to do:
- "Rank ligand-receptor interactions between my cell types with a consensus method"
- "Check whether my top interactions survive a different L-R database"
- "Which ligand from macrophages best explains the DE genes in T cells?"
- "Compare cell communication between control and disease"
- "Summarize WNT signaling at the pathway level with sender and receiver roles"

## Example Prompts

### Consensus Inference
> "Run LIANA rank_aggregate and report both magnitude and specificity ranks"
> "Find robust ligand-receptor pairs between fibroblasts and epithelial cells"

### Resource Sensitivity
> "Re-run the same method with the CellPhoneDB and CellChat resources and show which pairs are stable"
> "Is my VEGF-VEGFR hit robust to the choice of database?"

### Specificity Testing
> "Run CellPhoneDB with permutation p-values and proper complex handling"
> "Use the DEG-based CellPhoneDB method instead of one-vs-rest"

### Pathway-Level and Roles
> "Summarize signaling pathways between these cell types with CellChat"
> "Which cell types are the dominant senders and receivers of TGF-beta?"

### Downstream Mechanism (NicheNet)
> "Which ligands explain the activated-T-cell gene signature?"
> "Rank ligands by downstream target-gene activity in the receiver"

### Condition Comparison
> "Compare communication between control and treatment and show gained/lost interactions"

## What the Agent Will Do

1. Confirm cell-type annotations and check that ambient RNA and dissociation-stress genes were handled in preprocessing
2. Run LIANA `rank_aggregate` as the consensus default and report magnitude AND specificity ranks
3. Run a resource-sensitivity check by holding the method fixed and swapping the L-R database
4. Add CellPhoneDB for permutation specificity p-values, or CellChat for pathway-level summaries and sender/receiver roles, as the question requires
5. Use NicheNet only when the question is mechanism, building a clean receiver DE gene set first
6. For spatial data, switch to a proximity-aware method (Squidpy, COMMOT, CellChat v2 spatial, LIANA+ bivariate) and run a sensitivity analysis over the diffusion radius
7. Frame every interaction as a hypothesis and recommend orthogonal validation

## Tips

- **Consensus over single method** - LIANA `rank_aggregate` hedges against the structural discordance between methods; a single tool's ranking is one estimand, not the truth.
- **Resource can dominate** - The L-R database moves results as much as the method (Dimitrov 2022); report the resource and show the key claim survives >=2 resources.
- **Magnitude and specificity are different axes** - A pair can top one rank and bottom the other; report both rather than collapsing to one number.
- **Decontaminate before scoring** - Ambient RNA inflates the secreted-ligand half of pairs and manufactures "universal senders"; run SoupX/DecontX/CellBender first.
- **Abundance and depth are confounds** - Large clusters and deep sequencing produce more "significant" interactions independent of biology; never compare raw interaction counts across conditions without normalization.
- **Watch dissociation-stress genes** - FOS, JUN, JUNB, EGR1, HSPA, DUSP1 are bona fide ligands that fabricate AP-1 / heat-shock signaling; flag or regress them.
- **NicheNet needs a clean gene set** - Its ligand ranking is only as good as the receiver DE list, and its prior network is static and cell-type-agnostic; read top ligands as "consistent with the response", not proof.
- **Membrane-bound ligands need contact** - Notch-DLL/JAG and ephrins scored between non-adjacent types in dissociated data are almost certainly artifacts; restrict to secreted signaling or use spatial data.
- **Mouse data** - CellPhoneDB is human-only; prefer CellChatDB.mouse or LIANA `mouseconsensus` over ortholog mapping.
- **Validate orthogonally** - Confirm hits with downstream TF/pathway activity (decoupleR/PROGENy), receptor protein (CITE-seq), spatial co-localization, or perturbation.

## Related Skills

single-cell/cell-annotation - Cell-type labels define senders and receivers; annotation resolution is a hidden CCC hyperparameter
single-cell/clustering - Cluster granularity changes who is "specific"; fix it before running CCC
single-cell/doublet-detection - Doublets create fake co-expressing cells that masquerade as senders-receivers
single-cell/preprocessing - Ambient-RNA decontamination and stress-gene handling happen here, before CCC
single-cell/metabolite-communication - Metabolite-mediated CCC (enzyme-sensor) as the doubly-inferred counterpart to ligand-receptor
spatial-transcriptomics/spatial-communication - Proximity-constrained CCC when spatial coordinates are available
pathway-analysis/go-enrichment - Functional enrichment of NicheNet target genes or interacting receptors
differential-expression/deseq2-basics - Pseudobulk DE to build the receiver gene set NicheNet requires
