---
name: bio-immunoinformatics-tcr-epitope-binding
description: Infer or annotate TCR antigen specificity by unsupervised clustering (TCRdist/tcrdist3, GLIPH2, clusTCR, GIANA) and database lookup (VDJdb, IEDB, McPAS-TCR), and rank candidates with supervised predictors (ERGO-II, NetTCR-2.x, pMTnet) under explicit caveats. Encodes the central truth that general TCR-epitope prediction for UNSEEN epitopes essentially does not work (collapses to near-random; IMMREP22, Grazioli 2022) because labeled data is dominated by a few immunodominant epitopes and there is no true negative set — so clustering for discovery is the honest task and de-novo binding needs wet-lab validation. Use when annotating TCR specificity or grouping a repertoire. Epitope/MHC context lives in mhc-binding-prediction.
tool_type: python
primary_tool: tcrdist3
---

## Version Compatibility

Reference examples tested with: tcrdist3 0.2+, pandas 2.2+, scipy 1.12+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: tcrdist3 expects IMGT-style columns (`cdr3_b_aa`, `v_b_gene`, `j_b_gene`, and the `_a_` analogs, plus `count`). Tools disagree on whether CDR3 keeps the leading Cys / trailing Phe-Trp — a common silent input mismatch. Supervised predictors (ERGO-II, NetTCR, pMTnet) are separate repos with pretrained weights; their reported AUCs depend heavily on the train/test split and negative-sampling scheme. Re-verify before trusting any number.

# TCR-Epitope Binding

**"What antigen does this TCR recognize / which TCRs share specificity?"** -> Annotate specificity by clustering + database lookup; predict de-novo binding only as a validation-bound hypothesis.
- Python: `tcrdist3` (TCRrep distance + meta-clonotypes), GLIPH2, clusTCR, GIANA for clustering
- Python: ERGO-II / NetTCR-2.x / pMTnet for supervised scoring (caveated); VDJdb/IEDB/McPAS-TCR for lookup

## The Single Most Important Modern Insight -- general prediction for unseen epitopes does not work; clustering does

Every supervised TCR-epitope predictor performs respectably on epitopes seen in training and collapses to near-random on epitopes it has never seen (Grazioli 2022; IMMREP22, Meysman 2023 across 23 models). The cause is the data, not the architecture: the labeled TCR-pMHC universe is dominated by a few immunodominant epitopes (NLVPMVATV/CMV, GILGFVFTL/influenza M1, SARS-CoV-2 spike), so a model learns "is this an anti-CMV TCR" rather than the rules of TCR-peptide docking. Compounding this, there is no true negative set — experiments report binders, and absence of a measured non-binder is not non-binding — so every supervised model manufactures negatives, and that choice dominates the reported metric more than the architecture (Dens 2023). The honest, defensible task is unsupervised specificity clustering: "these TCRs are sequence-similar enough to likely share a specificity," a discovery statement used within one dataset and propagated by guilt-by-association to a known member. Clustering is honest because it never extrapolates into unseen-epitope space; per-pair prediction is dishonest when it pretends to. Route the user to the honest task and refuse to let a supervised per-pair probability substitute for a tetramer.

## Tool Taxonomy

| Tool | Citation | Task | Input | Note |
|------|----------|------|-------|------|
| TCRdist / tcrdist3 | Dash 2017; Mayer-Blackwell 2021 | Clustering (distance) | CDR3 + V/J, both chains | Multi-loop distance, 3x weight on CDR3; meta-clonotypes |
| GLIPH2 | Huang 2020 | Clustering (global + motif) | CDR3β + V/J + HLA | Predicts restricting allele; background-repertoire dependent |
| clusTCR | Valkiers 2021 | Clustering (Faiss+MCL) | CDR3β | Scales to millions; speed for specificity |
| GIANA / iSMART | Zhang 2021; Zhang 2020 | Clustering (fast) | CDR3β | Small high-specificity clusters |
| ERGO-II | Springer 2021 | Supervised prediction | CDR3β(+α,V,J,MHC) | Degrades gracefully; seen-epitope only |
| NetTCR-2.x | Montemurro 2021 | Supervised prediction | paired CDR3α+β | Paired beats single-chain; ~150 pos/epitope needed |
| pMTnet / PanPep | Lu 2021; Gao 2023 | Supervised, neoantigen-aimed | CDR3β + peptide + MHC | Zero-shot claims need skepticism |

## Reference Databases (training set AND lookup table)

| Database | Citation | Content | Caveat |
|----------|----------|---------|--------|
| VDJdb | Shugay 2018; Bagaev 2020 | Curated TCR-pMHC with confidence 0-3 | Filter on confidence; skewed to HLA-A*02:01 |
| IEDB | Vita 2019 | TCR + pMHC assays | The corpus most predictors draw on |
| McPAS-TCR | Tickotsky 2017 | Pathology-organized (infection/cancer/autoimmune) | Human + mouse |
| 10x dextramer | Zhang 2021 (Sci Adv) | Largest paired-chain set, 4 donors | Labels are threshold calls, not gold; multiplets/background |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Have known specificities (tetramer sort / DB hits) | Cluster (tcrdist3/GLIPH2) + lookup, propagate labels | The honest, bounded question |
| Group a repertoire by likely shared specificity | tcrdist3 or clusTCR within one cohort | Discovery within dataset; keep HLA as covariate |
| Truly de-novo novel epitope (e.g. neoantigen) | Rank with pMTnet/PanPep, label as hypothesis, validate | Prediction does not generalize; tetramer/functional assay decides |
| Millions of CDR3s | clusTCR (Faiss+MCL) | Speed at modest specificity cost |
| Predict restricting HLA from sequence | GLIPH2 | Infers allele from cross-donor co-occurrence |
| "Does this TCR bind peptide X?" for unseen X | No reliable computational answer | State plainly; there is no third branch |

## Cluster TCRs by Specificity (tcrdist3)

**Goal:** Group TCRs likely to share an antigen, for discovery within one cohort.

**Approach:** Build a TCRrep (which computes the position-weighted multi-loop distance, 3x on CDR3), then cluster the pairwise matrix and annotate clusters containing a known-specificity member. Keep HLA as an explicit covariate — the same CDR3 on a different allele is a different specificity.

```python
from tcrdist.repertoire import TCRrep
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

def cluster_tcrs(df, max_dist=50):
    '''df needs IMGT columns: cdr3_b_aa, v_b_gene, j_b_gene (+ _a_ analogs), count.
    Returns cluster labels; annotate clusters that contain a database/tetramer hit.'''
    tr = TCRrep(cell_df=df, organism='human', chains=['beta'])
    condensed = squareform(tr.pw_beta, checks=False)
    return fcluster(linkage(condensed, method='average'), t=max_dist, criterion='distance')
```

## Annotate by Database Lookup

**Goal:** Assign specificity to TCRs that match known TCR-pMHC pairs.

**Approach:** Match exactly or near-exactly against VDJdb/IEDB/McPAS, filtering VDJdb on its confidence score, and report the database hit and HLA restriction driving each annotation — not a per-pair probability dressed as certainty.

```python
import pandas as pd

def lookup_vdjdb(query_cdr3b, vdjdb, min_confidence=1):
    '''Exact CDR3b match against a confidence-filtered VDJdb. Near-matches (edit
    distance 1) belong to the clustering route, not a binding claim.'''
    db = vdjdb[vdjdb['vdjdb.score'] >= min_confidence]
    hits = db[db['cdr3'].isin(set(query_cdr3b))]
    return hits[['cdr3', 'antigen.epitope', 'antigen.species', 'mhc.a']]
```

## Per-Method Failure Modes

### Unseen-epitope collapse
**Trigger:** using a supervised model on an epitope absent from training. **Mechanism:** models learn a few well-sampled specificities, not docking rules. **Symptom:** great benchmark AUC, near-random on novel epitopes. **Fix:** route de-novo questions to ranking-plus-validation; never report a confident per-pair call.

### Negative-sampling artifact
**Trigger:** trusting a headline AUC. **Mechanism:** manufactured negatives (shuffled or random-TCR) create artificial separability; repeated-negative leakage lets the model count TCR frequency. **Symptom:** AUC > 0.85 with no discussion of negatives/splits. **Fix:** read the negative-sampling sentence first; require epitope-disjoint evaluation.

### CDR3β-only ceiling
**Trigger:** strong claims from a β-only model. **Mechanism:** alpha chain and V/J carry heavy signal; bulk β-only is information-poor. **Symptom:** big AUC from the least informative input (i.e. from artifacts). **Fix:** prefer paired-chain data; add V/J; discount β-only headline numbers.

### Clustering confounds (HLA + background)
**Trigger:** pooling multi-donor repertoires and clustering naively. **Mechanism:** same CDR3 on different HLA is a different specificity; motif enrichment depends on the reference background. **Symptom:** merged unrelated TCRs; spurious "enriched" motifs. **Fix:** cluster within a coherent cohort, keep HLA as a covariate, match the background repertoire.

### Metrics that lie
**Trigger:** a single global or per-epitope-averaged AUC. **Mechanism:** averaging over seen epitopes hides the novel-epitope collapse. **Symptom:** one trustworthy-looking number, no split description. **Fix:** demand epitope-disjoint (TPP-II/III) splits reported per-epitope with a peptide-distance decay analysis.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| ~150 distinct binders per epitope | Montemurro 2021 | Below this a per-epitope supervised model is unreliable |
| VDJdb confidence >= 1 (use 2-3 for high) | Shugay 2018 | Low-confidence records are weakly supported |
| 10x call: UMI > 10 and > 5x top negative-control | Zhang 2021 | Dextramer labels are threshold calls, not gold |
| Evaluate on epitope-disjoint split | IMMREP22; Grazioli 2022 | Seen-epitope/shuffled splits hide the collapse |
| tcrdist CDR3 weight 3x other loops | Dash 2017 | CDR3 is the chief specificity determinant |
| Paired α+β > single chain | Montemurro 2021; IMMREP22 | β-only caps achievable accuracy |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Confident de-novo binding call | Supervised model on unseen epitope | Reframe as hypothesis; validate by tetramer/assay |
| Irreproducible published AUC | Leaky split / negative-sampling bias | Re-evaluate on epitope-disjoint clean split |
| Merged unrelated TCR clusters | Mixed-HLA pooling | Cluster within cohort; HLA covariate |
| Input not recognized by tool | CDR3 Cys/Phe-Trp convention mismatch | Match the tool's IMGT trimming convention |
| Over-trusted 10x labels | Treated threshold calls as gold | Require replicate/donor concordance |
| Structure "solves" it | AlphaFold hype on a hard interface | Use TCRdock to rank/rationalize candidates, not screen |

## References

- Dash P, Fiore-Gartland AJ, Hertz T, et al. 2017. Quantifiable predictive features define epitope-specific T cell receptor repertoires (TCRdist). *Nature* 547(7661):89-93.
- Mayer-Blackwell K, Schattgen S, Cohen-Lavi L, et al. 2021. TCR meta-clonotypes for biomarker discovery with tcrdist3. *eLife* 10:e68605.
- Glanville J, Huang H, Nau A, et al. 2017. Identifying specificity groups in the T cell receptor repertoire (GLIPH). *Nature* 547(7661):94-98.
- Huang H, Wang C, Rubelt F, Scriba TJ, Davis MM. 2020. Analyzing the M. tuberculosis immune response by T-cell receptor clustering with GLIPH2. *Nature Biotechnology* 38:1194-1202.
- Valkiers S, Van Houcke M, Laukens K, Meysman P. 2021. clusTCR: a Python interface for rapid clustering of large sets of CDR3 sequences. *Bioinformatics* 37(24):4865-4867.
- Zhang H, Zhan X, Li B. 2021. GIANA allows computationally-efficient TCR clustering and multi-disease repertoire classification by isometric transformation. *Nature Communications* 12:4699.
- Springer I, Tickotsky N, Louzoun Y. 2021. Contribution of T cell receptor alpha and beta CDR3, MHC typing, V and J genes to peptide binding prediction (ERGO-II). *Frontiers in Immunology* 12:664514.
- Montemurro A, Schuster V, Povlsen HR, et al. 2021. NetTCR-2.0 enables accurate prediction of TCR-peptide binding using paired TCRα and β sequence data. *Communications Biology* 4:1060.
- Lu T, Zhang Z, Zhu J, et al. 2021. Deep learning-based prediction of the T cell receptor-antigen binding specificity (pMTnet). *Nature Machine Intelligence* 3:864-875.
- Meysman P, Barton J, Bravi B, et al. 2023. Benchmarking solutions to the T-cell receptor epitope prediction problem: IMMREP22 workshop report. *ImmunoInformatics* 9:100024.
- Grazioli F, Mösch A, Machart P, et al. 2022. On TCR binding predictors failing to generalize to unseen peptides. *Frontiers in Immunology* 13:1014256.
- Dens C, Bittremieux W, Affaticati F, Laukens K, Meysman P. 2023. The pitfalls of negative data bias for the T-cell epitope specificity challenge. *Nature Machine Intelligence* 5:1063-1065.
- Shugay M, Bagaev DV, Zvyagin IV, et al. 2018. VDJdb: a curated database of T-cell receptor sequences with known antigen specificity. *Nucleic Acids Research* 46(D1):D419-D427.
- Bradley P. 2023. Structure-based prediction of T cell receptor:peptide-MHC interactions (TCRdock). *eLife* 12:e82813.

## Related Skills

- immunoinformatics/mhc-binding-prediction - the pMHC context a TCR recognizes (HLA restriction)
- immunoinformatics/neoantigen-prediction - de-novo neoantigen TCRs are the unseen-epitope case where prediction fails
- tcr-bcr-analysis/mixcr-analysis - upstream TCR repertoire extraction from sequencing
- single-cell/clustering - paired single-cell TCR data and embedding-based grouping
- workflows/tcr-pipeline - end-to-end repertoire processing
