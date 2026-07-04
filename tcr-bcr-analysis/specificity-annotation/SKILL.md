---
name: bio-tcr-bcr-analysis-specificity-annotation
description: Maps TCR/BCR receptor sequences toward candidate antigen specificity and clusters repertoires by shared-specificity signal, while enforcing that a database match or a cluster label is a HYPOTHESIS, not a specificity call. Use when deciding among database annotation (VDJdb/McPAS/IEDB+TCRMatch, requiring V-gene and HLA concordance plus a confidence score) versus sequence clustering (tcrdist3 meta-clonotypes, GLIPH2, GIANA, clusTCR, which find enrichment not per-receptor labels) versus generation-probability nulls (OLGA Pgen, IGoR, SONIA Ppost) for testing public/convergent/shared claims; and when guarding against overclaiming specificity, base-rate false positives from bare CDR3 matches, unpaired beta-only annotation, ML predictor failure on unseen epitopes, and ignored MHC restriction. TCR-focused with a BCR/antibody note (SHM, conformational epitopes, IGHV3-53/3-66 public clonotypes). Keywords CDR3, pMHC, HLA restriction, cross-reactivity, meta-clonotype, Pgen, public clonotype, convergent recombination.
tool_type: mixed
primary_tool: tcrdist3
---

## Version Compatibility

Reference examples tested with: tcrdist3 0.2+, olga 1.2+, pandas 2.0+, numpy 1.24+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: OLGA's CLI entry point is `olga-compute_pgen` (underscore) and its bundled human beta model lives in `default_models/human_T_beta/`. GLIPH2 results are reference-repertoire- and parameter-sensitive; verify the reference set used and rerun with a second method before trusting cluster labels.

# Antigen-Specificity Annotation and Clustering for TCR/BCR Repertoires

**"Which antigens might my receptors recognize?"** -> Annotate CDR3s against curated TCR:pMHC records, keeping only defensible matches.
- Python: `pandas.merge` on VDJdb/McPAS export; IEDB `TCRMatch` for k-mer similarity to characterized receptors

**"Cluster my repertoire into shared-specificity groups."** -> Find sequence neighborhoods enriched for a common specificity.
- Python: `tcrdist.repertoire.TCRrep` (TCRdist neighborhoods / meta-clonotypes), GLIPH2, GIANA, clusTCR

**"Is this shared/public clonotype antigen-driven or just easy to generate?"** -> Test the sharing against a generation-probability null.
- Python/CLI: `olga` (Pgen), IGoR (learn the model), SONIA/soNNia (Ppost = Pgen x Q)

## The governing principle: a match or a cluster is a hypothesis, not a specificity call

Mapping receptor sequence to antigen specificity is largely UNSOLVED in the general case. Every method here does one of three things and none of them proves specificity: (a) it annotates against a biased database, (b) it finds a sequence neighborhood ENRICHED for shared specificity, or (c) it predicts binding only for epitopes seen in training. A VDJdb hit or a GLIPH2 cluster label is therefore a hypothesis with a confidence level, never "this receptor is specific for antigen Y". The entire job of this skill is to stop that overclaim: report annotations and clusters with their confidence, their required concordances (V-gene, HLA), and a generation-probability baseline, and reserve the word "specific" for tetramer/dextramer or functional confirmation.

Why the problem resists a clean solution:
- A TCR recognizes a peptide-MHC complex, not a peptide. Specificity is a property of the TCR:pMHC triple, so the same CDR3 can be specific for different peptides under different HLA. MHC restriction cannot be dropped.
- Massive cross-reactivity: a single TCR can recognize up to ~10^6 peptides (Sewell 2012 *Nat Rev Immunol* 12:669). One-receptor-one-antigen is false by design.
- Specificity is encoded jointly by the paired alpha and beta chains; bulk sequencing gives beta-only and discards the pairing that carries much of the signal.
- Training and database records are dominated by a few immunodominant epitopes (influenza GILGFVFTL, CMV NLVPMVATV, EBV GLCTLVAML, SARS-CoV-2), with heavy HLA-A*02:01 and CD8/MHC-I skew; CD4/MHC-II, gamma-delta TCR, and BCR data are sparse to absent. Any accuracy averaged over epitopes is inflated by the few easy ones, and a gamma-delta or CD4 repertoire will return almost nothing from these databases and models (absence of a hit is uninformative, not evidence of non-specificity).

## Three approaches: pick by the question, then corroborate

| Approach | Tools | What it answers | What it PROVES / does NOT | Best when | Fails when |
|----------|-------|-----------------|---------------------------|-----------|------------|
| Database annotation | VDJdb (score 0-3), McPAS-TCR, IEDB + TCRMatch | Does a curated TCR:pMHC record match this receptor? | A record matches; NOT that the receptor is specific (base-rate false positives) | Donor HLA known, V-gene available, high-confidence entries wanted | Bare CDR3 match, no HLA/V concordance, high-Pgen sequences match by chance |
| Sequence clustering | tcrdist3 meta-clonotypes, GLIPH2, GIANA, clusTCR, iSMART | Which receptors form a shared-specificity neighborhood? | A group is enriched for shared specificity; NOT a per-receptor antigen label | Discovering specificity groups, building reusable features from many receptors | Treating a cluster label as an antigen call; single-tool trust; no reference-null |
| Generation-probability null | OLGA (Pgen), IGoR, SONIA/soNNia (Ppost) | Is this sharing/convergence more than chance generation? | Whether a sequence is expected by recombination; the null for every sharing claim | Any "public"/convergent/shared/expanded-beyond-chance claim | Omitted entirely (the most common gap) -> publicity mistaken for antigen selection |

Default workflow: annotate with confidence and concordance, treat clusters as hypotheses, and attach a Pgen null to any sharing or convergence claim. Run at least two clustering methods and report agreement; benchmarks disagree on tool ranking by dataset and epitope, and there is no accepted gold standard (Meysman 2023 *ImmunoInformatics* 9:100024). Verify current best practice against each tool's latest docs before committing to one.

## Database annotation with the base-rate guardrail

**Goal:** Annotate a bulk beta repertoire against curated TCR:pMHC records without generating a flood of chance matches.

**Approach:** Restrict the database to confidence >= 1, join on CDR3 AND V-gene (never CDR3 alone), then drop matches whose restricting HLA the donor does not carry. A bare CDR3-beta match to a high-Pgen sequence is a base-rate false positive: the number of spurious hits scales with repertoire size x database size x match permissiveness.

```python
import pandas as pd

def annotate_by_db(repertoire, vdjdb, donor_hla, min_confidence=1):
    # repertoire, vdjdb: cdr3_b_aa + v_b_gene (IMGT, e.g. TRBV19*01); vdjdb also antigen_epitope, mhc_a, vdjdb_score
    # vdjdb_score 0-3: 0 = critical info missing, 3 = independently validated; >=1 drops single-observation noise
    db = vdjdb[vdjdb['vdjdb_score'] >= min_confidence]
    hits = repertoire.merge(db, on=['cdr3_b_aa', 'v_b_gene'], how='inner', suffixes=('', '_db'))  # V concordance, not CDR3 alone
    carries_hla = hits['mhc_a'].apply(lambda a: any(a.startswith(h) for h in donor_hla))  # restricting HLA must be present in donor
    hits = hits[carries_hla].copy()
    hits['annotation_confidence'] = 'hypothesis'  # a curated match, not a specificity call
    return hits
```

IEDB ships `TCRMatch` (Chronister 2021 *Front Immunol* 12:640725) for k-mer similarity of a CDR3-beta to characterized receptors: it returns a similarity score, which proves sequence similarity to a known receptor, not binding. Report match counts and an enrichment statistic against a size-matched synthetic/unexposed repertoire, not a binary "specific".

## Sequence clustering as neighborhood discovery

**Goal:** Group receptors into shared-specificity neighborhoods that can be quantified and reused, without mistaking a cluster for an antigen label.

**Approach:** TCRdist scores position-weighted CDR distances using the germline-encoded CDR1/CDR2/CDR2.5 loops (V-gene identity is baked in) plus the CDR3 up-weighted ~3x (Dash 2017 *Nature* 547:89). Build the pairwise beta matrix, then take fixed-radius neighborhoods in TCRdist units. A neighborhood is a meta-clonotype candidate (centroid + radius + optional motif), a testable feature, not a specificity assignment (Mayer-Blackwell 2021 *eLife* 10:e68605).

```python
import numpy as np
from tcrdist.repertoire import TCRrep

def beta_neighborhoods(clone_df, radius=50):
    # clone_df: cdr3_b_aa, v_b_gene, j_b_gene, count (IMGT gene names). organism/chains fix the germline loops used.
    tr = TCRrep(cell_df=clone_df, organism='human', chains=['beta'], db_file='alphabeta_gammadelta_db.tsv')
    # radius in TCRdist units; ~50 is a common meta-clonotype inclusion radius (Mayer-Blackwell 2021), tune per centroid
    neighbors = [set(np.where(row <= radius)[0]) for row in tr.pw_beta]
    return tr, neighbors
```

CDR3-only tools (GLIPH2 Huang 2020 *Nat Biotechnol* 38:1194; GIANA Zhang 2021 *Nat Commun* 12:4699; clusTCR Valkiers 2021 *Bioinformatics* 37:4865; iSMART Zhang 2020 *Clin Cancer Res* 26:1359) scale to millions of CDR3s but ignore the paired chain and often HLA, so they can merge receptors sharing a motif but differing in true restriction. GLIPH2 tends to produce large, low-specificity clusters and its enrichment is sensitive to the reference repertoire and input size; report the reference used and a null-corrected p-value, and prefer tcrdist3 meta-clonotypes as reusable features. For paired single-cell data, CoNGA (Schattgen 2022 *Nat Biotechnol* 40:54) links TCR neighborhoods to gene-expression state rather than to an antigen.

## Generation-probability nulls: the rigorous test for "public"/convergent

**Goal:** Decide whether a shared, convergent, or database-matched clonotype is antigen-driven or merely easy to generate.

**Approach:** High-Pgen sequences (few insertions, germline-like, short CDR3) recur across donors and match databases BY CHANCE, so publicity alone is not evidence of selection (Quigley 2010 *PNAS* 107:19414). Compute Pgen with OLGA (Sethna 2019 *Bioinformatics* 35:2974) for every match or shared clonotype and down-weight the high-Pgen ones; a clonotype is evidence of convergent selection only if observed more than its generation probability predicts. SONIA/soNNia add selection (Ppost = Pgen x Q) for a post-selection null.

```python
import os
import olga.load_model as load_model
import olga.generation_probability as generation_probability

def load_beta_pgen_model(model_dir):
    # model_dir = OLGA default_models/human_T_beta with model_params.txt, model_marginals.txt, V/J anchors CSV
    gen_data = load_model.GenomicDataVDJ()
    gen_data.load_igor_genomic_data(os.path.join(model_dir, 'model_params.txt'),
                                    os.path.join(model_dir, 'V_gene_CDR3_anchors.csv'),
                                    os.path.join(model_dir, 'J_gene_CDR3_anchors.csv'))
    gen_model = load_model.GenerativeModelVDJ()
    gen_model.load_and_process_igor_model(os.path.join(model_dir, 'model_marginals.txt'))
    return generation_probability.GenerationProbabilityVDJ(gen_model, gen_data)

def pgen(model, cdr3_b_aa, v_b_gene, j_b_gene):
    return model.compute_aa_CDR3_pgen(cdr3_b_aa, v_b_gene, j_b_gene)  # ~ms/seq; high value => expected by chance, down-weight
```

The CLI equivalent is `olga-compute_pgen --humanTRB CASSLGQAYEQYF` or `olga-compute_pgen --humanTRB -i seqs.tsv -o pgens.tsv`.

## ML binding predictors: interpolate within trained epitopes only

DeepTCR, ERGO-II, NetTCR-2.0, pMTnet, and TITAN predict TCR:epitope binding but interpolate WITHIN epitopes seen in training and collapse toward random on unseen epitopes (Moris 2021 *Brief Bioinform* 22:bbaa318; Grazioli 2022 *Front Immunol* 13:1014256). Published AUCs are inflated by data leakage (same epitope in train and test) and by negative-sampling artifacts, where models separate the negative-generation process rather than binding (Dens 2023 *Nat Mach Intell* 5:1060-1062). Gate any predictor to epitopes well-represented in its training set, never present a per-TCR score for a novel neoantigen as validated, and require leave-epitope-out evaluation with reference negatives. Cross-reference immunoinformatics/tcr-epitope-binding for the predictor details rather than duplicating them here.

## BCR/antibody note

Antibody specificity is harder than TCR: somatic hypermutation makes the functional sequence a moving target away from germline, and most antibody epitopes are conformational, so sequence-only epitope prediction is fundamentally limited and often needs structure. The sequence-level analogue of a public TCR is a convergent public antibody clonotype to a pathogen, e.g. the IGHV3-53/IGHV3-66 clonotype against the SARS-CoV-2 RBD (Robbiani 2020 *Nature* 584:437; Tan 2021 *Nat Commun* 12:4210). As with TCRs, a shared V-gene plus CDRH3 motif across donors is a hypothesis of convergent selection that still needs a Pgen/expected-sharing null and binding confirmation. Cluster BCR clones (shared V, J, junction length + within-partition distance) before any such analysis; see immcantation-analysis.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| A repertoire reported "specific for antigen Y" from a DB hit or cluster | A curated match or a cluster label treated as ground-truth specificity | Report it as an annotation/hypothesis with confidence; confirm with tetramer/dextramer or functional assay before saying "specific" |
| Many DB matches to many epitopes in a large repertoire | Base-rate false positives from bare CDR3-beta matching to high-Pgen sequences | Require V-gene concordance, require donor HLA carriage, filter to VDJdb score >= 1-2, and attach a Pgen null |
| "Public"/convergent response claimed from sharing across donors | Sharing tested without a generation-probability baseline | Compute OLGA Pgen (or SONIA Ppost); a clonotype is convergent only if observed above its generation expectation |
| Beta-only annotation trusted as full specificity | Bulk data is unpaired; specificity is a paired alpha/beta + HLA property | Down-weight beta-only annotations; use single-cell paired alpha/beta (scirpy) or paired TCRdist where possible |
| ML predictor score believed for a new neoantigen | Predictors fail out-of-distribution on unseen epitopes; benchmark AUCs leak | Gate to well-trained epitopes; use leave-epitope-out validation and reference negatives; treat novel-epitope scores as unreliable |
| Cluster from one tool taken as truth | No gold standard; GLIPH2 clusters are reference- and parameter-sensitive | Run >= 2 methods, report agreement and the reference repertoire, check HLA-concordance and cluster tightness |
| HLA ignored in annotation | Specificity is a TCR:pMHC triple; the same CDR3 differs by restricting HLA | Restrict DB entries to alleles the donor carries; report the restricting HLA with every annotation |

## Related Skills

- mixcr-analysis - Produce clonotype tables to annotate
- scirpy-analysis - Paired single-cell clonotypes for specificity work
- vdjtools-analysis - Public-clonotype context and overlap
- immunoinformatics/tcr-epitope-binding - ML epitope-binding prediction details
- immunoinformatics/mhc-binding-prediction - Upstream pMHC restriction
- immunoinformatics/neoantigen-prediction - Neoantigen-directed specificity

## References

- Dash P, et al. Quantifiable predictive features define epitope-specific T cell receptor repertoires. *Nature* 2017; 547(7661):89-93.
- Glanville J, et al. Identifying specificity groups in the T cell receptor repertoire. *Nature* 2017; 547(7661):94-98.
- Huang H, et al. Analyzing the Mycobacterium tuberculosis immune response by T-cell receptor clustering with GLIPH2 and genome-wide antigen screening. *Nat Biotechnol* 2020; 38:1194-1202.
- Mayer-Blackwell K, et al. TCR meta-clonotypes for biomarker discovery with tcrdist3. *eLife* 2021; 10:e68605.
- Schattgen SA, et al. Integrating T cell receptor sequences and transcriptional profiles by clonotype neighbor graph analysis (CoNGA). *Nat Biotechnol* 2022; 40:54-63.
- Zhang H, et al. Investigation of Antigen-Specific T-Cell Receptor Clusters in Human Cancers (iSMART). *Clin Cancer Res* 2020; 26(6):1359-1371.
- Zhang H, et al. GIANA allows computationally-efficient TCR clustering and multi-disease repertoire classification by isometric transformation. *Nat Commun* 2021; 12:4699.
- Valkiers S, et al. clusTCR: a Python interface for rapid clustering of large sets of CDR3 sequences with unknown antigen specificity. *Bioinformatics* 2021; 37(24):4865-4867.
- Shugay M, et al. VDJdb: a curated database of T-cell receptor sequences with known antigen specificity. *Nucleic Acids Res* 2018; 46(D1):D419-D427.
- Tickotsky N, et al. McPAS-TCR: a manually curated catalogue of pathology-associated T cell receptor sequences. *Bioinformatics* 2017; 33(18):2924-2929.
- Chronister WD, et al. TCRMatch: predicting T-cell receptor specificity based on sequence similarity to previously characterized receptors. *Front Immunol* 2021; 12:640725.
- Marcou Q, Mora T, Walczak AM. High-throughput immune repertoire analysis with IGoR. *Nat Commun* 2018; 9:561.
- Sethna Z, et al. OLGA: fast computation of generation probabilities of B- and T-cell receptor amino acid sequences and motifs. *Bioinformatics* 2019; 35(17):2974-2981.
- Quigley MF, et al. Convergent recombination shapes the clonotypic landscape of the naive T-cell repertoire. *PNAS* 2010; 107:19414-19419.
- Moris P, et al. Current challenges for unseen-epitope TCR interaction prediction and a new perspective derived from image classification. *Brief Bioinform* 2021; 22(4):bbaa318.
- Meysman P, et al. Benchmarking solutions to the T-cell receptor epitope prediction problem (IMMREP22). *ImmunoInformatics* 2023; 9:100024.
- Robbiani DF, et al. Convergent antibody responses to SARS-CoV-2 in convalescent individuals. *Nature* 2020; 584:437-442.
- Tan TJC, et al. Sequence signatures of two public antibody clonotypes that bind SARS-CoV-2 RBD. *Nat Commun* 2021; 12:4210.
- Sewell AK. Why must T cells be cross-reactive? *Nat Rev Immunol* 2012; 12:669-677.
