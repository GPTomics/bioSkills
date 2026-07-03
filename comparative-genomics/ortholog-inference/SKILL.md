---
name: bio-comparative-genomics-ortholog-inference
description: Infer orthologous genes and gene families across species using OrthoFinder3 (HOG-based phylogenetic orthology), SonicParanoid2, Broccoli, ProteinOrtho, OMA / FastOMA hierarchical orthologous groups, eggNOG-mapper, JustOrthologs, and TOGA whole-genome-alignment orthology. Use when building single-copy ortholog sets for phylogenomics, classifying co-orthologs and in/out-paralogs after gene duplication, propagating functional annotation via orthology with awareness of the ortholog conjecture, distinguishing speciation from duplication via gene-tree species-tree reconciliation, computing Quest-for-Orthologs benchmark performance, or running synteny-aware ortholog detection in WGD-affected lineages.
tool_type: mixed
primary_tool: OrthoFinder
---

## Version Compatibility

Reference examples tested with: OrthoFinder 3.0+ (bioRxiv 2025.07.15.664860), SonicParanoid 2.0.8+ (Cosentino 2024), Broccoli 1.2+ (Derelle 2020), ProteinOrtho 6.3.0+ (Lechner 2011 + recent), OMA standalone 2.6.0+, FastOMA 0.3.5+ (Majidian 2025), eggNOG-mapper 2.1.12+, JustOrthologs 2.0+, DIAMOND 2.1.10+, MMseqs2 17-b804f+, IQ-TREE 2.3.6+, BUSCO 5.7+, Compleasm 0.2.7+, BioPython 1.84+, R 4.4+ for downstream tree-based reconciliation.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `orthofinder --help`; `sonicparanoid --help`; `oma --help`
- Python: `pip show eggnog-mapper`; `which fastoma`

If code throws `Diamond requires N more sequences than provided`, `KeyError on species tree taxa`, `STAG branch length 0`, or `HOG file format mismatch`, the OrthoFinder v2 -> v3 file layout changed (Orthogroups/ -> Phylogenetic_Hierarchical_Orthogroups/; rooted gene trees are now per-HOG); update parsing accordingly.

# Ortholog Inference

**"Find the orthologs of my gene(s) across these species"** -> Choose between graph-based (RBH / similarity-clustering: fast, lower recall) and tree-based (gene-tree reconciliation: higher accuracy, slower) frameworks; recognize that "orthology" splits into 1:1, 1:many, many:many, and the practical unit for most pipelines is the **HOG (Hierarchical Orthologous Group)** -- a maximal cluster of genes descended from a single ancestral gene at a defined taxonomic level (Altenhoff 2013 PLoS ONE 8:e53786). The "ortholog conjecture" (orthologs more functionally similar than paralogs) is supported but weakly (Altenhoff 2012 PLoS Comp Biol 8:e1002514); don't treat 1:1 ortholog labeling as automatic functional equivalence.

- CLI: `orthofinder -f proteomes/ -t 16 -M msa` -- HOG output in v3 layout
- CLI: `sonicparanoid -i proteomes/ -o output --mode default` -- ML predictor + protein language model
- CLI: `broccoli.py -dir proteomes/ -threads 16` -- direct OG with chimeric handling
- CLI: `oma standalone` HOG inference at every taxonomic level
- CLI: `proteinortho6.pl --project=run proteomes/*.faa` -- graph clustering with optional synteny
- CLI: `emapper.py -i proteins.faa --output project --cpu 16` -- eggNOG annotation transfer

## Algorithmic Taxonomy

| Tool | Approach | Output | Strength | Fails when |
|------|----------|--------|----------|------------|
| OrthoFinder3 (Emms & Kelly 2025 bioRxiv) | DIAMOND2-ultra search -> gene trees -> rooted via STRIDE -> HOG inference at every node | HOGs at every taxonomic level, orthologs, species tree, gene-duplication events | Best Quest-for-Orthologs benchmark (Altenhoff 2024); HOGs are 12% more accurate than v2 orthogroups, 20% with outgroup | Slow with > 200 species without `-c 1` clustering pre-step; new file layout breaks old parsers |
| SonicParanoid2 (Cosentino 2024 GB 25:195) | RBH + ML classifier with protein language model embeddings | Pairwise orthologs + orthogroups | Best accuracy among graph methods at competitive speed; ML correction reduces InParanoid errors | InParanoid lineage still has paralogy confusion in WGD clades; close-relative duplications hard to resolve |
| Broccoli (Derelle 2020 MBE 37:3389) | k-mer similarity -> directed graph -> chimera-aware OG inference | OGs + chimera flagging | Robust to chimeric assemblies; runs without species tree | Less accurate than OrthoFinder on benchmark; no HOG output |
| ProteinOrtho 6 (Lechner 2011 BMC Bioinf 12:124) | Pairwise BLAST/DIAMOND + connectivity graph; optional synteny module | Orthogroups + synteny option | Fast; scales to 1000+ genomes; `-synteny` enables co-linear-anchor filtering | Lower recall than tree-based; synteny module slow and requires GFFs |
| OMA standalone (Altenhoff 2019 GR 29:1152) | Strict RBH + verification + HOG inference | HOG database; orthologs at each taxonomic level | Conservative; highest precision in QfO benchmarks; "Fast" mode for prefiltering | Lowest recall among methods (Altenhoff 2024); slow for large datasets |
| FastOMA (Majidian 2025 Nat Methods 22:269) | OMA HOG inference with GPU-accelerated DIAMOND + Roothap | Same HOG output as OMA, 10-100x faster | Scales OMA to 1000+ genomes | Newer; less benchmarked in production |
| eggNOG-mapper 2 (Cantalapiedra 2021 MBE 38:5825) | DIAMOND/MMseqs2 against eggNOG 5 -> map to precomputed orthogroups | OGs + functional annotation (GO/KEGG/COG) | Standard for functional annotation propagation; phylogeny-aware | Pre-computed OGs; cannot add novel species coherently; only as fresh as eggNOG release |
| JustOrthologs 2 (Miller 2019 Bioinformatics 35:546) | DNA-based; exon-aware; close-species RBH | Pairwise orthologs | Extremely fast for closely related species (same family); preserves splice variants | Only suitable for closely related species |
| TOGA (Kirilenko 2023 Science 380:eabn3107) | Whole-genome-alignment chain -> ML projection + intactness classification | Per-query orthologs with intact/lost/missing call | Modern paradigm for vertebrate-scale orthology; handles gene loss explicitly; integrates with CESAR 2.0 | Requires WGA (Cactus); not designed for prokaryotes or fungi |
| HOGENOM / HOGsuite (Penel 2009 BMC Bioinf 10:S3) | Tree-based HOGs in databases | Pre-computed HOG database | Legacy; for downstream use of stored HOGs | Not for new computation; outdated taxon sampling |

Methodology evolves; Quest-for-Orthologs benchmarks (Altenhoff 2024 NAR Genom Bioinform 6:lqae167) refresh annually. Verify the current QfO benchmark results before locking on a single tool for novel benchmarking-grade work.

## Decision Tree by Experimental Scenario

| Scenario | Recommended approach | Why |
|----------|------------------------|-----|
| 5-50 vertebrate / animal genomes for phylogenomics | OrthoFinder3 `-M msa -A mafft -T iqtree` | HOG quality at moderate scale; species tree included |
| 50-200 genomes, any clade | OrthoFinder3 default; SonicParanoid2 as second method | Cross-validate consensus HOGs |
| 200-1000 genomes, scaling required | FastOMA or SonicParanoid2 | OrthoFinder3 slow; cluster-first option (`-c 1`) helps |
| 1000+ bacterial genomes (pangenome scope) | Use [[pangenome-analysis]] pipeline; OrthoFinder unsuitable | Pangenome-specific tools needed |
| Closely related strains / same family | JustOrthologs or ProteinOrtho with synteny | Fast; preserves splice variants |
| WGD-affected lineage (plants, salmonids, teleosts) | OrthoFinder3 + synteny verification (see [[synteny-analysis]]); or GENESPACE | WGD inflates paralogy; synteny anchoring required |
| Functional annotation transfer | eggNOG-mapper 2 | Pre-computed; phylogeny-aware; integrates GO/KEGG |
| Single-copy orthologs for concatenation phylogenomics | OrthoFinder3 `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` filter for 1:1 | Most rigorous HOG layer; standard concatenation input |
| Gene-loss detection across mammals/birds at scale | TOGA + CESAR 2.0 | Explicit intact/lost classification; handles assembly-gap noise |
| Functional annotation transfer in poorly characterized genome | eggNOG-mapper + OrthoFinder3 single-copy orthologs cross-validated | Functional priors + phylogenetic confidence |
| Ortholog detection in highly fragmented assembly | BUSCO/Compleasm first to assess completeness; flag affected OGs | Assembly fragmentation creates false absence -> spurious lineage-specific losses |
| Distant homolog detection (50-200 Myr divergence) | MMseqs2 sensitive (`-s 7.5`) + OrthoFinder3 | Default DIAMOND misses distant homologs; sensitive search needed |
| Phylogenetic orthology with paralog-tolerant species trees | ASTRAL-Pro2 on OrthoFinder gene trees | Coalescent species tree from gene trees; handles paralogy explicitly |
| HOG-level functional propagation across many species | OMA standalone HOG database | Strict hierarchical orthology; supports functional propagation per taxonomic level |
| Synteny-anchored orthology in repeat-heavy genome | ProteinOrtho `-synteny` or GENESPACE | Filters tandem duplicates; uses gene order to disambiguate paralogs |

## Per-Method Failure Modes

### Hidden paralogy from missing outgroup

**Trigger:** OrthoFinder run without a sufficiently distant outgroup; only ingroup taxa included.

**Mechanism:** OrthoFinder roots gene trees via STRIDE / STAG using outgroup-based duplication signals. Without an outgroup, the root inferred for each gene tree may be internal, mistakenly classifying a duplication-then-loss pattern as orthology. The "1:1 orthologs" returned may actually be hidden paralogs (Emms & Kelly 2017 MBE 34:3267 STRIDE; benchmark shows 8% improvement in HOG accuracy with outgroup).

**Symptom:** Phylogenomic concatenation produces low-bootstrap species tree; per-gene trees show inconsistent rooting; same-clade species show longer-than-expected branches in single-copy ortholog trees.

**Fix:** Include >= 1 outgroup taxon at the next-higher taxonomic level (sister phylum / class / order). For vertebrates, use cyclostomes (lamprey) as outgroup for jawed vertebrates; for plants, use a non-flowering plant for angiosperms. Re-run with outgroup; cross-check HOG file `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` for 1:1 stability.

### Splice isoforms inflating copy number

**Trigger:** Proteome FASTA contains multiple isoforms per gene (Ensembl, RefSeq with `-NR` flag).

**Mechanism:** OrthoFinder / SonicParanoid / OMA treat each protein sequence as a gene; isoforms are clustered together within orthogroups but inflate the apparent gene count per species, producing artifactual "co-orthologs" that are just isoforms of one gene.

**Symptom:** Per-species gene count > 2x what gene-annotation pipeline reported; orthogroups contain multiple proteins from same gene; gene names contain isoform suffixes (`.1`, `-iso1`, etc.).

**Fix:** Pre-filter proteomes to longest isoform per gene. Use `agat_sp_keep_longest_isoform.pl` (AGAT toolkit), Biopython snippet on Ensembl gene-isoform mapping, or `Trinotate` longest-ORF picker. OrthoFinder3 ships a wrapper `tools/primary_transcript.py` for Ensembl/UniProt format.

### Annotation heterogeneity inflating lineage-specific OGs

**Trigger:** Genomes annotated by different pipelines (Augustus, MAKER, Funannotate, NCBI RefSeq) mixed in one analysis.

**Mechanism:** Pipelines differ in handling of intron predictions, gene boundaries, and small-gene filtering. A gene predicted by Augustus but missed by MAKER appears as an apparent lineage-specific gene in the MAKER-annotated species, producing spurious "species-specific" orthogroups.

**Symptom:** CAFE (gene family evolution) shows extreme expansions / contractions for species with different annotation pipelines; per-species "unique" gene count varies 5-10x between technically similar genomes.

**Fix:** Re-annotate all genomes with one pipeline (e.g. BRAKER3 or Funannotate) before orthology. When re-annotation isn't possible, normalize via BUSCO/Compleasm completeness (filter OGs absent from species with < 95% BUSCO complete). Document annotation pipeline per species in methods.

### Ortholog conjecture violations

**Trigger:** Transferring GO/KEGG functional annotation from 1:1 ortholog without testing functional divergence.

**Mechanism:** Studies show orthologs are *weakly* more functionally similar than paralogs (Altenhoff 2012 PLoS Comp Biol 8:e1002514; Stamboulian 2020 Bioinformatics 36:i219), but effect size is small and dependent on evolutionary distance. Subfunctionalization (Force 1999 Genetics 151:1531), neofunctionalization, and dosage subfunctionalization can rapidly differentiate 1:1 orthologs.

**Symptom:** Transferred annotation contradicts species-specific experimental data; ortholog has fold-change different expression patterns; rapid evolution (dN/dS > 0.3) on one branch only.

**Fix:** For high-confidence annotation transfer, require: (1) 1:1 orthology AND (2) low branch dN/dS (< 0.2 for both branches) AND (3) GO evidence code ISO with experimental support upstream. Tag annotation as "predicted from ortholog" rather than equating function. eggNOG ortholog-conjecture-aware mode provides confidence scoring.

### RBH symmetric-but-wrong errors

**Trigger:** Reciprocal best hits (RBH) method on a gene that has been replaced by a paralog in one lineage.

**Mechanism:** Lineage A retains the original ortholog; Lineage B lost the ortholog and replaced its function with a paralog (xenologous replacement). RBH between A and B identifies the paralog as the "ortholog" because it's the best hit.

**Symptom:** dN/dS on RBH ortholog pair > 1 on the B branch (suggesting positive selection but actually paralog substitution); gene tree shows the B sequence is sister to other paralogs of the A gene, not to A.

**Fix:** Use tree-based methods (OrthoFinder3, OMA HOGs) for organisms with extensive paralog history. RBH is appropriate for orthogroup pre-screening only. JustOrthologs is explicitly an RBH method, so it's affected by this; reserve for closely related species.

### Synteny ignored in WGD lineages

**Trigger:** Plant, yeast, fish, or salmonid analysis where one or more ancient WGD events occurred.

**Mechanism:** WGD doubles all genes; subsequent gene loss is biased toward certain functional categories (Dosage Balance Hypothesis; Birchler & Veitia 2007 Plant Cell 19:395). Sequence-only orthology cannot distinguish "ortholog" (single ancestral gene) from "homeolog" (paralog from WGD). In hexaploids (wheat, Brassica), this triples the confusion.

**Symptom:** Multiple co-orthologs per species in orthogroups for WGD-affected lineages; rate variation across "co-orthologs" suggests one is the true ortholog and the others are recent duplicates.

**Fix:** Use synteny-aware methods (GENESPACE -- Lovell 2022 eLife 78526; ProteinOrtho `-synteny`; OrthoFinder3 + post-hoc synteny verification). Restrict 1:1 ortholog phylogenomics to genes in single-copy syntenic regions. For deep WGD ancestry (e.g. 2R vertebrate WGD), modern orthologs of post-2R paralogs (ohnologs) can no longer be unambiguously identified by sequence; use synteny + duplication-dating.

### MAFFT-only alignment in OrthoFinder MSA mode

**Trigger:** `orthofinder -M msa -A mafft` default; very divergent orthogroups (deep taxonomy).

**Mechanism:** MAFFT default settings choose accurate-but-fast iterative refinement. For sequences with >50% divergence, MAFFT's auto choice may downgrade to FFT-NS-2, producing alignments with > 30% poorly aligned columns that distort gene-tree inference and downstream HOG construction.

**Symptom:** Per-OG MAFFT logs show "FFT-NS-2 selected"; gene trees have unstable rooting; bootstrap support < 60% for many internal branches.

**Fix:** Force `mafft-linsi` or specify `-A mafft --thread -1 --localpair --maxiterate 1000`. For deeply divergent OGs, run PRANK or MUSCLE5 post-hoc on critical OGs and re-build gene trees with IQ-TREE. PREQUAL or HmmCleaner segment-filtering improves resulting trees (Di Franco 2019 BMC Eco Evo 19:21).

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|-------------------|
| QfO benchmark adoption | 100% within-species ortholog pairs identified | Altenhoff 2024 NAR Genom Bioinform 6:lqae167; minimum competence threshold |
| BUSCO/Compleasm completeness for inclusion | >= 90% complete (single + duplicated) before orthology run | Below this, expect inflated lineage-specific OGs |
| OrthoFinder3 outgroup distance | >= one sister taxonomic level (sister phylum / class / order) | Emms & Kelly 2017 MBE 34:3267 STRIDE rooting; +20% HOG accuracy |
| Single-copy ortholog filter for phylogenomics | Present in >= 90% of species, exactly 1 copy each | Standard convention; below this, missing data biases tree inference |
| MMseqs2 sensitivity for divergent homologs | `-s 7.5` for > 50% divergence | mmseqs2 documentation; default `-s 4.0` misses distant homologs |
| ProteinOrtho `--conn` (connectivity) | >= 0.1 default; 0.2 stricter for less paralogy confusion | Lechner 2011 |
| OMA HOG inclusion criterion | RBH + Smith-Waterman score and pairwise stability | OMA convention; sub-clade HOGs nested in supergroup HOGs |
| Ortholog age threshold for functional transfer | divergence < 200 Myr OR dN/dS < 0.2 | Stamboulian 2020 Bioinformatics 36:i219; ortholog conjecture stronger at lower divergence |
| Annotation pipeline normalization | All species annotated with same pipeline OR BUSCO-completeness within 5% | Avoid annotation-heterogeneity bias on CAFE |
| eggNOG-mapper minimum score | bit score / e-value defaults; check `--seed_ortholog_score` | eggNOG-mapper docs |
| TOGA intactness classes (loss_summ_data.tsv) | I (intact), PI (partial intact), UL (uncertain loss), L (lost), M (missing/assembly gap), PM (partial missing) | Kirilenko 2023 + TOGA repo |
| TOGA orthology relationships (orthology_classification.tsv) | one2one, one2many, many2one, many2many, PG (paralogous projection) | Kirilenko 2023 |
| SonicParanoid2 confidence | >= 0.9 high-confidence orthologs; 0.5-0.9 moderate | Cosentino 2024 docs |
| Reasonable expected runtime (200 species, 16 cores) | OrthoFinder3 default ~8-24h; SonicParanoid2 ~2-6h; FastOMA ~3-8h | Hardware-dependent benchmarks |

## OrthoFinder3 Standard Workflow

**Goal:** Produce HOG-based orthology with species tree and gene-duplication events for any clade.

**Approach:** Prepare cleaned per-species proteomes (longest isoforms) -> include outgroup -> run with MSA + IQ-TREE option -> parse HOG output at appropriate taxonomic level.

```bash
# Pre-clean proteomes: longest isoform per gene
for f in raw_proteomes/*.faa; do
    python tools/primary_transcript.py $f > cleaned/$(basename $f)
done

# OrthoFinder v3 with MSA + IQ-TREE for tree-based HOG inference
orthofinder \
    -f cleaned/ \
    -t 16 \
    -a 4 \
    -M msa \
    -A mafft \
    -T iqtree \
    -S diamond_ultra_sens \
    -y \
    -o orthofinder_run

# Output of interest (v3 layout):
# orthofinder_run/Results_<date>/Phylogenetic_Hierarchical_Orthogroups/N0.tsv  (root-level HOGs)
# orthofinder_run/Results_<date>/Single_Copy_Orthologue_Sequences/             (single-copy MSA-ready)
# orthofinder_run/Results_<date>/Species_Tree/SpeciesTree_rooted.txt
# orthofinder_run/Results_<date>/Gene_Duplication_Events/                       (per-branch dup counts)
```

```python
'''Parse OrthoFinder v3 HOG output for downstream analysis.'''

import pandas as pd
from pathlib import Path


def load_hogs(results_dir, level='N0'):
    '''HOG file path differs from v2: Phylogenetic_Hierarchical_Orthogroups/{level}.tsv.'''
    p = Path(results_dir) / 'Phylogenetic_Hierarchical_Orthogroups' / f'{level}.tsv'
    df = pd.read_csv(p, sep='\t')
    return df.set_index('HOG')


def _is_present(v):
    '''OrthoFinder v3 HOG cells are either NaN (older versions), '' (newer), or a comma-separated gene list.'''
    return not (pd.isna(v) or v == '' or str(v).strip() == '')

def single_copy_hogs(hog_df, min_species_fraction=0.9):
    '''Return HOG IDs where every present species has exactly 1 gene and >= fraction of species present.'''
    sp_cols = [c for c in hog_df.columns if c not in ('OG', 'Gene Tree Parent Clade')]
    is_single = hog_df[sp_cols].apply(
        lambda r: all((not _is_present(v)) or ',' not in str(v) for v in r), axis=1
    )
    n_present = hog_df[sp_cols].apply(lambda r: sum(_is_present(v) for v in r), axis=1)
    keep = is_single & (n_present >= min_species_fraction * len(sp_cols))
    return hog_df.index[keep].tolist()


def classify_orthology(hog_df, sp_pair):
    '''Per-HOG classify pairwise orthology between two species.

    Returns one of: '1-1', '1-many', 'many-1', 'many-many', 'absent', 'sp1-only', 'sp2-only'.
    '''
    sp1, sp2 = sp_pair
    types = {}
    for hog, row in hog_df.iterrows():
        v1, v2 = row.get(sp1), row.get(sp2)
        n1 = 0 if not _is_present(v1) else len(str(v1).split(','))
        n2 = 0 if not _is_present(v2) else len(str(v2).split(','))
        if n1 == 0 and n2 == 0: types[hog] = 'absent'
        elif n1 == 0:           types[hog] = 'sp2-only'
        elif n2 == 0:           types[hog] = 'sp1-only'
        elif n1 == 1 and n2 == 1: types[hog] = '1-1'
        elif n1 == 1 and n2 > 1: types[hog] = '1-many'
        elif n1 > 1 and n2 == 1: types[hog] = 'many-1'
        else:                     types[hog] = 'many-many'
    return types
```

## SonicParanoid2 + Cross-Validation

**Goal:** Run a second orthology method to cross-validate OrthoFinder3 calls (consensus increases QfO benchmark precision).

**Approach:** SonicParanoid2 with default ML predictor -> intersect orthogroups with OrthoFinder HOGs at the root level -> compute Jaccard agreement.

```bash
sonicparanoid -i cleaned/ -o sp2_run --mode default --threads 16 --pfam pre-computed
# Output: sp2_run/runs/<date>/orthogroups/flat.ortholog_groups.tsv
```

```python
from itertools import combinations
import pandas as pd

def jaccard_ogs(ogs_a, ogs_b):
    '''Compute mean Jaccard between two orthogroup sets keyed by gene IDs.
    Each og is a frozenset of gene IDs.'''
    set_a = {gene: og for og in ogs_a for gene in og}
    set_b = {gene: og for og in ogs_b for gene in og}
    common = set(set_a) & set(set_b)
    jacs = []
    for g in common:
        a = set_a[g]
        b = set_b[g]
        j = len(a & b) / len(a | b) if a | b else 0
        jacs.append(j)
    return sum(jacs) / len(jacs) if jacs else 0
```

Consensus single-copy HOGs (1:1 in both methods) are the highest-confidence input for downstream phylogenomic concatenation or selection analyses.

## TOGA for Vertebrate-Scale Orthology with Gene-Loss

**Goal:** Identify orthologs and classify gene-loss/intactness across hundreds of mammal or bird genomes.

**Approach:** Run Progressive Cactus whole-genome alignment (see [[whole-genome-alignment]]) -> TOGA projects reference genes through chains -> classifies each query gene as I/PI/UL/L/M/PM.

```bash
# After Cactus alignment producing reference-query chain files
toga.py \
    --chain cactus_chain.bb \
    --bed reference_genes.bed \
    --tDB target.2bit \
    --qDB query.2bit \
    --nextflow_dir nf_pipeline_dir \
    --pn project_name \
    --cpu 64 \
    --quiet

# Output:
#   project_name/loss_summ_data.tsv          per-gene intactness call
#   project_name/orthology_classification.tsv  one-to-one / one-to-many / many-to-many
```

TOGA emits eight classes per gene in `loss_summ_data.tsv`: I (intact), PI (partial intact), UL (uncertain loss), L (lost), M (missing under assembly gap), PM (partial missing), PG (paralogous projection / no orthologous chain found), N (no data). For evolutionary analysis, treat I as functional; PI as fragmented (often functional but caveat); L + UL as candidates for true gene loss after manual review of read coverage; M / PM are assembly-quality issues, not biology; PG indicates the algorithm found only paralogous chains. CESAR 2.0 provides exon-aware coding annotation projection used internally.

## Functional Annotation Transfer with eggNOG-mapper

**Goal:** Propagate GO / KEGG / EC annotations from curated orthologs to a novel proteome.

**Approach:** Run eggNOG-mapper with appropriate evolutionary level -> integrate with OrthoFinder HOGs for confidence stratification.

```bash
emapper.py \
    --output project \
    -i novel_proteome.faa \
    --cpu 16 \
    --decorate_gff genome.gff \
    --tax_scope auto \
    --target_orthologs all \
    --evidence_type experimental \
    --pfam_realign denovo
```

Output `project.emapper.annotations` includes per-gene: best ortholog, taxonomic scope of OG, GO terms, KEGG pathways, COG category. Confidence stratification: experimental-evidence GO (EXP, IDA, IPI, IMP, IGI, IEP) > non-experimental (ISO from 1:1 close ortholog) > IEA (electronic). Tag annotation transfer evidence in output.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| OrthoFinder 1:1; SonicParanoid 1:many | Recent duplication in one lineage; OrthoFinder collapsed | Inspect gene tree; if duplication post-speciation, both are co-orthologs (correct in SP2) |
| OrthoFinder 1:1; ProteinOrtho merges into single OG | Different connectivity threshold | OrthoFinder HOG is more granular; both correct at different levels |
| OMA "1:1 across all species"; OrthoFinder "split orthogroups" | OMA more conservative (strict RBH) | Use OMA for stringent comparative analyses; OrthoFinder for broader recall |
| OrthoFinder OG contains tandem duplicates | Tandem duplications not collapsed | Use ProteinOrtho `-synteny` or GENESPACE post-filter; or apply MCScanX tandem detection (default 5-gene window) |
| eggNOG ortholog disagrees with OrthoFinder | eggNOG uses fixed reference set; OrthoFinder uses the sampled species | Trust OrthoFinder for sampled-species 1:1; trust eggNOG for functional annotation from curated references |
| TOGA reports "Lost" but OrthoFinder finds ortholog | TOGA chain-projection failed; assembly gap | Re-check assembly contig at locus; if gap, treat as missing (M) not lost (L) |
| All methods disagree on a single OG | Likely tandem duplicates, chimeric assembly, or hidden paralogy | Manual gene-tree inspection; treat OG as low-confidence |

**Operational rule for publication:** Cross-validate single-copy orthologs across 2+ methods (OrthoFinder3 + SonicParanoid2 or OrthoFinder3 + OMA); use consensus HOG set for phylogenomics. Functional annotation transfer requires either 1:1 orthology + low branch dN/dS, or eggNOG ISO/EXP evidence. Report annotation pipeline normalization in methods.

## Cohort Gotchas

- **WGD lineages:** Salmonids (recent Ss4R WGD), teleosts (Ts3R), plants (1-4 rounds), yeast (2 rounds), opisthokonts (1-2 ancient) -- require synteny-aware orthology
- **Highly fragmented assemblies:** N50 < 100 kb produces extensive false absence; flag affected species and exclude from CAFE
- **Polyploids:** modern polyploids are multi-genome individuals; assign subgenomes before orthology (see [[whole-genome-duplication]])
- **Genome size variation:** Drosophila ~150 Mb vs Locust ~5 Gb -- repeat-dominated large genomes have inflated false positive rates from TE-derived proteins; filter TE proteins before orthology
- **Domain rearrangements:** chimeric gene-fusion proteins (e.g. Jumonji-domain) place into multiple OGs; OrthoFinder3 + Broccoli chimera detection helpful
- **HGT-affected genes (prokaryotes):** standard orthology returns "vertical orthologs" that may not exist; use ALE-aware approaches (see [[gene-tree-species-tree-reconciliation]])

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why this orthology method?" | OrthoFinder3 / OMA / SonicParanoid2 chosen based on QfO benchmark (Altenhoff 2024) and clade-appropriate scale; consensus single-copy HOGs reported |
| "Outgroup?" | Outgroup taxon X from sister taxonomic level included; STRIDE-rooted gene trees |
| "Isoforms?" | Longest-isoform-per-gene pre-filter applied via AGAT / OrthoFinder primary_transcript.py |
| "Annotation pipeline heterogeneity?" | All species annotated with Y pipeline OR BUSCO completeness within Z% across species |
| "WGD?" | Acknowledged; synteny-aware verification via GENESPACE / ProteinOrtho synteny |
| "Functional transfer evidence?" | 1:1 ortholog + dN/dS < 0.2 + eggNOG EXP/IDA evidence + tagged as predicted |
| "Cross-validation?" | OrthoFinder3 + second method (SonicParanoid2 / OMA); consensus single-copy HOGs used |
| "Ortholog conjecture caveat?" | Acknowledged (Altenhoff 2012; Stamboulian 2020); high-confidence single-copy 1:1 within short divergence; functional divergence flagged |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| OrthoFinder "no orthogroups found" | Proteomes empty or wrong path | Check FASTA non-empty; default file extension `.fa` `.faa` `.fasta` recognized |
| HOG file Phylogenetic_Hierarchical_Orthogroups missing | v3 layout new; parser expects v2 path | Use `Phylogenetic_Hierarchical_Orthogroups/N0.tsv` for v3; OrthoFinder v2 used `Orthogroups/Orthogroups.tsv` |
| All single-copy HOGs are tiny (1-2 genes) | Outgroup too distant; only paralog-free genes survive | Add intermediate outgroup; relax single-copy fraction to 0.7-0.8 |
| SonicParanoid2 ML model fails to download | Network or model registry issue | Pre-download with `sonicparanoid-get-test-data`; use `--mode legacy` as fallback |
| OMA "infinity recursion in HOG" | Self-referential HOG (rare bug) | Update to latest OMA; bug reports on GitHub for known cases |
| eggNOG-mapper says "no hits" for half the proteome | Default DIAMOND sensitivity too low | Add `--sensmode more-sensitive` or use `--sensmode ultra-sensitive` |
| TOGA "no chain file" | Cactus output not properly converted | Use `halSynteny` and `chainNet` UCSC pipeline; verify chain file format |
| Gene IDs in HOG have unexpected suffix `__rev` | Synteny strand convention | Strip the suffix; check OrthoFinder log for the convention applied |
| Per-species gene count radically different from annotation | Isoforms not filtered | Re-filter to longest isoform per gene |
| Annotation transfer says "Hypothetical" everywhere | eggNOG version mismatch with proteome | Update eggNOG database to current; v5 vs v6 has different OG layout |
| ProteinOrtho synteny option times out | Slow synteny module on huge genomes | Skip synteny for screening; run only on candidate OGs |

## Tool Installation Notes

```bash
# OrthoFinder3
conda install -c bioconda orthofinder=3.0
# Or via Python wrapper: pip install orthofinder3

# SonicParanoid2
conda install -c bioconda sonicparanoid

# Broccoli
git clone https://github.com/rderelle/Broccoli && cd Broccoli && pip install .

# ProteinOrtho 6
conda install -c bioconda proteinortho

# OMA standalone (Linux)
wget https://omabrowser.org/standalone/OMA.tgz && tar xf OMA.tgz && cd OMA && ./install.sh
# FastOMA
pip install fastoma

# eggNOG-mapper 2
conda install -c bioconda eggnog-mapper

# JustOrthologs
git clone https://github.com/ridgelab/justOrthologs && cd justOrthologs && python setup.py install

# TOGA
conda env create -f https://raw.githubusercontent.com/hillerlab/TOGA/master/toga_env.yml
# Requires nf-core and Nextflow

# Quality control
conda install -c bioconda busco compleasm
```

For Quest-for-Orthologs benchmark submission, follow https://orthology.benchmarkservice.org/ instructions; the benchmark refreshes annually.

## References

- Fitch WM 1970 Syst Zool 19:99 (ortholog / paralog definitions)
- Sonnhammer EL & Koonin EV 2002 Trends Genet 18:619 (in-paralog / out-paralog)
- Altenhoff AM & Dessimoz C 2009 PLoS Comp Biol 5:e1000262 (orthology-method benchmark/assessment)
- Altenhoff AM et al 2012 PLoS Comp Biol 8:e1002514 (ortholog conjecture quantification)
- Altenhoff AM et al 2013 PLoS ONE 8:e53786 (HOGs from orthologous gene pairs, GETHOGs)
- Altenhoff AM et al 2019 Genome Res 29:1152 (OMA standalone)
- Altenhoff AM et al 2024 NAR Genom Bioinform 6:lqae167 (Quest for Orthologs benchmark)
- Emms DM & Kelly S 2019 Genome Biol 20:238 (OrthoFinder 2)
- Emms DM & Kelly S 2017 MBE 34:3267 (STRIDE rooting)
- Emms DM & Kelly S 2025 bioRxiv 2025.07.15.664860 (OrthoFinder 3)
- Cosentino S, Sriswasdi S & Iwasaki W 2024 Genome Biol 25:195 (SonicParanoid2)
- Derelle R et al 2020 MBE 37:3389 (Broccoli)
- Lechner M et al 2011 BMC Bioinf 12:124 (ProteinOrtho)
- Majidian S et al 2025 Nat Methods 22:269 (FastOMA)
- Cantalapiedra CP et al 2021 MBE 38:5825 (eggNOG-mapper 2)
- Miller JB et al 2019 Bioinformatics 35:546 (JustOrthologs 2)
- Kirilenko BM et al 2023 Science 380:eabn3107 (TOGA + CESAR)
- Sharma V, Schwede P & Hiller M 2017 Bioinformatics 33:3985 (CESAR 2.0)
- Stamboulian M et al 2020 Bioinformatics 36:i219 (ortholog conjecture revisited)
- Force A et al 1999 Genetics 151:1531 (subfunctionalization)
- Birchler JA & Veitia RA 2007 Plant Cell 19:395 (gene balance hypothesis)
- Nehrt NL et al 2011 PLoS Comp Biol 7:e1002073 (ortholog conjecture challenge)
- Studer RA & Robinson-Rechavi M 2009 Trends Genet 25:210 (ortholog conjecture critique)
- Penel S et al 2009 BMC Bioinf 10:S3 (HOGENOM)
- Zhang C & Mirarab S 2022 Bioinformatics 38:4949 (ASTRAL-Pro2)
- Di Franco A et al 2019 BMC Eco Evo 19:21 (alignment filtering improves trees)

## Related Skills

- comparative-genomics/synteny-analysis - Synteny-anchored ortholog disambiguation in WGD lineages
- comparative-genomics/whole-genome-duplication - Distinguishing ohnologs (WGD paralogs) from orthologs
- comparative-genomics/gene-family-evolution - Birth-death modeling on OG counts requires consistent orthology
- comparative-genomics/gene-tree-species-tree-reconciliation - DTL-aware orthology for prokaryotes / HGT-affected lineages
- comparative-genomics/positive-selection - Selection analysis on single-copy ortholog alignments
- comparative-genomics/comparative-annotation-projection - TOGA / CESAR projects orthology with gene-loss classification
- phylogenetics/modern-tree-inference - Single-copy concatenation phylogenomics input
- phylogenetics/species-trees - Coalescent species-tree from OrthoFinder gene trees
- alignment/multiple-alignment - MSA quality affects OrthoFinder MSA-mode HOG output
- genome-annotation/functional-annotation - Functional annotation propagation via eggNOG / orthology
- read-qc/quality-reports - BUSCO / Compleasm completeness affects orthology reliability
