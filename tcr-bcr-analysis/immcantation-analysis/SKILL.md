---
name: bio-tcr-bcr-analysis-immcantation-analysis
description: Reconstructs B-cell clonal families, quantifies somatic hypermutation and selection, and builds antibody lineage trees with the Immcantation R suite (alakazam, shazam, scoper, dowser, tigger) on AIRR-format BCR data. Use when deriving the clonal-clustering threshold from the distToNearest bimodal valley (never a hardcoded 0.15); choosing hierarchicalClones vs spectralClones (vj vs novj) for SHM-diverged repertoires; personalizing the germline with TIGGER before mutation counting; reconstructing D-masked germlines with createGermlines; measuring R/S mutation frequency by CDR and FWR region; testing antigen-driven selection with BASELINe; comparing Hill-number diversity at equal sampling depth; and inferring IgPhyML lineage trees for affinity maturation, class-switch, and ancestral-antibody analysis.
tool_type: r
primary_tool: alakazam
---

## Version Compatibility

Reference examples tested with: alakazam 1.3+, shazam 1.2+, scoper 1.3+, dowser 2.x, tigger 1.1+ (Immcantation R suite), plus IgBLAST, Change-O, and PHYLIP/IgPhyML as external dependencies.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: `createGermlines` now lives in dowser (not shazam); BASELINe selection uses `calcBaseline`/`groupBaseline` (the old `estimateBaseline` name is gone); mutation R/S classification is set by `regionDefinition`, not a fake `mutationDefinition=MUTATION_SCHEMES$S5F` (that has no `S5F` member); the clonal threshold must come from `findThreshold`, never a literature constant.

# Immcantation Analysis

**"Find the B-cell clones and measure their affinity maturation"** -> partition SHM-diverged sequences into clonal families, quantify somatic hypermutation and selection against a reconstructed germline, and build antibody lineage trees.
- R: `shazam::distToNearest()` + `shazam::findThreshold()` (threshold), `scoper::hierarchicalClones()`/`scoper::spectralClones()` (clones), `dowser::createGermlines()` + `shazam::observedMutations()` (SHM), `shazam::calcBaseline()` (selection), `dowser::getTrees()` (lineage trees)

## The governing principle: the clonal threshold is derived, not assumed

Every downstream number in a BCR analysis -- clone counts, diversity, selection strength, tree topology -- inherits its error from one quantity: the nucleotide-distance cutoff used to group sequences into clonal families. That cutoff is NOT a literature constant. `distToNearest` computes each sequence's Hamming distance to its nearest neighbor within the same V gene, J gene, and junction length; because unrelated rearrangements almost never share V/J plus a near-identical junction by chance while clonally related sequences differ only by SHM, the resulting `dist_nearest` distribution is bimodal. `findThreshold` locates the VALLEY between the clonally-related mode (small distances) and the unrelated mode (large distances). That valley is the per-dataset threshold. A hardcoded `threshold = 0.15` is the exact anti-pattern to avoid: the valley shifts with subject, locus, sequencing depth, and chemistry, and a wrong threshold silently merges independent lineages or shatters one clone into many (Gupta 2015 *Bioinformatics* 31:3356; Nouri 2018 *Bioinformatics* 34:i341).

If the `dist_nearest` histogram is UNIMODAL (no clear valley), a fixed threshold is undefined -- switch to `spectralClones(method="novj")`, whose adaptive local threshold does not require `findThreshold`.

## Why BCR needs a different clonotype definition than TCR

TCR does not hypermutate, so all progeny of a founding T cell share the exact CDR3 nucleotide sequence and exact-CDR3 matching is correct. BCR hypermutates: members of one lineage are NOT identical, so exact-CDR3 shatters a single clone into hundreds of fragments. The field-standard BCR clone groups sequences sharing the same V gene, same J gene, and same junction LENGTH, then clusters within that partition by junction nucleotide distance at the derived threshold. Use nucleotide (not amino-acid) junction distance -- SHM is a nucleotide process and codon degeneracy would blur it.

| Method | How it clusters | Best when | Fails when |
|--------|-----------------|-----------|------------|
| `hierarchicalClones` | Single-linkage on junction Hamming distance within V/J/length partitions, cut at the `findThreshold` value | `dist_nearest` is clearly bimodal; a defensible fixed threshold exists | Unimodal distance histogram (threshold undefined); heavily diverged clones fragment |
| `spectralClones(method="novj")` | Spectral clustering with an adaptive local junction-similarity threshold; no fixed cutoff needed | Unimodal repertoires where no `findThreshold` valley exists | Very small groups (spectral needs several sequences) |
| `spectralClones(method="vj")` | Adds shared V/J SHM (targeting model) to junction homology | SHM-driven within-clone divergence pulls junctions apart; a mutated clone would otherwise be split | Needs `germline_alignment`/`sequence_alignment` and is slower |

Verify current best practice against the SCOPer vignette before committing to a method; the spectral `vj` model is the reason spectral clustering holds diverged clones together where a fixed threshold fragments them.

## Pipeline order (load-bearing)

This order is not interchangeable; getting it wrong silently corrupts mutation and selection counts.

0. TIGGER genotype FIRST. An unrecorded personal germline polymorphism otherwise reads as recurrent SHM at a fixed position -- it inflates mutation and selection counts AND adds spurious junction distance that corrupts `distToNearest`.
1. `createGermlines` (per-sequence) to reconstruct the D-masked germline BEFORE any mutation counting (mutation = observed vs inferred germline).
2. `distToNearest` -> `findThreshold` to derive the threshold.
3. Clonal clustering (`hierarchicalClones`/`spectralClones`).
4. `createGermlines` again per-clone (clone consensus germline), then `observedMutations` with the CDR3/junction MASKED (the D-masked germline handles this; junctional N/P bases have no template).
5. BASELINe selection (`calcBaseline` -> `groupBaseline`) with a codon+motif-aware null -- raw R/S is biased by germline codon structure and SHM hotspot/transition bias, so naive R/S is not selection.
6. Dowser lineage trees.

Immcantation reads and writes one AIRR TSV. Expected columns: `sequence_id`, `v_call`, `j_call`, `junction`, `junction_length`, `sequence_alignment`, `germline_alignment_d_mask`, `clone_id` (plus `locus` and `cell_id` for single-cell). These are lowercase snake_case; legacy UPPERCASE Change-O names (`V_CALL`, `JUNCTION`, `CLONE`) are deprecated and mixing schemas is a silent failure.

## Personalize the germline with TIGGER

**Goal:** Build the subject's own V-gene genotype so germline polymorphisms are not miscounted as somatic mutations.

**Approach:** Detect novel alleles from the mutation-frequency-vs-position signature, infer the personal genotype, and re-call V alleles against it before anything downstream.

```r
library(tigger)

ighv <- readIgFasta('IMGT_Human_IGHV.fasta')             # named vector of germline V alleles
novel <- findNovelAlleles(db, germline_db = ighv, v_call = 'v_call', nproc = 1)
genotype <- inferGenotypeBayesian(db, germline_db = ighv, novel = novel, find_unmutated = TRUE)
gt_seqs <- genotypeFasta(genotype, germline_db = ighv, novel = novel)
db <- reassignAlleles(db, genotype_db = gt_seqs)         # collapse ambiguous calls to alleles the subject carries
```

## Derive the clonal threshold

**Goal:** Obtain the per-dataset nucleotide-distance cutoff that separates clonally related from unrelated sequences.

**Approach:** Compute each sequence's distance to its nearest same-V/J/length neighbor, then find the valley of the bimodal distribution. Inspect the histogram before trusting the value.

```r
library(shazam)

db <- distToNearest(db, sequenceColumn = 'junction', vCallColumn = 'v_call',
                    jCallColumn = 'j_call', model = 'ham', normalize = 'len', nproc = 1)
# Single-cell: add cellIdColumn='cell_id', locusColumn='locus', onlyHeavy=TRUE
#   (light chains lack the junction diversity to define clones alone)

thr_obj <- findThreshold(db$dist_nearest, method = 'density')   # 'gmm' makes the FP/FN tradeoff explicit
threshold <- thr_obj@threshold                                   # S4 slot; NA/unimodal -> use spectralClones('novj')
plot(thr_obj)                                                    # confirm bimodality before proceeding
```

## Cluster sequences into clonal families

**Goal:** Group SHM-diverged sequences descended from one naive B cell into clones.

**Approach:** Cluster within V/J/junction-length partitions at the derived threshold; for single-cell paired data, cluster on heavy chains, then resolve light chains as a separate step.

```r
library(scoper)

results <- hierarchicalClones(db, threshold = threshold, method = 'nt', linkage = 'single')
db <- as.data.frame(results)                       # adds clone_id

# Single-cell paired BCR: cluster on heavy only, then split clones by light-chain V/J.
# The scoper only_heavy/split_light args are DEPRECATED; use dowser::resolveLightChains:
# db <- dowser::resolveLightChains(db)

# Unimodal repertoire (no clear threshold): adaptive, SHM-aware alternative
# db <- as.data.frame(spectralClones(db, method = 'vj',
#     germline = 'germline_alignment', sequence = 'sequence_alignment'))
```

## Reconstruct germline and quantify SHM

**Goal:** Measure somatic hypermutation as replacement (R) and silent (S) frequency by region, the signal of affinity maturation.

**Approach:** Rebuild the D-masked clonal germline, then compare each observed V-region to it. Use frequency (not raw counts) when coverage varies, and restrict to the V segment so the untemplated junction is excluded.

```r
library(dowser)

references <- readIMGT('imgt/human/vdj')           # IMGT-gapped V/D/J reference dir
db <- createGermlines(db, references)              # per-clone germline; adds germline_alignment_d_mask

db <- observedMutations(db, sequenceColumn = 'sequence_alignment',
                        germlineColumn = 'germline_alignment_d_mask',
                        regionDefinition = IMGT_V,             # V only; stops before CDR3/junction
                        frequency = TRUE, nproc = 1)
# Adds mu_freq_cdr_r, mu_freq_cdr_s, mu_freq_fwr_r, mu_freq_fwr_s
# For property-based R/S use mutationDefinition = CHARGE_MUTATIONS (or HYDROPATHY/POLARITY/VOLUME).
# S5F is a TARGETING model (HH_S5F) for selection, NOT a mutationDefinition.
```

## Test for selection (BASELINe)

**Goal:** Decide whether replacement mutations are enriched (positive selection, typically CDR) or depleted (purifying, typically FWR) beyond what SHM alone produces.

**Approach:** Compute the expected R/S per region from the germline under an SHM targeting model, form a posterior over selection strength per sequence, then convolve posteriors within groups. Analyze one representative per clone so shared ancestral mutations are not double-counted.

```r
baseline <- calcBaseline(db, testStatistic = 'focused', regionDefinition = IMGT_V, nproc = 1)
grouped <- groupBaseline(baseline, groupBy = 'sample_id')   # convolves per-sequence PDFs
# testBaseline(grouped, groupBy='sample_id') for significance; sigma>0 = positive selection
```

## Compare diversity at equal depth

**Goal:** Compare clonal diversity across samples without confounding by sequencing depth.

**Approach:** Report a Hill-number profile with uniform resampling to equal N and bootstrap CIs; comparing raw diversity across unequal-depth libraries measures depth, not biology.

```r
library(alakazam)

div <- alphaDiversity(db, group = 'sample_id', clone = 'clone_id',
                      min_q = 0, max_q = 2, step_q = 0.1,      # q=0 richness, q=1 Shannon, q=2 Simpson
                      ci = 0.95, nboot = 200)                  # uniform=TRUE (default) resamples to equal N
plot(div)
```

## Build lineage trees

**Goal:** Reconstruct each clone's antibody lineage to trace affinity maturation, class switching, and ancestral (intermediate) antibodies.

**Approach:** Build clonally-collapsed, germline-rooted trees under IgPhyML's HLP codon model, which encodes SHM's context-dependence, non-reversibility, and known germline root -- assumptions that standard phylogenetics violates.

```r
clones <- formatClones(db, traits = 'c_call', minseq = 3)     # collapse duplicates, attach clonal germline
trees <- getTrees(clones, build = 'igphyml',
                  igphyml = '/usr/local/share/igphyml/src/igphyml', nproc = 1)
plots <- plotTrees(trees)                                     # ggtree, germline-rooted; color tips by trait
# findSwitches(clones, ...) + testSP/testSC reconstruct isotype/tissue switching across bootstrap trees.
# Legacy: alakazam::buildPhylipLineage() (PHYLIP dnapars max-parsimony) still exists but is superseded.
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Clone counts differ wildly from a published study | Hardcoded `threshold = 0.15` instead of the data's valley | Run `distToNearest` -> `findThreshold`; read `@threshold`; inspect the histogram |
| `observedMutations` gives near-zero or nonsensical mutations | Counted before `createGermlines` (no reconstructed germline) | Run `createGermlines` first; compare against `germline_alignment_d_mask` |
| Inflated R mutations concentrated in CDR3 | Junction/CDR3 not masked; junctional N/P bases have no template | Use the D-masked germline and `regionDefinition = IMGT_V` (V only) |
| `MUTATION_SCHEMES$S5F` errors or gives odd R/S | No `S5F` member exists; S5F is a targeting model, not a mutation definition | Drop it (default R/S by AA identity) or use `CHARGE_MUTATIONS`; use `HH_S5F` only as a targeting model |
| `estimateBaseline` not found | Renamed | Use `calcBaseline` then `groupBaseline`/`testBaseline` |
| Recurrent "mutation" at the same position across many sequences | Unrecorded personal germline allele scored as SHM | Run TIGGER (`findNovelAlleles`/`inferGenotypeBayesian`/`reassignAlleles`) before germline reconstruction |
| Diversity differences vanish or invert after resequencing | Compared raw diversity across unequal-depth samples | Use `alphaDiversity` with uniform resampling (default) and bootstrap CIs |
| Same clone appears in two individuals | Pooled clones across subjects with private genotypes | Cluster clones within each subject; treat cross-subject sharing as a separate convergence question |
| Unimodal `dist_nearest` histogram, `findThreshold` returns NA | No clear valley (e.g. low-SHM or shallow repertoire) | Use `spectralClones(method = 'novj')` (adaptive threshold) |

## Related Skills

- mixcr-analysis - Produce AIRR/clonotype input for BCR
- scirpy-analysis - Single-cell BCR integration and handoff
- specificity-annotation - Convergent/public antibody signatures
- phylogenetics/tree-visualization - General lineage-tree plotting concepts
- phylogenetics/modern-tree-inference - Phylogenetic inference background
- workflows/tcr-pipeline - End-to-end orchestration

## References

- Gupta NT, Vander Heiden JA, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH. Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. *Bioinformatics* 2015, 31(20):3356-3358.
- Vander Heiden JA, Yaari G, Uduman M, Stern JNH, O'Connor KC, Hafler DA, Vigneault F, Kleinstein SH. pRESTO: a toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires. *Bioinformatics* 2014, 30(13):1930-1932.
- Yaari G, Uduman M, Kleinstein SH. Quantifying selection in high-throughput immunoglobulin sequencing data sets (BASELINe). *Nucleic Acids Research* 2012, 40(17):e134.
- Yaari G, Vander Heiden JA, Uduman M, et al. Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data (S5F). *Frontiers in Immunology* 2013, 4:358.
- Gadala-Maria D, Yaari G, Uduman M, Kleinstein SH. Automated analysis of high-throughput B-cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles (TIGGER). *PNAS* 2015, 112(8):E862-E870.
- Nouri N, Kleinstein SH. A spectral clustering-based method for identifying clones from high-throughput B cell repertoire sequencing data (SCOPer). *Bioinformatics* 2018, 34(13):i341-i349.
- Hoehn KB, Pybus OG, Kleinstein SH. Phylogenetic analysis of migration, differentiation, and class switching in B cells (Dowser). *PLoS Computational Biology* 2022, 18(4):e1009885.
- Hoehn KB, Lunter G, Pybus OG. A phylogenetic codon substitution model for antibody lineages (IgPhyML). *Genetics* 2017, 206(1):417-427.
- Stern JNH, Yaari G, Vander Heiden JA, et al. B cells populating the multiple sclerosis brain mature in the draining cervical lymph nodes. *Science Translational Medicine* 2014, 6(248):248ra107.
