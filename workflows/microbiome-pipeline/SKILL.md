---
name: bio-workflows-microbiome-pipeline
description: End-to-end 16S/ITS amplicon workflow from demultiplexed FASTQ to a consensus differential-abundance result, orchestrating cutadapt primer removal, per-run DADA2 ASV inference (learnErrors/mergeSequenceTables/removeBimeraDenovo), region-matched taxonomy assignment, a SEPP/Greengenes2 tree, alpha/beta diversity at a declared sampling depth (phyloseq/vegan, adonis2 paired with betadisper), compositional DA as a consensus of >=2 tools (ALDEx2/ANCOM-BC2) on unrarefied counts, and optional PICRUSt2 functional prediction gated on NSTI. Covers the stage-ordering decisions (primers before truncation, per-run error model, rarefy for diversity not DA, predicted potential not activity) and defers each per-step choice to the six microbiome skills. Use when staging an amplicon study end to end or chaining ASV inference, taxonomy, diversity, and differential abundance. For shotgun reads see workflows/metagenomics-pipeline.
tool_type: mixed
primary_tool: DADA2
workflow: true
depends_on:
  - read-qc/adapter-trimming
  - microbiome/amplicon-processing
  - microbiome/taxonomy-assignment
  - microbiome/diversity-analysis
  - microbiome/differential-abundance
  - microbiome/functional-prediction
qc_checkpoints:
  - after_denoising: "Per-sample reads tracked through filter/denoise/merge/nonchim; no merge cliff"
  - after_diversity: "Sampling depth declared; dropped-sample list reported"
  - after_da: "Consensus of >=2 CoDA tools on unrarefied counts; tools named"
---

## Version Compatibility

Reference examples tested with: DADA2 1.30+, cutadapt 4.6+, phyloseq 1.46+, vegan 2.6+, ALDEx2 1.34+, ANCOMBC 2.4+, QIIME2 2024.2+, PICRUSt2 2.5+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The error model is a PER-RUN artifact and the reference database is a versioned dependency, not just a tool: run `learnErrors` once per sequencing run, never pooled; record the SILVA/GTDB/UNITE release and the classifier training region, the tree-build method (SEPP/Greengenes2 vs de novo), the rarefaction depth, and the PICRUSt2 reference release alongside every result.

# Microbiome Pipeline

**"Run my 16S amplicon study end to end from FASTQ"** -> Stage primer removal, per-run denoising, region-matched taxonomy, a placed tree, declared-depth diversity, and a consensus differential-abundance result - deferring every per-step decision to the owning microbiome skill, because each stage's parameters silently set what the next stage can find.

Complete workflow from demultiplexed amplicon FASTQ to a confidence-graded differential-abundance result. This is an ORCHESTRATION skill: it sequences the stages and routes the science to the six microbiome category skills; it does not re-teach the per-step decisions.

Scope: amplicon (16S/ITS) reads -> ASV table -> taxonomy -> diversity + consensus DA + optional predicted function. Shotgun/WGS reads -> workflows/metagenomics-pipeline. Each method choice -> its owning microbiome skill (links below). QIIME2 artifact/provenance route -> microbiome/qiime2-workflow.

## The Single Most Important Modern Insight -- The Pipeline Is a Chain of Modeling Choices, Not a Conveyor Belt

The feature table, the diversity number, and the differential-taxa list are not observations of the community - each is the residue of a knob turned at an EARLIER stage, before any result existed. Turn the knobs differently and the answer changes. The orchestration discipline is therefore to declare and defend each stage's choice, in order, because every stage silently constrains the next:

1. **What is stripped and where reads are truncated decides what is detectable.** Primers left on corrupt the error model and inflate chimeras; truncating below the merge-overlap budget erases taxa by arithmetic. These are set before any ASV exists (microbiome/amplicon-processing).
2. **The error model is per-run, the tree is a model, and the depth is a sample-deletion knob.** Pooling runs into one `learnErrors` denoises wrong; a de novo tree from short reads injects topology noise into UniFrac; a sampling depth above some samples' totals silently drops the lowest-biomass samples (microbiome/diversity-analysis).
3. **Rarefy for diversity, never for DA; and which taxa are "significant" depends more on the tool than the biology.** Keep the raw counts, rarefy only into the diversity branch, and report a CONSENSUS of >=2 compositionally-aware tools rather than one tool's hit list (Nearing 2022 *Nat Commun* 13:342; microbiome/differential-abundance).
4. **PICRUSt2 predicts potential, not activity, and is circular with taxonomy.** Predicted function is the ASV table re-encoded through a fixed lookup - hypothesis-generating, gated on NSTI, never "more active" (microbiome/functional-prediction).

## Pipeline Stages

Each stage hands its decisions to the owning skill; this table is the routing map, not a parameter sheet.

| Stage | Operation | Owning skill (decisions) |
|-------|-----------|--------------------------|
| 0 | Read QC; confirm primers/region known | read-qc/quality-reports, read-qc/adapter-trimming |
| 1 | Remove primers (cutadapt, BEFORE truncation) | microbiome/amplicon-processing, read-qc/adapter-trimming |
| 2 | Per-RUN DADA2: filterAndTrim -> learnErrors -> dada -> mergePairs | microbiome/amplicon-processing |
| 3 | mergeSequenceTables across runs, then ONE removeBimeraDenovo | microbiome/amplicon-processing |
| 4 | Region-matched taxonomy (genus, not species, for 16S) | microbiome/taxonomy-assignment |
| 5 | Build phyloseq + a SEPP/Greengenes2 tree (NOT de novo for short reads) | microbiome/diversity-analysis, phylogenetics/tree-io |
| 6 | Alpha/beta diversity at a DECLARED depth; adonis2 + betadisper | microbiome/diversity-analysis |
| 7 | Consensus DA of >=2 CoDA tools on UNrarefied counts | microbiome/differential-abundance |
| 8 | Optional PICRUSt2 functional PREDICTION, gated on NSTI | microbiome/functional-prediction |

QIIME2 artifact/provenance alternative for the whole chain: microbiome/qiime2-workflow (`qiime cutadapt` -> `dada2 denoise-paired` -> `feature-classifier classify-sklearn` -> `fragment-insertion sepp` -> `diversity core-metrics-phylogenetic` -> `composition ancombc` -> `q2-picrust2`). Same decisions, `.qza`/`.qzv` provenance ecosystem.

## Workflow Overview

```
Demultiplexed amplicon FASTQ (per sample, per run; primers/region known)
    |
    v
[0. Read QC]                         FastQC/MultiQC -> read-qc/quality-reports
    |
    v
[1. Remove primers]  cutadapt -g FWD -G REV --discard-untrimmed   (BEFORE truncation)
    |
    v
[2. Per-RUN DADA2]   filterAndTrim -> learnErrors (one run) -> dada -> mergePairs
    |                (repeat per sequencing run; never pool FASTQs into one learnErrors)
    v
[3. Combine + chimeras]  mergeSequenceTables(run1, run2, ...) -> ONE removeBimeraDenovo
    |
    v
[4. Taxonomy]        region-matched SILVA/GTDB/UNITE classifier -> genus (16S)
    |
    v
[5. phyloseq + tree] otu+tax+sample_data + SEPP/Greengenes2 placed tree (NOT de novo)
    |
    +--> [6. Diversity]  rarefy_even_depth(declared depth) -> alpha; UniFrac/Bray + adonis2 + betadisper
    |        (report the depth AND the dropped-sample list)
    |
    +--> [7. Differential abundance]  UNRAREFIED counts -> ALDEx2 AND ANCOM-BC2/LinDA -> CONSENSUS
    |
    +--> [8. Functional prediction]  PICRUSt2 (optional) -> KO/MetaCyc POTENTIAL, report NSTI
    |
    v
ASV table + taxonomy + diversity + consensus DA hits (+ predicted potential)
```

## Stage 1 to 3: Primers, Per-Run Denoising, Chimeras

**Goal:** Turn demultiplexed reads into one chimera-free ASV table, with primers off first and the error model fit per run.

**Approach:** Strip primers with cutadapt before any quality step; run the DADA2 block (filter -> learnErrors -> dada -> mergePairs) ONCE PER sequencing run; merge the run-level tables by exact sequence string; remove chimeras once on the combined table. Decisions (truncLen budget, maxEE, pooling mode, ITS handling) live in microbiome/amplicon-processing.

```bash
# Primers come OFF before truncation - leftover primer bases corrupt the error model and look chimeric.
# -g = forward primer (5' on R1), -G = reverse primer (5' on R2); --discard-untrimmed drops primerless pairs.
cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT --discard-untrimmed \
    -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz sample_R1.fastq.gz sample_R2.fastq.gz
```

```r
library(dada2)
# Run this block ONCE PER sequencing run. truncLen is a merge-overlap budget, not just a quality cut:
# truncLen_F + truncLen_R >= amplicon_length + ~12 (mergePairs minOverlap). See amplicon-processing.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160),
                     maxEE = c(2, 2), truncQ = 2, maxN = 0, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
errF <- learnErrors(filtFs, multithread = TRUE)   # fit THIS run only - never pool runs
errR <- learnErrors(filtRs, multithread = TRUE)
mergers <- mergePairs(dada(filtFs, err = errF, multithread = TRUE), filtFs,
                      dada(filtRs, err = errR, multithread = TRUE), filtRs)
seqtab_run <- makeSequenceTable(mergers)
# Combine per-run tables (exact-sequence string is the join key), THEN one chimera removal:
seqtab <- mergeSequenceTables(seqtab_run1, seqtab_run2)   # single run: skip; pass seqtab_run
seqtab_nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = TRUE)
```

A merge cliff (near-zero `merged` column) or a large READ fraction lost as chimeric is a stage-1/2 knob problem, not bad data - see microbiome/amplicon-processing. For low-biomass studies, sequence negative/positive controls and run decontam on the ASV table before downstream analysis (microbiome/amplicon-processing; metagenomics/contamination-controls).

## Stage 4 to 6: Taxonomy, Tree, Diversity

**Goal:** Label ASVs against a region-matched reference, assemble a phyloseq object with a PLACED tree, and summarize diversity at a depth that retains samples.

**Approach:** Assign taxonomy with a classifier trained on the amplicon region (report genus for 16S, not species); filter host mitochondria/chloroplast features (universal 16S primers amplify them); build the phyloseq object and attach a SEPP/Greengenes2 placed tree rather than a de novo tree from short reads; rarefy ONLY into the diversity branch at a declared depth, report the dropped samples, and pair adonis2 with betadisper. Decisions live in microbiome/taxonomy-assignment and microbiome/diversity-analysis.

```r
library(phyloseq); library(vegan)
# minBoot 50 = DADA2/RDP default for reads <=250 nt; ranks below it return NA, not a guess.
taxa <- assignTaxonomy(seqtab_nochim, 'silva_nr99_v138.1_train_set.fa.gz', minBoot = 50, multithread = TRUE)
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), tax_table(taxa),
               sample_data(metadata), phy_tree(placed_tree))   # placed_tree from SEPP/GG2, NOT de novo
ps <- subset_taxa(ps, is.na(Order) | Order != 'Chloroplast')      # drop host organelle 16S before diversity/DA
ps <- subset_taxa(ps, is.na(Family) | Family != 'Mitochondria')

# Diversity branch ONLY: rarefy to a DECLARED depth (not min(sample_sums)) and report who was dropped.
depth <- 10000   # choose on the alpha-rarefaction plateau; below this samples are dropped - report them
ps_rare <- rarefy_even_depth(ps, sample.size = depth, rngseed = 42, replace = FALSE)
adonis2(UniFrac(ps_rare, weighted = TRUE) ~ Group, data = data.frame(sample_data(ps_rare)), permutations = 999)
permutest(betadisper(UniFrac(ps_rare, weighted = TRUE), sample_data(ps_rare)$Group))   # location vs dispersion
```

SEPP placement is `qiime fragment-insertion sepp`; a de novo `align-to-tree-mafft-fasttree` is acceptable only when no reference package fits the marker, and unweighted UniFrac on it must be treated as suspect (microbiome/diversity-analysis).

## Stage 7: Consensus Differential Abundance

**Goal:** Identify differentially abundant taxa as a confidence-graded consensus, not one tool's volcano plot.

**Approach:** On the UNrarefied counts (rarefying discards information DA needs), filter rare features, run ALDEx2 plus a second compositionally-aware tool (ANCOM-BC2 or LinDA), gate ALDEx2 on q AND effect size, and report the intersection as high-confidence. Tool mechanics live in microbiome/differential-abundance.

```r
library(ALDEx2); library(ANCOMBC)
ps_filt <- filter_taxa(ps, function(x) sum(x > 0) >= 0.10 * nsamples(ps), TRUE)   # >=10% prevalence; declared knob
counts <- as.matrix(otu_table(ps_filt)); if (!taxa_are_rows(ps_filt)) counts <- t(counts)   # taxa in ROWS, integer counts
groups <- as.character(sample_data(ps_filt)$Group)

ax <- aldex(counts, groups, mc.samples = 128, test = 't', effect = TRUE, denom = 'all')
sig_aldex <- rownames(ax)[ax$we.eBH < 0.05 & abs(ax$effect) > 1]   # q AND |effect|>1 (between-group diff exceeds within-condition dispersion); NOT p alone

# Multi-run study: carry run as a batch covariate -> fix_formula = 'run_id + Group'. (ALDEx2 cannot
# take a covariate via aldex(); use aldex.clr() + aldex.glm() on a model matrix for the run-adjusted arm.)
ab <- ancombc2(data = ps_filt, fix_formula = 'Group', p_adj_method = 'BH',   # default is 'holm' - set BH deliberately
               prv_cut = 0.10, group = 'Group', struc_zero = TRUE, pseudo_sens = TRUE)$res
dcol <- grep('^diff_Group', names(ab), value = TRUE)[1]   # coefficient = variable+factor level, verbatim case (e.g. 'Grouptreated')
sig_ancombc <- ab$taxon[ab[[dcol]] & ab[[sub('^diff_', 'passed_ss_', dcol)]]]    # significant AND pseudo-count-robust

confident <- intersect(sig_aldex, sig_ancombc)   # high-confidence; union = exploratory; name both tools
```

## Stage 8: Functional Prediction (optional)

**Goal:** Summarize predicted community functional POTENTIAL, framed as hypothesis-generating and gated on NSTI.

**Approach:** Run PICRUSt2 on the rep-seqs and ASV table, report the NSTI distribution and the read fraction dropped at `--max_nsti 2`, and restrict every claim to "potential" - never "activity". Decisions live in microbiome/functional-prediction; for MEASURED function use shotgun (metagenomics/functional-profiling).

```bash
# Predicts KO/MetaCyc POTENTIAL from who-is-there - never measured genes, never activity.
picrust2_pipeline.py -s asv_seqs.fna -i asv_table.biom -o picrust2_out -p 8 --max_nsti 2 --hsp_method mp
# Report mean/median NSTI and the ASV+read fraction dropped at NSTI>2 (marker_predicted_and_nsti.tsv.gz).
```

## Per-Stage Failure Modes

### Primers left on before truncation (stage 1)
**Trigger:** running filterAndTrim/learnErrors on reads still carrying primers. **Mechanism:** degenerate primer bases read as sequencing error, corrupting the error model and masquerading as chimeras. **Symptom:** wrong error fit, a large READ fraction lost as chimeric, inflated ASV count. **Fix:** cutadapt `--discard-untrimmed` first; order is primers -> filter -> learnErrors (microbiome/amplicon-processing).

### Pooling runs into one error model (stage 2)
**Trigger:** concatenating multiple runs' FASTQs before `learnErrors`. **Mechanism:** one error model is fit to a mixture of run-specific error structures. **Symptom:** ASVs that vanish or appear when runs are split. **Fix:** per-run inference, then `mergeSequenceTables`, then one chimera removal; carry run as a batch covariate into DA.

### Merge cliff from over-truncation (stage 2)
**Trigger:** truncLen_F + truncLen_R below amplicon length + ~12. **Mechanism:** denoised pairs no longer overlap enough to merge. **Symptom:** near-zero `merged` column, misread as "low diversity". **Fix:** compute the overlap budget first; keep length and loosen `maxEE` on the reverse for long amplicons.

### De novo tree on short reads (stage 5)
**Trigger:** UniFrac/Faith PD on a MAFFT+FastTree tree from ~250 bp reads. **Mechanism:** short reads give an unstable topology and arbitrary midpoint root. **Symptom:** unweighted-UniFrac separation that vanishes under SEPP. **Fix:** use SEPP fragment-insertion or Greengenes2 placement (microbiome/diversity-analysis).

### Sampling-depth sample massacre (stage 6)
**Trigger:** a rarefaction depth above some samples' totals (or `min(sample_sums)`). **Mechanism:** samples below the depth are silently dropped; the lost ones skew low-biomass. **Symptom:** fewer points in the PCoA than samples in the metadata. **Fix:** pick the depth from the alpha-rarefaction plateau, report the dropped-sample list, confirm at a nearby depth.

### Rarefied table reused for DA (stage 6 -> 7)
**Trigger:** feeding `ps_rare` into the DA tools. **Mechanism:** rarefaction discards count information the compositional model needs. **Symptom:** underpowered or distorted DA. **Fix:** keep raw counts; rarefy only into the diversity branch; run DA on the unrarefied `ps`.

### Single-tool DA hit list (stage 7)
**Trigger:** reporting only ALDEx2 (or only the tool that flagged the favored taxon). **Mechanism:** the significant-taxa list depends more on the tool than the biology (Nearing 2022). **Symptom:** "the method found X" with no mention of disagreeing tools. **Fix:** run >=2 CoDA tools, report the intersection as confident and the union as exploratory, name every tool.

### PICRUSt2 reported as activity (stage 8)
**Trigger:** "increased butyrate production" / "upregulated" from predicted KOs. **Mechanism:** PICRUSt2 measured no genes and no transcripts - it re-encodes taxonomy through a fixed lookup. **Symptom:** an activity verb on a predicted pathway, or predicted+taxonomic differences claimed as two independent findings. **Fix:** restrict claims to "potential"; report NSTI; for activity use metatranscriptomics (microbiome/functional-prediction).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| truncLen budget: truncLen_F + truncLen_R >= amplicon_len + ~12 | DADA2 `mergePairs` `minOverlap` default | below this, denoised pairs cannot merge; the merge cliff is arithmetic, not data |
| `maxEE` c(2,2) (loosen R for long amplicons) | Callahan 2016 *Nat Methods* 13:581 | expected-errors filter beats a hard Q cut; computed on the truncated read |
| assignTaxonomy `minBoot` 50 (default) | Wang 2007 *Appl Environ Microbiol* 73:5261 | RDP floor for reads <=250 nt; ranks below it return NA, not a guess |
| Sampling depth on the alpha-rarefaction plateau (not `min(sample_sums)`) | McMurdie 2014 *PLoS Comput Biol* 10:e1003531 | depth must saturate richness while retaining samples; report the dropped list |
| PERMANOVA permutations >= 999, paired with betadisper | vegan docs; Anderson & Walsh 2013 *Ecol Monogr* 83:557 | resolution floor; separates a location shift from a dispersion difference |
| Rarefy for diversity, NOT for DA | McMurdie 2014; Schloss 2024 *mSphere* 9:e00354-23 | per-analysis decision; DA needs the raw counts |
| Prevalence cut 10-25% before DA | Nearing 2022 *Nat Commun* 13:342; tool defaults | rare features crush the BH denominator; declare and test sensitivity |
| ALDEx2 `we.eBH` <= 0.05 AND `|effect|` > 1 | Gloor 2016 *J Comput Graph Stat* 25:971 | gate on the standardized effect (|effect|>1 = between-group diff exceeds within-condition dispersion), not p alone; large n makes trivial diffs "significant" |
| Consensus of >=2 CoDA tools | Nearing 2022 *Nat Commun* 13:342 | tool choice drives the hit list more than biology; intersection = confident |
| PICRUSt2 `--max_nsti` 2.0 (default) | Douglas 2020 *Nat Biotechnol* 38:685 | ASVs >2 substitutions/site from a reference genome are too extrapolated; report the dropped read fraction |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Near-zero merge rate | truncLen below the overlap budget | recompute the budget; keep length, loosen `maxEE` R |
| Large read fraction "chimeric" | primers not trimmed | cutadapt `--discard-untrimmed` before filtering |
| ASVs vanish/appear when runs split | one error model fit across runs | per-run `learnErrors`, then `mergeSequenceTables` |
| PCoA has fewer points than samples | rarefaction depth dropped low-count samples | lower the depth or report the loss; never assume zero drops |
| adonis2 p<0.05 but groups overlap | dispersion difference, not location | run betadisper/permutest; report both |
| ALDEx2 returns NA effects / errors | proportions or non-integer matrix passed | feed integer COUNTS with taxa in rows |
| Far fewer ANCOM-BC2 hits than expected | `p_adj_method` left at `holm` | set `p_adj_method = 'BH'` deliberately if FDR is wanted |
| Tools disagree on the DA hit list | normal - tool choice drives results | report the consensus and the disagreement, do not cherry-pick |
| PICRUSt2 result with no NSTI numbers | NSTI distribution not reported | summarize `metadata_NSTI`; report ASV+read fraction dropped at NSTI>2 |

## References

- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: high-resolution sample inference from Illumina amplicon data. *Nat Methods* 13:581-583.
- Martin M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet J* 17:10-12.
- Wang Q, Garrity GM, Tiedje JM, Cole JR. 2007. Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Appl Environ Microbiol* 73:5261-5267.
- Janssen S, McDonald D, Gonzalez A, et al. 2018. Phylogenetic placement of exact amplicon sequences improves associations with clinical information. *mSystems* 3:e00021-18.
- McDonald D, Jiang Y, Balaban M, et al. 2024. Greengenes2 unifies microbial data in a single reference tree. *Nat Biotechnol* 42:715-718.
- McMurdie PJ, Holmes S. 2014. Waste not, want not: why rarefying microbiome data is inadmissible. *PLoS Comput Biol* 10:e1003531.
- Schloss PD. 2024. Rarefaction is currently the best approach to control for uneven sequencing effort in amplicon sequence analyses. *mSphere* 9:e00354-23.
- Anderson MJ, Walsh DCI. 2013. PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: what null hypothesis are you testing? *Ecol Monogr* 83:557-574.
- Fernandes AD, Reid JNS, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB. 2014. Unifying the analysis of high-throughput sequencing datasets by compositional data analysis. *Microbiome* 2:15.
- Gloor GB, Macklaim JM, Fernandes AD. 2016. Displaying Variation in Large Datasets: Plotting a Visual Summary of Effect Sizes. *J Comput Graph Stat* 25:971-979.
- Lin H, Peddada SD. 2020. Analysis of compositions of microbiomes with bias correction. *Nat Commun* 11:3514.
- Nearing JT, Douglas GM, Hayes MG, et al. 2022. Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun* 13:342.
- Douglas GM, Maffei VJ, Zaneveld JR, et al. 2020. PICRUSt2 for prediction of metagenome functions. *Nat Biotechnol* 38:685-688.

## Related Skills

- microbiome/amplicon-processing - Primer removal, per-run error model, truncLen budget, chimeras, ITS
- reporting/automated-qc-reports - Aggregate FastQC/MultiQC across samples (sample-name resolution; the report is a snapshot, not a gate)
- microbiome/taxonomy-assignment - Region-matched classifier and reference-database choice
- microbiome/diversity-analysis - Sampling depth, tree choice, metric choice, adonis2 + betadisper
- microbiome/differential-abundance - Compositional DA tools and the consensus deliverable
- microbiome/functional-prediction - PICRUSt2 predicted potential gated on NSTI
- microbiome/qiime2-workflow - The QIIME2 artifact/provenance route for the whole chain
- read-qc/adapter-trimming - cutadapt primer removal mechanics before DADA2
- metagenomics/kraken-classification - Shotgun (not amplicon) read classification
- metagenomics/abundance-estimation - Shared compositional/normalization/rarefaction theory
- workflows/metagenomics-pipeline - The shotgun (WGS) equivalent of this amplicon pipeline
