---
name: bio-workflows-crispr-screen-pipeline
description: End-to-end pooled and single-cell CRISPR screen analysis from FASTQ to hit genes. Orchestrates library design QC, guide counting, six-stage screen QC (plasmid Gini, replicate Pearson, CEGv2 PR-AUC, copy-number artifact), method-appropriate hit calling across MAGeCK RRA/MLE, BAGEL2, drugZ, JACKS, and Chronos, cancer-cell-line copy-number correction (CRISPRcleanR / Chronos), batch correction for multi-batch screens, and the specialized branches for combinatorial paralog screens, single-cell Perturb-seq, base-editor variant-function screens, prime-editor screens, and in vivo bottleneck-aware screens. Use when analyzing any pooled CRISPR screen end-to-end, matching the hit-calling method to the experimental design, integrating copy-number correction into the pipeline, or branching the workflow for single-cell, combinatorial, base-editor, prime-editor, or in vivo variants.
tool_type: mixed
primary_tool: MAGeCK
workflow: true
depends_on:
  - crispr-screens/library-design
  - crispr-screens/screen-qc
  - crispr-screens/mageck-analysis
  - crispr-screens/bagel-essentiality
  - crispr-screens/drugz-chemogenomic
  - crispr-screens/jacks-analysis
  - crispr-screens/hit-calling
  - crispr-screens/copy-number-correction
  - crispr-screens/batch-correction
  - crispr-screens/crispresso-editing
  - crispr-screens/base-editing-analysis
  - crispr-screens/prime-editing-screens
  - crispr-screens/perturb-seq-analysis
  - crispr-screens/combinatorial-screens
  - crispr-screens/in-vivo-screens
qc_checkpoints:
  - after_counting: ">65% mapping rate; <0.5% zero-count in plasmid; Gini <0.1 on plasmid"
  - after_qc: "Replicate Pearson on log-counts >0.85; Spearman >0.7; CEGv2 PR-AUC >0.7"
  - after_cn_correction: "Spearman ρ between CN and gene LFC abs <0.10 post-correction (literature 'significant bias' band; <0.05 is a stricter target). Requires a matched CN profile, which the unsupervised CRISPRcleanR path never loads -- compute in crispr-screens/copy-number-correction, or use Chronos, which takes CN as input"
  - after_hit_calling: "Tier-1 hits = 3-method consensus; Tier-2 = 2 of 3; Tier-3 = single-method exploratory"
---

## Version Compatibility

Reference examples tested with: MAGeCK 0.5.9+, BAGEL2 1.0.5+, drugZ Aug 2019+, JACKS 0.2.0+, Chronos 2.0+, CRISPRcleanR 3.0+ (R), Pertpy 0.6+, PRIDICT2, CRISPResso2 2.2.14+, MAGeCKFlute 2.0+, pandas 2.2+, numpy 1.26+, matplotlib 3.8+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mageck --version`, `BAGEL.py fc --help`, `drugz -h`, `CRISPResso --version`
- Python: `pip show pertpy scanpy anndata` (mageck-vispr via conda `mageck --version`; JACKS/Chronos are GitHub installs — check their repos)
- R: `packageVersion('CRISPRcleanR')`, `packageVersion('MAGeCKFlute')`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

## CRISPR Screen Pipeline

**"Analyze my pooled or single-cell CRISPR screen end-to-end"** -> Pick the screen design branch, run guide counting, audit six QC stages, apply copy-number and batch correction as needed, run the design-matched hit-calling method, and consolidate across methods for high-confidence hits.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

Every LFC, QC gate, and hit call is computed against a reference that is committed once at library-order time; a wrong-but-silent commitment invalidates the endpoint with no error thrown.

1. **The guide LIBRARY definition (guide->gene map + control classes) is the denominator, the calibrator, and the training reference — committed once.** The library must carry non-targeting controls (NTCs, ~1%, the null distribution) AND CEGv2 reference essentials + NEGv1 non-essentials (the positive/negative calibrators for PR-AUC and BAGEL2/Chronos priors). NTCs calibrate the null/FDR; CEGv2/NEGv1 calibrate PR-AUC — swapping or dropping a class silently breaks FDR or QC.
2. **The baseline choice has a right answer and rescales every hit.** Dropout/enrichment LFC is against a baseline: plasmid pool for the cloning-bottleneck baseline, Day-0/T0 for the biology baseline, and vehicle (NOT Day-0) for drug screens — drug-vs-Day-0 conflates drug effect with normal proliferation.
3. **Copy-number correction MUST precede hit calling in cancer cell lines.** Multiple simultaneous Cas9 cuts at amplified loci trigger a gene-independent DNA-damage/G2 arrest (Aguirre 2016; Munoz 2016 shows the amplicon effect is p53-INDEPENDENT -- separately, Ihry 2018 / Haapaniemi 2018 report p53-dependent toxicity of Cas9 cutting generally), so amplified regions look essential regardless of gene function; calling hits first yields false essentials at ERBB2/MYC/FGFR1. Run CRISPRcleanR/Chronos BEFORE hit calling, or use CRISPRi to bypass the DSB. This is a pipeline step, not a post-hoc interpretation. Verifying `abs(rho(LFC,CN)) < 0.1` afterwards needs a matched CN profile: CRISPRcleanR corrects unsupervised without one, so the check happens in crispr-screens/copy-number-correction (or use Chronos, which consumes CN directly).
4. **CEGv2 essential-gene depletion is the screen's built-in positive control.** If known essentials do not deplete (CEGv2 PR-AUC below ~0.7), the screen failed selection and NO novel hit is trustworthy regardless of its p-value — the seam analog of a spike-in. Normalize -> QC -> (CN correct) -> hit-call, never hit-call first; add batch as an MLE covariate, never pre-corrected with ComBat on counts (distorts the NB mean-variance the caller assumes).

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Guide library (guide->gene map + NTC/CEGv2/NEGv1 control classes) | The counting denominator, the QC calibrator, the hit-calling priors; a missing/misassigned class breaks FDR or PR-AUC |
| Baseline (plasmid pool / Day-0 / vehicle) | Every LFC; drug-vs-Day-0 conflates drug effect with proliferation |
| Screen type (dropout / enrichment / FACS / drug-modifier) | Which hit-calling method is even valid |
| Copy-number profile (cancer lines) | Whether amplicon artifacts are removed before hit calling; residual rho(LFC,CN) is the tell |

## Pipeline Branches by Screen Design

```
                    Library Design ([[library-design]])
                              |
                              v
                FASTQ Files -> mageck count -> count matrix
                              |
                              v
                Six-Stage QC ([[screen-qc]])
                              |
        +---------------------+---------------------+
        |                                            |
        v                                            v
  Cancer cell line?                          Non-cancer?
  Apply CN correction                        No CN correction needed
  ([[copy-number-correction]])
        |                                            |
        +---------------------+---------------------+
                              v
                  Multi-batch? Apply batch covariate
                  ([[batch-correction]])
                              |
                              v
                 Pick hit-calling method by design ([[hit-calling]])
                              |
        +-----------+---------+---------+-----------+-----------+
        |           |         |         |           |           |
        v           v         v         v           v           v
    2-cond       Time      Drug      Essential   Multi-       Specialized
    MAGeCK RRA   MAGeCK    drugZ     BAGEL2      screen       (PE/BE/SC/
                 MLE                              JACKS or     in vivo/
                                                  Chronos      combinat)
        |           |         |         |           |           |
        +-----------+---------+---------+-----------+-----------+
                              v
                   Tier-based consensus
                              v
                Orthogonal validation
```

## Step 1: Library Design and Pre-Screen Validation

Reference [[library-design]] for full library composition. Verify before sequencing:
- Plasmid pool Gini <0.1 (Joung 2017 Nat Protoc 12:828)
- >=99% guides detected at >25 reads/guide
- Skew (p90/p10) <2
- NTCs comprise ~1% of library; CEGv2 reference essentials + NEGv1 non-essentials included

## Step 2: Guide Counting

**Goal:** Turn raw FASTQ into a per-guide count matrix with consistent sample labels.

**Approach:** Run mageck count with the library CSV, sample labels in column order, the vector adapter trimmed off the 5' end, and median normalization.

```bash
mageck count \
    --list-seq library.csv \
    --sample-label Plasmid,Day0,Veh_r1,Veh_r2,Drug_r1,Drug_r2 \
    --fastq Plasmid.fq.gz Day0.fq.gz Veh_r1.fq.gz Veh_r2.fq.gz Drug_r1.fq.gz Drug_r2.fq.gz \
    --norm-method median \
    --output-prefix experiment \
    --trim-5 5   # integer base-count (or AUTO), NOT an adapter sequence; 5 trims the CACCG scaffold
```

For Cas12a libraries (Inzolia, in4mer): see [[combinatorial-screens]]. For 10X single-cell direct capture: use cellranger-arc or pertpy-aware counting; see [[perturb-seq-analysis]].

## Step 3: Six-Stage Quality Control

**Goal:** Decide whether the screen is analyzable before calling any hits, using six orthogonal QC stages.

**Approach:** Load the count matrix, compute per-sample Gini, zero-fraction, and depth plus replicate correlation against the hard gates below. Essential-gene recovery (CEGv2 PR-AUC) is a separate check computed once endpoint-vs-baseline LFCs exist (it needs CEGv2/NEGv1 labels) -- see screen-qc.

```python
import pandas as pd
import numpy as np

counts = pd.read_csv('experiment.count.txt', sep='\t', index_col=0)
genes = counts['Gene']
count_matrix = counts.drop('Gene', axis=1)

def gini(x):
    x = np.sort(x[x > 0].astype(float))
    if x.size == 0:
        return np.nan
    n = x.size
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

per_sample = pd.DataFrame({
    'pct_zero': (count_matrix == 0).sum() / len(count_matrix) * 100,
    'gini': count_matrix.apply(gini),
    'reads_per_sgrna': count_matrix.sum() / len(count_matrix),
})

log_counts = np.log10(count_matrix + 1)
pearson = log_counts.corr()
print(per_sample)
print('Replicate Pearson:', pearson.values[pearson.values < 1].mean())
```

Hard gates from [[screen-qc]]:
- Plasmid Gini <0.1; endpoint <0.3 (or <0.55 for heavy drug screens)
- Replicate Pearson on log-counts >0.85
- CEGv2 PR-AUC >0.7 against Hart 2017 reference essential gene set
- Reads per sgRNA per sample >=300 (DepMap convention)

## Step 4: Copy-Number Correction (Cancer Cell Lines Only)

If screening in a cancer cell line, apply CRISPRcleanR (unsupervised, no CN profile needed) or Chronos (joint with CN profile). Required to remove Aguirre 2016 / Munoz 2016 amplicon artifact.

**Goal:** Strip the copy-number amplicon artifact that makes amplified regions look essential in cancer lines.

**Approach:** Run CRISPRcleanR unsupervised genome-wide LFC correction (no CN profile needed), then feed the corrected counts downstream; for DepMap-scale panels with matched CN, use Chronos instead.

```r
library(CRISPRcleanR)
data(KY_Library_v1.0)
norm <- ccr.NormfoldChanges('experiment.count.txt', min_reads = 30, EXPname = 'screen',
                              libraryAnnotation = KY_Library_v1.0)   # arg 1 is the file PATH
gw_lfc <- ccr.logFCs2chromPos(norm$logFCs, KY_Library_v1.0)          # $logFCs, not $norm_fold_changes
cleaned <- ccr.GWclean(gw_lfc, display = TRUE, label = 'screen')
corrected_counts <- ccr.correctCounts('screen', norm$norm_counts, cleaned,
                                        KY_Library_v1.0)              # (CL, normalised_counts, correctedFCs, libraryAnnotation)
# ccr.correctCounts returns an in-memory frame; it does NOT write this file. Persist it, because the
# hit callers below read a count TABLE from disk -- the CN-correction commitment in rule 3 is only
# honored if that file, not the raw experiment.count.txt, is what MAGeCK / BAGEL2 / drugZ consume.
write.table(corrected_counts, 'screen_cleanr_corrected_counts.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)
```

For DepMap-scale panels with longitudinal data + matched CN, use Chronos. See [[copy-number-correction]].

## Step 5: Batch Correction (Multi-Batch Screens)

For multi-batch screens, add batch as a covariate in MAGeCK MLE rather than pre-correcting with ComBat. See [[batch-correction]] for full decision tree.

## Step 6: Method-Matched Hit Calling

**Goal:** Call hits with the method that matches the experimental design, plus at least one orthogonal method for consensus.

**Approach:** Pick by design - RRA or BAGEL2 for two-condition essentiality, MLE for time course, drugZ for drug-modifier, JACKS for multi-screen, Chronos for cancer panels - and run two methods so the consensus step has something to reconcile.

### 6a. Two-condition essentiality (MAGeCK RRA or BAGEL2)

Cancer cell lines: pass the CRISPRcleanR-corrected count file (`screen_cleanr_corrected_counts.txt` from Step 4) as `--count-table`/`-i` below, NOT the raw `experiment.count.txt` — the CN-correction commitment is only honored if the corrected counts are what the hit caller reads.

The `Day0`/`Day14_r*` columns below illustrate a time-course dropout design; they must match the count step's `--sample-label` (the drug-screen count above uses `Plasmid,Day0,Veh_r*,Drug_r*`).

```bash
mageck test \
    --count-table experiment.count.txt \
    --treatment-id Day14_r1,Day14_r2,Day14_r3 \
    --control-id Day0 \
    --norm-method median \
    --output-prefix essentiality_rra
```

```bash
BAGEL.py fc -i experiment.count.txt -o experiment -c Day0 --min-reads 30   # -o is an output LABEL; fc writes experiment.foldchange
BAGEL.py bf -i experiment.foldchange -o bayes_factor.txt -e CEGv2.txt -n NEGv1.txt \
    -c Day14_r1,Day14_r2,Day14_r3   # add -b -NB 1000 for bootstrapping; -k is not a bf option
```

### 6b. Time-course / multi-condition (MAGeCK MLE)

```bash
mageck mle --count-table experiment.count.txt --design-matrix design.txt \
    --output-prefix timecourse_mle --norm-method median
```

### 6c. Drug-modifier (drugZ)

```bash
python drugz.py \
    -i experiment.count.txt \
    -o drugz_output.txt \
    -c Veh_r1,Veh_r2 \
    -x Drug_r1,Drug_r2 \
    -p 5
```

drugZ requires vehicle as control, not Day-0. See [[drugz-chemogenomic]].

### 6d. Multi-screen joint analysis (JACKS)

```bash
python run_JACKS.py experiment.count.txt replicatemap.txt guidemap.txt \
    --rep_hdr Replicate --sample_hdr Sample --ctrl_sample_hdr Control \
    --sgrna_hdr sgRNA --gene_hdr Gene --outprefix jacks_out --apply_w_hp
```

### 6e. Cancer cell-line panels (Chronos)

```python
import chronos
# All three inputs must be dicts of DataFrame keyed by library name, not bare DataFrames.
model = chronos.Chronos(sequence_map={'screen': sequence_map},
                          guide_gene_map={'screen': guide_gene_map},
                          readcounts={'screen': counts_df})   # readcounts=, not reads=
model.train(nepochs=301)                                 # nepochs (default 301), not n_steps
gene_effects = model.gene_effect                         # attribute, not a method call
# Copy-number correction is a separate post-hoc step (chronos.alternate_CN(gene_effect, copy_number) / a CN matrix), not a constructor arg
```

DepMap quarterly standard; handles CN bias + screen quality + longitudinal jointly.

## Step 7: Tier-Based Consensus

**Goal:** Consolidate the per-method calls into confidence tiers.

**Approach:** Merge each method's gene-level result, threshold each to a per-method hit flag, and tier by how many methods agree (Tier 1 = all three, Tier 2 = two of three).

```python
mageck = pd.read_csv('essentiality_rra.gene_summary.txt', sep='\t')[['id', 'neg|fdr']].rename(
    columns={'id': 'gene', 'neg|fdr': 'mageck_neg_fdr'})
bagel = pd.read_csv('bayes_factor.txt', sep='\t')[['GENE', 'BF']].rename(
    columns={'GENE': 'gene', 'BF': 'bagel_bf'})
drugz_df = pd.read_csv('drugz_output.txt', sep='\t')[['GENE', 'fdr_synth']].rename(
    columns={'GENE': 'gene', 'fdr_synth': 'drugz_synth_fdr'})

merged = mageck.merge(bagel, on='gene', how='outer').merge(drugz_df, on='gene', how='outer')
merged['mageck_hit'] = merged['mageck_neg_fdr'] < 0.05
merged['bagel_hit'] = merged['bagel_bf'] > 6
merged['drugz_hit'] = merged['drugz_synth_fdr'] < 0.05
merged['tier'] = merged[['mageck_hit', 'bagel_hit', 'drugz_hit']].astype(int).sum(axis=1)
tier1 = merged[merged['tier'] >= 3]
tier2 = merged[merged['tier'] == 2]
merged.to_csv('tier_consensus.csv', index=False)   # the documented deliverable; the frame above is otherwise in-memory only
```

## Specialized Branches

| Screen design | Specialized workflow |
|---------------|----------------------|
| Single-cell Perturb-seq / CROP-seq / Multiome | [[perturb-seq-analysis]] -- Pertpy + Mixscape + SCEPTRE |
| Combinatorial paralog (Cas12a Inzolia / Big Papi) | [[combinatorial-screens]] -- GI scoring; synthetic-lethal identification |
| Base-editor variant-function (Hanna 2021 style) | [[base-editing-analysis]] + [[crispresso-editing]] |
| Prime-editor variant installation | [[prime-editing-screens]] -- PRIDICT2 pegRNA design |
| In vivo tumor / immune screens | [[in-vivo-screens]] -- focused library; per-animal meta-analysis |

## Visualization

**Goal:** Show the hit landscape as a volcano of effect size against significance.

**Approach:** Plot log2 fold change against -log10(FDR), highlight genes past the FDR gate, and save the figure to file.

```python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 8))
gene_summary = pd.read_csv('essentiality_rra.gene_summary.txt', sep='\t')
sig = gene_summary['neg|fdr'] < 0.05
ax.scatter(gene_summary.loc[~sig, 'neg|lfc'],
            -np.log10(gene_summary.loc[~sig, 'neg|fdr'].clip(lower=1e-10)),
            c='lightgray', alpha=0.5, s=10)
ax.scatter(gene_summary.loc[sig, 'neg|lfc'],
            -np.log10(gene_summary.loc[sig, 'neg|fdr'].clip(lower=1e-10)),
            c='red', alpha=0.7, s=18)
ax.axhline(-np.log10(0.05), ls='--', c='black', lw=0.5)
ax.set_xlabel('Log2 Fold Change')
ax.set_ylabel('-Log10(FDR)')
plt.savefig('volcano.png', dpi=150)
```

MAGeCKFlute R package provides one-shot FluteRRA / FluteMLE dashboards with KEGG/Reactome enrichment.

## Output Files

| File | Source step | Description |
|------|-------------|-------------|
| experiment.count.txt | mageck count | Raw count matrix |
| experiment.countsummary.txt | mageck count | Per-sample Gini, mapping, % zero |
| screen_cleanr_corrected_counts.txt | CRISPRcleanR | CN-corrected counts (cancer lines) |
| essentiality_rra.gene_summary.txt | mageck test | Gene-level RRA scores |
| bayes_factor.txt | BAGEL2 | Per-gene Bayes factors |
| drugz_output.txt | drugZ | sumZ, normZ, per-direction FDR |
| jacks_out_gene_JACKS_results.txt | JACKS | Gene effect + sgRNA efficacy |
| tier_consensus.csv | Custom aggregation | Tier-1/2/3 hits across methods |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| False essentials at ERBB2/MYC/FGFR1 | Hit calling before copy-number correction (gene-independent, p53-independent DNA-damage arrest at amplicons) | Run CRISPRcleanR/Chronos BEFORE hit calling; verify abs(rho(LFC,CN)) < 0.1; or use CRISPRi to bypass the DSB |
| FDR broken or PR-AUC uncomputable | NTC (null) and CEGv2 essential (positive control) classes swapped or one absent | Keep both classes; NTCs calibrate the null/FDR, CEGv2/NEGv1 calibrate PR-AUC and BAGEL2/Chronos priors |
| Every hit rescaled / drug effect confounded | Wrong baseline (Day-0 for a drug screen) | Drug screen -> vehicle control; plasmid pool for the cloning-bottleneck baseline |
| "Everything significant at FDR<0.01" | Heavy selection breaks median normalization (>40% guides change) | Switch to `--norm-method control` on NTCs, or BAGEL2 |
| Underpowered / method mismatch | RRA on a time course; single-line Chronos | Pick method by design (fork table); RRA fails multi-condition, Chronos is overkill single-line |
| Distorted NB mean-variance | Batch pre-corrected with ComBat on counts | Add batch as a MAGeCK MLE covariate instead |
| Novel hits from a failed screen | CEGv2 essentials did not deplete (PR-AUC < 0.7) | The screen failed selection; no hit is trustworthy regardless of p-value |

## References

- Li W, Xu H, Xiao T, et al (2014) MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. *Genome Biology* 15:554. DOI 10.1186/s13059-014-0554-4.
- Aguirre AJ, Meyers RM, Weir BA, et al (2016) Genomic copy number dictates a gene-independent cell response to CRISPR/Cas9 targeting. *Cancer Discovery* 6:914-929. DOI 10.1158/2159-8290.CD-16-0154. (the amplicon artifact.)
- Hart T, Moffat J (2016) BAGEL: a computational framework for identifying essential genes from pooled library screens. *BMC Bioinformatics* 17:164. DOI 10.1186/s12859-016-1015-8.
- Iorio F, Behan FM, Goncalves E, et al (2018) Unsupervised correction of gene-independent cell responses to CRISPR-Cas9 targeting (CRISPRcleanR). *BMC Genomics* 19:604. DOI 10.1186/s12864-018-4989-y.
- Joung J, Konermann S, Gootenberg JS, et al (2017) Genome-scale CRISPR-Cas9 knockout and transcriptional activation screening. *Nature Protocols* 12:828-863. DOI 10.1038/nprot.2017.016. (library QC / Gini.)

## Related Skills

- crispr-screens/library-design - Library composition and design rules
- crispr-screens/screen-qc - Six-stage QC + CEGv2 PR-AUC
- crispr-screens/mageck-analysis - MAGeCK RRA + MLE detail
- crispr-screens/bagel-essentiality - BAGEL2 Bayes factor essentiality
- crispr-screens/drugz-chemogenomic - drugZ for drug-modifier screens
- crispr-screens/jacks-analysis - Joint multi-screen analysis with shared efficacy
- crispr-screens/hit-calling - Cross-method decision tree + reconciliation
- crispr-screens/copy-number-correction - CRISPRcleanR / CERES / Chronos
- crispr-screens/batch-correction - Multi-batch design matrix
- crispr-screens/crispresso-editing - CRISPResso2 editing quantification
- crispr-screens/base-editing-analysis - Variant-function BE screens
- crispr-screens/prime-editing-screens - PRIDICT2 pegRNA design
- crispr-screens/perturb-seq-analysis - Single-cell screen analysis
- crispr-screens/combinatorial-screens - Cas12a multiplex + GI scoring
- crispr-screens/in-vivo-screens - Bottleneck-aware in vivo design
- pathway-analysis/go-enrichment - Functional enrichment of hits
- pathway-analysis/gsea - Pre-ranked GSEA on hit lists
