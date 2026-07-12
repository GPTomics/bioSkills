---
name: bio-isoform-switching
description: Analyzes differential transcript usage (DTU) and isoform switches with functional consequence prediction (NMD via 50nt rule, ORF disruption, protein domain loss/gain, signal peptide changes, IDR alterations, coding-potential shifts). Tools include IsoformSwitchAnalyzeR v2 (auto-selects satuRn for >5 reps else DEXSeq), the manual DRIMSeq -> DEXSeq/satuRn -> stageR DTU pipeline, and fishpond/swish for inferential-uncertainty-aware DTE. Distinguishes DTU from DGE and DTE; integrates external annotators (CPC2, Pfam, SignalP, IUPred2A or DeepTMHMM). Use when investigating how splicing differences alter protein function or trigger NMD-mediated degradation.
tool_type: r
primary_tool: IsoformSwitchAnalyzeR
---

## Version Compatibility

Reference examples tested with: IsoformSwitchAnalyzeR 2.11+, DRIMSeq 1.34+, DEXSeq 1.52+, satuRn 1.14+, stageR 1.28+, fishpond 2.14+, tximport 1.34+, tximeta 1.24+, Salmon 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Isoform Switching and Differential Transcript Usage

Identify shifts in *which* transcript a gene predominantly uses between conditions, and predict functional consequences. Statistically distinct from DGE and DTE; biologically distinct because the same gene-level expression can hide a complete isoform switch with major protein-level consequences.

## DGE vs DTE vs DTU: Which Question Is Being Asked?

| Question | Statistic | Tool | Example claim |
|----------|-----------|------|----------------|
| **DGE** Does the gene total change? | Sum of transcript counts | DESeq2, edgeR, limma-voom | "Gene X is upregulated 2-fold" |
| **DTE** Does this transcript change in absolute abundance? | Per-transcript count | swish (fishpond), DESeq2 on transcripts, sleuth | "Transcript X-201 is upregulated 2-fold" |
| **DTU** Do proportions of transcripts within the gene shift? | Vector of per-transcript proportions | DRIMSeq, DEXSeq, satuRn (+ stageR) | "Gene X switches from isoform 201 (50% -> 10%) to 202 (50% -> 90%)" |

DTU is statistically harder than DGE because:
1. The null is **compositional** (proportions sum to 1; one transcript up means another down).
2. **Multi-stage testing** is required: gene-level "any DTU" + transcript-level "which transcript" -> stageR formalizes this.
3. **Quantification uncertainty propagates** when transcripts are similar (Salmon EM ambiguity).

DTU and event-level differential splicing answer related but distinct questions: rMATS' `IncLevelDifference` is essentially a 1-D projection of a DTU shift onto a single event coordinate. The pragmatic 2026 default: run both an event-level tool (rMATS or leafcutter) and a DTU pipeline; reconcile.

## Tool Selection for DTU

| Tool | Model | When to use | Fails when |
|------|-------|-------------|------------|
| IsoformSwitchAnalyzeR v2 | Wraps DEXSeq or satuRn + functional consequence annotation | Standard interpretation workflow with NMD/domain output | Manual DTU control needed; very large cohorts (>200) |
| DRIMSeq | Dirichlet-multinomial on transcript counts; gene-level DTU | Pre-filter step before DEXSeq/satuRn | Cannot annotate functional consequences alone |
| DEXSeq | Negative-binomial GLM on exon-bin or transcript counts | Classic DTU; conservative; <=5 replicates per condition | Slow at scale; uses bins not transcripts in default mode |
| satuRn | Quasi-binomial GLM with empirical-Bayes shrinkage | DTU at scale (single-cell, large bulk cohorts) | Newer; less battle-tested than DEXSeq |
| swish (fishpond) | Non-parametric SAMseq across Salmon Gibbs samples | DTE/DGE incorporating quantification uncertainty | Requires Gibbs samples; not strictly DTU |
| stageR | Two-stage testing framework | Required for proper OFDR control on top of DRIMSeq/DEXSeq/satuRn | Standalone — wraps another tool's output |
| sleuth | Bootstrap-based DTE on kallisto | When committed to kallisto pipeline | Less active development; superseded by fishpond+swish |

The IsoformSwitchAnalyzeR v2 default rule is: **satuRn if any condition has >5 replicates; else DEXSeq.** For exactly 5 replicates per condition (boundary), explicitly choose; results may differ.

## Decision Tree by Research Question

| Question | Recommended approach |
|----------|----------------------|
| Functional consequences of switches (domains, NMD, signal peptide) | IsoformSwitchAnalyzeR v2 with full external annotator pipeline |
| Pure statistical DTU (gene-level + transcript-level OFDR) | DRIMSeq (filter) -> DEXSeq -> stageR; or -> satuRn -> stageR for n>5 |
| DTU with proper quantification uncertainty | Salmon `--numGibbsSamples 20` -> tximeta -> swish for DTE; concurrent DTU |
| Single-cell DTU | satuRn (DEXSeq doesn't scale to scRNA-seq) |
| Long-read DTU (PacBio Iso-Seq, ONT) | IsoformSwitchAnalyzeR v2 long-read input mode (no Salmon EM uncertainty) |
| Time-course DTU | DEXSeq with time as factor + interaction; or limma::lmFit on logit-prop matrix |
| Cancer / disease — switch hits -> mechanism | Standard pipeline + cross-reference with eCLIP, ClinVar, COSMIC |
| Therapeutic ASO target identification | Standard pipeline + sashimi visualization + SpliceAI design |

## IsoformSwitchAnalyzeR v2 Workflow

**Goal:** Identify isoform switches with functional consequences in one integrated workflow.

**Approach:** Import Salmon, pre-filter, run statistical test (satuRn auto-selected if any condition has >5 replicates, else DEXSeq), annotate switches with external tools (CPC2, Pfam, SignalP, IUPred2A or DeepTMHMM), then summarize consequences.

```r
library(IsoformSwitchAnalyzeR)

salmonQuant <- importIsoformExpression(
    parentDir = 'salmon_quant/',
    addIsofomIdAsColumn = TRUE
)

design <- data.frame(
    sampleID = colnames(salmonQuant$counts)[-1],
    condition = c('control', 'control', 'control', 'treatment', 'treatment', 'treatment')
)

aSwitchList <- importRdata(
    isoformCountMatrix = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix = design,
    isoformExonAnnoation = 'annotation.gtf',
    isoformNtFasta = 'transcripts.fa',
    showProgress = TRUE
)

aSwitchList <- preFilter(
    aSwitchList,
    geneExpressionCutoff = 1,
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,
    removeSingleIsoformGenes = TRUE,
    keepIsoformInAllConditions = TRUE
)

aSwitchList <- isoformSwitchTestSatuRn(
    aSwitchList,
    reduceToSwitchingGenes = TRUE,
    alpha = 0.05,
    dIFcutoff = 0.1
)
```

`preFilter` parameters:
- `geneExpressionCutoff = 1` — minimum TPM for gene to be tested (raise for stricter)
- `isoformExpressionCutoff = 0` — minimum TPM per isoform (set to 1 for stricter)
- `IFcutoff = 0.01` — minimum isoform fraction; below = noise
- `removeSingleIsoformGenes = TRUE` — drop genes with only one detectable isoform (cannot have DTU)
- `keepIsoformInAllConditions = TRUE` — require expression across all conditions

For long-read input, use `importRdata` with long-read transcript counts directly — bypasses Salmon EM uncertainty entirely.

## Functional Consequence Annotation

**Goal:** Predict how each switch alters protein structure, function, and stability.

**Approach:** Extract sequences, run external annotators outside R, then re-import results into the switchAnalyzeRlist.

```r
aSwitchList <- extractSequence(
    aSwitchList,
    pathToOutput = 'sequences/',
    writeToFile = TRUE
)

# Run external tools on sequences/isoformSwitchAnalyzeR_isoform_*.fasta
# Then import:

# IMPORTANT: ORF analysis must run BEFORE analyzeSwitchConsequences for NMD_status,
# ORF_seq_similarity, and coding_potential consequences to be computed.
aSwitchList <- analyzeORF(aSwitchList, orfMethod = 'longest', genomeObject = NULL)

aSwitchList <- analyzeCPC2(aSwitchList, pathToCPC2resultFile = 'cpc2_results.txt', removeNoncodinORFs = TRUE)
aSwitchList <- analyzePFAM(aSwitchList, pathToPFAMresultFile = 'pfam_results.txt')
aSwitchList <- analyzeSignalP(aSwitchList, pathToSignalPresultFile = 'signalp_results.txt')
aSwitchList <- analyzeIUPred2A(aSwitchList, pathToIUPred2AresultFile = 'iupred2_results.txt')
aSwitchList <- analyzeAlternativeSplicing(aSwitchList, onlySwitchingGenes = TRUE)

aSwitchList <- analyzeSwitchConsequences(
    aSwitchList,
    consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
        'IDR_identified',
        'IDR_type',
        'signal_peptide_identified'
    ),
    dIFcutoff = 0.1
)
```

| External tool | Purpose | Required for |
|----------------|---------|--------------|
| CPC2 (Coding Potential Calculator 2) | Coding vs non-coding classification | `coding_potential` consequence |
| Pfam (HMMER hmmscan against Pfam-A) | Protein domain identification | `domains_identified` |
| SignalP 6.0+ | Signal peptide prediction | `signal_peptide_identified` |
| IUPred2A or DeepTMHMM | Intrinsic disorder regions / TM domains | `IDR_identified`, `IDR_type` |
| NetSurfP-3 | Surface accessibility (optional) | Extended IDR analysis |

The external tools must be run *outside* R; IsoformSwitchAnalyzeR provides FASTA outputs and re-imports the parsed results. Plan for ~30-60 minutes of external compute on typical mammalian transcriptomes.

## NMD Prediction (The 50-nt Rule)

A transcript is predicted NMD-sensitive if its premature termination codon (PTC) lies **>50-55 nt upstream of the last exon-exon junction** (Maquat 2004 *Nat Rev Mol Cell Biol*; Lykke-Andersen & Jensen 2015 *Nat Rev Mol Cell Biol*).

**Mechanism:** Spliceosome deposits the Exon Junction Complex (EJC) ~20-24 nt upstream of every exon-exon junction. During the pioneer round of translation, ribosome reading through removes EJCs upstream of the stop codon. If a stop codon precedes the last EJC by >50 nt, the EJC remains, recruits UPF1 -> SMG1 phosphorylation -> SMG6/SMG7 -> mRNA decay.

**Caveats and exceptions:**
- **Last-exon PTCs escape NMD** — can be dominant-negative or gain-of-function (e.g. MYH7 truncating variants).
- **3'UTR length matters**: very long 3' UTRs (>1 kb past stop) trigger NMD via UPF1 binding even without EJCs (faux-3'UTR rule).
- **Tissue-specific NMD**: SMG6 vs SMG5/7 ratios vary; UPF1 stress conditions modulate.
- **PTC distance must be measured on the spliced transcript**, not the genomic distance.
- **~10-20% of "predicted NMD" transcripts escape NMD per orthogonal RNA-seq** (Lindeboom 2016 *Nat Genet*; ~22% of canonical PTC-bearing transcripts escape in some tissues). Treat NMD prediction as probabilistic, not certain.

IsoformSwitchAnalyzeR's `analyzeSwitchConsequences` with `'NMD_status'` evaluates this from the predicted ORF + transcript model.

## AS-NMD as a Regulatory Layer

A large class of conserved alternative splicing events is **deliberately PTC-introducing** to titrate functional protein levels:

- **All major SR proteins** (SRSF1-12) autoregulate via poison exons (Lareau 2007 *Nature*; Ni 2007 *Genes Dev*)
- **All major hnRNPs** likewise
- **Ribosomal protein genes** use AS-NMD autoregulation (e.g. rpL3, rpL12; Cuccurese 2005 *NAR*)
- **SCN1A** poison exon -> Stoke STK-001 ASO in Phase 1/2 for Dravet syndrome (Han 2020 *Sci Transl Med*)

**Functional implication:** an *increase* in PSI of a poison exon *decreases* functional protein. Sign-of-effect in DTU output is opposite from intuition for these genes. Always check whether the alternative form is PTC-bearing before interpreting direction.

**Disease examples:**
- TDP-43 cryptic exons (UNC13A, STMN2) introduce PTCs -> NMD on disease-relevant transcript (Brown 2022 *Nature*)
- Last-exon truncating variants in TTN (and MYH7): escape NMD -> stable poison/dominant-negative protein

## Manual DTU Pipeline (DRIMSeq + DEXSeq + stageR)

The canonical reference is the *F1000Research* "Swimming downstream" workflow (Love, Soneson, Patro 2018; Bioconductor `rnaseqDTU`).

```r
library(tximeta); library(DRIMSeq); library(DEXSeq); library(stageR)

se <- tximeta(coldata)
counts <- assays(se)$counts

samples <- data.frame(
    sample_id = colnames(counts),
    condition = c('control', 'control', 'control', 'treatment', 'treatment', 'treatment')
)

txdf <- data.frame(
    gene_id = rowData(se)$gene_id,
    feature_id = rowData(se)$tx_id,
    counts
)

d <- dmDSdata(counts = txdf, samples = samples)
d <- dmFilter(d, min_samps_feature_expr = 3, min_feature_expr = 10,
              min_samps_feature_prop = 3, min_feature_prop = 0.1,
              min_samps_gene_expr = 6, min_gene_expr = 10)

design_full <- model.matrix(~ condition, data = samples(d))

dxd <- DEXSeqDataSet(
    countData = round(as.matrix(counts(d)[, -c(1, 2)])),
    sampleData = samples(d),
    design = ~ sample + exon + condition:exon,
    featureID = counts(d)$feature_id,
    groupID = counts(d)$gene_id
)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE)
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
qval <- perGeneQValue(DEXSeqResults(dxd))

dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
pConfirmation <- matrix(dxr$pvalue, ncol = 1)
rownames(pConfirmation) <- dxr$featureID
tx2gene <- as.data.frame(dxr[, c('featureID', 'groupID')])

stageRObj <- stageRTx(
    pScreen = qval,
    pConfirmation = pConfirmation,
    pScreenAdjusted = TRUE,
    tx2gene = tx2gene
)
stageRObj <- stageWiseAdjustment(stageRObj, method = 'dtu', alpha = 0.05)

results <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)
```

**stageR semantics:**
- **Stage 1 (screening)**: gene-level p-value (`perGeneQValue` from DEXSeq, or DRIMSeq's gene-level p) is filtered at the desired Overall FDR.
- **Stage 2 (confirmation)**: only within significant genes, individual transcripts are tested at a within-gene FWER computed to maintain global OFDR.
- **Net effect**: gene-level FDR is properly controlled, AND the transcript that drove the call is known.
- Without stageR: naive transcript-level BH overcounts because the gene-level multiple-testing burden is ignored.

## fishpond/swish for Inferential-Uncertainty-Aware Testing

**Goal:** Test DTE while propagating quantification uncertainty from Salmon's Gibbs samples.

**Approach:** Run Salmon with `--numGibbsSamples 20`, import with tximeta, then use swish to average a non-parametric SAMseq-style test across inferential replicates.

```r
library(fishpond); library(tximeta)

se <- tximeta(coldata)
y <- scaleInfReps(se)
y <- labelKeep(y)
y <- y[mcols(y)$keep, ]

set.seed(1)
y <- swish(y, x = 'condition')

dte_results <- as.data.frame(mcols(y))
sig <- subset(dte_results, qvalue < 0.05)
```

**`infRV`** (inferential relative variance) is a per-feature uncertainty diagnostic; high-infRV transcripts are unreliable and can be filtered before testing. Critical for genes with many similar isoforms (TTN, MAPT, NEFM) where Salmon's EM is uncertain.

## Per-Tool Failure Modes

### DEXSeq: Slowness at Scale

**Trigger:** Bulk cohort with >50 samples or single-cell DTU.

**Mechanism:** DEXSeq fits a NB GLM per exon-bin per gene; computational cost scales linearly with samples × bins.

**Symptom:** `estimateDispersions` takes hours; `testForDEU` exhausts memory.

**Fix:** Switch to satuRn (designed for scale, including scRNA-seq); run with parallelization (`BPPARAM = MulticoreParam(8)`).

### DRIMSeq: Filtering Sensitivity

**Trigger:** Default `dmFilter` parameters too strict for low-expression cohort.

**Mechanism:** `min_samps_feature_expr = 3, min_feature_expr = 10` drops transcripts seen in <=2 samples or with <10 counts.

**Symptom:** Most candidate genes filtered out; few testable genes.

**Fix:** Tune to dataset: lower thresholds for low-coverage data, raise for high-coverage. Document choice.

### satuRn: Empirical-Bayes Shrinkage Limits

**Trigger:** Very small cohort (n=2 vs n=2) or very heterogeneous.

**Mechanism:** Empirical-Bayes shrinkage assumes shared dispersion across genes; collapses with too-few or too-heterogeneous samples.

**Symptom:** Inflated p-values; few discoveries despite real effects.

**Fix:** Aggregate replicates (pseudobulk), or switch to DEXSeq for small cohorts; use larger cohorts when possible.

### swish: Salmon Gibbs Requirements

**Trigger:** Running swish on Salmon output without Gibbs samples.

**Mechanism:** swish averages over inferential replicates from Salmon's Gibbs sampler; requires `--numGibbsSamples 20` (or bootstrap with `--numBootstraps`) at Salmon time.

**Symptom:** `scaleInfReps` errors about missing inferential replicates.

**Fix:** Re-run Salmon with `--numGibbsSamples 20`; this triples Salmon runtime but enables uncertainty-aware testing.

### IsoformSwitchAnalyzeR: External-Annotator Failure

**Trigger:** Forgetting to run all 4 external annotators (CPC2, Pfam, SignalP, IUPred2A).

**Mechanism:** `analyzeSwitchConsequences` silently drops consequence types for which annotation wasn't imported.

**Symptom:** `extractConsequenceSummary` shows fewer types than requested; specific consequence reports missing.

**Fix:** Verify all 4 result files exist before `analyzeSwitchConsequences`; check `aSwitchList$AlternativeSplicingAnalysis` slot for completeness.

## Reconciliation: When DTU and Event-Level Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Significant DTU, no rMATS hit | DTU shift across many transcripts; no single canonical event captures it | Examine isoform structure in switchPlot; report at gene level |
| rMATS sig, no significant DTU | Single event in single isoform; not a gene-level DTU | Report as event-level result; DTU not the right framing |
| Both sig, same gene, different "main" isoforms | Annotation differs (rMATS uses GENCODE basic; ISA uses comprehensive) | Standardize annotation; re-run |
| DTU shows poison-exon switch, gene-level DGE shows decrease | NMD-coupled regulation: AS-NMD reducing protein on top of transcription | Mechanism: AS-NMD; report direction carefully |

For high-confidence reporting: concordant DTU + event-level + sashimi visualization.

## Single-Cell DTU

For scRNA-seq, **satuRn** scales where DEXSeq does not. IsoformSwitchAnalyzeR v2 supports single-cell input via `importRdata` with single-cell count matrices, and the underlying satuRn test has explicit single-cell calibration (Gilis 2022 *F1000Research*).

**Strong recommendation:** pseudobulk by cell type first; per-cell DTU is rarely powered with droplet 3' chemistry. See `single-cell-splicing` for chemistry-specific limitations.

## Visualization

```r
extractTopSwitches(
    aSwitchList,
    filterForConsequences = TRUE,
    n = 25,
    sortByQvals = TRUE
)

switchPlot(
    aSwitchList,
    gene = 'TARGET_GENE',
    condition1 = 'control',
    condition2 = 'treatment',
    localTheme = theme_bw(base_size = 12)
)

extractConsequenceSummary(aSwitchList, consequencesToAnalyze = 'all', plotGenes = FALSE)
extractConsequenceEnrichment(aSwitchList, consequencesToAnalyze = 'all')
extractSplicingSummary(aSwitchList, asFractionTotal = FALSE)
```

## Significance Thresholds

| Parameter | Default | Notes |
|-----------|---------|-------|
| isoform_switch_q_value | < 0.05 | Switch significance |
| dIF (delta isoform fraction) | > 0.1 | Minimum biological effect |
| Consequence q-value | < 0.05 | Significance per consequence type |
| Gene-level OFDR (stageR) | < 0.05 | Gene-level screening FDR |
| satuRn alpha | 0.05 | Empirical-Bayes alpha |
| swish qvalue | < 0.05 | Local FDR from qvalue package |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `Error in importRdata: ... transcript_ids do not match` | Salmon index built from different annotation than provided GTF | Rebuild Salmon index with matching transcripts.fa |
| `analyzeSwitchConsequences: not enough switching genes` | preFilter too strict; few candidate switches | Lower `geneExpressionCutoff`, `dIFcutoff` in test step |
| `dmFilter: empty result` | Filter parameters too strict | Reduce `min_feature_expr`, `min_samps_feature_expr` |
| `satuRn: rank-deficient design matrix` | Confounder perfectly correlated with condition | Drop the confounder or stratify analysis |
| `swish: no inferential replicates found` | Salmon run without `--numGibbsSamples` | Re-run Salmon with `--numGibbsSamples 20` |
| `analyzeIUPred2A: file not found` | External annotator output missing or wrong path | Verify CPC2/Pfam/SignalP/IUPred2A all completed and paths match |

## Common Pitfalls

- **Skipping stageR** -> inflated transcript-level FDR; gene-level multiple-testing burden ignored.
- **Forgetting NMD direction** -> sign-of-effect on protein opposite to sign-of-effect on transcript when alternative form is a PTC-bearer. Always check.
- **Treating short-read-derived isoform calls as ground truth** -> Salmon EM is uncertain; use Gibbs samples + swish if quantification uncertainty matters.
- **Comparing across annotations** -> GENCODE basic vs comprehensive, RefSeq, Ensembl all have different transcript catalogs; switches "appear" or "disappear" with annotation choice. Document version.
- **Not running long-read where possible** -> Iso-Seq / ONT removes ambiguity for genes with many similar isoforms (TTN, MAPT, NEFM, DSCAM).
- **Choosing satuRn or DEXSeq blindly at the n=5 boundary** -> IsoformSwitchAnalyzeR v2 auto-selects based on >5 vs <=5; results may differ. Document choice.
- **Reporting a "switch" without a sashimi plot** -> reviewers will demand it; do it upfront.
- **Forgetting stageR also corrects gene-level p when starting from DRIMSeq** -> DRIMSeq's `gene_p` should be passed as `pScreen`, not raw transcript p-values.

## Related Skills

- differential-splicing - Event-level (rMATS, leafcutter, MAJIQ) complementary to DTU
- splicing-quantification - PSI is a 1D projection of DTU shifts
- splicing-qc - Verify upstream library, depth, alignment before DTU
- sashimi-plots - Required visualization for switch validation and reporting
- splice-variant-prediction - Connects SpliceAI variant predictions to specific isoforms
- long-read-splicing - Full-isoform DTU bypasses transcript-quant uncertainty; preferred for many-isoform genes
- pathway-analysis/go-enrichment - Pathway enrichment of switching genes
- rna-quantification/alignment-free-quant - Salmon with `--numGibbsSamples` is upstream

## References

- Han et al 2025 *bioRxiv* 10.64898/2025.12.08.693027 - IsoformSwitchAnalyzeR v2
- Vitting-Seerup & Sandelin 2019 *Bioinformatics* 35:4469-4471 - IsoformSwitchAnalyzeR original
- Anders et al 2012 *Genome Res* - DEXSeq
- Nowicka & Robinson 2016 *F1000Research* - DRIMSeq
- Gilis et al 2022 *F1000Research* - satuRn
- Zhu et al 2019 *NAR* - swish / fishpond
- Van den Berge et al 2017 *Genome Biol* - stageR
- Love, Soneson, Patro 2018 *F1000Research* - Swimming downstream DTU workflow
- Maquat 2004 *Nat Rev Mol Cell Biol* - NMD review
- Lykke-Andersen & Jensen 2015 *Nat Rev Mol Cell Biol* - NMD update
- Lindeboom et al 2016 *Nat Genet* - NMD escape rates from RNA-seq
- Lareau et al 2007 *Nature* - SR protein AS-NMD autoregulation
- Ni et al 2007 *Genes Dev* - ultraconserved-element AS-NMD in splicing regulators
- Cuccurese et al 2005 *NAR* 33:5965-5977 - ribosomal protein AS-NMD autoregulation
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43)
- Han et al 2020 *Sci Transl Med* - SCN1A poison exon ASO
