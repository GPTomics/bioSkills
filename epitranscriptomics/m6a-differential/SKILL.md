---
name: bio-epitranscriptomics-m6a-differential
description: Identifies differential m6A methylation between conditions from MeRIP-seq paired IP/input data using exomePeak2 (GC-bias-aware differential via its bam_ip/bam_input control + bam_treated_ip/bam_treated_input treatment arms), QNB beta-binomial, MeTDiff HMM, and RADAR, plus the paired-symmetric edgeR/DESeq2-on-peak-counts route when batch/lot covariates need fixed-effect handling that exomePeak2's API does not accept. Covers paired vs unpaired vs interaction designs, batch confounding and per-lot meta-analysis, the stoichiometry-vs-expression-vs-IP-efficiency confound, and effect-size filtering against under-powered N=2 designs. Use when comparing m6A across two or more conditions, choosing between exomePeak2/QNB/RADAR/MeTDiff for a design, handling batch confounding when exomePeak2's API is too rigid, distinguishing real hyper/hypo-methylation from expression shifts, applying effect-size thresholds, or planning orthogonal stoichiometry validation (GLORI/SAC-seq/m6Anet mod_ratio).
tool_type: r
primary_tool: exomePeak2
---

## Version Compatibility

Reference examples tested with: exomePeak2 1.14+ (Bioconductor 3.18+), QNB 1.1.11 (GitHub `lzcyzm/QNB`), MeTDiff (GitHub, bundled with MeTPeak), RADAR 0.2.4+ (GitHub `scottzijiezhang/RADAR`), DESeq2 1.42+, edgeR 4.0+, GenomicFeatures 1.54+, ggplot2 3.5+, GenomicAlignments 1.38+, Rsubread 2.16+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('exomePeak2')` then `?exomePeak2` to verify parameters
- R: `packageVersion('QNB')` then `?qnbtest` to confirm argument signature

If R throws `unused argument` or `argument is missing`, the API moved between Bioconductor minor releases; consult `?exomePeak2` directly. QNB function signature has been stable since 2017 but the GitHub source has occasional changes; pin the commit SHA.

exomePeak2's differential interface is NOT a `mode=` argument — populate the `bam_treated_ip=` and `bam_treated_input=` arguments alongside the standard `bam_ip=` / `bam_input=` (control arm) to trigger paired differential calling. `peak_calling_mode` is a separate argument controlling locus scope (`exon` | `full_transcript` | `whole_genome`). Verify against `?exomePeak2` and Bioconductor 3.20+ release notes. QNB is on GitHub only; install via `devtools::install_github('lzcyzm/QNB')`. RADAR is on GitHub only; install via `devtools::install_github('scottzijiezhang/RADAR')`; its differential workflow is `countReads -> normalizeLibrary -> adjustExprLevel -> filterBins -> diffIP -> reportResult`.

# Differential m6A Analysis

**"Compare m6A methylation between my conditions"** -> Quantify how much each m6A peak's IP/input enrichment shifts between conditions, after normalising for transcript-abundance changes (which all show up in input) and for IP-efficiency drift (which the design matrix and within-run replicates control for). Then apply effect-size filtering to distinguish real biology from the technical noise floor that all MeRIP differential methods inherit (McIntyre 2020 *Sci Rep* 10:6590: between-study m6A peak overlap is ~45% median; differential calls within that noise envelope routinely fail to replicate). Critically: a higher MeRIP signal in condition A vs B can mean (1) more transcripts of the peak-bearing gene, (2) more methylation per transcript, OR (3) higher IP efficiency in batch A — distinguishing requires careful normalisation or an orthogonal absolute-stoichiometry method.

- R: `exomePeak2::exomePeak2(bam_ip, bam_input, bam_treated_ip, bam_treated_input, ...)` -- integrated peak + differential GLM (modern default)
- R: `QNB::qnbtest(control_ip, treated_ip, control_input, treated_input)` -- beta-binomial for small N (Liu 2017 *BMC Bioinformatics* 18:387)
- R: `RADAR::diffIP()` then `reportResult()` -- Poisson-NB on peak windows with TMM normalisation (Zhang 2019 *Genome Biol* 20:294)
- R: `MeTDiff::metdiff()` -- HMM-based differential paired with MeTPeak
- R: edgeR / DESeq2 on featureCounts-on-peaks matrix -- defensible for paired symmetric designs with strong input normalisation; also the route for arbitrary batch / lot covariates exomePeak2's API does not accept

## The Single Most Important Modern Insight -- IP fold-change between conditions conflates stoichiometry change with expression change and IP efficiency drift

A higher m6A peak signal in condition A vs B can mean ANY of: (1) more transcripts of the peak-bearing gene (expression up; the input increases proportionally, so the ratio should not change — but residual normalisation noise leaks through), (2) more methylation per transcript (stoichiometry up; the real biology of interest), (3) higher IP efficiency in batch A (technical; antibody lot, IP day, RNA prep). Differential MeRIP WITHOUT per-window input-normalisation OR an orthogonal stoichiometry-aware method (GLORI Liu 2023 *Nat Biotechnol* 41:355; SAC-seq Hu 2022 *Nat Biotechnol* 40:1210; MAZTER-seq Garcia-Campos 2019 *Cell* 178:731; m6Anet per-read modification rate) cannot separate the three. exomePeak2 differential mode, QNB, RADAR, and MeTDiff all implement per-window IP/input ratio modelling, but each has different default normalisations and the choice matters. Equally critical: McIntyre 2020 *Sci Rep* 10:6590 showed that "differential" m6A peaks from MeRIP-seq routinely do not replicate between independent studies in nominally identical conditions. The empirical noise floor is high; effect-size filtering (|log2FC| >= 0.5 minimum, often >= 1) AND replicate-direction concordance (the change is consistent in direction across replicates) AND minimum N=3 per condition are needed for differential calls to survive replication. For any absolute stoichiometry claim ("this peak is 80% methylated in tumour vs 20% in normal"), require an orthogonal stoichiometry method, not MeRIP alone.

## Algorithmic Taxonomy

| Tool / mode | Mechanism | Inputs | Output | Strength | Fails when |
|-------------|-----------|--------|--------|----------|------------|
| exomePeak2 differential (Liu 2022) | Transcript-windowed Poisson GLM with GC-bias correction; integrated peak + differential | control IP/input BAM vectors + treated IP/input BAM vectors + TxDb + BSgenome | Differential peaks with log2FC + FDR per peak | Most use cases; integrates with peak calling; modern default | Small-N (<=2) overdispersion poorly estimated; top-level API does NOT accept arbitrary covariates (use DESeq2/edgeR for batch-aware designs) |
| QNB (Liu 2017 *BMC Bioinformatics* 18:387) | Quad-negative-binomial joint model of IP, input, condition | per-peak count matrices (4 matrices: ip1, ip2, input1, input2) | Differential peaks with p-value + log2 RR | Designed for small N (2-3 per group); handles overdispersion explicitly | Requires pre-computed count matrices; not integrated with peak calling |
| MeTDiff (bundled with MeTPeak) | HMM + Beta-binomial differential paired with MeTPeak | paired IP/input BAM + GTF + condition factor | Differential peaks per window | Pairs naturally with MeTPeak output; HMM smoothing helps low-coverage | GitHub-only; less benchmarked than exomePeak2 / QNB |
| RADAR (Zhang 2019 *Genome Biol* 20:294) | Poisson-NB with TMM normalisation; reproducibility-aware | paired IP/input BAM + condition factor | Differential peaks with logFC + FDR | Explicit replicate-variance modeling; reproducibility-aware framework | GitHub-only; slower than exomePeak2 |
| DRME (Liu 2016 *Anal Biochem* 499:15) | Count-based small-sample alternative | per-peak count matrices | Differential peaks | Sister to QNB from same group; small-N alternative | Less benchmarked than QNB / exomePeak2 |
| edgeR / DESeq2 on peak counts | Generic RNA-seq differential framework applied to featureCounts-on-peaks | peak count matrix + sample sheet | Differential peaks with log2FC + FDR | Familiar; flexible designs; well-tested in RNA-seq | Treats peak counts as RNA counts; loses IP/input pairing structure; defensible only for paired symmetric designs with strong input normalisation |
| Ratio-of-ratios (heuristic) | Compute per-peak log2 (IP_A / Input_A) - log2 (IP_B / Input_B) per sample, then t-test | per-peak count matrix | per-peak t-test | Transparent; no model assumptions | No multiple-testing correction; ignores overdispersion; not recommended for primary analysis |

## Decision Tree by Scenario

| Scenario | Recommended | Why wrong choices fail |
|----------|-------------|------------------------|
| Standard 3-vs-3 paired-design MeRIP differential | exomePeak2 with `bam_ip` (ctrl) + `bam_treated_ip` (treat) and matching inputs | QNB usable but designed for smaller N; edgeR/DESeq2 loses IP/input pairing |
| Very small N (2 vs 2) | QNB (designed for small-sample overdispersion); supplement with exomePeak2 if possible | edgeR / DESeq2 dispersion estimation collapses; exomePeak2 GLM also struggles at N=2 |
| Paired design (patient as blocking factor) | QNB per-pair then aggregate; OR featureCounts-on-peaks -> DESeq2 with `~patient + condition` design; exomePeak2 cannot encode patient blocking via its top-level API | Unpaired analysis inflates within-group variance |
| Interaction design (genotype × treatment) | featureCounts-on-peaks -> DESeq2 / edgeR with interaction term; exomePeak2 top-level API is two-group only | QNB pairwise only; build interaction model from pairwise contrasts manually |
| Batch confounding (antibody lot, sequencing run, IP day) | Include batch as fixed effect in DESeq2 / edgeR model on featureCounts-on-peaks matrix; OR run exomePeak2 per-batch and meta-analyse | Pooling cross-batch counts without batch term attributes lot-effect to condition |
| Time-course differential | featureCounts-on-peaks -> DESeq2 / limma with time as numeric covariate; OR pairwise time-point exomePeak2 contrasts | Naive group-vs-group ignores time structure |
| Stoichiometry claims (not just enrichment) | NOT MeRIP differential -- orthogonal GLORI / SAC-seq / m6Anet per-read | MeRIP IP fold-change is relative; cannot give per-molecule stoichiometry |
| Cross-batch differential (different antibody lots) | featureCounts-on-peaks -> DESeq2 / edgeR with `lot` in design; OR run exomePeak2 per-lot then meta-analyse; ideally avoid confounding lot with condition at the experimental-design stage | Lot effect inflates false positives; exomePeak2 top-level API cannot encode lot |
| Validation of differential calls | Run >=2 differential methods; require concordant direction across replicates AND |log2FC| > 0.5 minimum (>= 1 stringent); orthogonal validation at top hits | Single-method differential calls within technical noise floor (per McIntyre 2020) |
| Visualising differential peaks | Volcano plot with |log2FC| + FDR thresholds; MA plot to inspect normalisation; per-peak boxplot for top hits | Single number summaries hide stoichiometry vs expression confound |
| Wanting to test a single gene / locus | Targeted: per-peak boxplot across replicates with condition factor; manual t-test or Wilcoxon at high-coverage peak | Whole-transcriptome differential testing wastes multiple-testing budget for single-locus questions |

Methodology evolves; before any high-stakes differential analysis, web-search "exomePeak2 differential mode Bioconductor 3.20" and "MeRIP differential benchmark McIntyre" for current consensus parameters.

## exomePeak2 Differential Workflow

**Goal:** Identify m6A peaks that differ in methylation level between conditions, controlling for transcript-abundance differences (via input normalisation) and GC bias (via internal correction), with an integrated peak-calling + differential pipeline.

**Approach:** Build TxDb from the matched GTF; pass control IP/input BAM vectors via `bam_ip` and `bam_input` AND treatment IP/input BAM vectors via `bam_treated_ip` and `bam_treated_input` — populating the treated arms triggers differential mode (there is no separate `mode=` argument). Output is per-peak log2FC + FDR.

```r
library(exomePeak2)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb <- makeTxDbFromGFF('refs/annotation.gtf', format='gtf')

ctrl_ip      <- c('aligned/ctrl_IP1.bam', 'aligned/ctrl_IP2.bam', 'aligned/ctrl_IP3.bam')
ctrl_input   <- c('aligned/ctrl_Input1.bam', 'aligned/ctrl_Input2.bam', 'aligned/ctrl_Input3.bam')
treat_ip     <- c('aligned/treat_IP1.bam', 'aligned/treat_IP2.bam', 'aligned/treat_IP3.bam')
treat_input  <- c('aligned/treat_Input1.bam', 'aligned/treat_Input2.bam', 'aligned/treat_Input3.bam')

result <- exomePeak2(
    bam_ip             = ctrl_ip,
    bam_input          = ctrl_input,
    bam_treated_ip     = treat_ip,
    bam_treated_input  = treat_input,
    txdb               = txdb,
    genome             = BSgenome.Hsapiens.UCSC.hg38,
    paired_end         = TRUE,
    library_type       = 'unstranded',
    peak_calling_mode  = 'exon',
    save_dir           = 'exomepeak2_diff_output',
    experiment_name    = 'ctrl_vs_treat'
)

diff_table <- as.data.frame(result)
nrow(diff_table)
head(diff_table[, c('seqnames', 'start', 'end', 'log2FC', 'pvalue', 'padj')])
```

`peak_calling_mode` accepts `'exon'` (transcript-aware, default), `'full_transcript'`, or `'whole_genome'`; the meaning is locus scope, NOT differential-vs-non-differential. For arbitrary covariate adjustment (batch, antibody lot, patient blocking), the exomePeak2 top-level API is insufficient — move counts into DESeq2 / edgeR via the featureCounts-on-peaks route below.

## QNB Beta-Binomial for Small-N Designs

**Goal:** Test differential m6A at pre-called peaks using a quad-negative-binomial model that handles small-N overdispersion better than generic GLM frameworks.

**Approach:** Count reads in each IP and Input BAM at each peak using featureCounts or summarizeOverlaps; pass the four count matrices (ip1, ip2, input1, input2 — ip/input per group) to `qnbtest()`.

```r
library(QNB)
library(Rsubread)
library(rtracklayer)

peaks <- import('exomepeak2_output/m6a_run1/peaks.bed')
peak_saf <- data.frame(
    GeneID = paste0('peak_', seq_along(peaks)),
    Chr    = as.character(seqnames(peaks)),
    Start  = start(peaks),
    End    = end(peaks),
    Strand = as.character(strand(peaks))
)

count_matrix <- function(bam_paths, peak_saf) {
    fc <- featureCounts(
        files       = bam_paths,
        annot.ext   = peak_saf,
        isPairedEnd = TRUE,
        nthreads    = 8,
        allowMultiOverlap = TRUE
    )
    fc$counts
}

ip_ctrl   <- count_matrix(c('aligned/ctrl_IP1.bam', 'aligned/ctrl_IP2.bam', 'aligned/ctrl_IP3.bam'), peak_saf)
ip_treat  <- count_matrix(c('aligned/treat_IP1.bam', 'aligned/treat_IP2.bam', 'aligned/treat_IP3.bam'), peak_saf)
in_ctrl   <- count_matrix(c('aligned/ctrl_Input1.bam', 'aligned/ctrl_Input2.bam', 'aligned/ctrl_Input3.bam'), peak_saf)
in_treat  <- count_matrix(c('aligned/treat_Input1.bam', 'aligned/treat_Input2.bam', 'aligned/treat_Input3.bam'), peak_saf)

qnb_result <- qnbtest(
    control_ip    = ip_ctrl,
    treated_ip    = ip_treat,
    control_input = in_ctrl,
    treated_input = in_treat,
    mode          = 'per-condition'
)

head(qnb_result)
sig <- qnb_result[qnb_result$padj < 0.05 & abs(qnb_result$log2.RR) > 0.5, ]
nrow(sig)
```

Verify QNB argument names against `?qnbtest` for the installed version; older tutorials may show different signatures.

## RADAR Reproducibility-Aware Differential

**Goal:** Test differential m6A using a Poisson-NB framework with TMM normalisation and explicit replicate-variance modeling; useful when replicate variability is a known issue.

**Approach:** RADAR's documented workflow is `countReads -> normalizeLibrary -> adjustExprLevel -> filterBins -> diffIP -> reportResult`. CRITICAL: RADAR expects matched BAMs in `bamFolder` named `<sample>.input.bam` and `<sample>.m6A.bam` per replicate; the generic IP / Input naming used elsewhere must be re-conformed or symlinked. `variable()` is set with a data.frame, NOT a bare factor.

```r
library(RADAR)

radar <- countReads(
    samplenames  = c('ctrl_rep1', 'ctrl_rep2', 'ctrl_rep3', 'treat_rep1', 'treat_rep2', 'treat_rep3'),
    gtf          = 'refs/annotation.gtf',
    bamFolder    = 'aligned_radar/',
    modification = 'm6A',
    strandToKeep = 'opposite',
    threads      = 8
)

radar <- normalizeLibrary(radar)
radar <- adjustExprLevel(radar)

variable(radar) <- data.frame(group = c('ctrl', 'ctrl', 'ctrl', 'treat', 'treat', 'treat'))

radar <- filterBins(radar, minCountsCutOff = 15)
radar <- diffIP(radar)

result <- reportResult(radar, cutoff = 0.1, Beta_cutoff = 0.5)
sig <- result[result$padj < 0.05 & abs(result$logFC) > 0.5, ]
nrow(sig)
```

`reportResult` thresholds (`cutoff` = p-value cutoff; `Beta_cutoff` = effect-size cutoff in beta units) are RADAR-specific; convert to the project's standard reporting thresholds downstream. The `aligned_radar/` directory should contain BAMs with RADAR's expected naming (`<sample>.input.bam`, `<sample>.m6A.bam`).

## Volcano Plot of Differential Peaks

**Goal:** Visualise the differential peak set with effect size on the x-axis and statistical significance on the y-axis; flag peaks passing |log2FC| and FDR thresholds.

**Approach:** Standard ggplot2 volcano with colour-coded significance and threshold lines.

```r
library(ggplot2)

diff_table$significance <- with(diff_table,
    ifelse(padj < 0.05 & abs(log2FC) > 0.5, 'differential', 'not_sig'))

ggplot(diff_table, aes(x=log2FC, y=-log10(padj), colour=significance)) +
    geom_point(alpha=0.5, size=0.8) +
    geom_vline(xintercept=c(-0.5, 0.5), linetype='dashed') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') +
    scale_colour_manual(values=c(differential='red', not_sig='grey60')) +
    labs(x='log2 (treat / ctrl) MeRIP enrichment ratio',
         y='-log10 (FDR)',
         title='Differential m6A peaks: ctrl vs treat',
         caption='Per-window IP/input ratio normalised; not absolute stoichiometry') +
    theme_minimal()
```

The caption is intentional: MeRIP differential reports CHANGES IN ENRICHMENT RATIO, NOT changes in absolute stoichiometry. For stoichiometry claims, cross-validate with GLORI / SAC-seq / m6Anet.

## Per-Method Failure Modes

### Reporting "hyper-methylation" without orthogonal calibration

**Trigger:** "Peak X shows hyper-methylation in treatment" inferred from MeRIP IP fold-change alone.

**Mechanism:** MeRIP IP fold-change conflates per-molecule methylation stoichiometry, transcript abundance, and IP efficiency variation between libraries. An IP fold-change increase can reflect any or all of these.

**Symptom:** Reported m6A "hyper-methylation" tracks RNA-seq expression changes between conditions; reverse-direction effects when properly normalised against input.

**Fix:** Use "increased / decreased enrichment" terminology for MeRIP-only studies. Reserve "hyper- / hypo-methylated" for studies with absolute quantification orthogonal validation (GLORI, SAC-seq, MAZTER-seq, m6Anet per-read). For high-stakes claims at named loci, run GLORI on a subset of conditions.

### Effect-size threshold absence

**Trigger:** Reporting "1,500 differential m6A peaks" with FDR < 0.05 (uncorrected p-value or naive multiple testing) and no effect-size filter.

**Mechanism:** With sufficient sequencing depth, MeRIP-seq has high statistical power to detect very small (~1.1-1.2x) IP-ratio changes that lie within antibody / technical noise. McIntyre 2020 *Sci Rep* 10:6590 showed these changes do not replicate.

**Symptom:** Differential peak set has many peaks with small effect sizes; replication in an independent study recovers <30% of original calls.

**Fix:** Apply effect-size filter (|log2FC| >= 0.5 minimum, often >= 1) AND adjusted p-value (FDR < 0.05) AND replicate-direction concordance. Differential peaks should be reported with effect size, not just p-value. Report effect-size distribution alongside peak count.

### Underpowered N=2 design

**Trigger:** Differential m6A study with N=2 IP and N=2 input per condition.

**Mechanism:** Per McIntyre 2020 and many subsequent benchmarks, MeRIP replicate variance is high; N=2 estimates of dispersion are unreliable; differential calls are unstable.

**Symptom:** Many "differential" peaks; volcano plot dense; small fraction replicates in held-out replicate.

**Fix:** Minimum N=3 per condition (per condition per IP/input arm = 12 BAMs for a 2-condition study); N=4-5 preferred for high-stakes claims. Underpowered studies should report effect-size-only filtered subsets (the most extreme peaks) and acknowledge the noise floor explicitly.

### Batch confounded with condition (antibody lot, IP day)

**Trigger:** Control samples processed in batch 1 with antibody lot A; treatment samples processed in batch 2 with antibody lot B.

**Mechanism:** Anti-m6A antibody lots have batch-to-batch variability in pulldown efficiency and m6A-vs-m6Am cross-reactivity. Pooling cross-batch counts in a differential model attributes batch-effect to condition.

**Symptom:** "Differential" peaks at high-abundance transcripts; effect sizes track batch rather than condition; reanalysis with `batch` in the model removes most differential peaks.

**Fix:** Include `antibody_lot` / `batch` / `prep_day` as a fixed effect in a DESeq2 / edgeR model on featureCounts-on-peaks counts; OR run exomePeak2 separately per lot and meta-analyse the per-lot differential peak sets; OR re-design the experiment to avoid lot-condition confounding. exomePeak2's top-level API does NOT accept arbitrary covariates — DESeq2 / edgeR is the route when covariate handling is required.

### exomePeak2 invoked with a fabricated `mode=` argument

**Trigger:** `exomePeak2(..., mode='differential')` OR `exomePeak2(..., mode='diff_peak')` returns "unused argument" error.

**Mechanism:** exomePeak2 has NO `mode=` argument. Differential is triggered by populating `bam_treated_ip` and `bam_treated_input` alongside the standard `bam_ip` and `bam_input` (control arm). The `peak_calling_mode` argument is unrelated — it accepts `'exon' | 'full_transcript' | 'whole_genome'` and controls locus scope, not differential-vs-non-differential.

**Fix:** Populate the four BAM-vector arguments (`bam_ip`, `bam_input`, `bam_treated_ip`, `bam_treated_input`); drop any `mode=` reference; consult `?exomePeak2` for the authoritative signature in the installed version.

### Treating peak count matrices like RNA count matrices in edgeR/DESeq2

**Trigger:** Compute featureCounts at peaks, build a count matrix, run edgeR / DESeq2 on the matrix as if peaks were genes.

**Mechanism:** Peak counts reflect both IP enrichment AND transcript abundance. Generic RNA-seq DE on peak counts mixes the two; size-factor normalisation on IP-only counts loses the input-pair information.

**Fix:** edgeR / DESeq2 on peak counts is defensible ONLY for paired symmetric designs where input is modelled as an offset (per-sample size factor on input counts AND per-sample size factor on IP counts, then differential on the ratio). For most uses, exomePeak2 / QNB / RADAR's purpose-built models are more appropriate.

### Counting reads at peaks WITHOUT featureCounts strand-awareness

**Trigger:** `featureCounts(...)` invoked without `strandSpecific=` flag; or with wrong strand setting.

**Mechanism:** Strand-aware counting matters when the protocol is stranded (most modern MeRIP is unstranded; some are reverse-stranded). Wrong strand counts include antisense reads as if they were sense.

**Fix:** Verify protocol strandedness from sequencing-core notes or by inspecting featureCounts summary at a few transcripts. Pass `strandSpecific=0` (unstranded), `1` (forward), or `2` (reverse) explicitly.

## Reconciliation: When Differential Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| exomePeak2 calls a differential peak; QNB does not | Different dispersion estimates; QNB more conservative at low coverage | Trust intersection; report concordant set as high-confidence |
| Most "differential" peaks fall at high-abundance housekeeping transcripts | Expression / IP-efficiency confound | Check log2FC vs input log2FC; if correlated, batch effect or expression-driven |
| RADAR vs exomePeak2 disagree on direction | Different normalisation (TMM vs internal) | Inspect normalisation diagnostic; choose method aligned with experimental design |
| Differential peaks anti-correlated with RNA-seq DE | Expression conflated with methylation in the differential model | Re-run with stronger input normalisation; consider per-peak ratio normalisation |
| Single differential peak survives across all methods | High-confidence call | Orthogonally validate (GLORI / SAC-seq) at the named locus |
| Differential calls scatter randomly across genome | Underpowered; technical noise dominates | Increase N; apply stricter effect-size filter; report null result honestly |
| Cross-condition peak overlap < 50% before differential | Conditions are biologically very different; OR antibody lot effect | Inspect cross-replicate concordance; check antibody lot metadata |
| Volcano shows extreme outliers at low-coverage peaks | Per-peak variance dominated by Poisson sampling | Filter peaks by minimum coverage (>=30 reads in IP AND input) before differential |

## Quantitative Thresholds

| Quantity | Threshold | Source / rationale |
|----------|-----------|--------------------|
| Minimum biological replicates per condition | 3 (4-5 preferred) | McIntyre 2020 *Sci Rep* 10:6590 — N=2 routinely under-powered |
| FDR threshold | 0.05 | Standard convention |
| Effect-size threshold (|log2FC|) | >= 0.5 minimum; >= 1 stringent | Below 0.5, calls within technical noise floor |
| Minimum coverage per peak window (each IP AND input) | 30 reads | Standard convention; below this, statistical calls noisy |
| Replicate-direction concordance | Same direction in >=2/N replicates | Guardrail against single-replicate artifacts |
| Antibody lot tracking | Mandatory in design matrix when lots differ | Lot confounding is a known false-positive source |
| MeRIP per-window IP/input ratio inflation | log2(IP/input) >= 1 for "enriched"; >=2 for "strongly enriched" | Convention |
| Cross-study peak overlap baseline | ~45% median between labs | McIntyre 2020 *Sci Rep* 10:6590 — bounds inter-study reproducibility |
| Effect size for biological validation | |log2FC| >= 1 AND padj < 0.05 typically used for downstream wet-lab follow-up | Field convention; tighten for low-N studies |
| GLORI orthogonal validation threshold | Stoichiometry change >= 10% at named site | Liu C 2023 *Nat Biotechnol* 41:355 calibration |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| exomePeak2 `mode='differential'` rejected | exomePeak2 has no `mode=` arg; differential is via `bam_treated_ip` + `bam_treated_input` | Populate the four BAM-vector args; consult `?exomePeak2` |
| QNB error "unused argument" | Function signature changed between versions | Pin commit SHA; verify against installed `?qnbtest` |
| RADAR install fails | GitHub-only; requires devtools and specific Bioconductor deps | `devtools::install_github('scottzijiezhang/RADAR')`; check Bioconductor requirements |
| featureCounts returns zero counts | Strand setting wrong; OR peak BED has different chromosome naming than BAM | Verify `strandSpecific=`; reconcile chromosome names |
| Volcano plot all peaks near origin | Low effect sizes; technical noise dominant | Check input normalisation; increase N; report null result if appropriate |
| FDR-significant but small effect size | Large dataset finds small differences with high power | Apply effect-size filter; report effect-size distribution |
| Many differential peaks track expression changes | Input normalisation insufficient | Re-run with stronger input adjustment; OR use ratio-of-ratios |
| edgeR `estimateDisp()` fails at N=2 | Small-N dispersion estimation collapses | Use QNB instead; designed for this case |
| MeTDiff install fails | GitHub-only, bundled with MeTPeak | `devtools::install_github('compgenomics/MeTPeak')` |
| exomePeak2 result has no log2FC / padj columns | `bam_treated_ip` / `bam_treated_input` were not populated; only ran peak calling on the control arm | Pass all four BAM-vector arguments (`bam_ip`, `bam_input`, `bam_treated_ip`, `bam_treated_input`) to trigger differential output |
| Per-peak boxplot shows huge within-condition variance | High biological noise; OR one replicate is an outlier | Inspect plotCorrelation in merip-preprocessing for outlier; consider exclusion |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "How many replicates per condition?" | N=3 minimum; N=4-5 preferred; rationale McIntyre 2020 |
| "What's the effect-size threshold?" | |log2FC| >= 0.5 minimum; results table reports effect size alongside FDR |
| "Was batch / antibody lot controlled for?" | Yes — `antibody_lot` included as fixed effect in DESeq2 model on featureCounts-on-peaks counts (exomePeak2's top-level API does not accept covariates); per-lot exomePeak2 + meta-analysis when DESeq2 path infeasible; lot-condition confounding assessed |
| "Does the differential signal track expression changes?" | Cross-checked log2FC vs input log2FC per peak; only peaks with strong IP/input ratio shift reported as differential |
| "Was orthogonal validation done?" | Top hits orthogonally validated via GLORI / m6Anet / per-locus assay |
| "Why exomePeak2 over QNB?" | exomePeak2 default for standard 3-vs-3; QNB used for small-N sensitivity analyses; both reported as concordance check |
| "What about absolute stoichiometry?" | MeRIP differential reports relative enrichment changes; absolute stoichiometry requires GLORI / SAC-seq |
| "Was the right normalisation used?" | Per-window IP/input ratio normalisation; exomePeak2 internal; alternative TMM (RADAR) compared |
| "How does this replicate in independent studies?" | Cross-study peak overlap reported; differential subset checked against published m6A-Atlas |
| "Why not edgeR / DESeq2?" | Generic RNA-seq DE on peak counts loses IP/input pairing structure; used only for sensitivity analysis with paired symmetric design |

## References

- Dominissini D, Moshitch-Moshkovitz S, Schwartz S et al (2012) Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. *Nature* 485(7397):201-206. doi:10.1038/nature11112
- Meyer KD, Saletore Y, Zumbo P, Elemento O, Mason CE, Jaffrey SR (2012) Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. *Cell* 149(7):1635-1646. doi:10.1016/j.cell.2012.05.003
- Liu L, Zhang SW, Huang Y, Meng J (2017) QNB: differential RNA methylation analysis for count-based small-sample sequencing data with a quad-negative binomial model. *BMC Bioinformatics* 18(1):387. doi:10.1186/s12859-017-1808-4
- Cui X, Meng J, Zhang S, Chen Y, Huang Y (2016) A novel algorithm for calling mRNA m6A peaks by modeling biological variances in MeRIP-seq data. *Bioinformatics* 32(12):i378-i385. doi:10.1093/bioinformatics/btw281
- Liu L, Zhang SW, Gao F et al (2016) DRME: count-based differential RNA methylation analysis at small sample size scenario. *Anal Biochem* 499:15-23. doi:10.1016/j.ab.2016.01.014
- Zhang Z, Zhan Q, Eckert M et al (2019) RADAR: differential analysis of MeRIP-seq data with a random effect model. *Genome Biol* 20(1):294. doi:10.1186/s13059-019-1915-9
- Meng J, Lu Z, Liu H et al (2014) A protocol for RNA methylation differential analysis with MeRIP-Seq data and exomePeak R/Bioconductor package. *Methods* 69(3):274-281. doi:10.1016/j.ymeth.2014.06.008
- Liu J, Zhang Z, Meng J et al (2022) exomePeak2: a peak calling and differential analysis tool for MeRIP-Seq with bias awareness. *NAR Genom Bioinform* 4(3):lqac046. doi:10.1093/nargab/lqac046
- McIntyre ABR, Gokhale NS, Cerchietti L, Jaffrey SR, Horner SM, Mason CE (2020) Limits in the detection of m6A changes using MeRIP/m6A-seq. *Sci Rep* 10(1):6590. doi:10.1038/s41598-020-63355-3
- Liu C, Sun H, Yi Y et al (2023) Absolute quantification of single-base m6A methylation in the mammalian transcriptome using GLORI. *Nat Biotechnol* 41(3):355-366. doi:10.1038/s41587-022-01487-9
- Hu L, Liu S, Peng Y et al (2022) m6A RNA modifications are measured at single-base resolution across the mammalian transcriptome. *Nat Biotechnol* 40(8):1210-1219. doi:10.1038/s41587-022-01243-z
- Garcia-Campos MA, Edelheit S, Toth U et al (2019) Deciphering the m6A code via antibody-independent quantitative profiling. *Cell* 178(3):731-747.e16. doi:10.1016/j.cell.2019.06.013
- Robinson MD, McCarthy DJ, Smyth GK (2010) edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 26(1):139-140. doi:10.1093/bioinformatics/btp616
- Love MI, Huber W, Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Liao Y, Smyth GK, Shi W (2014) featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics* 30(7):923-930. doi:10.1093/bioinformatics/btt656

## Related Skills

- merip-preprocessing - Upstream IP/input BAM preparation; design-matrix metadata (antibody lot, batch) originates here
- m6a-peak-calling - Peak calling step that produces input to differential analysis
- m6anet-analysis - Orthogonal ONT-direct-RNA validation for stoichiometry claims at high-stakes loci
- modification-visualization - Volcano / MA / per-peak boxplot rendering of differential results
- differential-expression/deseq2-basics - Canonical DE design philosophy; m6a-differential defers to this for general design-matrix patterns
- differential-expression/de-results - Post-DE interpretation, ranking, gene-list extraction
- differential-expression/edger-basics - edgeR fundamentals for the paired-symmetric sensitivity case
- chip-seq/differential-binding - Closest sibling for IP-vs-input differential binding (general framework)
- rna-quantification/featurecounts-counting - Peak count matrix construction
- data-visualization/volcano-and-ma-plots - Volcano + MA plot recipes
- data-visualization/multipanel-figures - Figure assembly
- pathway-analysis/go-enrichment - GO enrichment on differential-peak-bearing gene lists
- workflows/rnaseq-to-de - End-to-end pipeline orchestration patterns
