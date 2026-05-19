---
name: bio-clip-seq-differential-clip
description: Identify differentially bound regions across CLIP-seq conditions (knockdown vs control, treatment vs vehicle, disease vs healthy) using DEWSeq (sliding-window DESeq2), Flipper (Skipper-downstream), ASpeak, edgeR, or limma-voom. Use when computing condition-level changes in RBP binding intensity, choosing peak-level vs window-level vs crosslink-level testing, designing replicate experiments, or distinguishing biological binding shifts from technical confounders.
tool_type: r
primary_tool: DEWSeq
---

## Version Compatibility

Reference examples tested with: DEWSeq 1.18+, htseq-clip 2.0+, DESeq2 1.44+, edgeR 4.2+, limma 3.60+, Flipper (commit 2024.04+), Skipper (commit 2023.05+), pybedtools 0.10+, pyranges 0.0.129+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Differential CLIP-seq Analysis

**"Identify regions with changed RBP binding across conditions"** -> Test for condition-level differences in IP enrichment relative to SMInput, accounting for replicate variance and (where available) sequencing depth normalization. Three statistical scales are possible: peak-level (test each peak as a unit), window-level (test fixed transcriptome windows; DEWSeq, Flipper), or crosslink-site level (test single-nt positions). The choice depends on the biology (narrow regulatory shift vs broad binding-mode change) and on which upstream peak caller was used (CLIPper -> peak-level; Skipper -> window-level; PureCLIP -> CL-level).

- R (window-level, DEWSeq + htseq-clip): `library(DEWSeq); dds <- DESeqDataSetFromSlidingWindows(counts, colData, design=~condition); dds <- DESeq(dds); res <- results(dds)`
- R (Skipper downstream, Flipper): `flipper differential -i skipper_out/ --design design.tsv --contrast treatment vs control -o flipper_out/`
- R (peak-level, edgeR): `dge <- DGEList(counts=peak_counts, group=condition); dge <- calcNormFactors(dge); design <- model.matrix(~condition); dge <- estimateDisp(dge, design); fit <- glmQLFit(dge, design); res <- glmQLFTest(fit)`
- R (peak-level, limma-voom): `v <- voom(dge, design); fit <- lmFit(v, design); fit <- eBayes(fit); res <- topTable(fit, coef=2, number=Inf)`
- CLI (htseq-clip preprocessing): `htseq-clip extract -i annotation.gff -o annotation_windows.bed -w 50 -s 20 && htseq-clip count -i sample.bam -w annotation_windows.bed -o sample.counts.txt`

DEWSeq is the EMBL/Hentze-group windowed-NB approach designed specifically for differential CLIP. Flipper is the Skipper-companion tool (Tu 2024) for the modern Skipper workflow. Peak-level edgeR/limma-voom work when the upstream peak caller produced a comparable peak BED across conditions (e.g., consensus peaks from CLIPper).

## Algorithmic Taxonomy

| Tool | Scale | Statistical model | Replicate requirement | Strength | Fails when |
|------|-------|-------------------|----------------------|----------|------------|
| DEWSeq (Schwarzl 2020) | 50-100 nt sliding window | Negative binomial GLM (DESeq2 internals) | >= 2 reps per condition | Designed specifically for CLIP; integrates SMInput | Slow on dense libraries; output window-resolution |
| Flipper (Tu 2024) | Skipper window (100 nt) | Negative binomial; designed for Skipper output | >= 2 reps per condition | Modern; pairs with Skipper peak caller | Only useful if upstream is Skipper |
| ASpeak | Peak | Negative binomial | >= 2 reps | Peak-level for CLIPper output | Less popular; legacy |
| edgeR (general) | Peak or window | Quasi-likelihood F-test on NB | >= 2 reps | Mature, widely cited | Generic; not CLIP-aware; needs careful normalization |
| limma-voom | Peak or window | Linear model with mean-variance trend | >= 2 reps | Fast; well-validated; handles small samples | Generic; treats counts as continuous after voom |
| DESeq2 (direct) | Peak or window | Negative binomial GLM | >= 2 reps | Mature; same engine as DEWSeq | Same as edgeR caveats |
| MASTR-seq (Banks 2023) | Cell-resolved (scCLIP) | Cell-mixture model | Cells | Single-cell CLIP differential | Specialized; few published applications |
| diffbind (CLIP adaptation) | Peak | DESeq2 or edgeR backend | >= 2 reps | Familiar from ChIP/ATAC | Designed for ChIP; needs CLIP-specific normalization |
| MAnorm2 | Peak | Hierarchical empirical Bayes | >= 2 reps | Tested on ChIP; less on CLIP | Less CLIP-specific |

Methodology evolves; verify DEWSeq vignette and Flipper paper for current best practice. The DEWSeq + htseq-clip pipeline (Schwarzl 2020) is the most-published CLIP-specific differential framework; Flipper (Tu 2024) is the modern Skipper-coupled alternative.

## Critical Decision: The Interaction-Term Design

**Use `~ type + condition + type:condition` and test the interaction-term coefficient.** This is the single most consequential statistical choice in differential CLIP.

- `type` = `ip` vs `sminput` (whether the library is IP or size-matched input)
- `condition` = `treated` vs `untreated` (or KD vs control)
- `type:condition` interaction = "Does the IP-vs-input ratio shift across conditions?"

A naive `~ condition` design tests whether read counts differ regardless of whether they come from IP or SMInput - this confounds binding changes with expression changes. The interaction term explicitly tests for **differential binding** (the biologically meaningful signal) rather than differential expression at peak loci.

When testing, extract the interaction-term coefficient: `results(dds, name = 'typeip.conditiontreat')`. The log2FoldChange returned is the change in IP/input ratio in `treated` vs `untreated` - this is what "differential binding" means.

## Critical Choice: Peak-Level vs Window-Level vs Crosslink-Level

Three scales differ in resolution and statistical power:

**Peak-level (edgeR, limma-voom, DESeq2 on CLIPper peaks):** Test each consensus peak's IP/SMInput log2 FC across conditions. Pro: peak boundaries are biologically meaningful; output interpretable. Con: peak set changes across conditions (a new peak in treatment but missing in control complicates testing); SMInput normalization must be applied consistently.

**Window-level (DEWSeq, Flipper):** Tile transcriptome into fixed 50-100 nt windows; test each window. Pro: comparable across conditions (windows are pre-defined); handles binding-mode shifts within a peak; high statistical power. Con: window-resolution; multiple-testing burden (millions of windows); biological meaning of a window needs translation.

**Crosslink-level (custom):** Test each single-nt CL position. Pro: nucleotide resolution; captures motif-level shifts. Con: very low coverage per position; massive multiple-testing burden; rarely used in published differential CLIP.

| Goal | Scale | Tool |
|------|-------|------|
| ENCODE-style peak-level differential (CLIPper upstream) | Peak | DESeq2 / edgeR / limma-voom on CLIPper consensus peaks |
| Maximum sensitivity windowed differential | Window | DEWSeq (with htseq-clip) or Flipper (with Skipper) |
| Modern Skipper-coupled workflow | Window | Flipper |
| Single-cell scCLIP differential | Cell | MASTR-seq or scCLAM (sparse) |
| Compare binding-mode shifts within a peak | Window | DEWSeq |
| Allele-specific differential | CL site | BEAPR + custom logistic |
| RBP-KD effect on binding profile | Peak/window | DEWSeq (handles KD-effect on RBP itself) |

## DEWSeq Workflow (Window-Level Differential)

DEWSeq is the EMBL/Hentze-group differential framework specifically designed for CLIP. The pipeline is:

**Goal:** Identify transcriptome windows where the IP-vs-SMInput ratio shifts across conditions, accounting for replicate variance with the negative-binomial GLM.

**Approach:** Use htseq-clip to extract sliding 50 nt windows across annotated features, count reads per window per sample, build a DESeqDataSetFromSlidingWindows object with the `~ type + condition + type:condition` interaction design, extract the `typeip.conditiontreat` interaction coefficient as the differential-binding effect size, and aggregate adjacent significant windows with `bedtools merge -d 100`.

```bash
# Step 1: htseq-clip generates sliding-window count matrices
htseq-clip extract \
    -i gencode.v38.annotation.gff \
    -o annotation_windows.bed \
    --window-size 50 \
    --window-step 20 \
    --feature-type CDS,UTR

# Step 2: count IP and SMInput reads per window per sample
for sample in ip_rep1 ip_rep2 sminput_rep1 sminput_rep2; do
    htseq-clip count \
        -i ${sample}.dedup.bam \
        -a annotation_windows.bed \
        -o ${sample}.counts.txt \
        --mate 2
done
# (For eCLIP, --mate 2 because R2 5' is the truncation site; for iCLIP single-end use --mate 1)

# Step 3: DEWSeq differential testing
htseq-clip mergeCounts \
    -i ip_rep1.counts.txt ip_rep2.counts.txt sminput_rep1.counts.txt sminput_rep2.counts.txt \
    -o merged_counts.tsv
```

```r
library(DEWSeq)

counts <- read.table('merged_counts.tsv', sep='\t', header=TRUE, row.names=1)
colData <- data.frame(
    sample = c('ip_rep1','ip_rep2','sminput_rep1','sminput_rep2'),
    type = c('ip','ip','sminput','sminput'),
    condition = c('treated','treated','untreated','untreated')
)

dds <- DESeqDataSetFromSlidingWindows(
    countData = counts,
    colData = colData,
    annotObj = 'annotation_windows.bed',
    design = ~ type + condition
)

dds <- DESeq(dds)
res <- results(dds, name='condition_treated_vs_untreated')

# Window-level FDR adjustment
res_filtered <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]

# Aggregate adjacent significant windows into differential regions
sig_windows <- as.data.frame(res_filtered)
sig_windows$chr <- gsub('_.*', '', rownames(sig_windows))
# Custom reduce: combine adjacent windows within 100 nt
```

## Flipper Workflow (Skipper-Coupled)

Flipper (Tu 2024) is the differential companion to Skipper, using the same 100 nt feature-respecting windows and the same beta-binomial framework.

```bash
# Assume Skipper has been run on all samples; Skipper output is at skipper_out/
flipper differential \
    -i skipper_out/ \
    --design design.tsv \
    --contrast treatment vs control \
    -o flipper_out/

# design.tsv format:
# sample_id   condition   replicate   ip_or_input
# ip_treat_r1 treatment   1           ip
# ip_treat_r2 treatment   2           ip
# in_treat_r1 treatment   1           input
# ...
```

Output: differential window BED with log2 FC, p, padj per window.

## Peak-Level Differential (CLIPper Upstream)

```r
library(DESeq2)
library(GenomicRanges)
library(Rsubread)

# Step 1: union of CLIPper peaks across conditions
# (See bedtools merge upstream)
peaks <- read.table('consensus_peaks.bed', sep='\t', col.names=c('chr','start','end','name','score','strand'))

# Step 2: count reads per peak per sample with featureCounts
saf <- data.frame(
    GeneID = peaks$name,
    Chr = peaks$chr,
    Start = peaks$start + 1,  # 1-based for featureCounts
    End = peaks$end,
    Strand = peaks$strand
)
counts_ip <- featureCounts(
    files = c('ip_treat_r1.bam','ip_treat_r2.bam','ip_ctrl_r1.bam','ip_ctrl_r2.bam'),
    annot.ext = saf,
    isGTFAnnotationFile = FALSE,
    strandSpecific = 1,
    isPairedEnd = TRUE
)$counts

counts_in <- featureCounts(
    files = c('in_treat_r1.bam','in_treat_r2.bam','in_ctrl_r1.bam','in_ctrl_r2.bam'),
    annot.ext = saf,
    isGTFAnnotationFile = FALSE,
    strandSpecific = 1,
    isPairedEnd = TRUE
)$counts

# Step 3: DESeq2 with IP vs SMInput interaction
all_counts <- cbind(counts_ip, counts_in)
colData <- data.frame(
    type = rep(c('ip','input'), each=4),
    condition = rep(c('treat','treat','ctrl','ctrl'), 2),
    replicate = rep(c('r1','r2','r1','r2'), 2)
)

dds <- DESeqDataSetFromMatrix(countData = all_counts, colData = colData,
                               design = ~ type + condition + type:condition)
dds <- DESeq(dds)

# The interaction term `typeip.conditiontreat` tests:
# does the IP/input ratio differ in treatment vs control?
res <- results(dds, name = 'typeip.conditiontreat')

# Filter
res_sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
```

The interaction-term design (`type:condition`) is the correct statistical model for differential CLIP: it tests whether the IP-vs-input ratio differs across conditions, which is what "differential binding" means. Naive testing of just `condition` (ignoring SMInput) confounds binding changes with expression changes.

## RBP Knockdown Experiment Design

The canonical differential CLIP design is to knock down the RBP and observe what binding sites are lost. Caveats:

| Issue | Implication | Mitigation |
|-------|-------------|------------|
| RBP KD also depletes the RBP protein in cells | siRNA/shRNA reduces RBP -> reduces IP yield -> reduces unique fragments | Normalize against SMInput WITHIN each condition; the relative IP/SMInput captures binding, not protein level |
| RBP KD changes transcript abundance | mRNA stability regulators (HuR, PUM2) when knocked down change target abundance | Both IP and SMInput see the change; ratio still works |
| Off-target effects of siRNA | Multiple binding profiles change | Use multiple independent siRNAs; require concordance |
| KD efficiency varies | Lower KD -> smaller binding-loss signal | Validate KD by WB on the same IP lysate; > 70% protein loss target |
| Rescue requires re-introducing RBP | siRNA-resistant RBP cDNA for rescue | The standard differential validation experiment |

## Per-Tool Failure Modes

### DEWSeq -- Slow on dense libraries

**Trigger:** Whole-genome window tiling at 20 nt step; dense library (50M unique fragments); 4+ samples.

**Mechanism:** DEWSeq runs DESeq2 internals on millions of windows; the dispersion fit on this many features is slow.

**Symptom:** Runtime > 6 h; out-of-memory; "size of object exceeds vector limit".

**Fix:** Increase window step size to 50 nt; pre-filter windows with low counts; or restrict to expressed transcripts only. DEWSeq vignette suggests `keep <- rowSums(counts(dds)) >= 30; dds <- dds[keep,]` before testing.

### DEWSeq -- Custom adjacency aggregation needed

**Trigger:** Windows are 50 nt; biological binding sites are 50-500 nt; user expects DEWSeq to output continuous "differential regions" but gets individual windows.

**Mechanism:** DEWSeq outputs per-window results; aggregating adjacent significant windows into regions is a separate step.

**Symptom:** Output has 10,000 individual windows; user expects 1,000 biological regions.

**Fix:** Use the DEWSeq utility `resultsDEWSeq()` then `bedtools merge -d 100` on the significant-window BED. Or use the `top_hits_to_bed.R` script from DEWSeq examples.

### Peak-level differential -- Peak set differs between conditions

**Trigger:** CLIPper called peaks separately per condition; treatment has peaks at sites missing in control (and vice versa).

**Mechanism:** Peak unification requires a consensus peakset; testing on a "treatment-only" peak underestimates evidence in control (zero reads) and produces spurious DE.

**Symptom:** "Treatment-specific" peaks dominate DE results; biologically implausible.

**Fix:** Generate consensus peakset across all conditions (bedtools merge of all per-condition peak BEDs); count reads per consensus peak across all samples; THEN run differential. The Yeo lab convention is consensus peakset across all samples.

### Interaction term forgotten

**Trigger:** DESeq2 design `~ condition` instead of `~ type + condition + type:condition`.

**Mechanism:** Simple `~ condition` tests whether read counts differ between treatment and control regardless of whether reads are from IP or SMInput. A condition-driven expression change in SMInput is detected as DE binding.

**Symptom:** DE results dominated by transcripts with global expression changes (housekeeping shifts).

**Fix:** Always use the interaction-term design. The biologically meaningful test is the interaction `type:condition` p-value.

### Normalization assumptions

**Trigger:** edgeR `calcNormFactors(method='TMM')` on CLIP-seq data.

**Mechanism:** TMM assumes most features (genes) are not differentially expressed. CLIP-seq peak counts can be globally shifted if the RBP itself is knocked down; TMM normalization would force the shift to be invisible.

**Symptom:** Knockdown experiment shows ~0 DE peaks; expected hundreds.

**Fix:** Use SMInput as the control library; spike-in normalization if available; or skip TMM and use library-size normalization only. Some CLIP-specific tools (DEWSeq, Flipper) handle this internally.

### Flipper requires Skipper upstream

**Trigger:** Flipper called on CLIPper output.

**Mechanism:** Flipper expects Skipper's window-level output format with beta-binomial estimates.

**Symptom:** Flipper crashes or produces nonsense.

**Fix:** Use DEWSeq with htseq-clip for CLIPper-upstream workflows; use Flipper only with Skipper.

## Decision Tree by Scenario

| Scenario | Tool + design | Why |
|----------|---------------|-----|
| KD vs control eCLIP, CLIPper upstream | DEWSeq + htseq-clip | CLIP-specific NB GLM with interaction term |
| KD vs control eCLIP, Skipper upstream | Flipper | Skipper companion; matches windowing |
| Treatment vs vehicle (small effect) | DEWSeq (window-level higher power) | Sliding windows capture small shifts |
| Multiple time points | DEWSeq with time as covariate | Continuous design with time vector |
| Allele-specific differential | BEAPR per-allele + custom logistic | See clip-seq/clip-alignment for WASP |
| Single-cell CLIP differential | MASTR-seq or scCLAM | Specialized; few options |
| Differential motif occupancy | Window-level + DEWSeq + motif overlap | Combine differential windows with motif BED |
| RBP overexpression vs control | Same as KD reversed | Same statistical framework |
| Compare two RBPs | NOT differential CLIP; use SPIDR or separate CLIPs | Different RBPs need separate IPs |
| Spike-in normalization needed | DEWSeq + spike-in size factors | For global occupancy shifts |
| chimeric eCLIP differential miRNA targets | Custom; treat each miRNA-target chimera as feature | Specialized; see clip-seq/ago-clip-mirna-targets |

## Reconciliation: When Differential Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| DEWSeq finds many DE windows; edgeR peak-level finds few | Window-level higher power for narrow shifts | Aggregate DEWSeq windows; cross-check |
| edgeR many DE; DEWSeq few | edgeR not accounting for SMInput | Re-run edgeR with interaction term |
| DE peaks dominated by expression changes | No interaction term | Use `~ type + condition + type:condition` |
| KD experiment shows ~0 DE | TMM over-corrects global shift | Switch to library-size norm only; use SMInput |
| siRNA replicates discordant | Off-target effects vary | Use multiple independent siRNAs; require concordance |
| Treatment-only peaks dominate DE | No consensus peakset | Generate consensus first; then test on unified set |
| Significant windows scattered | Window aggregation step skipped | `bedtools merge -d 100` on significant-window BED |
| Same gene appears in many DE windows | Multiple binding sites per gene differential | Report at gene-level too; not just window-level |
| Flipper fails with non-Skipper input | Upstream mismatch | Use DEWSeq for non-Skipper workflows |
| DESeq2 dispersion fit fails | Too few replicates (n=2 per condition); too few features after filtering | Increase replicates; or relax filtering |

**Operational rule for high-confidence differential reporting:** (a) Use SMInput-aware design (`~ type + condition + type:condition`); (b) generate consensus peakset across conditions; (c) require padj < 0.05 AND |log2FC| > 1; (d) require concordance with at least one orthogonal method (e.g., DEWSeq + edgeR peak-level on same data); (e) for KD experiments, validate KD efficiency by WB and require multiple independent siRNAs.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| DESeq2 `~ condition` finds many DE | No interaction with type | Use `~ type + condition + type:condition` and test interaction |
| DEWSeq output is gene-level not window-level | `keep <- ...` filter too aggressive | Loosen pre-filter |
| Adjacent significant windows not aggregated | Forgot bedtools merge | `bedtools merge -d 100` on sig windows |
| edgeR TMM fits all libraries to one value | Most-features-not-DE assumption violated | Use library-size norm; or DEWSeq for CLIP-specific |
| Flipper crashes on CLIPper input | Tool mismatch | Switch to DEWSeq |
| Few replicates -> unstable estimates | n=2 not enough for dispersion | Increase n; or use limma-voom (more tolerant) |
| Peaks differ across conditions | Per-condition peak calls | Unify with consensus peakset |
| KD experiment yields 0 DE | Normalization over-corrected; OR KD efficiency too low | Validate KD WB; check normalization |
| Global shift in IP relative to SMInput | RBP itself knocked down so IP yield lower | Normalize WITHIN each condition |
| Lots of "treatment-only" peaks | Caller stringency higher in one condition | Use consensus peakset for fairness |

## References

- Schwarzl T et al 2020 Bioconductor pkg DEWSeq (windowed CLIP differential)
- Sahadevan S et al 2022 Methods Mol Biol 2404:189 (DEWSeq + htseq-clip pipeline)
- Tu Y et al 2024 PMC PMC13060198 (Flipper, Skipper-companion differential)
- Boyle EA et al 2023 Cell Genomics 3:100317 (Skipper, parent of Flipper)
- Love MI et al 2014 Genome Biol 15:550 (DESeq2)
- McCarthy DJ et al 2012 Nucleic Acids Res 40:4288 (edgeR)
- Ritchie ME et al 2015 Nucleic Acids Res 43:e47 (limma)
- Wu B et al 2018 Nat Commun 9:5117 (BEAPR allele-specific differential)
- Van Nostrand EL et al 2020 Nature 583:711 (ENCODE 150 RBP shRNA + eCLIP comparison)

## Related Skills

- clip-seq/clip-peak-calling - CLIPper / Skipper outputs feed differential
- clip-seq/clip-qc - Replicate QC required for valid differential
- clip-seq/binding-site-annotation - Annotate differential regions
- clip-seq/clip-motif-analysis - Motif analysis on differential windows
- clip-seq/ago-clip-mirna-targets - Differential miRNA targeting from chimeric eCLIP
- differential-expression/deseq2-basics - Underlying NB GLM model
- differential-expression/de-results - DE results interpretation
- differential-expression/edger-basics - edgeR for peak counts
- chip-seq/differential-binding - DNA-protein analogue
