---
name: bio-atac-seq-footprinting
description: Detect transcription factor binding footprints in ATAC-seq using TOBIAS, HINT-ATAC, Wellington, or scprinter. Use when identifying bound TF sites within accessible regions, correcting Tn5 insertion bias before footprinting, choosing between cleavage-based and aggregate-based footprinters, or comparing differential TF activity between conditions.
tool_type: cli
primary_tool: tobias
---

## Version Compatibility

Reference examples tested with: TOBIAS 0.16+, RGT HINT-ATAC 1.0.2+, Wellington (pyDNase) 0.3+, scprinter 0.1+, samtools 1.19+, bedtools 2.31+, pyBigWig 0.3+, MEME suite 5.5+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# TF Footprinting

**"Identify TF binding footprints in my ATAC-seq data"** -> Detect short DNA stretches (typically 6-20 bp) of reduced Tn5 cleavage within accessible regions, where a bound TF physically protects DNA. Requires (1) Tn5 sequence-bias correction, (2) per-base footprint scoring, (3) motif-anchored detection.

- CLI: `TOBIAS ATACorrect` -> `TOBIAS ScoreBigwig` (formerly `FootprintScores`) -> `TOBIAS BINDetect`
- CLI: `rgt-hint footprinting --atac-seq` (HINT-ATAC, single-step)
- CLI: `wellington_footprints.py` (legacy DNase, adapted for ATAC)
- Python: `scprinter` (multi-scale, single-cell aware; Bao Yu 2024 bioRxiv)

Tn5 has a strong sequence preference (Lazarovici 2013, Calviello 2019), reading approximately +/- 4 bp around the insertion site. Without bias correction, "footprints" reflect Tn5 sequence preference rather than TF binding. This is the single most important step.

## Algorithmic Taxonomy

| Tool | Bias model | Scoring | Min depth | Strength | Fails when |
|------|-----------|---------|-----------|----------|------------|
| TOBIAS (BINDetect) | k-mer (default 9bp) PWM-style; per-base correction | Two-step: continuous footprint score then motif-anchored bound/unbound classification | >= 50M nuclear reads | Mature, peer-reviewed (Bentsen 2020), differential support, modular pipeline | Below 50M reads; sequencing errors near motif inflate background |
| HINT-ATAC | Hidden-Markov + dinucleotide bias correction | HMM emits open/footprint/closed states; calls ranked footprints | >= 50M | Single-step; integrates motif matching; handles DNase too | Less control over individual stages; HMM occasionally over-segments |
| Wellington (pyDNase) | DNase-developed; ATAC adaptation by post-shift | Cleavage-rate Poisson Z-score | >= 50M (DNase >= 80M) | Original footprinting framework; well-validated for DNase | Designed for DNase II; ATAC-specific bias not corrected as carefully |
| PIQ | Bayesian latent variable on cut sites | Genome-wide PWM scan + cleavage profile | >= 30M (lower because of model) | Per-TF posterior probabilities; works on lower depth | Outdated; not actively maintained; harder to install |
| scprinter | Multi-scale CNN-based footprint and TF activity | Resolves footprints at multiple TF size scales (CTCF vs nuclear receptors) | >= 1M cells (sc) or 50M (bulk) | Modern ML approach; single-cell; multi-scale resolves problematic TF families | Newer tool; benchmarks evolving; GPU recommended |
| TOBIAS + scprinter combination | TOBIAS bias correction + scprinter scoring | Two-step bridging | >= 50M | Combines the best bias model with multi-scale scoring | Manual pipeline, no single CLI |

Methodology evolves; verify against the current Bentsen 2020, Calviello 2019, and scprinter 2024 benchmarks. ATAC footprinting power saturates above 100M nuclear reads; below 50M, weak-binding TFs (transient occupancy) cannot be reliably called.

## Tn5 Bias and Why Correction Matters

**Trigger:** Tn5 inserts preferentially at certain k-mers (Lazarovici 2013 reported strong A/T preference at -3 to -1 and +1 to +3; Calviello 2019 measured the full PWM). The preference is reproducible and biologically uninteresting.

**Mechanism:** Without correction, every "TF footprint" near a high-bias k-mer reads as occupancy. Conversely, regions with low-bias flanks but real binding may show no footprint dip relative to corrected expectation.

**Symptom:** Aggregate footprint at random GC-rich motifs shows V-shape; aggregate at AT-rich motifs shows inverse-V (peak instead of dip).

**Fix:** Apply ATACorrect (TOBIAS), seqOutBias (Martins 2018), or HINT's dinucleotide model. Bias correction subtracts the Tn5-expected per-base profile from observed cleavage. After correction, V-shape is preserved only at TF-bound sites.

**Goal:** Subtract Tn5 sequence-bias from the per-base cleavage signal so residual footprints reflect TF binding rather than enzyme preference.

**Approach:** Run TOBIAS ATACorrect over the deduplicated BAM with the reference genome, consensus peaks, and ENCODE blacklist; it emits per-condition uncorrected, bias, expected, and corrected bigWigs for downstream scoring.

```bash
# TOBIAS ATACorrect: produces uncorrected, bias, expected, and corrected bigWigs
TOBIAS ATACorrect \
    --bam sample.dedup.bam \
    --genome hg38.fa \
    --peaks consensus_peaks.bed \
    --blacklist hg38-blacklist.v2.bed \
    --outdir corrected/ \
    --cores 8
# Output: sample_uncorrected.bw, sample_bias.bw, sample_expected.bw, sample_corrected.bw
```

## Tn5 Cut Geometry: +4 / -5 Dual-Cut

Tn5 dimers cut both strands of DNA but with a 9 bp staggered offset. The cleavage event creates two free 5' ends: one shifted +4 bp from the binding center on the forward strand and -5 bp on the reverse strand. Footprinting tools must apply this shift before per-base counting:

| Strand | Read 5' end correction |
|--------|------------------------|
| + strand | shift +4 bp downstream |
| - strand | shift -5 bp upstream |

**Trigger:** Computing per-base Tn5 cut signal manually.

**Symptom:** Footprint aggregates show ~9 bp asymmetry (apex shifted from motif center).

**Fix:** Apply +4/-5 shift before counting; TOBIAS, HINT-ATAC, and scprinter handle this internally. Custom analyses must apply explicitly. deepTools provides this via `alignmentSieve --ATACshift` (which applies the canonical +4 / -5 shift in one step).

## Bias Correction Alternatives

| Method | Approach | When to use |
|--------|---------|-------------|
| TOBIAS ATACorrect | 9-bp k-mer PWM | Default for most ATAC; fast |
| chromBPNet bias model (Pampari 2024) | CNN trained on naked-DNA control or k-mer baseline | Best when sequence context complex; handles low-complexity flanks |
| seqOutBias (Martins 2018) | Genome-wide naive Bayes on k-mer counts | Independent of footprinting tool; works upstream |
| HINT-ATAC dinucleotide | HMM-integrated dinucleotide bias | Built into HINT pipeline; less control |
| Naked-DNA empirical | Sequence Tn5 on protein-free DNA | Gold standard for non-model organisms; expensive wet-lab |

For non-model organisms (no published Tn5 bias model), naked-DNA control is required. chromBPNet's bias model is the modern standard for human/mouse and outperforms TOBIAS at low-complexity sequence contexts (Pampari 2024). See atac-seq/deep-learning-atac.

## In Silico Variant Effect at Footprinted TF Motifs

**Trigger:** A GWAS-fine-mapped or rare variant falls inside a TOBIAS-bound motif site.

**Mechanism:** Sequence-based DL models (chromBPNet, EnFormer) predict per-base accessibility at ref vs alt allele; combined with footprint evidence (TOBIAS bound site overlap), this produces a mechanistic hypothesis: "variant disrupts binding of TF X at enhancer Y."

**Workflow:** Run TOBIAS BINDetect to identify bound motif sites; for variants in bound sites, score with chromBPNet (atac-seq/deep-learning-atac) for ref/alt log2FC; |log2FC| > 1 supports functional disruption. Cross-reference with allele-specific accessibility (atac-seq/allele-specific-accessibility) for observed evidence.

## Per-TF Footprinting Failure Modes

Different TF families produce different footprint signatures. The same tool can report a clean V-shape for one family and noise for another.

### CTCF -- The gold standard

**Trigger:** Footprinting CTCF.

**Mechanism:** CTCF has high ChIP-seq concordance, deep V-shaped footprint (~19-20 bp protected), strong sequence specificity. ChIP-seq overlap is typically >70%.

**Symptom:** Aggregate corrected footprint shows clean ~20 bp dip with bilateral cleavage shoulders.

**Verification:** Always validate footprinting output by checking CTCF first; if CTCF footprint is shallow, the bias correction or depth is the problem, not the biology.

### Nuclear receptors (ER, AR, GR) -- Transient binding

**Trigger:** Glucocorticoid response, hormone-stimulated systems.

**Mechanism:** Steroid receptors bind transiently (residence time minutes vs hours for CTCF); average ATAC sample captures binding probability < 30% per allele.

**Symptom:** Aggregate footprint is shallow or absent despite ChIP-seq peaks at the same sites.

**Fix:** Use scprinter's multi-scale model OR limit to ChIP-validated sites OR pool replicates for higher effective depth. Do not interpret absence of footprint as absence of binding.

### Pioneer TFs (FOXA1, GATA, OCT4) -- Half-site footprint

**Trigger:** Pioneer-factor binding to nucleosomal DNA.

**Mechanism:** Pioneer factors bind one DNA face; the back face is on the histone octamer. Footprint is asymmetric (one side protected, other side accessible).

**Symptom:** Aggregate plot shows asymmetric V; one shoulder is taller than the other.

**Fix:** Use single-stranded scoring; HINT-ATAC has stranded mode. Treat asymmetric footprints as biologically meaningful, not artefactual.

### AP-1 family (FOS, JUN) -- Heterodimer composite footprint

**Trigger:** AP-1 enrichment; the JASPAR motif is composite of multiple heterodimer combinations.

**Mechanism:** Different AP-1 dimers (FOS+JUN, FOS+JUNB, JUNB+JUNB) bind slightly different motifs. JASPAR entries are degenerate; scoring averages over all.

**Fix:** Use specific HOCOMOCO motifs per heterodimer when distinguishing matters. Otherwise accept the composite call.

### ZBTB family / BTB-zinc finger -- Dynamic / unfootprintable

**Trigger:** Footprinting ZBTB16, BCL6, others.

**Mechanism:** Dynamic binding kinetics + cofactor-mediated stabilization mean steady-state occupancy is highly variable. Some ZBTBs simply do not produce reliable ATAC footprints despite genuine ChIP-seq binding.

**Fix:** Document the failure; use ChIP-seq for these TFs. Footprinting cannot rescue everything.

### Forkhead / homeobox (FOX, HOX) -- Short footprint < 8 bp

**Trigger:** Short-motif TFs.

**Mechanism:** Footprint extent matches motif length; <8 bp footprints are at the resolution limit of Tn5 (which has ~4 bp positional uncertainty).

**Fix:** Multi-scale scoring (scprinter); aggregate over thousands of sites; do not rely on per-site calls for short motifs.

## Decision Tree by Goal

| Goal | Recommended pipeline |
|------|---------------------|
| Identify all TFs differentially bound between two conditions | TOBIAS ATACorrect (per condition) -> ScoreBigwig -> BINDetect with `--cond_names` |
| Find the strongest single-TF binding (e.g., CTCF) | TOBIAS PlotAggregate over JASPAR CTCF motif sites; verify V-shape |
| Per-cell footprinting (scATAC) | scprinter (single-cell mode); avoid TOBIAS unless pseudobulking by cluster |
| Multi-scale TF activity (handle short and long simultaneously) | scprinter or TOBIAS + custom multi-scale | 
| Differential nuclear-receptor binding | TOBIAS pooled-replicate footprints + ChIP cross-validation; raw ATAC alone often misses transient binding |
| Plant / non-model organism | TOBIAS or HINT-ATAC with custom motifs; bias model retrained from genomic background |

## TOBIAS Three-Step Pipeline

**Goal:** Call bound/unbound TF motif sites per condition and detect differential occupancy across two conditions.

**Approach:** Run ATACorrect to subtract Tn5 bias from cleavage counts, ScoreBigwig to compute a continuous per-base footprint score, then BINDetect to anchor footprints to motif positions and produce per-TF differential bound calls with p-values.

```bash
# Step 1: Bias correction
TOBIAS ATACorrect \
    --bam cond1.bam --genome hg38.fa \
    --peaks consensus.bed --blacklist hg38-blacklist.v2.bed \
    --outdir cond1_corrected/ --cores 16

# Step 2: Per-base footprint scoring (continuous)
TOBIAS ScoreBigwig \
    --signal cond1_corrected/cond1_corrected.bw \
    --regions consensus.bed \
    --output cond1_footprints.bw \
    --cores 16

# Step 3: Motif-anchored bound/unbound calls + differential
TOBIAS BINDetect \
    --motifs JASPAR2024_CORE_vertebrates.pfm \
    --signals cond1_footprints.bw cond2_footprints.bw \
    --genome hg38.fa --peaks consensus.bed \
    --outdir bindetect/ \
    --cond_names cond1 cond2 \
    --cores 16
```

BINDetect output columns: `output_prefix`, motif info, condition counts (`cond1_bound`, `cond2_bound`), `cond1_mean_score`, `cond2_mean_score`, `cond1_cond2_change` (differential), `cond1_cond2_pvalue` (one-sided per direction).

## Differential Reading

| BINDetect output | Interpretation |
|------------------|----------------|
| `cond1_cond2_change` > 0, low pvalue | TF more bound in cond1 |
| `cond1_cond2_change` < 0, low pvalue | TF more bound in cond2 |
| Both `cond1_bound` and `cond2_bound` near 0 | Motif present but no footprint either condition; TF likely not active |
| `cond1_bound` >> `cond2_bound` but change small | High dynamic range; differential per-site rather than aggregate |

The differential score is the difference in mean footprint score across motif sites, not a fold-change. Magnitudes around 0.1-0.5 are typical for biologically relevant changes (TOBIAS BINDetect tutorials / Bentsen 2020 examples; no formally published cutoff — calibrate against positive controls in the current dataset).

## Reconciling TOBIAS vs HINT-ATAC

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| TOBIAS calls binding, HINT does not | TOBIAS more sensitive; HINT's HMM filters edge calls | Trust if motif is canonical; suspect for novel/weak motifs |
| HINT calls binding, TOBIAS does not | HINT's HMM occasionally over-segments and reports spurious | Verify by aggregate footprint at the called sites |
| Both call same TF as differential but opposite directions | Different bias correction model; different bound/unbound thresholds | Re-check ATACorrect output; one bias model may be miscalibrated |
| Both flat | Library too shallow; chromatin too closed at motif sites | Pool replicates; consider scprinter multi-scale |

**Operational rule:** For high-confidence reporting, require two-tool concordance (TOBIAS + HINT-ATAC OR TOBIAS + ChIP-seq overlap > 50%). Single-tool calls should be reported as exploratory.

## NFR-Filtering Before Footprinting

**Goal:** Restrict footprinting input to nucleosome-free (sub-100 bp) fragments where TF binding signal lives.

**Approach:** Stream the BAM through awk, keep header lines and fragments whose insert size is between -100 and 100 bp, then re-index.

```bash
# Filter to fragments < 100 bp (NFR) -- TF binding lives here, not on nucleosomes
samtools view -h sample.bam | \
    awk 'substr($0,1,1)=="@" || ($9 > 0 && $9 < 100) || ($9 < 0 && $9 > -100)' | \
    samtools view -b > sample.nfr.bam
samtools index sample.nfr.bam
# Use this NFR BAM as input to TOBIAS ATACorrect
```

Filtering NFR strengthens footprint signal but discards di-nucleosome-borne information. Keep the unfiltered BAM for nucleosome-positioning analysis.

## Motif Database Choice

| Database | Coverage | Format | Notes |
|----------|---------|--------|-------|
| JASPAR 2024 CORE vertebrates | ~1900 motifs (curated, ChIP-validated) | JASPAR PFM, MEME, etc. | Default for vertebrate ATAC |
| HOCOMOCO v12 | ~1400 high-confidence + 600 derived | JASPAR PFM | Best for resolving paralogues; provides per-cell-type variants |
| CIS-BP 2.0 | ~80,000 motifs across 1000+ species | PWM, .meme | Broadest coverage including non-model species |
| MEME-CHIP / homer | Custom from peaks | .meme, .motif | When de novo motif needed |

JASPAR motifs are conservatively curated; HOCOMOCO is comprehensive for human/mouse with quality scores per motif (A/B/C/D); CIS-BP excels for non-model organisms.

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Aggregate footprint inverted (peak instead of dip) | No bias correction; or wrong genome FASTA | Run ATACorrect; verify FASTA matches BAM build |
| BINDetect reports zero bound sites | Default cutoff too stringent; or peakset too narrow | Lower `--prefix` and inspect; verify peaks include where binding expected |
| TOBIAS ATACorrect out of memory | Genome FASTA huge or many cores | Reduce `--cores`; use `samtools faidx` to confirm FASTA index exists |
| Differential score noisy / random | Per-condition bias correction inconsistent | Re-run ATACorrect with identical peakset and blacklist for each |
| Empty motif file warning | JASPAR PFM format mismatch | Use `MEME suite` to convert; TOBIAS expects JASPAR format |
| HINT-ATAC reports many tiny footprints | Default HMM over-segments | Use `--organism` flag explicitly; check `--region-file` is consensus, not raw peaks |
| Wellington crashes on paired-end ATAC | Wellington was DNase-targeted, single-end model | Use TOBIAS instead, or convert paired-end to cuts-only BED |
| Per-site footprint is V-shape but aggregate is flat | Mixing strands; some motifs on - strand | Aggregate function should handle strand; verify input motif strand column |

## References

- Buenrostro JD et al 2013 Nat Methods 10:1213 (ATAC-seq protocol)
- Lazarovici A et al 2013 PNAS 110:6376 (Tn5 insertion bias in ATAC)
- Calviello AK et al 2019 BMC Genomics (Tn5 bias PWM characterization; verify exact volume/pages in current literature search)
- Bentsen M et al 2020 Nat Commun 11:4267 (TOBIAS framework, benchmark)
- Li Z et al 2019 Genome Biol 20:45 (HINT-ATAC)
- Piper J et al 2013 NAR 41:e201 (Wellington / pyDNase)
- Sherwood RI et al 2014 Nat Biotechnol 32:171 (PIQ)
- Bao Y et al 2024 bioRxiv (scprinter; multi-scale single-cell footprinting)
- Martins AL et al 2018 NAR 46:e9 (seqOutBias bias correction alternative)
- Castro-Mondragon JA et al 2022 NAR 50:D165 (JASPAR 2022, recently 2024)
- Vorontsov IE et al 2024 NAR 52:D116 (HOCOMOCO v12)

## Related Skills

- atac-seq/atac-peak-calling - Generate input peakset (NFR-only optional)
- atac-seq/atac-qc - Confirm depth >= 50M before footprinting
- atac-seq/single-cell-atac - scprinter for single-cell footprinting
- atac-seq/motif-deviation - Complementary chromVAR for accessibility variability
- atac-seq/deep-learning-atac - chromBPNet bias correction alternative; in silico variant effect at footprints
- atac-seq/allele-specific-accessibility - Observed allelic imbalance at TF-bound sites
- chip-seq/peak-annotation - Cross-validate footprints with ChIP peaks
- sequence-manipulation/motif-search - Underlying motif scanning patterns
- gene-regulatory-networks/scenic-regulons - Downstream regulatory network inference
