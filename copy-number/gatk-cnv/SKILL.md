---
name: bio-copy-number-gatk-cnv
description: Call copy number variants with the GATK best-practices workflows — the somatic CNV pipeline (CollectReadCounts, DenoiseReadCounts with tangent normalization, ModelSegments, CallCopyRatioSegments) and the germline GATK-gCNV pipeline (DetermineGermlineContigPloidy, GermlineCNVCaller cohort/case mode, PostprocessGermlineCNVCalls). Covers panel-of-normals construction, AnnotateIntervals/FilterIntervals, allelic-count integration, and QS-based filtering. Use when integrating CNV calling into a GATK variant pipeline, calling rare germline CNVs from an exome cohort, deciding between the somatic and germline GATK workflows, or diagnosing why tangent normalization removed a real event or why gCNV output has low precision.
tool_type: cli
primary_tool: gatk
---

## Version Compatibility

Reference examples tested with: GATK 4.5+ (gatk4), Python 3.10+ (gcnv conda env), R 4.3+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `gatk --version` then `gatk <ToolName> --help` to confirm arguments
- gCNV requires a working `gatkcondaenv` (theano/tensorflow stack) — `gatk` will report if the Python environment is missing

GATK 4.5+ gCNV inference defaults are tuned for whole-exome data; whole-genome runs generally need parameter changes. If a tool reports an unrecognized argument, check the help for that exact GATK version rather than retrying.

# GATK CNV Workflows

**"Call CNVs the GATK way"** -> GATK has two *separate* CNV workflows that share almost no tools. Picking the wrong one is the most common mistake.

- Somatic CNV: `CollectReadCounts` -> `DenoiseReadCounts` -> `ModelSegments` -> `CallCopyRatioSegments`. Tumor copy-ratio segments, optionally allele-aware.
- Germline gCNV: `DetermineGermlineContigPloidy` -> `GermlineCNVCaller` -> `PostprocessGermlineCNVCalls`. Per-sample germline CN genotypes (VCF).

## Critical: What GATK Somatic CNV Does NOT Provide

`ModelSegments` + `CallCopyRatioSegments` produce **copy-ratio segments** and a **minor-allele fraction** per segment, and the "call" is a simple t-test emitting `+` / `-` / `0`. This is **not** integer allele-specific copy number, **not** tumor purity, and **not** ploidy. Practitioners routinely assume parity with ASCAT/FACETS and there is none. For integer allele-specific CN, purity, ploidy, LOH state, or whole-genome-doubling status, use allele-specific-copy-number (ASCAT, Sequenza, FACETS, or PureCN — PureCN can even reuse the GATK `ModelSegments` segmentation as input).

## Somatic vs Germline — Choosing the Workflow

| Question | Somatic CNV | Germline gCNV |
|----------|-------------|---------------|
| Input | One tumor (+ optional matched normal) | A cohort of constitutional samples |
| Output | Copy-ratio segments, +/-/0 call, minor-allele fraction | Integer germline CN genotype VCF per sample |
| Normalization | Tangent (projection onto PoN subspace) | PCA batching + Bayesian read-depth model |
| Cohort needed | PoN of normals for denoising | >= ~100 technically matched samples (cohort mode) |
| Use for | Tumor SCNAs, focal amplifications/deletions | Rare/de novo germline CNVs, NDD/Mendelian cohorts |

## Decision Tree by Scenario

| Scenario | Workflow | Key parameters |
|----------|----------|----------------|
| Tumor-normal WGS/WES, want SCNAs | Somatic, with matched-normal allelic counts | `PreprocessIntervals --bin-length 1000` (WGS) or `0` (WES) |
| Tumor-only somatic CNV | Somatic, no matched-normal allelic counts | Genotype hets in the case sample; expect more no-calls |
| Rare germline CNV, exome cohort >= 100 | gCNV cohort mode | Run `DetermineGermlineContigPloidy` cohort first |
| New sample vs an existing gCNV model | gCNV case mode | Must reuse identical scatter count and interval list |
| Need integer ASCN / purity / ploidy | Neither — escalate | Use allele-specific-copy-number |
| Targeted panel (< few hundred genes) | Prefer CNVkit | GATK interval models are unstable on tiny panels |

## Somatic CNV Pipeline

```bash
# 1. Preprocess and annotate intervals (WES: bin-length 0 = use exome targets as-is)
gatk PreprocessIntervals -R ref.fa -L targets.interval_list \
    --bin-length 0 --interval-merging-rule OVERLAPPING_ONLY -O preprocessed.interval_list
gatk AnnotateIntervals -R ref.fa -L preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY -O annotated.tsv     # GC content for FilterIntervals

# 2. Collect read counts (each BAM)
gatk CollectReadCounts -R ref.fa -I sample.bam -L preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY -O sample.counts.hdf5

# 3. Build the panel of normals (tangent-normalization basis)
# --minimum-interval-median-percentile 5.0 is the GATK CNV tutorial value (tool default 10.0)
gatk CreateReadCountPanelOfNormals \
    -I normal1.counts.hdf5 -I normal2.counts.hdf5 -I normalN.counts.hdf5 \
    --annotated-intervals annotated.tsv \
    --minimum-interval-median-percentile 5.0 -O cnv.pon.hdf5

# 4. Denoise tumor against the PoN (tangent normalization)
gatk DenoiseReadCounts -I tumor.counts.hdf5 --count-panel-of-normals cnv.pon.hdf5 \
    --standardized-copy-ratios tumor.standardizedCR.tsv \
    --denoised-copy-ratios tumor.denoisedCR.tsv

# 5. Allelic counts at common biallelic SNPs (tumor and matched normal)
gatk CollectAllelicCounts -R ref.fa -I tumor.bam -L common_snps.interval_list \
    -O tumor.allelicCounts.tsv
gatk CollectAllelicCounts -R ref.fa -I normal.bam -L common_snps.interval_list \
    -O normal.allelicCounts.tsv

# 6. Joint segmentation of copy ratio and allele fraction
gatk ModelSegments --denoised-copy-ratios tumor.denoisedCR.tsv \
    --allelic-counts tumor.allelicCounts.tsv \
    --normal-allelic-counts normal.allelicCounts.tsv \
    --output-prefix tumor -O segments/

# 7. Call each segment +/-/0 (simple t-test against the copy-ratio baseline)
gatk CallCopyRatioSegments -I segments/tumor.cr.seg -O segments/tumor.called.seg
```

`AnnotateIntervals` (step 1) and supplying `--annotated-intervals` to the PoN are frequently skipped — they enable explicit GC-bias correction and are recommended.

## Germline gCNV Pipeline

```bash
# 1. Determine contig ploidy across the cohort (karyotype + global depth)
gatk DetermineGermlineContigPloidy -L preprocessed.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -I sample1.counts.hdf5 -I sampleN.counts.hdf5 \
    --contig-ploidy-priors ploidy_priors.tsv --output-prefix cohort -O ploidy-calls/

# 2. FilterIntervals — remove low-mappability / extreme-GC / low-count intervals
gatk FilterIntervals -L preprocessed.interval_list --annotated-intervals annotated.tsv \
    -I sample1.counts.hdf5 -I sampleN.counts.hdf5 \
    --interval-merging-rule OVERLAPPING_ONLY -O filtered.interval_list

# 3. GermlineCNVCaller, cohort mode (builds the model AND calls the cohort)
gatk GermlineCNVCaller --run-mode COHORT -L filtered.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    --contig-ploidy-calls ploidy-calls/cohort-calls \
    -I sample1.counts.hdf5 -I sampleN.counts.hdf5 \
    --output-prefix cohort -O gcnv-calls/

# 4. Post-process per sample into a genotyped VCF
gatk PostprocessGermlineCNVCalls \
    --calls-shard-path gcnv-calls/cohort-calls \
    --model-shard-path gcnv-calls/cohort-model \
    --contig-ploidy-calls ploidy-calls/cohort-calls \
    --sample-index 0 \
    --output-genotyped-intervals sample0.intervals.vcf.gz \
    --output-genotyped-segments sample0.segments.vcf.gz \
    --output-denoised-copy-ratios sample0.denoisedCR.tsv
```

Case mode (`--run-mode CASE`) scores a new sample against the cohort `*-model` shards; it must use the **identical** `filtered.interval_list` and the **same scatter count** as the cohort run, or it fails or produces incomparable calls.

## Failure Modes

### Tangent normalization removes a real CNV

**Trigger:** PoN is small (< ~20 normals) or contains samples that share a recurrent CNV (e.g. a common germline CNV, or a PoN accidentally built from tumors).

**Mechanism:** `DenoiseReadCounts` projects the tumor coverage profile onto the subspace spanned by the PoN's principal components. Any copy-number pattern present in that subspace is treated as "systematic noise" and subtracted. A CNV shared by PoN members is therefore normalized out of the tumor.

**Symptom:** A known event (recurrent amplification/deletion, or a common germline CNV) is absent from `denoisedCR.tsv`; denoised profile is suspiciously flat at that locus.

**Fix:** Build the PoN from >= 20-40 unrelated, tumor-free, process-matched normals. Never put tumors in the PoN. Cross-check against the `standardizedCR.tsv` (pre-tangent) profile — if the event is there but gone after denoising, the PoN ate it.

### Mistaking ModelSegments output for allele-specific integer CN

**Trigger:** Treating `tumor.modelFinal.seg` minor-allele fraction as integer minor copy number, or expecting a purity/ploidy field.

**Mechanism:** GATK somatic CNV models copy ratio and allele fraction but never fits the purity/ploidy grid that converts log-ratio to integer absolute CN.

**Symptom:** No purity/ploidy in any output; "copy number" is continuous log2; LOH is a low minor-allele fraction, not an explicit CN-LOH state.

**Fix:** Accept GATK somatic CNV as a relative caller. For integer ASCN, feed the data to PureCN (`segmentationGATK4` reuses GATK segments), FACETS, ASCAT, or Sequenza — see allele-specific-copy-number.

### FilterIntervals silently drops intervals containing real variants

**Trigger:** Aggressive mappability or segmental-duplication cutoffs in `FilterIntervals`.

**Mechanism:** A minimum mappability > 0 or a maximum segmental-duplication content < 1 excludes intervals overlapping segdups and low-mappability regions — exactly where many disease-relevant CNVs (e.g. recurrent genomic-disorder loci flanked by segdups) live.

**Symptom:** Known recurrent CNVs at segdup-mediated loci are never called; the gene of interest has no intervals in `filtered.interval_list`.

**Fix:** Inspect `filtered.interval_list` for genes of interest before calling. Relax mappability/segdup cutoffs for targeted analyses; for genomic-disorder loci, depth-based callers are inherently limited near segdups — confirm with an orthogonal assay.

### Raw gCNV output has ~22% precision

**Trigger:** Using unfiltered `GermlineCNVCaller` / `PostprocessGermlineCNVCalls` output for association or de novo analysis.

**Mechanism:** gCNV is tuned for high recall (~95% of rare coding CNVs >= 2 exons) at the cost of precision; raw calls are dominated by false positives.

**Symptom:** Implausibly many rare CNVs per sample; de novo CNV rate far above the expected ~0.01-0.02/genome.

**Fix:** Apply the QS (quality score) filter. QS > 100 is a common starting threshold; QS > 1000 reaches ~96% precision. Also apply sample-level filters (call rate, number of CNVs per sample) per Babadi 2023.

### gCNV cohort too small or mismatched

**Trigger:** Cohort mode with < ~100 samples, or a cohort spanning multiple capture kits / library protocols.

**Mechanism:** The Bayesian model needs enough technically similar samples to learn coverage bias; mixed protocols are not separable and the model misattributes batch effects to copy number.

**Symptom:** Unstable calls; many CNVs tracking sequencing batch; model fails to converge.

**Fix:** Use >= 100 process-matched samples per cohort model; split heterogeneous cohorts by capture kit. For < 100 samples, gCNV case mode against an external compatible model, or a cohort caller like ExomeDepth, is more appropriate — see germline-cnv-interpretation.

## Reconciliation: GATK vs Other Callers

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| GATK somatic flat where CNVkit calls a focal event | Tangent normalization absorbed it | Check `standardizedCR.tsv`; rebuild PoN without the event |
| GATK and ASCAT disagree on a "deletion" | GATK has no purity model; the event is subclonal or impure | Trust ASCAT/FACETS integer ASCN |
| gCNV calls a CNV ExomeDepth misses | Different sensitivity profiles; both have poor inter-tool concordance | Require QS filtering + a second caller for rare-CNV claims |
| gCNV CNV count tracks batch | Cohort mixes protocols | Re-batch by capture kit |

**Operational rule:** GATK somatic CNV output is relative copy ratio — report it as gain/loss/neutral, not absolute CN, unless downstream-fit by an allele-specific tool. gCNV calls are reportable only after QS and sample-level filtering, and rare-CNV or de novo claims need orthogonal confirmation.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| PoN size (somatic) | >= 20-40 normals | Larger PoN = stabler tangent subspace; small PoN over-fits |
| gCNV cohort size | >= ~100 technically matched | Babadi 2023 Nat Genet; model needs coverage-bias signal |
| gCNV QS for high precision | QS > 1000 -> ~96% precision | Babadi 2023; raw output ~22% precision, ~95% recall |
| `minimum-interval-median-percentile` | 10.0 default; 5.0 in the GATK CNV tutorial | Drops the lowest-coverage intervals from the PoN |
| WGS bin length | ~1000 bp | `PreprocessIntervals --bin-length 1000`; WES uses 0 (targets as-is) |
| Het sites for somatic ModelSegments | >= ~10,000 (WGS) | Sparse hets give noisy minor-allele-fraction segmentation |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `gatkcondaenv` / theano error in gCNV | gCNV Python env not installed | Install the GATK conda env; gCNV cannot run without it |
| "At least one interval must remain" in FilterIntervals | All intervals filtered out | Relax mappability/GC/count cutoffs; check annotated.tsv |
| Case-mode gCNV fails or gives odd calls | Different scatter count or interval list vs cohort model | Reuse the exact cohort scatter count and `filtered.interval_list` |
| Denoised profile flat at a known event | Tangent normalization removed it | Rebuild PoN larger, tumor-free; inspect standardizedCR |
| ModelSegments has few het sites | SNP interval list misses captured regions | Use a common-SNP list intersected with the capture targets |
| Expecting purity/ploidy in output | Somatic CNV does not estimate them | Use allele-specific-copy-number |

## References

- GATK Best Practices: Somatic copy number variant discovery (CNV). Broad Institute documentation.
- Babadi M et al 2023. GATK-gCNV enables the discovery of rare copy number variants from exome sequencing data. Nat Genet 55:1589
- Gao GF, Oh S, Saksena G, Tabak B, Beroukhim R, Getz G et al 2022. Tangent normalization for somatic copy-number inference in cancer genome analysis. Bioinformatics 38:4677 (Tabak is a middle, not first, author).

## Related Skills

- copy-number/allele-specific-copy-number - Integer ASCN, purity, ploidy (ASCAT/Sequenza/FACETS/PureCN)
- copy-number/copy-ratio-segmentation - Segmentation algorithms and depth normalization theory
- copy-number/cnvkit-analysis - Read-depth CNV calling for panels and exomes
- copy-number/germline-cnv-interpretation - ACMG/ClinGen classification of germline CNV calls
- copy-number/cnv-visualization - Plotting GATK denoised ratios and modeled segments
- copy-number/recurrent-cnv - Cohort-level recurrent and driver CNV
- variant-calling/gatk-variant-calling - GATK SNV/indel pipeline for allelic-count SNP sites
