# Spike-In Normalization - Usage Guide

## Overview

Normalize ChIP-seq data using exogenous spike-in chromatin: Drosophila chromatin (ChIP-Rx; Orlando 2014 / Egan 2016) for human/mouse ChIP, or E. coli carryover from bacterial pA-MNase/pA-Tn5 for CUT&RUN/CUT&Tag. Required when global signal shifts are expected (HDACi, BETi, EZH2i, target knockdown, dosage). Covers RRPM vs Rx-Input scaling formulas, integration with DiffBind/DESeq2/edgeR/csaw, ChIPseqSpikeInFree for post-hoc detection without spike-in, and the Hammond Norris 2024 review's failure-mode framework. Emphasizes the internal-control sanity check (blacklist regions should NOT shift post-scaling).

## Prerequisites

```bash
conda install -c bioconda samtools bowtie2 bedtools deeptools
```

```r
BiocManager::install(c('DiffBind', 'DESeq2', 'edgeR', 'csaw', 'ChIPseqSpikeInFree'))
```

```bash
# SpikeFlow (Snakemake wrapper, optional)
# https://github.com/sebastian-gregoricchio/SpikeFlow
```

## Quick Start

Tell the agent what to do:
- "I have H3K27ac ChIP-seq from DMSO vs JQ1 (BET inhibitor). Compute Drosophila spike-in scaling factors and apply via DiffBind."
- "Align reads to combined hg38 + dm6 genome and count Drosophila reads per sample"
- "Apply spike-in scaling at the read level via DESeq2 sizeFactors (inverse convention)"
- "Generate spike-in-scaled bigWigs without `--normalizeUsing` (mutually exclusive with `--scaleFactor`)"
- "Verify internal-control sanity: blacklist regions should show no signal change post-scaling"
- "I didn't add spike-in. Run ChIPseqSpikeInFree for post-hoc detection of global shifts (diagnostic only)."
- "Compute Rx-Input scaling factor incorporating spike-in fractions in BOTH ChIP and input"

## Example Prompts

### Standard ChIP-Rx workflow
> "I have 4 samples (2 DMSO, 2 JQ1) for H3K27ac. Build a combined hg38+dm6 bowtie2 index, align reads, filter with -F 1804 -q 30, count Drosophila reads per sample, compute scale factors as min/each, and apply to DESeq2 as `sizeFactors(dds) <- 1/scale_factors`."

### DiffBind integration
> "Use DiffBind 3.20+ with `dba.normalize(obj, spikein = TRUE)`. Sample sheet has `Spikein` column with Drosophila read counts."

### Spike-in scaled tracks
> "Generate spike-in-scaled bigWigs with `bamCoverage --scaleFactor X`. Do NOT also pass `--normalizeUsing`."

### Internal-control validation
> "After applying spike-in scaling, verify ENCODE blacklist regions show no signal change between conditions. Also check signal at U6 snRNA promoter (constitutively bound by Pol III; should be stable across H3K27ac conditions)."

### Post-hoc detection
> "I didn't run spike-in but suspect a global shift in my BET-inhibitor experiment. Run ChIPseqSpikeInFree to estimate scaling factors from signal-distribution shape."

### Rx-Input variant
> "Apply Rx-Input scaling (Fursova 2019) which scales by both ChIP and input spike-in fractions to correct for IP efficiency variation."

### Failure-mode debugging
> "My results have wrong signs after scaling. Check inverse convention; check whether scaling was applied to peak counts instead of read counts; verify internal-control regions."

## What the Agent Will Do

1. **Determine spike-in protocol**: Drosophila ChIP-Rx (deliberate) or E. coli carryover (CUT&RUN/CUT&Tag)
2. **Build combined alignment index**: host + spike genome (avoid chromosome naming collisions; prefix Drosophila chromosomes if needed)
3. **Align reads**: bowtie2 or bwa-mem to combined genome
4. **Filter and deduplicate**: ENCODE filter `-F 1804 -q 30` applied to both target and spike alignments
5. **Count spike reads**: at high mapq, post-dedup; verify counts are in linear range (0.5-5% of total reads)
6. **Compute scaling factors**: RRPM: `min(spike) / per_sample_spike`; Rx-Input: ratio of ChIP/input spike fractions
7. **Apply at correct layer:**
   - bigWig tracks: `bamCoverage --scaleFactor X` (NOT with `--normalizeUsing`)
   - DiffBind: `dba.normalize(obj, spikein = TRUE)`
   - DESeq2: `sizeFactors(dds) <- 1 / scale_factors` (inverse convention)
   - edgeR: `normFactors` via DGEList
   - csaw: pass to `windowCounts` library factors
8. **Internal-control sanity check**: verify blacklist + housekeeping signal stability post-scaling
9. **Differential testing**: run DESeq2 / edgeR / csaw with spike-in-based size factors
10. **Document**: spike-in organism + cell count + protocol, alignment %, mapq filter, scaling formula, internal-control validation

## Tips

- **Spike-in is required for global-shift biology.** HDACi, BETi, EZH2i, dosage, target knockdown all confound reads-in-peaks normalization.
- **Apply at the read level, never to peak counts.** Hammond Norris 2024: ~25% of published papers violate this.
- **DESeq2 / edgeR `sizeFactors` use inverse convention.** `sizeFactors(dds) <- 1 / scale_factors`, not direct.
- **`--normalizeUsing` and `--scaleFactor` are mutually exclusive in bamCoverage.** Use one only.
- **Verify titration linearity.** Spike-in counts should be 0.5-5% of total reads; outside this range, saturation or noise makes scaling unreliable.
- **Filter spike reads to mapq ≥ 30 before counting.** Low-mapq reads at low-complexity regions inflate counts spuriously.
- **Deduplicate spike reads before counting.** PCR duplicates of spike-in don't reflect chromatin amount.
- **Use a single enzyme lot for CUT&RUN/CUT&Tag cross-condition comparison.** E. coli carryover varies between bacterial batches.
- **For high-stakes CUT&Tag cross-condition claims, supplement with deliberate Drosophila spike-in.**
- **ChIPseqSpikeInFree is diagnostic only, not publication-grade.** Re-do experiment with deliberate spike-in if publication is goal.
- **The internal-control sanity check is mandatory.** Blacklist regions (no biology) should show NO signal change after spike-in scaling. If they do, scaling is broken.

## Troubleshooting

### Scaling factors all near 1.0

Spike-in not added at fixed amount, OR samples not properly randomized. Re-check Egan 2016 protocol (50k Drosophila nuclei per 5M target cells).

### Spike-in fraction <0.1%

Spike-in lost during library prep, OR too few cells used. Scaling unreliable; re-do experiment.

### Spike-in fraction >5%

Spike-in too concentrated; saturating. Titration not in linear range; scaling factors unreliable.

### Replicate-to-replicate scaling factors vary >2×

Library prep batch effects on spike-in fraction. Re-extract spike reads after deduplication and mapq filter; if still variable, library prep needs improvement.

### Results sign-flipped after scaling

Inverse convention bug. `sizeFactors(dds) <- 1 / scale_factors`, not `sizeFactors(dds) <- scale_factors`.

### Internal-control regions shift after scaling

Scaling broken. Common causes:
1. Scaling applied to peak counts → fix application layer
2. Spike-in reads not dedup'd → re-extract from filtered BAM
3. mapq filter too loose → use `-q 30`
4. Spike-in saturated → verify titration linearity

### Cross-condition comparison still shows no effect

If spike-in is properly applied and internal controls are stable, the biology may not have the global shift. Consider:
1. Mark biology; H3K27ac may shift less than H3K27me3 with EZH2i
2. Time-course; acute (1 hr) vs chronic (24 hr) differs
3. Inhibitor concentration; sub-saturating doesn't produce global shift

### ChIPseqSpikeInFree predicts shift but spike-in says no

Trust spike-in measurement. ChIPseqSpikeInFree is a distribution-shape heuristic; it may misinterpret biology as a global shift.

## Related Skills

- chip-seq/peak-calling - Peak calling upstream
- chip-seq/chipseq-qc - Spike-in fraction QC + library complexity
- chip-seq/differential-binding - Apply spike-in scaling via DiffBind / DESeq2 / csaw
- chip-seq/cut-and-run-tag - E. coli spike-in carryover specifics
- chip-seq/super-enhancers - SE calling requires spike-in for cross-condition
- chip-seq/chipseq-visualization - Spike-in-scaled bigWig generation
- alignment-files/sam-bam-basics - Combined-genome alignment
- differential-expression/deseq2-basics - DESeq2 sizeFactors conventions
