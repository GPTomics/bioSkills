---
name: bio-atac-seq-consensus-peakset
description: Build a differential-ready consensus peakset from per-replicate ATAC-seq peaks using iterative overlap removal, fixed-width re-centering, and majority-rule overlap. Use when generating a stable peak coordinate system for downstream differential accessibility, ML feature engineering, cross-sample comparison, or fixed-width peak counts; covers Corces 2018 iterative overlap (501 bp), DiffBind summit re-centering, and ENCODE consistency rules.
tool_type: mixed
primary_tool: bedtools
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, samtools 1.19+, BEDOPS 2.4.41+, GenomicRanges 1.54+, DiffBind 3.12+, Subread 2.0+ (featureCounts), pybedtools 0.10+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Consensus Peakset Construction

**"Build a single peakset I can count reads against across all my samples"** -> Combine per-replicate or per-condition peak calls into a non-redundant, fixed-width set of regions. The strategy chosen drives FDR calibration, peak-width fairness, and reproducibility downstream.

- CLI: `bedtools merge` (simple union) and `bedtools multiinter -j` (per-sample membership)
- CLI: Corces 2018 iterative overlap removal (custom shell)
- R: `DiffBind::dba.count(summits=250)` (built-in fixed-width)
- Python: `pybedtools` for programmatic merging

The peakset choice is rarely default-correct. Wrong width or wrong overlap rule propagates to every downstream analysis (differential, motif, footprint, ML).

## Why a Consensus Peakset Matters

ATAC peaks vary in width across replicates: same regulatory element might be called 200 bp in rep1 and 800 bp in rep2 because of stochastic Tn5 cuts at edges. Counting reads in different-width intervals confounds peak width with biological signal. A fixed-width consensus avoids this.

For ENCODE-style differential analysis: ALL samples must be counted against the SAME peak coordinates; otherwise the count matrix is non-rectangular and statistical models are misspecified.

## Strategy Taxonomy

| Strategy | Implementation | Width | When to use | Fails when |
|----------|---------------|-------|-------------|------------|
| Naive union | `bedtools merge` of all peaks | Variable, tends wide | Quick exploratory; never for differential | Width inflation drives spurious differential |
| Naive intersection | `bedtools multiinter` requiring all samples | Variable | High-stringency reproducibility | Loses real condition-specific peaks |
| Majority-rule overlap | `multiinter` requiring >= n/2 samples | Variable | Balance; DiffBind default with `minOverlap` | Width still varies; counts are width-biased |
| Iterative overlap removal (Corces 2018) | Sort by significance, greedily keep non-overlapping at fixed width | 501 bp fixed | ML features; cross-study comparison; modern ATAC standard | Loses sub-501bp resolution; overweights high-significance peaks |
| Summit-centered fixed width (DiffBind) | `dba.count(summits=250)` re-centers all peaks on summit +/- 250 bp | 501 bp fixed | Default for DiffBind workflows; integrates with replicate counts | Requires summit info (MACS narrowPeak); broad peaks lose width info |
| IDR-filtered union | Union of IDR-passed peaks across rep pairs | Variable | ENCODE pipeline-compliant; reproducibility-aware | Requires running IDR per pair; computationally heavier |
| Per-condition union, then global union | Each group consensus separately, then merge | Variable | Different cell types / strong condition shift | Same width issues as naive union |
| Width-controlled extension | Extend each peak to median width centered on midpoint | User-set | Quick fixed-width without summit info | Midpoint != summit; can shift biology |

Methodology evolves; verify against current ENCODE 4 ATAC standards (encodeproject.org/atac-seq) and the Corces 2018 iterative overlap algorithm before locking pipelines.

## Iterative Overlap Removal (Corces 2018)

This is the modern standard for fixed-width consensus peaksets used in cross-study comparison and machine-learning feature matrices.

**Algorithm:**
1. Pool all peaks from all samples; re-center each on its summit; extend +/- 250 bp -> 501 bp fixed-width peaks.
2. Sort by significance (narrowPeak column 7, signalValue, descending).
3. Walk down the sorted list. Keep each peak if it does not overlap any previously kept peak. Drop if overlap.
4. Output the kept peaks as the consensus.

**Goal:** Produce a non-overlapping, fixed-width peakset weighted toward strongest evidence.

**Approach:** Greedy non-overlap on summit-centered fixed-width peaks ranked by signalValue.

```bash
#!/bin/bash
# Corces 2018 iterative overlap removal
PEAKS_DIR=peaks/per_sample
GENOME_SIZES=hg38.chrom.sizes
WIDTH_HALF=250                                       # 501 bp total width

# 1. Pool all peaks, re-center on summit, extend
awk -v w=$WIDTH_HALF 'BEGIN{OFS="\t"}
    {summit=$2+$10; print $1, summit-w, summit+w, $4, $7, $6}' \
    $PEAKS_DIR/*.narrowPeak | \
    awk '$2 >= 0' | \
    bedtools slop -i - -g $GENOME_SIZES -b 0 | \
    sort -k1,1 -k2,2n > pooled_recentered.bed

# 2. Sort by signalValue descending (column 5 in our BED -- which was column 7 of narrowPeak)
sort -k5,5gr pooled_recentered.bed > pooled_by_sig.bed

# 3. Iterative greedy non-overlap (process in significance order)
python3 - <<'EOF'
import sys
kept = []
with open('pooled_by_sig.bed') as f:
    for line in f:
        chrom, start, end, name, sig, strand = line.strip().split('\t')[:6]
        start, end = int(start), int(end)
        overlap = any(c == chrom and not (end <= s or start >= e) for c, s, e in kept)
        if not overlap:
            kept.append((chrom, start, end))

with open('consensus_iterative.bed', 'w') as f:
    for c, s, e in sorted(kept):
        f.write(f'{c}\t{s}\t{e}\n')
EOF

sort -k1,1 -k2,2n consensus_iterative.bed > consensus_final.bed
echo "Consensus peakset: $(wc -l < consensus_final.bed) fixed-width 501bp peaks"
```

For better performance at scale, replace the Python loop with `bedtools cluster` followed by per-cluster top-significance selection.

## Per-Strategy Failure Modes

### Naive union -- Width-driven differential

**Trigger:** Using `bedtools merge` of all per-rep peaks; counting reads in merged intervals.

**Mechanism:** A peak appearing as 200 bp in one rep but 800 bp in another merges to 800 bp in the union. Read count in 800 bp interval is biased high; differential analysis flags it as condition-specific even when it's just width difference.

**Symptom:** Top differential peaks track peak width (mean width different by 100+ bp between conditions).

**Fix:** Use fixed-width strategy (Corces iterative or DiffBind `summits=250`). Never use merged variable-width peaks for differential.

### Intersection (peak in all reps) -- Loses condition-specific biology

**Trigger:** Using `bedtools multiinter -i ... -j -intervals` requiring all samples.

**Mechanism:** Peaks present only in one condition fail intersection requirement and are excluded from the consensus, even though they are the biology of interest.

**Symptom:** Differential analysis returns near-zero significant peaks; the condition-specific peaks were filtered out before testing.

**Fix:** Use per-condition consensus, then union of consensus. Or majority rule with per-condition minOverlap.

### Majority rule (DiffBind default-ish) -- Borderline peaks dropped

**Trigger:** `dba.count(minOverlap=ceiling(N/2))` where N = total replicates and one condition has fewer reps.

**Mechanism:** A peak in 2/2 reps of cond1 but 0/3 reps of cond2 is in 2/5 = 40% < 50% -> dropped.

**Fix:** Compute consensus per condition first (e.g. peak in >= 2/3 reps), then union across conditions. This preserves condition-specific peaks.

### Width-controlled extension -- Shifts off summit

**Trigger:** Extending peak to fixed width using midpoint (`(start+end)/2`).

**Mechanism:** Midpoint of the called peak is rarely at the actual summit; particularly for asymmetric peaks (TSS-flanking) or wide broadPeaks.

**Fix:** Use summit position from narrowPeak column 10 (`start + summit_offset`). If summit unavailable (broadPeak), use peak start + half median width as approximation.

### IDR-filtered union -- Computational cost; threshold mismatch

**Trigger:** Running IDR for each rep pair, then unioning IDR-passed peaks.

**Mechanism:** IDR thresholds differ (true reps 0.05; pseudoreps 0.10). Mixing both into a union biases the consensus toward looser pseudorep peaks unless filtering carefully.

**Fix:** Decide upfront -- use only true-rep IDR (stricter, fewer peaks) OR pseudorep IDR (looser, more peaks). Don't mix.

## Decision Tree by Goal

| Goal | Strategy |
|------|----------|
| Standard 2-3 rep DA analysis with DiffBind | DiffBind `summits=250, minOverlap=2` |
| Modern ATAC publication / ML features | Corces 2018 iterative overlap (501 bp fixed) |
| ENCODE-compliant differential | IDR-passed per-rep-pair union (true-rep threshold 0.05) |
| Multi-condition with strong biology shift | Per-condition consensus then union |
| Single-cell ATAC pseudobulk | MACS3 per cluster -> iterative overlap across clusters |
| Cross-study peak comparison | Iterative overlap on shared genome build; lift over if needed |
| Quick exploratory | `bedtools merge` (acknowledge width bias; not for stats) |

## Reconciliation Across Methods

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Iterative overlap peakset much smaller than DiffBind summits=250 | Iterative removes overlapping summits within 500 bp; DiffBind keeps them | Both valid; iterative is sparser, DiffBind preserves resolution |
| Per-condition union then merged peak count >> within-condition consensus | Strong condition-specific peaks; biology is real | Use this strategy; do NOT collapse to global majority rule |
| IDR-pass peaks much fewer than majority-rule | IDR is stricter (reproducibility-aware) than overlap counting | Use IDR for high-confidence reporting; majority for exploratory |
| Same biology gives different peak counts on different genome builds | Lift-over noise; build-specific blacklist differences | Match builds; use ENCODE blacklist v2 specific to build |

**Operational rule:** Document the strategy explicitly in methods. The strategy choice can shift downstream results by 30-50% in peak count and FDR. Reproducibility requires explicit specification.

## Width-Re-centering Patterns

```python
# pybedtools: re-center on summit and fix width to 501 bp
import pybedtools as pbt

def fix_width_recenter(narrowpeak_path, half_width=250, out_path='consensus.bed'):
    """Re-center each narrowPeak on its summit and fix width to 2 * half_width + 1 bp."""
    rows = []
    for line in open(narrowpeak_path):
        f = line.strip().split('\t')
        chrom = f[0]; start = int(f[1]); summit_off = int(f[9])
        signal = float(f[6])
        summit = start + summit_off
        rows.append((chrom, max(0, summit - half_width), summit + half_width + 1,
                     f[3], signal, f[5]))
    rows.sort(key=lambda r: -r[4])               # Sort by signalValue desc
    bt = pbt.BedTool([list(map(str, r)) for r in rows])
    return bt.saveas(out_path)
```

## Per-Condition Consensus Then Union

```bash
# For each condition, build a within-condition consensus first
for cond in cond1 cond2; do
    bedtools multiinter -i ${cond}_rep1.narrowPeak ${cond}_rep2.narrowPeak ${cond}_rep3.narrowPeak \
        -names rep1 rep2 rep3 \
        -empty | \
        awk '$4 >= 2 {print $1"\t"$2"\t"$3}' > ${cond}_consensus.bed
done

# Union across condition consensuses (preserves condition-specific peaks)
cat cond1_consensus.bed cond2_consensus.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > all_conditions_consensus.bed
```

This is the right pattern when conditions differ enough that a global majority-rule would discard condition-specific peaks.

## R-native Fixed-Width Re-centering (GenomicRanges)

For Bioconductor pipelines that avoid intermediate BED files. The iterative-overlap step requires a loop; `!duplicated(findOverlaps(...))` does NOT implement greedy non-overlap (every peak self-overlaps).

```r
library(GenomicRanges); library(rtracklayer)

peaks <- import('peaks.narrowPeak')                    # GRanges
peaks_summit <- peaks                                  # column 10 of narrowPeak is summit offset
start(peaks_summit) <- start(peaks) + peaks$peak       # peak summit position
peaks_fixed <- resize(peaks_summit, width=501, fix='center')
peaks_fixed <- trim(peaks_fixed)                       # clamp to chromosome bounds

# Sort by signalValue (column 7) descending; greedy iterative non-overlap
peaks_sorted <- peaks_fixed[order(-peaks_fixed$signalValue)]
kept <- logical(length(peaks_sorted))
already_used <- GRanges()
for (i in seq_along(peaks_sorted)) {
    if (!any(overlapsAny(peaks_sorted[i], already_used))) {
        kept[i] <- TRUE
        already_used <- c(already_used, peaks_sorted[i])
    }
}
consensus <- peaks_sorted[kept]
export(consensus, 'consensus_iterative.bed')
```

For very large peaksets (>500k), use `bedtools` shell pipeline (faster) or precompute the per-rank reduce in chunks. The greedy loop is the canonical Corces 2018 algorithm; approximations via `reduce()` or `disjoin()` are NOT equivalent.

## Cross-Organism / Cross-Build Consensus

When merging peaks across genome builds (hg19 -> hg38; human -> mouse for cross-species analysis):

```bash
# Lift over hg19 peaks to hg38
liftOver published_peaks.hg19.bed hg19ToHg38.over.chain.gz \
    published_peaks.hg38.bed published_peaks.unmapped.bed

# Then merge with current hg38 peaks
cat published_peaks.hg38.bed my_peaks.hg38.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > merged_consensus.bed
```

**Trigger:** Cross-cohort meta-analysis or comparing to published datasets on older builds.

**Mechanism:** Genome assemblies differ in coordinates, gap regions, alt-haplotypes; chain files (UCSC) provide the position translation but ~3-5% of regions fail liftover (especially complex regions, sex chromosomes, MHC).

**Fix:** Always document liftOver chain version; report unmapped fraction; for high-stakes cross-build comparisons, re-align rather than liftover.

For human -> mouse comparison, use synteny-based mapping (`bnMapper` or `halLiftover`) rather than chain liftOver because conservation is non-trivial.

## Peak-Width Adaptive (When ATAC + Broad Histone Marks Co-analyzed)

**Trigger:** Combined ATAC narrow peaks + H3K27ac (broadPeak) or H3K4me1 in the same workflow.

**Mechanism:** Fixed 501 bp from ATAC under-counts H3K27ac which spreads across 2-5 kb domains. Counting H3K27ac reads in 501 bp windows misses most signal.

**Fix:** Build a dual-resolution consensus: 501 bp for ATAC narrow peaks; broadPeak union (2-5 kb) for H3K27ac. Run differential separately on each resolution. Or use a unified per-feature width derived from the broader signal.

## ENCODE-rE2G Candidate Element Lists

**Alternative to iterative-overlap:** ENCODE-rE2G (atac-seq/enhancer-gene-linking) publishes candidate cis-regulatory element (cCRE) lists for ENCODE cell types, ranked by predicted regulatory function rather than just ATAC signal strength. These are biologically prioritized and serve as an alternative consensus peakset for ML feature engineering. Trade-off: limited to ENCODE cell types; cell-type-specific cCREs may not align with the actual peakset of any specific dataset.

## Counting Reads in Consensus Peaks

```bash
# Convert BED to SAF for featureCounts
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
     {print $1"_"$2"_"$3, $1, $2, $3, "+"}' consensus.bed > consensus.saf

featureCounts -F SAF -a consensus.saf \
    -o consensus_counts.tsv \
    -p --countReadPairs \
    -T 8 \
    sample1.bam sample2.bam sample3.bam
```

`-p --countReadPairs` is required for paired-end ATAC; `-T 8` parallelizes across 8 threads. Output is the count matrix consumed by DESeq2 / edgeR / DiffBind.

## Blacklist Filtering

Blacklist filtering belongs after consensus construction, before or alongside differential testing.

```bash
# ENCODE blacklist v2 (Amemiya 2019)
bedtools intersect -v -a consensus_final.bed -b hg38-blacklist.v2.bed.gz \
    > consensus_no_blacklist.bed
echo "After blacklist: $(wc -l < consensus_no_blacklist.bed) peaks"
```

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Differential analysis flags hundreds of peaks; biology is implausible | Variable-width consensus driving width effects | Switch to fixed-width (iterative or summits=250) |
| `Negative coordinates` in fixed-width output | Peak summit too close to chromosome start | Clip with `awk '$2 >= 0'` or `bedtools slop -g chrom.sizes` |
| featureCounts SAF complaint about non-unique GeneID | Duplicate peak names in SAF | Use chr_start_end as GeneID; ensure peakset is non-redundant |
| Iterative overlap output much smaller than expected | Too many peaks within 500 bp; strict non-overlap | Try smaller half-width (e.g. 100) for higher resolution |
| `multiinter` output empty | Peak files not sorted; or wrong chromosome naming | Sort all inputs first; confirm chr name consistency |
| Per-condition union has 2x peaks of any per-condition | Strong condition-specific biology -- expected | Confirm with browser tracks; use this consensus |
| DiffBind takes hours on consensus | All-by-all counting is N samples x M peaks | Use `bParallel=TRUE`; provide `BPPARAM` |

## References

- Corces MR et al 2018 Science 362:eaav1898 (Iterative overlap, fixed-width 501 bp consensus standard)
- Stark R & Brown G 2011 DiffBind (Bioconductor; canonical reference; `summits=250` parameter)
- Quinlan AR & Hall IM 2010 Bioinformatics 26:841 (bedtools)
- Liao Y et al 2014 Bioinformatics 30:923 (featureCounts)
- Amemiya HM et al 2019 Sci Rep 9:9354 (ENCODE blacklist v2)
- Li Q et al 2011 Ann Appl Stat 5:1752 (IDR framework)
- ENCODE 4 ATAC-seq Standards (encodeproject.org/atac-seq)

## Related Skills

- atac-seq/atac-peak-calling - Generate per-replicate peaks
- atac-seq/differential-accessibility - Use consensus for DA testing
- atac-seq/atac-qc - Filter samples before consensus construction
- atac-seq/single-cell-atac - Per-cluster consensus for pseudobulk DA
- genome-intervals/bed-file-basics - bedtools operations on the consensus
- genome-intervals/interval-arithmetic - merge / intersect / subtract patterns
- chip-seq/peak-calling - Same consensus strategy applies to ChIP-seq
