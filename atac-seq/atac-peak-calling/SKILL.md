---
name: bio-atac-seq-atac-peak-calling
description: Call accessible chromatin regions from ATAC-seq BAM files using MACS3, MACS2, Genrich, or HMMRATAC. Use when identifying open chromatin from aligned ATAC-seq, choosing between point-source vs HMM peak callers, applying ENCODE-style pseudoreplicate IDR, removing blacklist regions, or fixing 501bp consensus peaks for downstream differential analysis.
tool_type: cli
primary_tool: macs3
---

## Version Compatibility

Reference examples tested with: MACS3 3.0.2+, MACS2 2.2.9+, Genrich 0.6.1+, HMMRATAC 1.2+ (now bundled in MACS3 as `macs3 hmmratac`), samtools 1.19+, bedtools 2.31+, IDR 2.0.4+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws unexpected errors, introspect the installed binary (`<tool> -h`) and adapt the example to match the actual CLI rather than retrying.

# ATAC-seq Peak Calling

**"Call accessible regions from my ATAC-seq BAM"** -> Identify Tn5-hypersensitive open chromatin, treating fragments as point insertion events (not protein-bound regions as in ChIP-seq) and accounting for the lack of input control.

- CLI (canonical, ENCODE 4): `macs2 callpeak -t atac.bam -f BAM -g hs -n sample --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01` (use `-f BAM`, not `-f BAMPE` -- BAMPE reads true fragment ends and ignores `--shift`/`--extsize`)
- CLI (HMM-based, single sample): `macs3 hmmratac -i atac.bam -n sample --outdir hmm_out`
- CLI (joint replicates): `Genrich -j -t rep1.bam,rep2.bam -o peaks.narrowPeak -e chrM -E blacklist.bed`

The `-p 0.01` (loose) plus IDR is the ENCODE pattern: low stringency increases peak overlap between replicates, and IDR rescues the reproducible set. Single-sample workflows usually swap to `-q 0.05` instead.

## Algorithmic Taxonomy

| Tool | Model | Treats fragments as | Min reps | Strength | Fails when |
|------|-------|---------------------|----------|----------|------------|
| MACS3/MACS2 | Local Poisson lambda + FDR | Point-source insertions (+/- shift) | 1 | Mature, ENCODE-default, fast, narrow + broad modes | Confounds NFR with broad accessible domains; no input means lambda from local genome only |
| Genrich (ATAC mode -j) | q-value on log-transformed p-value, joint replicate model | Whole fragments (paired-end intervals) | 1 (multi-rep optional) | Treats reps jointly; can exclude chrM via `-e chrM`; auto blacklist via `-E`; PCR-dup removal via `-r` | Less peer-reviewed than MACS; thin literature; slow on deep libraries |
| MACS3 hmmratac (was HMMRATAC) | 3-state HMM (open / nucleosomal / background) on fragment-size signal | Fragment-size classes | 1 | Models nucleosome periodicity directly; differentiates NFR and flanking nucleosomes | Needs >= 30M de-duplicated nuclear reads; memory-hungry; slow; flat fragment distribution -> garbage HMM |
| HOMER `findPeaks -style dnase` | Fixed window + fold-change cutoff | Tag positions | 1 | Convenient for downstream HOMER motif analysis | Less calibrated p-values than MACS; window-size sensitive |
| nf-core/atacseq | Wrapper (MACS2 by default) | Same as MACS2 | 1 | Reproducible Nextflow pipeline with QC built in | Only as good as the underlying caller |

Methodology evolves; verify the current ENCODE ATAC-seq Standards (encodeproject.org pipelines/atac-seq) before locking parameters. ENCODE 4 still defaults to MACS2 (not MACS3) at time of writing; `macs3 callpeak` is API-compatible for ATAC parameters but not yet the official ENCODE binary.

## Shift-Extend vs BAMPE: The Critical Choice

Two valid ways to feed paired-end ATAC into MACS:

**Pattern A (ENCODE / "single-end-ified"):** `-f BAMPE` actually IGNORES `--shift/--extsize`. To activate them, use `-f BAM` and treat each end independently. ENCODE's pipeline uses `-f BAM --shift -75 --extsize 150` to model each Tn5 cut as a 150 bp window centered on the insertion site, ignoring fragment lengths.

**Pattern B (paired-fragment):** `-f BAMPE` uses the full paired-end fragment span as the signal interval. Best when fragment lengths are biologically meaningful (e.g., NFR-only peak calling at 38/75 bp). In BAMPE mode, do NOT set `--shift/--extsize` (silently ignored, but confusing).

For most bulk ATAC, Pattern A matches ENCODE convention and is reproducible against published peak sets. Pattern B can be more sensitive at narrow regulatory elements but does not match ENCODE outputs.

## Effective Genome Size

`-g hs` and `-g mm` are MACS shorthands for old defaults. Modern values:

| Genome | MACS shorthand | Actual mappable size | Source |
|--------|---------------|----------------------|--------|
| hg38 | `-g hs` (2.7e9) | 2.701e9 (50bp), 2.748e9 (75bp), 2.806e9 (100bp), 2.862e9 (150bp) | deepTools `effectiveGenomeSize` |
| hg19 | `-g hs` (2.7e9) | 2.686e9 (50bp), 2.777e9 (100bp) | deepTools |
| mm10 | `-g mm` (1.87e9) | 2.308e9 (50bp), 2.408e9 (75bp), 2.467e9 (100bp) | deepTools |
| mm39 | none | 2.310e9 (50bp), 2.468e9 (100bp) | deepTools |

Wrong size shifts every q-value but rarely changes peak ranks. Use `unique-kmers.py` (khmer) or the deepTools tabulated values for exact sizes; the shorthand is a decade-old approximation.

## Effective Genome Size: When It Matters

**Trigger:** Comparing peaks across genome builds or species; reproducing published q-value cutoffs; hi-resolution lambda estimation.

**Mechanism:** MACS estimates genome-wide lambda as `total_reads / effective_size`. Wrong size -> wrong null -> shifted q-values, especially at the marginal cutoff.

**Symptom:** Peak counts diverge ~10-20% from published numbers when re-running an old dataset.

**Fix:** Pull the read-length-matched value from deepTools `effectiveGenomeSize` table. For pipelines, parameterize this; never inline the shorthand for cross-study comparisons.

## Per-Tool Failure Modes

### MACS2/MACS3 -- Confounded NFR + broad accessibility

**Trigger:** Cell type with extended open domains (e.g., active super-enhancers, MYOD1 regulons, locus-control regions).

**Mechanism:** Default narrow-peak mode segments wide accessible domains into multiple smaller peaks at local lambda spikes; `--broad --broad-cutoff 0.1` merges them but inflates total length and breaks IDR comparability.

**Symptom:** Peak count >> 200k for human bulk ATAC at ENCODE depth; mean peak width < 200 bp; visual inspection in IGV shows 3-5 calls under one continuous accessibility block.

**Fix:** Run both narrow and broad; use narrow for differential analysis, broad for domain-level enrichment (e.g., super-enhancer overlap). Do NOT use `--call-summits` for broad mode.

### Genrich -- Replicate weighting and chrM exclusion

**Trigger:** Replicates with very different library sizes; high-mitochondrial samples not pre-filtered.

**Mechanism:** Genrich's joint mode combines the per-replicate p-values at each position via Fisher's method. Library-size imbalance dominates the joint p-value; chrM reads inflate background unless `-e chrM` is set.

**Symptom:** Most-significant peaks cluster on chrM or on the largest-library replicate's high-coverage regions.

**Fix:** Always pass `-e chrM` (Genrich 0.6+) and `-E blacklist.bed`. Down-sample BAMs to common depth (`samtools view -s`) before joint calling if libraries differ >2x. Add `-r` to remove PCR duplicates inside Genrich, OR pre-deduplicate (do not do both).

### MACS3 hmmratac (HMMRATAC) -- Depth and fragment-size dependence

**Trigger:** Library < 25M nuclear reads, or libraries with degraded chromatin and flat fragment-size distribution.

**Mechanism:** The 3-state HMM is trained from fragment-size classes (NFR ~50 bp, mono ~200 bp, di ~400 bp peaks). Without periodicity the emission distributions collapse and the HMM cannot separate states.

**Symptom:** Output BED is empty, or all peaks are tiny (~150 bp) with no nucleosome flanks called; runtime explodes (>24h) on shallow data.

**Fix:** Verify fragment-size periodicity in QC first (atac-qc skill). If flat, fall back to MACS3 callpeak. HMMRATAC needs deep coverage (a practical minimum around 30M deduplicated nuclear reads).

### HOMER findPeaks -- Window-size sensitivity

**Trigger:** Default `-style dnase` uses 75 bp peaks; ATAC peaks are 250-500 bp typically.

**Mechanism:** HOMER's window-based caller does not auto-fit width to ATAC.

**Fix:** Use `-style factor -size 150` for narrow ATAC peaks, or skip HOMER for peak calling and use it only for downstream motif analysis on MACS peaks.

### Aligner choice -- chromap vs bwa-mem2 vs bowtie2 affects peak shape

**Trigger:** Switching aligners between datasets and expecting reproducible peaks.

**Mechanism:** chromap (Zhang 2021) applies its own ATAC-specific 4 bp / -5 bp Tn5 shift before fragment output; bwa-mem2 and bowtie2 do not. Downstream `--shift -75 --extsize 150` parameters are calibrated for unshifted bwa/bowtie BAMs; applying them to chromap output double-shifts the signal.

**Symptom:** Peaks called from chromap output are shifted by ~5-10 bp relative to bwa output at the same locus.

**Fix:** When using chromap, drop `--shift` and `--extsize` (chromap's pre-shift is sufficient) OR omit chromap's `--Tn5-shift` (the shift is opt-in, applied only when that flag or an ATAC preset is set) so it is not double-applied, then proceed with standard MACS parameters. Document the aligner version and any shift choices in methods. Within a project, pin the aligner.

### Single-sample (no replicate) -- Rotation / circular-shift permutation

**Trigger:** Single biological sample without any replicate for IDR.

**Mechanism:** IDR requires two replicates by construction. For n=1, statistical confidence per peak comes from local background (Poisson p-value) but reproducibility cannot be assessed.

**Fix:** Apply a stricter `-q 0.01` (vs ENCODE `-p 0.01` + IDR pattern) and additionally apply rotation/circular-shift permutation: shift the BAM cuts by a random distance modulo each chromosome and re-call peaks; the per-peak persistence rate across rotations is a non-parametric reproducibility proxy. Document this is a single-sample heuristic, not ENCODE-compliant.

## ENCODE 3 vs ENCODE 4 Differences

| Feature | ENCODE 3 (legacy) | ENCODE 4 (current) |
|---------|-------------------|---------------------|
| Per-rep significance threshold | `-q 0.05` directly | `-p 0.01` (loose) + IDR |
| Pseudoreplicate IDR cutoff | Not formalized | `--idr-threshold 0.10` self-consistency |
| TSS enrichment threshold | >= 6 (older) | >= 7 (hg38, GENCODE v29) |
| Mt fraction expectation | < 25% | < 20% (Omni-ATAC < 5%) |
| Blacklist | v1 | v2 (Amemiya 2019) |
| Default genome size | hardcoded `hs`/`mm` | encouraged: deepTools effectiveGenomeSize |

To reproduce a published ENCODE 3 dataset, pin the original pipeline and threshold exactly. ENCODE 4 results are not directly numerically comparable to ENCODE 3 even on the same input BAM.

## Super-Enhancer Detection

For active super-enhancer (SE) annotation alongside narrow-peak workflow, ROSE (Whyte 2013) and LILY (Boeva 2017) stitch ATAC or H3K27ac peaks separated by < 12.5 kb and rank by signal:

```bash
# ROSE expects H3K27ac BAM but works on ATAC narrowPeak with care
ROSE_main.py -g hg38 -i atac_peaks.gff -r atac.bam -o rose_out/ -t 2500
```

ROSE-style stitching is complementary to MACS3 narrow peaks: narrow peaks for differential analysis; SE annotation for biology interpretation. SE calls require H3K27ac input for definitive annotation; ATAC alone produces "stretch enhancers" that overlap but are not identical to H3K27ac SE.

## ENCODE 4 ATAC-seq Pipeline (Reference Implementation)

The exact ENCODE pattern produces the most-comparable peak sets:

```bash
# Per-replicate peak calling (loose threshold)
macs2 callpeak \
    -t rep1.filt.dedup.bam \
    -f BAM -g hs \
    -n rep1 --outdir peaks/rep1/ \
    --nomodel --shift -75 --extsize 150 \
    --keep-dup all \
    -B --SPMR \
    -p 0.01

# Pooled (all replicates)
macs2 callpeak \
    -t rep1.filt.dedup.bam rep2.filt.dedup.bam \
    -f BAM -g hs -n pooled --outdir peaks/pooled/ \
    --nomodel --shift -75 --extsize 150 --keep-dup all -B --SPMR -p 0.01

# Pseudoreplicates (two independent 50% subsamples; approximate, not a disjoint partition)
samtools view -b -h -s 1.5 rep1.filt.dedup.bam > rep1.psr1.bam   # seed.fraction
samtools view -b -h -s 2.5 rep1.filt.dedup.bam > rep1.psr2.bam   # different seed
# (call peaks on each pseudoreplicate the same way)
```

`--SPMR` writes signal as Signal Per Million Reads (normalized bedGraph). `-p 0.01` is intentionally loose; IDR will tighten to a reproducible set.

## IDR for Reproducible Peaks

**Goal:** Find peaks reproducible across biological replicates at controlled IDR.

**Approach:** Score paired peak lists by signalValue, fit IDR's two-component mixture (reproducible + noise), threshold at IDR <= 0.05 (true reps) or 0.10 (pseudoreplicates).

```bash
# Sort peaks by p-value (column 8) so IDR scores by significance
sort -k8,8nr rep1_peaks.narrowPeak > rep1.sorted.narrowPeak
sort -k8,8nr rep2_peaks.narrowPeak > rep2.sorted.narrowPeak

# True replicates -- threshold IDR <= 0.05
idr --samples rep1.sorted.narrowPeak rep2.sorted.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file true_reps.idr \
    --idr-threshold 0.05 --plot --log-output-file idr.log

# Pseudoreplicates -- threshold IDR <= 0.10 (looser, ENCODE Nself <= 2 rule)
idr --samples psr1_peaks.narrowPeak psr2_peaks.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file psr.idr --idr-threshold 0.10 --plot
```

**ENCODE consistency rules:** Nt = peaks passing IDR on true reps; Nself = peaks passing IDR on pseudoreps. Library passes if `max(Nt, Nself) / min(Nt, Nself) <= 2`. If both ratios > 2, the library is rejected.

**IDR fails when:** Ranking column choice matters. `--rank p.value` (column 8) is robust; `--rank signal.value` (column 7) breaks if MACS pile-up scaling differs between replicates.

## Decision Tree by Experimental Scenario

| Scenario | Recommended caller | Why |
|----------|-------------------|-----|
| Bulk ATAC, 2-3 reps, depth >= 25M | MACS2 ENCODE pipeline + IDR | Reproducible, comparable to published peaksets |
| Bulk ATAC, 1 sample (no rep) | MACS3 callpeak with `-q 0.05`; do not run IDR | IDR is meaningless without reps; tighter q-value substitutes |
| Bulk ATAC, depth >= 30M, want NFR + flanking nuc structure | MACS3 hmmratac | HMM separates NFR from nucleosome flanks |
| Multi-replicate joint analysis where rep weighting is symmetric | Genrich `-j` ATAC mode | Joint p-value across reps; built-in chrM and blacklist |
| Cell type with broad super-enhancer accessibility | MACS3 `--broad --broad-cutoff 0.1` for SE; narrow for differential | Domain-level inference vs site-level |
| FFPE / degraded chromatin (flat fragment dist) | MACS3 callpeak with stringent `-q 0.01`; never HMMRATAC | HMM needs fragment periodicity |
| scATAC pseudobulk per cluster | MACS3 callpeak per cluster + iterative overlap | See atac-seq/single-cell-atac |
| Want fixed-width consensus peaks for differential | Call broadly, then re-center to summits +/- 250 bp | See atac-seq/consensus-peakset |
| Plant / non-model organism | MACS3 with `-g <effective_size>`; verify size empirically | Default `-g hs/mm` invalid; compute via khmer |

## Reconciliation: When Callers Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MACS narrow peaks much fewer than Genrich | Genrich q-cutoff different default (`-q 0.05` log-scale, MACS `-q 0.05` linear) | Re-run with `-q 0.01` (Genrich) for parity |
| HMMRATAC misses peaks MACS finds | Library too shallow OR fragment-size periodicity weak | Trust MACS; HMMRATAC is depth-sensitive |
| HMMRATAC calls peaks MACS misses | HMM is sensitive to mid-strength accessibility flanked by phased nucleosomes | Inspect; often genuine but unconfirmed by short-fragment signal |
| Same peak called by all but width 2x different | Broad mode vs narrow mode mismatch | Standardize: re-center to summit +/- 250 bp for differential |
| Per-rep MACS calls peak; pooled MACS does not | One rep dominates; lambda smoothes it out in pooled | Trust pooled + IDR over per-rep counts |

**Operational rule for high-confidence reporting:** Require a peak to pass IDR <= 0.05 on true replicates AND survive blacklist/greylist filtering AND have mean signalValue >= 5 across reps. Two callers from different families (MACS + Genrich) agreeing within 250 bp is acceptable evidence when IDR is unavailable.

## Blacklist and Greylist

```bash
# ENCODE blacklist (Amemiya 2019) -- always remove
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
bedtools intersect -v -a peaks.narrowPeak -b hg38-blacklist.v2.bed > peaks.no_blacklist.narrowPeak

# Sample-specific greylist (input-derived high-signal regions; rarely available for ATAC)
# For ATAC, ENCODE recommends pooling all samples' top-percentile signal and removing
# regions exceeding 100x median coverage as a "soft greylist"
```

Blacklist is mandatory; greylist is optional and most useful when the same library prep produces consistent artifact regions across samples.

## NFR-Only Peak Calling

**Goal:** Call peaks using only sub-nucleosomal fragments (<100 bp) for sharper TF-binding-relevant accessibility.

**Approach:** Pre-filter BAM to short fragments, then call peaks with parameters scaled to the smaller fragment length.

```bash
samtools view -h sample.dedup.bam | \
    awk 'substr($0,1,1)=="@" || ($9 > 0 && $9 < 100) || ($9 < 0 && $9 > -100)' | \
    samtools view -b > nfr.bam
samtools index nfr.bam

macs2 callpeak -t nfr.bam -f BAM -g hs -n sample_nfr \
    --nomodel --shift -37 --extsize 75 \
    --keep-dup all -p 0.01
```

`--shift -37 --extsize 75` halves both parameters to match shorter fragments; this is a fragment-scaled convention for NFR-focused input, not a TOBIAS-specified setting (TOBIAS instead applies the +4/-5 Tn5 correction to the full BAM via ATACorrect).

## Output Files (narrowPeak)

| Column | Field | Notes |
|--------|-------|-------|
| 1-3 | chrom, start, end | 0-based, half-open |
| 4 | name | MACS auto-numbers |
| 5 | score | Min(int(-10*log10(qvalue)), 1000) |
| 6 | strand | `.` for ATAC |
| 7 | signalValue | Fold enrichment over local lambda |
| 8 | pValue | -log10 p |
| 9 | qValue | -log10 q (BH-FDR) |
| 10 | summit_offset | Peak summit relative to start |

Convert to bigWig for browsers: `sort -k1,1 -k2,2n sample_treat_pileup.bdg > sample.sorted.bdg && bedGraphToBigWig sample.sorted.bdg chrom.sizes sample.bw` (bedGraphToBigWig is multi-pass and cannot read from a pipe/stdin, so sort to a file first).

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `--shift/--extsize ignored` warning | Used `-f BAMPE` with these flags | Switch to `-f BAM` or remove the flags |
| 0 peaks called | Forgot `--nomodel`; MACS tries to build a shifting model and fails | Add `--nomodel --shift -75 --extsize 150` |
| Peak count >> 500k | Did not deduplicate; or did not remove chrM; or `-q` too loose | Pre-filter (samtools view -F 1804 -q 30; samtools idxstats); use `-q 0.01` |
| `Sequence chrM not found` (Genrich) | Wrong chromosome name in `-e` flag (chrM vs MT) | Match BAM header naming convention |
| HMMRATAC out of memory | Deep library / large genome; the current tool is `macs3 hmmratac` (Python, not Java) | Increase available RAM and use a scratch `--outdir`; the `-Xmx`/`HMMRATAC.jar` heap flags apply only to the deprecated standalone Java HMMRATAC |
| Peaks shifted by 75 bp from expected positions | Forgot `--shift -75` (cuts at one end of read) | Add the shift; positions are now centered on Tn5 cut site |
| IDR returns 0 reproducible peaks | Sorted by wrong column; ranks are random | Sort each peakset by `-k8,8nr` (p-value descending) |

## References

- Buenrostro JD et al 2013 Nat Methods 10:1213 (ATAC-seq protocol)
- Corces MR et al 2017 Nat Methods 14:959 (Omni-ATAC protocol)
- Corces MR et al 2018 Science 362:eaav1898 (iterative-overlap fixed-width 501 bp consensus peaks)
- Tarbell ED & Liu T 2019 Nucleic Acids Res 47:e91 (HMMRATAC)
- Gaspar JM, Genrich: detecting sites of genomic enrichment (github.com/jsh58/Genrich; no published paper)
- Li Q et al 2011 Ann Appl Stat 5:1752 (IDR framework)
- Landt SG et al 2012 Genome Res 22:1813 (ENCODE/modENCODE peak calling guidelines, IDR Nself rule)
- Amemiya HM et al 2019 Sci Rep 9:9354 (ENCODE blacklist v2)
- ENCODE ATAC-seq Standards (encodeproject.org/atac-seq) -- canonical pipeline parameters

## Related Skills

- atac-seq/atac-qc - Verify TSS enrichment, FRiP, and fragment periodicity before calling
- atac-seq/consensus-peakset - Combine per-sample peaks into a fixed-width differential-ready set
- atac-seq/single-cell-atac - Pseudobulk peak calling per cluster
- atac-seq/differential-accessibility - Downstream DiffBind/csaw/DESeq2 testing
- atac-seq/deep-learning-atac - chromBPNet bias-corrected per-base profiles as alternative input
- read-alignment/bowtie2-alignment - Upstream ATAC alignment
- alignment-files/duplicate-handling - Pre-call dedup with Picard MarkDuplicates
- chip-seq/peak-calling - ChIP-seq comparison (uses input control)
- chip-seq/super-enhancers - ROSE / LILY for super-enhancer annotation
- genome-intervals/bed-file-basics - Peak file manipulation
