---
name: bio-workflows-atacseq-pipeline
description: Orchestrates the end-to-end bulk ATAC-seq pipeline from FASTQ to differential accessibility and TF footprints, chaining Nextera-aware fastp QC, Bowtie2 alignment, chrM removal, dedup, a single Tn5 +4/-5 shift, MACS3 peak calling, Corces fixed-width consensus, DiffBind/csaw differential accessibility, and TOBIAS footprinting. Use when committing the reference build + blacklist once, recognizing ATAC has NO input control (the shift-extend model IS the background), applying the Tn5 shift exactly once (never combining -f BAMPE with --shift/--extsize), removing chrM before calling, building a fixed-width consensus so per-sample counts are comparable, or choosing MACS3 vs Genrich vs HMMRATAC. Hands mechanism to the atac-seq component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: MACS3
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - read-alignment/bowtie2-alignment
  - alignment-files/duplicate-handling
  - atac-seq/atac-peak-calling
  - atac-seq/atac-qc
  - atac-seq/consensus-peakset
  - atac-seq/differential-accessibility
  - atac-seq/footprinting
  - atac-seq/motif-deviation
  - atac-seq/nucleosome-positioning
qc_checkpoints:
  - after_qc: "Q30 >85%, adapter content <5% (Nextera)"
  - after_alignment: "Mapping rate >80%, mitochondrial <20% (Omni-ATAC lower)"
  - before_dedup: "NRF >0.8, PBC1 >0.8 (computed PRE-dedup)"
  - after_peaks: "FRiP >0.2, TSS enrichment >5 (ENCODE v3; v4 thresholds differ, do not mix)"
  - after_consensus: "Fixed-width (Corces 501 bp) consensus built before counting for differential"
---

## Version Compatibility

Reference examples tested with: Bowtie2 2.5.3+, MACS3 3.0+, Genrich 0.6+, bedtools 2.31+, deepTools 3.5+ (alignmentSieve), fastp 0.23+, samtools 1.19+, DiffBind 3.12+, TOBIAS 0.16+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: `macs3 callpeak -f BAMPE` uses real fragment lengths and IGNORES `--shift/--extsize/--nomodel`; the cut-site style needs `-f BAM`/`-f BED` on Tn5-shifted reads. `alignmentSieve --ATACshift` applies the +4/-5 shift once. ENCODE ATAC-seq v3 and v4 QC thresholds are not interchangeable. Confirm in-tool before quoting.

# ATAC-seq Pipeline

**"Run ATAC-seq from FASTQ to differential accessibility and footprints"** -> Chain QC/trim, alignment, chrM removal, dedup, a single Tn5 shift, peak calling, fixed-width consensus, differential accessibility, and footprinting.
- CLI + R: fastp -> bowtie2 -> drop chrM -> markdup -> Tn5 shift (once) -> macs3 -> Corces consensus -> DiffBind/csaw -> TOBIAS

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

ATAC-seq differs from ChIP-seq at four seams, and each is where the analysis goes wrong.

1. **There is NO input control -- the shift-extend cut-site model IS the background.** ATAC has no matched IP/input; peak significance comes from local lambda over the Tn5 insertion signal. Do not invent a "control"; commit instead to the build + ENCODE blacklist (removed before calling) as the coordinate frame.
2. **The Tn5 +4/-5 shift is applied EXACTLY ONCE, after dedup and chrM removal.** `alignmentSieve --ATACshift` (or one bedtools awk) applies it. Applying it twice, or combining `-f BAMPE` with `--shift/--extsize` (silently ignored), misplaces every cut site. Pick ONE calling mode: cut-site (`-f BAM`/`-f BED` + `--nomodel --shift -75 --extsize 150`) OR fragment (`-f BAMPE` on shifted reads, NO `--shift`).
3. **chrM is removed BEFORE peak calling.** Mitochondrial reads dominate ATAC libraries (often 20-50%, less with Omni-ATAC); leaving them in inflates depth and distorts FRiP and normalization.
4. **Differential accessibility requires a FIXED-WIDTH consensus peakset.** Variable-width MACS peaks make per-sample counts non-comparable. Build the Corces 501 bp iterative-overlap consensus (Corces 2018) so every region is the same width before counting; DiffBind/csaw then count into uniform intervals.

Reporting corollary: ENCODE ATAC v3 and v4 define TSS-enrichment/FRiP thresholds differently -- pick one standard and state which; do not mix rows across versions.

## Pipeline map

```
FASTQ (paired, Nextera)
  | [1] QC & trim -----------------> fastp (Nextera adapters)   (read-qc/fastp-workflow)
  v
  | [2] Align ---------------------> bowtie2 --very-sensitive -X 2000   (read-alignment/bowtie2-alignment)
  v     ^-- commitment: build + ENCODE blacklist (NO input control)
  | [3] Drop chrM (BEFORE dedup/peaks) -> mito can be 20-50% of reads
  v
  | [4] Dedup --------------------> markdup -r   (alignment-files/duplicate-handling)
  v
  | [5] Tn5 shift ONCE ------------> alignmentSieve --ATACshift (+4/-5)
  v     ^-- pick ONE calling mode; never BAMPE + --shift
  | [6] Peak calling -------------> macs3 (cut-site -f BAM --shift/--extsize | -f BAMPE)  (atac-seq/atac-peak-calling)
  v
  | [7] Fixed-width consensus -----> Corces 501 bp iterative overlap   (atac-seq/consensus-peakset)
  v
  | [8] QC + differential + footprints -> TSS/FRiP/fragment; DiffBind/csaw; TOBIAS  (atac-seq/atac-qc, differential-accessibility, footprinting)
  v
Accessibility peaks + differential regions + TF activity
```

## Made-once commitments

| Commitment | Choice | Consequence inherited downstream |
|------------|--------|----------------------------------|
| Build + blacklist | One build; ENCODE blacklist (removed before calling) | ATAC has no input, so the blacklist + shift-extend model ARE the background control |
| Tn5 shift | Applied ONCE (`--ATACshift`), then ONE calling mode | Double-shift or BAMPE+`--shift` misplaces cut sites |
| chrM handling | Removed before dedup/peaks | Mito reads (20-50%) inflate depth, FRiP, normalization |
| Differential interval | Fixed-width Corces 501 bp consensus | Variable-width peaks make per-sample counts non-comparable |

## The canonical order and why

1. **QC/trim** with Nextera adapters (`CTGTCTCTTATACACATCT`).
2. **Align** (bowtie2 `--very-sensitive -X 2000`) so the full nucleosome-spanning fragment distribution is captured.
3. **Remove chrM, then compute NRF/PBC, then dedup** -- order-trap on both ends: `markdup -r` physically removes duplicates, so NRF/PBC1 computed afterwards are identically 1.0; and mito reads are over-amplified, so computing them before chrM removal measures chrM chemistry, not nuclear-library complexity. The binding constraint is PRE-DEDUP. Mito must go before peak calling regardless.
4. **Dedup** (collate -> fixmate -m -> sort -> markdup -r).
5. **Tn5 shift ONCE** (`alignmentSieve --ATACshift`).
6. **Call peaks in ONE mode** -- order-trap: `-f BAMPE` + `--shift/--extsize` silently drops the flags.
7. **Build the fixed-width consensus** (Corces 501 bp) -- order-trap: differential on variable-width peaks is not comparable.
8. **QC, differential (DiffBind/csaw on the consensus), footprinting (TOBIAS)**.

## Choosing the caller and calling mode

Pipeline-level selection only; mechanism lives in the component skills.

| Fork | Lean toward | Hand off to |
|------|-------------|-------------|
| Caller | MACS3 (standard); Genrich (`-j` ATAC mode: handles replicates + chrM + blacklist in one pass); HMMRATAC (nucleosome-aware HMM) | atac-seq/atac-peak-calling |
| Calling mode | Cut-site `-f BAM`/`-f BED` + `--nomodel --shift -75 --extsize 150` (ENCODE smoothing window on shifted reads) vs fragment `-f BAMPE` on shifted reads (no `--shift`) | atac-seq/atac-peak-calling |
| Consensus | Corces 2018 iterative-overlap fixed-width 501 bp | atac-seq/consensus-peakset |
| Differential | DiffBind / csaw / DESeq2 on the fixed-width count matrix; spike-in for global shifts | atac-seq/differential-accessibility |

## Primary path: Bowtie2 + Tn5 shift + MACS3

**Goal:** turn Nextera FASTQ into shifted, chrM-free peaks ready for a fixed-width consensus.

**Approach:** align with a wide insert window, drop chrM, dedup, Tn5-shift once, then call in ONE mode. Full runnable script: `examples/atacseq_workflow.sh`; differential: `examples/differential_atac.R`.

```bash
bowtie2 -p 8 -x bt2_index/genome -1 trimmed/${s}_R1.fq.gz -2 trimmed/${s}_R2.fq.gz \
    --very-sensitive --no-mixed --no-discordant -X 2000 2> aligned/${s}.log \
  | samtools view -@4 -bS -q 30 -f 2 - | samtools sort -@4 -o aligned/${s}.sorted.bam
samtools index aligned/${s}.sorted.bam

# Drop chrM BEFORE dedup/peaks (mito dominates ATAC), then dedup
samtools idxstats aligned/${s}.sorted.bam | cut -f1 | grep -v -e '^chrM$' -e '^MT$' \
  | xargs samtools view -b aligned/${s}.sorted.bam > aligned/${s}.noMT.bam
samtools collate -@8 -O -u aligned/${s}.noMT.bam | samtools fixmate -m -u - - \
  | samtools sort -@8 -u - | samtools markdup -r -@8 - aligned/${s}.dedup.bam
samtools index aligned/${s}.dedup.bam            # alignmentSieve needs an indexed input BAM

# Tn5 +4/-5 shift ONCE
alignmentSieve -b aligned/${s}.dedup.bam -o aligned/${s}.shifted.bam --ATACshift -p 8
samtools index aligned/${s}.shifted.bam

# Remove ENCODE blacklist regions BEFORE calling (the made-once commitment above; see the example script)
# Everything downstream (peaks, counts, footprints) consumes ${s}.filt.bam, never ${s}.shifted.bam.
# NOTE: examples/atacseq_workflow.sh names its blacklist-FILTERED output `.shifted.bam`; same reads,
# different name. Match on the step, not the suffix.
bedtools intersect -v -a aligned/${s}.shifted.bam -b "$BLACKLIST" > aligned/${s}.filt.bam
samtools index aligned/${s}.filt.bam

# Cut-site calling on the shifted, blacklist-filtered reads (ONE mode; do NOT also use -f BAMPE with these flags)
macs3 callpeak -t aligned/${s}.filt.bam -f BAM -g hs -n ${s} --outdir peaks \
    --nomodel --shift -75 --extsize 150 --keep-dup all -q 0.01
```

For the ENCODE 4 IDR + pseudoreplicate pipeline and the Corces 501 bp iterative-overlap consensus, see atac-seq/atac-peak-calling and atac-seq/consensus-peakset.

## Differential accessibility and footprinting

**Goal:** compare accessibility across conditions on comparable intervals, then read TF activity.

**Approach:** count into the fixed-width consensus with DiffBind (or csaw), then run the TOBIAS three-step (ATACorrect -> ScoreBigwig -> BINDetect) for footprints.

```r
library(DiffBind)                                  # counts into the fixed-width consensus
dba <- dba(sampleSheet = samples)                  # bamReads = shifted BAMs, Peaks = per-sample narrowPeak
dba <- dba.count(dba)                              # use summits/consensus for uniform width
dba <- dba.normalize(dba); dba <- dba.contrast(dba, categories = DBA_CONDITION)
dba <- dba.analyze(dba); report <- dba.report(dba)
```

```bash
# peaks/consensus.bed is the Corces 501 bp FIXED-WIDTH consensus from atac-seq/consensus-peakset (step 7).
# It is NOT peaks/consensus_peaks.narrowPeak, which is the variable-width pooled MACS3 call; build the
# fixed-width set first or these three commands have no input.
# TOBIAS three-step: bias-correct -> score -> detect bound motifs (differential across two conditions).
# Footprint on the BLACKLIST-FILTERED reads (${s}.filt.bam), the same reads MACS3 called peaks from --
# blacklist regions are artifact pileups, and bias-correcting over them corrupts the footprint scores.
TOBIAS ATACorrect -b aligned/${s}.filt.bam -g genome.fa -p peaks/consensus.bed --outdir foot --cores 8
TOBIAS ScoreBigwig --signal foot/${s}_corrected.bw --regions peaks/consensus.bed --output foot/${s}.bw --cores 8
TOBIAS BINDetect --motifs motifs.jaspar --signals foot/ctrl.bw foot/treat.bw --genome genome.fa \
    --peaks peaks/consensus.bed --outdir foot/bindetect --cores 8
```

## QC checkpoints between steps

| After | Gate | Interpretation |
|-------|------|----------------|
| Alignment | Mapping >80%, mito <20% (Omni-ATAC lower) | High mito = suboptimal lysis; drop before calling |
| PRE-dedup | NRF >0.8, PBC1 >0.8 | Low complexity = over-amplification/low input; compute before dedup |
| Peaks | FRiP >0.2, TSS enrichment >5 (v3) | Low TSS/FRiP = over/under-digestion or degraded chromatin (atac-seq/atac-qc) |
| Fragment size | NFR <100 bp, mono ~200 bp, di ~400 bp periodicity | Loss of nucleosome periodicity = over-digestion (Tn5:DNA too high) |
| Consensus | Fixed-width (501 bp) built before counting | Variable-width peaks make counts non-comparable |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Depth/FRiP dominated by one contig; few real peaks | chrM not removed before calling | Drop chrM/MT before dedup and peak calling |
| Cut sites offset / footprints smeared | Tn5 shift applied twice, or `-f BAMPE` used with `--shift/--extsize` | Shift ONCE; pick ONE calling mode (cut-site `-f BAM` OR fragment `-f BAMPE`) |
| Differential counts not comparable across samples | Counted into variable-width MACS peaks | Build the Corces 501 bp fixed-width consensus first |
| Looked for an input/IgG track and found none | ATAC has no input control | Use the blacklist + shift-extend model as background; do not fabricate a control |
| QC numbers disagree with a reference | Mixed ENCODE v3 and v4 thresholds | Pick one ENCODE version and report which |

## Pipeline map (hand-offs)

- read-qc/fastp-workflow - Nextera adapter trimming
- read-alignment/bowtie2-alignment - the aligner, wide insert window
- alignment-files/duplicate-handling - collate/fixmate/sort/markdup order
- atac-seq/atac-peak-calling - MACS3/Genrich/HMMRATAC, ENCODE 4 IDR, calling modes
- atac-seq/atac-qc - TSS enrichment, FRiP, NRF/PBC, fragment periodicity
- atac-seq/consensus-peakset - Corces 2018 iterative-overlap fixed-width consensus
- atac-seq/differential-accessibility - DiffBind/csaw/DESeq2 on the consensus
- atac-seq/footprinting - TOBIAS three-step and per-TF failure modes
- atac-seq/nucleosome-positioning - V-plot, NucleoATAC, +1 nucleosome

The complete runnable scripts are in this skill's examples/ (`atacseq_workflow.sh`, `differential_atac.R`).

## Related Skills

- database-access/sra-data - Pull ATAC-seq FASTQ from SRA / ENA
- database-access/geo-data - Resolve GEO accessions for ATAC datasets
- read-qc/fastp-workflow - Nextera adapter trimming and quality filtering
- read-alignment/bowtie2-alignment - Standard ATAC-seq aligner
- alignment-files/duplicate-handling - MarkDuplicates pre-peak-calling
- atac-seq/atac-peak-calling - MACS3 / Genrich / HMMRATAC details, ENCODE 4 IDR
- atac-seq/atac-qc - TSS enrichment, FRiP, NRF/PBC1/PBC2 details
- atac-seq/consensus-peakset - Corces 2018 iterative-overlap fixed-width consensus
- atac-seq/differential-accessibility - DiffBind / csaw / DESeq2; spike-in normalization
- atac-seq/footprinting - TOBIAS three-step; per-TF failure modes
- atac-seq/motif-deviation - chromVAR for motif accessibility variability
- atac-seq/nucleosome-positioning - V-plot, NucleoATAC, +1 nucleosome
- atac-seq/single-cell-atac - For scATAC instead of bulk
- atac-seq/co-accessibility - Cicero cis-regulatory inference
- atac-seq/enhancer-gene-linking - ABC, ENCODE-rE2G enhancer-gene mapping
- atac-seq/deep-learning-atac - chromBPNet variant-effect prediction
- atac-seq/allele-specific-accessibility - WASP + caQTL mapping
- chip-seq/peak-annotation - Annotate ATAC peaks to genes

## References

- Buenrostro JD, Giresi PG, Zaba LC, Chang HY, Greenleaf WJ (2013) Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position. *Nature Methods* 10:1213-1218. DOI 10.1038/nmeth.2688. (original ATAC-seq.)
- Corces MR, Trevino AE, Hamilton EG, et al (2017) An improved ATAC-seq protocol reduces background and enables interrogation of frozen tissues. *Nature Methods* 14:959-962. DOI 10.1038/nmeth.4396. (Omni-ATAC.)
- Corces MR, Granja JM, Shams S, et al (2018) The chromatin accessibility landscape of primary human cancers. *Science* 362:eaav1898. DOI 10.1126/science.aav1898. (fixed-width iterative-overlap consensus peakset.)
- Bentsen M, Goymann P, Schultheis H, et al (2020) ATAC-seq footprinting unravels kinetics of transcription factor binding during zygotic genome activation. *Nature Communications* 11:4267. DOI 10.1038/s41467-020-18035-1. (TOBIAS.)
