---
name: bio-workflows-chipseq-pipeline
description: Orchestrates the end-to-end ChIP-seq pipeline from FASTQ to blacklist-filtered, annotated peaks, chaining fastp QC, Bowtie2 alignment, pre-dedup library-complexity QC (NRF/PBC), duplicate removal, chrM + ENCODE-blacklist filtering, MACS3 peak calling against a matched input, IDR/consensus reproducibility, deepTools signal tracks, and ChIPseeker annotation. Use when committing the reference build + blacklist version + effective genome size once, pairing each IP with its matched control, computing complexity metrics BEFORE dedup, choosing narrow vs broad and MACS3 vs SEACR/Genrich, keeping per-replicate peaks for IDR, or avoiding depth-normalization that erases a spike-in global shift. Hands mechanism to the chip-seq component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: MACS3
workflow: true
depends_on:
  - read-qc/fastp-workflow
  - read-alignment/bowtie2-alignment
  - alignment-files/duplicate-handling
  - chip-seq/chipseq-qc
  - chip-seq/peak-calling
  - chip-seq/peak-annotation
  - chip-seq/differential-binding
  - chip-seq/chipseq-visualization
  - chip-seq/motif-analysis
qc_checkpoints:
  - after_qc: "Q30 >85%, adapter content <5%"
  - after_alignment: "Mapping rate >80%, unique mapping >70%"
  - before_dedup: "NRF >0.8, PBC1 >0.8 (computed on the PRE-dedup BAM; after dedup the metric is meaningless)"
  - after_peaks: "FRiP >1% (TF) or >5% (sharp histone; broad marks run lower); NSC >1.05; RSC >0.8; fingerprint separates IP from input"
  - after_idr: "IDR rescue ratio max(Np,Nt)/min and self-consistency ratio max(N1,N2)/min both <=2 (ENCODE); IDR run on PER-REPLICATE peaks"
---

## Version Compatibility

Reference examples tested with: Bowtie2 2.5.3+, MACS3 3.0+, HOMER 4.11+, bedtools 2.31+, deepTools 3.5+, fastp 0.23+, samtools 1.19+, ChIPseeker 1.38+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: `macs3 callpeak -f BAMPE` uses real fragment lengths and IGNORES `--shift/--extsize/--nomodel` (those apply to single-end `-f BAM`); the `-g` shortcut (`hs/mm`) sets the effective genome size and must match the build/read length. Confirm in-tool before quoting.

# ChIP-seq Pipeline

**"Process my ChIP-seq data from FASTQ to annotated peaks"** -> Chain QC/trim, alignment, pre-dedup complexity QC, dedup + blacklist filtering, control-matched peak calling, reproducibility, signal tracks, and annotation.
- CLI + R: fastp -> bowtie2 -> (NRF/PBC pre-dedup) -> samtools markdup -> chrM/blacklist filter -> macs3 callpeak (IP vs input) -> IDR -> bamCoverage -> ChIPseeker

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step. Every step below cross-references the component skill that teaches its mechanism.

## The governing principle

A ChIP-seq peakset is decided at four seams, not inside the caller.

1. **The reference build + blacklist version + effective genome size is one coordinate commitment made once and inherited by everything downstream** — peak coordinates, signal-track scaling, and every overlap. Blacklist filtering is a committed pipeline step, not optional cleanup: ENCODE-blacklisted regions (satellite/rDNA/high-signal artifacts) produce reproducible false peaks in every dataset regardless of biology, so they are removed before calling (Amemiya 2019).
2. **An IP is only interpretable against its matched control — the control IS the enrichment background.** Calling peaks without the right input/IgG fabricates peaks at open/accessible and copy-number-amplified regions. Pair each IP with its control at the calling step.
3. **Library-complexity QC (NRF/PBC1/PBC2) is computed on the PRE-dedup BAM.** After `markdup -r` the duplicates are gone, so computing complexity afterward reads ~1.0 and is meaningless. Compute it on the filtered, position-sorted BAM before removing duplicates.
4. **Normalization must not silently undo the experiment.** deepTools RPKM/CPM rescales every library to the same depth, which ERASES a spike-in global-shift signal (the whole point of ChIP-Rx). For spike-in experiments use `--scaleFactor` + `--normalizeUsing None`; for standard experiments RPKM/CPM is fine (chip-seq/spike-in-normalization).

Reproducibility corollary: pool replicates for a consensus peakset, but keep PER-REPLICATE peaks — IDR needs individual replicates plus pooled pseudo-replicates; running IDR on an already-pooled peakset is not IDR.

## Pipeline map

```
FASTQ (IP + matched Input, replicates)
  | [1] QC & trim -----------------> fastp                (read-qc/fastp-workflow)
  v
  | [2] Align ---------------------> bowtie2 (-q30 unique) (read-alignment/bowtie2-alignment)
  v     ^-- commitment: build + blacklist version + effective genome size
  | [3] Complexity QC (PRE-dedup) -> NRF/PBC1/PBC2         (chip-seq/chipseq-qc)
  v
  | [4] Dedup + filter ------------> markdup -r; drop chrM; SUBTRACT ENCODE blacklist  (alignment-files/duplicate-handling)
  v
  | [5] Peak calling (IP vs input)-> macs3 callpeak (narrow | --broad)  (chip-seq/peak-calling)
  v     ^-- keep PER-REPLICATE peaks for IDR
  | [6] Reproducibility -----------> IDR (per-rep + pooled pseudo-reps)  (chip-seq/peak-calling)
  v
  | [7] Signal tracks -------------> bamCoverage (RPKM | spike-in scaleFactor)  (chip-seq/chipseq-visualization)
  v
  | [8] QC + Annotate -------------> FRiP/NSC/RSC/fingerprint; ChIPseeker  (chip-seq/chipseq-qc, peak-annotation)
  v
Blacklist-filtered, annotated, reproducible peaks
```

## Made-once commitments

| Commitment | Choice | Consequence inherited downstream |
|------------|--------|----------------------------------|
| Build + blacklist + effective genome size | One genome build; the matching ENCODE blacklist BED; `-g hs/mm`/numeric | Mixed builds mis-place peaks; skipping the blacklist plants reproducible false peaks; wrong `-g` mis-scales p-values |
| Control pairing | Each IP has its input/IgG | No control => peaks at open chromatin / CN-amplified loci |
| Peak shape | Narrow (TF, H3K4me3, H3K27ac) vs broad (H3K27me3, H3K36me3, H3K9me3) | Broad marks called with narrow settings fragment into many small peaks |
| Fragment model | PE: `-f BAMPE` (real fragments); SE: `-f BAM` + `--nomodel --extsize` from predictd/xcorr | BAMPE silently ignores `--shift/--extsize` |

## The canonical order and why

1. **QC/trim** (fastp) both IP and input.
2. **Align** (bowtie2), keep uniquely-mapped (`samtools view -q 30`), coordinate-sort.
3. **Compute NRF/PBC1/PBC2 on the PRE-dedup BAM** — order-trap: after dedup they are meaningless.
4. **Mark/remove duplicates** (collate -> fixmate -m -> sort -> markdup -r), then **drop chrM** and **subtract the ENCODE blacklist** — order-trap: skipping the blacklist leaves reproducible artifact peaks.
5. **Call peaks against the matched control** (narrow or `--broad`).
6. **IDR on per-replicate peaks** (+ pooled pseudo-replicates) — order-trap: IDR on a pooled peakset is not IDR.
7. **Signal tracks** — RPKM/CPM for standard; `--scaleFactor` + `--normalizeUsing None` for spike-in (order-trap: RPKM erases the spike-in global shift).
8. **QC (FRiP/NSC/RSC/fingerprint) and annotate** (ChIPseeker).

## Choosing the caller and peak shape

Pipeline-level selection only; mechanism lives in the component skills.

| Fork | Lean toward | Hand off to |
|------|-------------|-------------|
| Caller | MACS3 (standard IP+input); SEACR (CUT&RUN/CUT&Tag, low background); Genrich (some ChIP/ATAC, built-in blacklist/replicate handling) | chip-seq/peak-calling, chip-seq/cut-and-run-tag |
| Narrow vs broad | Narrow: TFs, H3K4me3, H3K27ac. Broad (`--broad --broad-cutoff 0.1`): H3K27me3, H3K36me3, H3K9me3 | chip-seq/peak-calling |
| Reproducibility | ENCODE IDR (per-rep + pooled pseudo-reps) for TFs; naive overlap acceptable for exploratory histone | chip-seq/peak-calling |
| Consensus set | Pool for a union/consensus set AFTER IDR selects the reproducible threshold | chip-seq/differential-binding |

## Primary path: Bowtie2 + MACS3 + ChIPseeker

**Goal:** turn IP+input FASTQ into a blacklist-filtered, control-matched, annotated peakset.

**Approach:** align and keep unique reads, measure complexity before dedup, dedup + drop chrM + subtract the blacklist, call against the control, then annotate. Full runnable script: `examples/narrow_peak_workflow.sh`; annotation: `examples/peak_annotation.R`.

```bash
bowtie2 -p 8 -x bt2_index/genome -1 trimmed/${s}_R1.fq.gz -2 trimmed/${s}_R2.fq.gz \
    --no-mixed --no-discordant --maxins 1000 2> aligned/${s}.log \
  | samtools view -@4 -bS -q 30 - | samtools sort -@4 -o aligned/${s}.sorted.bam
samtools index aligned/${s}.sorted.bam

# Complexity QC on the PRE-dedup BAM (NRF = distinct positions / total; PBC1 = singletons / distinct).
# Counted per-mate here (close to ENCODE fragment-level values); use `bamtobed -bedpe` for exact parity.
bedtools bamtobed -i aligned/${s}.sorted.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | sort | uniq -c \
  | awk '{tot+=$1; dist++; if($1==1) one++} END{printf "NRF=%.3f PBC1=%.3f\n", dist/tot, one/dist}'

# Dedup, drop chrM, then SUBTRACT the ENCODE blacklist (committed step, not optional)
samtools collate -@8 -O -u aligned/${s}.sorted.bam | samtools fixmate -m -u - - \
  | samtools sort -@8 -u - | samtools markdup -r -@8 - aligned/${s}.dedup.bam
samtools index aligned/${s}.dedup.bam
samtools idxstats aligned/${s}.dedup.bam | cut -f1 | grep -v -e '^chrM$' -e '^MT$' \
  | xargs samtools view -b aligned/${s}.dedup.bam > aligned/${s}.nochrM.bam
bedtools intersect -v -a aligned/${s}.nochrM.bam -b ENCODE_blacklist.bed > aligned/${s}.final.bam
samtools index aligned/${s}.final.bam
```

```bash
# Narrow (TFs, sharp marks) vs broad (spreading marks). -f BAMPE uses real fragment sizes.
macs3 callpeak -t aligned/IP_rep1.final.bam aligned/IP_rep2.final.bam \
    -c aligned/Input_rep1.final.bam aligned/Input_rep2.final.bam \
    -f BAMPE -g hs -n experiment --outdir peaks -q 0.01 --keep-dup all   # dedup done upstream (markdup -r); tell MACS3 to keep all
# Broad marks: add  --broad --broad-cutoff 0.1  (do NOT call H3K27me3 with narrow settings)
```

For IDR, call peaks PER REPLICATE (and on pooled pseudo-replicates) with a relaxed `-q`, then run `idr` across them (chip-seq/peak-calling). For higher confidence, intersect a second caller (HOMER `-style histone` for all histone marks).

## Signal tracks and annotation

```bash
# Standard experiment: RPKM/CPM is fine. SPIKE-IN experiment: this would ERASE the global shift.
bamCoverage -b aligned/IP_rep1.final.bam -o bigwig/IP_rep1.bw --normalizeUsing RPKM -p 8
# Spike-in (ChIP-Rx): bamCoverage --scaleFactor <spike-in factor> --normalizeUsing None  (chip-seq/spike-in-normalization)
```

Annotation uses a project GTF via `makeTxDbFromGFF()` when provided, else a pre-built TxDb. `overlap='all'` couples gene assignment with feature overlap (host-gene convention); default `overlap='TSS'` assigns the nearest-TSS gene independently. Full code: `examples/peak_annotation.R`.

## QC checkpoints between steps

| After | Gate | Interpretation |
|-------|------|----------------|
| QC/trim | Q30 >85%, adapter <5% | DNA higher quality than RNA |
| Alignment | Mapping >80%, unique >70% | Low unique = repeats/contamination/wrong build |
| PRE-dedup | NRF >0.8, PBC1 >0.8 | Low complexity = over-amplification/low input; MUST be computed before dedup |
| Peaks | FRiP >1% (TF) / >5% (sharp histone; broad marks run lower); NSC >1.05; RSC >0.8; fingerprint separates IP/input | Low FRiP/flat fingerprint = weak antibody or failed enrichment (chip-seq/chipseq-qc) |
| IDR | rescue ratio and self-consistency ratio both <=2 | Poor replicate consistency; run IDR on PER-replicate peaks |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Reproducible peaks over satellite/rDNA/high-signal regions | ENCODE blacklist never subtracted | `bedtools intersect -v` the blacklist BED before calling (committed step) |
| NRF/PBC ~1.0 and uninformative | Computed after `markdup -r` | Compute complexity on the PRE-dedup, filtered BAM |
| Peaks at open chromatin / CN-amplified loci | Called without a matched control | Pair each IP with its input/IgG in `callpeak -c` |
| H3K27me3/H3K9me3 fragmented into many tiny peaks | Broad mark called with narrow settings | Add `--broad --broad-cutoff 0.1` |
| `--shift/--extsize` had no effect | Used with `-f BAMPE` (ignored for PE) | Use `-f BAM` + `--nomodel` for SE; BAMPE derives fragments |
| Spike-in global shift disappears in tracks | bamCoverage RPKM/CPM re-equalized depth | `--scaleFactor` + `--normalizeUsing None` (chip-seq/spike-in-normalization) |
| "IDR" numbers look too good | IDR run on a pooled peakset | Run IDR on per-replicate peaks + pooled pseudo-replicates |

## Pipeline map (hand-offs)

- read-qc/fastp-workflow - adapter/quality trimming
- read-alignment/bowtie2-alignment - the standard ChIP-seq aligner, build/index
- alignment-files/duplicate-handling - collate/fixmate/sort/markdup order
- chip-seq/chipseq-qc - NRF/PBC, FRiP, NSC/RSC, fingerprint, hyper-ChIPable detection
- chip-seq/peak-calling - MACS3/SEACR/Genrich/HOMER, IDR vs naive overlap
- chip-seq/peak-annotation - ChIPseeker/HOMER/GREAT
- chip-seq/differential-binding - DiffBind/csaw and the normalization-problem framing
- chip-seq/chipseq-visualization - deepTools tracks and normalization choices
- chip-seq/spike-in-normalization - ChIP-Rx global-shift experiments
- chip-seq/motif-analysis - HOMER/MEME-ChIP/monaLisa

The complete runnable scripts are in this skill's examples/ (`narrow_peak_workflow.sh`, `peak_annotation.R`).

## Related Skills

- database-access/sra-data - Pull ChIP-seq FASTQ from SRA / ENA for re-analysis
- database-access/geo-data - Resolve ENCODE / Roadmap GSE accessions to SRA
- read-qc/fastp-workflow - Upstream adapter trimming and quality filtering
- read-alignment/bowtie2-alignment - Standard ChIP-seq aligner
- alignment-files/duplicate-handling - MarkDuplicates pre-peak-calling
- chip-seq/chipseq-qc - FRiP, NSC/RSC, library complexity, antibody validation
- chip-seq/peak-calling - MACS3/MACS2/HOMER/SPP, IDR vs naive overlap, per-tool failure modes
- chip-seq/peak-annotation - ChIPseeker, HOMER, ENCODE cCRE classification, GREAT regulatory domains
- chip-seq/differential-binding - DiffBind, DESeq2, csaw with the three-normalization-problems framing
- chip-seq/chipseq-visualization - deepTools, pyGenomeTracks, heatmaps with bigWig normalization choices
- chip-seq/motif-analysis - HOMER, MEME-ChIP (STREME), monaLisa with background-selection theory
- chip-seq/super-enhancers - ROSE/ROSE2/LILY for SE calling (H3K27ac vs MED1 vs BRD4)
- chip-seq/cut-and-run-tag - SEACR + MACS2 consensus for CUT&RUN/CUT&Tag (different protocol)
- chip-seq/spike-in-normalization - ChIP-Rx Drosophila spike-in for global-shift experiments
- chip-seq/chromatin-state-segmentation - ChromHMM multi-mark integration into chromatin states
- chip-seq/chip-deep-learning - BPNet/chromBPNet/Enformer for variant-effect prediction
- chip-seq/allele-specific-binding - WASP/BaalChIP/RASQUAL for allele-specific TF binding

## References

- Zhang Y, Liu T, Meyer CA, et al (2008) Model-based analysis of ChIP-Seq (MACS). *Genome Biology* 9:R137. DOI 10.1186/gb-2008-9-9-r137.
- Landt SG, Marinov GK, Kundaje A, et al (2012) ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. *Genome Research* 22:1813-1831. DOI 10.1101/gr.136184.111. (NSC/RSC, FRiP, IDR practice.)
- Li Q, Brown JB, Huang H, Bickel PJ (2011) Measuring reproducibility of high-throughput experiments. *Annals of Applied Statistics* 5:1752-1779. DOI 10.1214/11-AOAS466. (the IDR framework.)
- Amemiya HM, Kundaje A, Boyle AP (2019) The ENCODE blacklist: identification of problematic regions of the genome. *Scientific Reports* 9:9354. DOI 10.1038/s41598-019-45839-z.
