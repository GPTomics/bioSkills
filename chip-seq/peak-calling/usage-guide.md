# Peak Calling - Usage Guide

## Overview

Calls ChIP-seq peaks with MACS3/MACS2, HOMER, SPP, or Genrich. Handles narrow (TF) vs broad (histone) mode selection, fragment-size modeling vs `--nomodel`, effective genome size, ENCODE-style IDR (TFs) or naive overlap (histones), aligner-specific shift caveats, and hyper-ChIPable region awareness. Cross-references chipseq-qc for pre-call validation (antibody, fragment size, FRiP) and chip-seq/cut-and-run-tag for CUT&RUN/CUT&Tag protocols (which use SEACR + MACS2 rather than MACS3 alone).

## Prerequisites

```bash
# Core peak callers
conda install -c bioconda macs3 macs2 homer samtools bedtools idr

# Phantompeakqualtools for SPP / cross-correlation fragment length
conda install -c bioconda phantompeakqualtools

# Optional: Genrich for joint replicate analysis
conda install -c bioconda genrich
```

## Quick Start

Tell the agent what to do:
- "Call peaks from chip.bam with input.bam control using ENCODE-compliant MACS2 parameters"
- "Call broad peaks for H3K27me3 with input control and naive overlap across 3 replicates"
- "Estimate fragment length via cross-correlation before peak calling"
- "Run IDR on two MACS narrowPeak replicates and apply ENCODE Nself/Nt consistency rules"
- "Call peaks on chr21-only data; adjust genome size accordingly"
- "Find peaks at the level of high-confidence consensus across MACS3 and HOMER"
- "My MACS model building failed; diagnose and rescue"

## Example Prompts

### ENCODE-compliant TF ChIP-seq
> "Run the ENCODE TF peak-calling pipeline: per-replicate MACS2 at `-p 0.01 --nomodel --shift 0 --extsize {fraglen}` with fraglen from cross-correlation, then IDR at 0.05 across true replicates with the Nself/Nt consistency check."

### Histone marks (broad)
> "Call H3K27me3 broad peaks with `--broad --broad-cutoff 0.1` for three replicates and combine via naive overlap with ≥40% reciprocal overlap."

### Narrow histone marks
> "Call H3K4me3 peaks. Use narrow mode in MACS3 and `-style histone` in HOMER (not `-style factor` despite the sharp signal; Omnipeak 2025 benchmark)."

### Subset data
> "Call peaks from chr21-only ChIP-seq data with ~400k reads. Use numeric genome size, skip MACS model building, and use a narrow `--extsize 147` for the H3K4me3 mark."

### Diagnosing peak quality
> "Why are my peaks shifted ~75bp from known motif positions? Check fragment-size modeling and shift parameters."

### Hyper-ChIPable artifact suspicion
> "I'm seeing strong peaks at rRNA loci and histone gene clusters. Filter against a custom blacklist of top-1% input signal regions."

### Reconciliation
> "MACS3 and HOMER disagree on which peaks pass. Walk through the reconciliation rules from the skill and decide what to trust."

## What the Agent Will Do

1. **Pre-call validation**: confirm antibody is validated (KO/KD or peptide-array); inspect fragment-size distribution via `samtools view -f 0x2 sample.bam | awk '{print $9}'`; verify input control matches library prep / fragmentation; check FRiP, NSC, RSC if available (see chipseq-qc)
2. **Choose narrow vs broad mode** from target biology (TFs narrow; H3K27me3/H3K9me3/H3K36me3 broad; H3K4me3/H3K27ac narrow)
3. **Choose modeling strategy**: full model only for whole-genome, ≥1M treatment reads; `--nomodel --extsize` for subsets, low-depth, or ENCODE consistency
4. **Determine effective genome size**: read-length-matched value from deepTools table; numeric for subset data
5. **Call peaks per replicate** at loose threshold (`-p 1e-2` ENCODE pattern); generate signal tracks with `-B --SPMR`
6. **Replicate handling**: IDR @ 0.05 for TFs with Nself/Nt consistency rule; naive overlap with ≥40% reciprocal for histones
7. **Quality sanity check**: peak counts in expected range for target (e.g., 20-50k for H3K4me3); FRiP threshold met
8. **Hyper-ChIPable filter** if claims are made at rRNA / tRNA / histone clusters / mitochondrial DNA
9. **Output conversion**: narrowPeak/broadPeak to BED if downstream requires
10. **Document parameters** in methods (genome version + size, fragment length source, IDR threshold, blacklist version)

## Pipeline Quick Reference

### ENCODE TF pipeline

```bash
# 1. Estimate fragment length via cross-correlation
Rscript run_spp.R -c=chip.bam -savp=chip_cc.pdf -out=chip_cc.txt
FRAGLEN=$(awk -F'\t' '{print $3}' chip_cc.txt | cut -d',' -f1)

# 2. Per-replicate (loose threshold for IDR)
for i in rep1 rep2; do
    macs2 callpeak -t ${i}.bam -c input.bam \
        -f BAM -g 2.701e9 -n ${i} \
        --nomodel --shift 0 --extsize $FRAGLEN \
        --keep-dup all -B --SPMR -p 1e-2
done

# 3. IDR across true replicates
sort -k8,8nr rep1_peaks.narrowPeak > rep1.sorted
sort -k8,8nr rep2_peaks.narrowPeak > rep2.sorted
idr --samples rep1.sorted rep2.sorted \
    --input-file-type narrowPeak --rank p.value \
    --idr-threshold 0.05 --output-file true_reps.idr --plot
```

### ENCODE histone pipeline (broad)

```bash
# Per-replicate broad peaks
for i in rep1 rep2 rep3; do
    macs2 callpeak -t ${i}.bam -c input.bam \
        -f BAM -g 2.701e9 -n ${i} \
        --broad --broad-cutoff 0.1 \
        --nomodel --shift 0 --extsize $FRAGLEN \
        --keep-dup all -B --SPMR -p 1e-2
done

# Naive overlap (≥40% reciprocal, present in ≥2 reps)
bedtools intersect -a rep1_peaks.broadPeak -b rep2_peaks.broadPeak -f 0.40 -r -u > tmp12.bed
bedtools intersect -a tmp12.bed -b rep3_peaks.broadPeak -f 0.40 -r -u > naive_overlap.bed
```

## Troubleshooting

### Model building fails

MACS2/3 needs ≥100 paired plus/minus enrichment regions within `--mfold` (default `[5, 50]`). Causes:
- Single-chromosome or targeted data → use `--nomodel`
- Low read count (<500k) → use `--nomodel`
- Sparse enrichment → widen with `--mfold 3 50` first, then fall back to `--nomodel`

Always inspect `<sample>_model.r` plot; silent model failure produces wrong fragment size.

### Zero peaks called

- Wrong genome size on subset data (`hs` on chr21 inflates lambda 60×) → use numeric `-g`
- Wrong `-f` flag (BAM vs BAMPE vs BED) → match input file type
- Swapped treatment / control → verify `-t` is enriched sample
- Library too shallow → check sequenced read count; minimum ~10M unique mapped for TFs, 20M+ for broad marks

### Peak count >> 500k

- Did not deduplicate → `samtools view -F 1804 -q 30` filter pre-call
- chrM reads dominate → remove chrM from BAM before calling
- `-q` too loose → tighten to `-q 0.01`
- Hyper-ChIPable artifacts → build custom blacklist from top-1% input signal

### Peaks shifted from motif by ~75 bp

- `--shift` not set with `-f BAM` → add `--shift 0 --extsize {fraglen}`
- Cross-correlation fraglen wrong → inspect `_cc.txt` and verify against `predictd` output
- Aligner pre-applied shift (chromap) → drop downstream shift OR use `--no-correction`

### IDR returns 0 reproducible peaks

- Sorted by wrong column → use `sort -k8,8nr` (p-value descending)
- Library size imbalance >2× → pseudoreplicate IDR will fail Nself/Nt rule; rebalance via downsampling or re-do shallow rep
- One replicate is bad → check per-rep FRiP, NSC, RSC; do not average

### FRiP < 1%

This is a hard fail. Re-validate the antibody (KO/KD control), check fragment-size distribution, confirm input control matches. Proceeding with low-FRiP data wastes downstream effort.

## Output Files

| File | Description |
|------|-------------|
| `*_peaks.narrowPeak` | BED6+4: chr, start, end, name, score, strand, signalValue, -log10(p), -log10(q), summit offset |
| `*_peaks.broadPeak` | BED6+3: broad region coordinates without summit |
| `*_summits.bed` | Per-peak summit positions (narrow only) |
| `*_model.r` | R script for fragment-size model visualization (inspect this!) |
| `*_treat_pileup.bdg` | Treatment signal (with `-B`); pileup of fragments |
| `*_control_lambda.bdg` | Control signal (with `-B`); local lambda estimate |

Convert narrowPeak to BED: `cut -f1-5 peaks.narrowPeak > peaks.bed`. For browser viewing: `sort -k1,1 -k2,2n peaks_treat_pileup.bdg | bedGraphToBigWig - chrom.sizes peaks.bw`.

## Tips

- ENCODE pattern uses `--keep-dup all` because deduplication happens upstream during BAM filtering. If skipping the ENCODE filter, set `--keep-dup auto` and accept MACS's binomial p-value estimate.
- `-p 1e-2` is intentionally loose for ENCODE; IDR (TF) or naive overlap (histone) tightens downstream. Single-replicate workflows should use `-q 0.05` or `-q 0.01` instead.
- For HOMER, ALWAYS use `-style histone` for histone marks (including narrow H3K4me3, H3K27ac per Omnipeak 2025 benchmark). Reserve `-style factor` for TFs.
- Effective genome size matters: use deepTools `effectiveGenomeSize` table read-length-matched value, not `hs`/`mm` shorthand.
- For paired-end ChIP, `-f BAMPE` uses actual fragment spans (good for histones); `-f BAM --shift 0 --extsize {fraglen}` is the ENCODE TF pattern.
- TFs and most histone marks need 20-25M unique mapped reads per replicate (ENCODE 2012 standard); H3K27me3/H3K9me3 need 40-60M for adequate broad-domain detection.
- Always document: aligner version, genome assembly + effective size, fragment length source, peak caller version + parameters, IDR/overlap threshold, blacklist version. Reviewers ask.

## Related Skills

- chip-seq/chipseq-qc - FRiP, NSC, RSC, library complexity, antibody validation
- chip-seq/cut-and-run-tag - CUT&RUN/CUT&Tag protocols (different caller choice, lower depth)
- chip-seq/spike-in-normalization - When global signal shifts expected (HDACi, BETi, EZH2i)
- chip-seq/differential-binding - DiffBind/csaw downstream
- chip-seq/peak-annotation - Annotate to genes and ENCODE cCREs
- chip-seq/motif-analysis - Motif discovery on peak sequences
- chip-seq/super-enhancers - Stitch H3K27ac peaks for SE calls
- atac-seq/atac-peak-calling - ATAC-specific Tn5 shift; no input control
- alignment-files/sam-bam-basics - Pre-call BAM filtering
- alignment-files/duplicate-handling - MarkDuplicates before peak calling
