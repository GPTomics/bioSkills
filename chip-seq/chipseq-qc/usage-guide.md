# ChIP-seq QC - Usage Guide

## Overview

Compute the ENCODE-compliant ChIP-seq QC battery (FRiP, NSC/RSC, library complexity, fingerprint metrics, replicate concordance) and detect hyper-ChIPable artifacts before committing to downstream analysis. Embeds antibody validation logic and fragment-size diagnostics; both of which are upstream of any numeric metric and determine whether the metrics themselves are interpretable.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda samtools bedtools phantompeakqualtools idr deeptools featurecounts

# Python utilities for custom metrics
pip install pysam pybedtools pandas
```

```r
# Bioconductor ChIPQC for one-call full report
BiocManager::install('ChIPQC')
```

## Quick Start

Tell the agent what to do:
- "Compute the full ENCODE QC battery (FRiP, NSC, RSC, NRF, PBC1, PBC2, fingerprint) for chip.bam vs input.bam"
- "Run phantompeakqualtools and interpret NSC/RSC against ENCODE thresholds for an H3K4me3 sample"
- "Generate a fragment-size distribution to check for over-sonication"
- "Check whether my peaks are concentrated at hyper-ChIPable regions (rRNA, tRNA, mtDNA, histone clusters)"
- "Run IDR on rep1 and rep2 narrowPeak files with ENCODE Nself/Nt consistency rule"
- "Compute replicate Spearman correlation across 3 H3K27me3 replicates"
- "Generate a ChIPQC R report for 6 samples in a sample sheet"

## Example Prompts

### Full QC sweep
> "Run the full ENCODE QC battery: FRiP, NSC/RSC via phantompeakqualtools, NRF/PBC1/PBC2, deepTools plotFingerprint JS distance and AUC, and replicate Spearman correlation. Grade each metric against ENCODE thresholds."

### Antibody validation interpretation
> "I have ChIP-seq from WT and KO cells with the same antibody. Compute FRiP and signal at known target loci to verify antibody specificity (KO should drop to background)."

### Fragment-size diagnostic
> "Extract the fragment-size distribution from properly-paired reads and tell me whether the chromatin was over-sonicated."

### Hyper-ChIPable detection
> "Build a custom cell-type-specific blacklist from the top 1% input signal and intersect my peaks against it to flag suspicious calls at hyper-ChIPable loci."

### IDR with Nself/Nt
> "Run IDR at threshold 0.05 on true replicates and at threshold 0.10 on pseudoreplicates. Apply the ENCODE Nself/Nt rule (both ratios must be ≤ 2) and tell me whether the library passes."

### Fingerprint interpretation
> "Generate a plotFingerprint curve for ChIP vs Input and interpret AUC, JS distance, and synthetic JS against the appropriate threshold for an H3K27me3 broad mark."

### Failure diagnosis
> "My NSC is 1.02 but FRiP is 8%. What's going on and should I trust these peaks?"

## What the Agent Will Do

1. **Antibody status check**: confirm validation (KO/KD or peptide-array); request the catalog number + lot if missing from metadata
2. **Fragment-size diagnostic**: extract from properly-paired BAM; classify as sub-nucleosomal (TF), nucleosomal (histone), or over-sonicated (rescue impossible)
3. **Compute the QC battery:**
   - FRiP: `bedtools intersect` reads-in-peaks / total mapped
   - NSC, RSC, QualityTag: `Rscript run_spp.R` (phantompeakqualtools)
   - NRF, PBC1, PBC2: pre-deduplication BAM
   - JS distance, AUC, synthetic JS: `plotFingerprint --outQualityMetrics`
   - Replicate Spearman: `multiBamSummary` + `plotCorrelation`
4. **Grade against ENCODE thresholds**: TF-specific vs histone-specific cutoffs (broad histones tolerate lower NSC, JS distance; FRiP threshold higher for histones)
5. **Hyper-ChIPable artifact check**: build top-1% input signal blacklist; flag peaks in suspect regions; verify by motif enrichment + KO/KD if claims are made
6. **Replicate concordance**: IDR with Nself/Nt rule (TFs) or naive overlap (histones); replicate Spearman correlation; identify failing replicates
7. **Recommendation**: pass / caution / reject with the specific failing metric(s) and remediation path
8. **Output**: qc_report.txt + diagnostic plots (fingerprint, correlation heatmap, fragment-size histogram, cross-correlation plot)

## QC Decision Logic

The agent applies the following operational rules:

| Condition | Action |
|-----------|--------|
| Antibody unvalidated AND no KO/KD control AND no peptide-array data | Flag; request validation before trusting metrics |
| Fragment-size distribution flat or 100-1000 bp continuum | Reject; over-sonication unrecoverable |
| FRiP < 1% (TF) or < 5% (histone) | Reject or re-do ChIP; metric implies antibody / IP failure |
| NSC < 1.05 OR RSC < 0.8 (TF) | Caution; check fragment length consistency; for broad histones may be acceptable with strong FRiP |
| NRF < 0.5 OR PBC1 < 0.5 | Severe PCR bottleneck; re-do library or accept reduced statistical power |
| Replicate Spearman < 0.6 (narrow), < 0.4 (broad) | Reject pair; check sample swap or batch effect |
| Nself/Nt ratio > 2 in either direction | Reject library per ENCODE rule |
| FRiP excellent but signal concentrated at rRNA/mtDNA/HIST clusters | Build custom blacklist; recompute on filtered peaks |
| All metrics pass + replicates concordant | Proceed to peak calling and downstream analysis |

## Tips

- **Always compute QC before peak calling.** Failing QC samples produce false-positive peaks that consume downstream resources.
- **NRF / PBC1 / PBC2 must be computed BEFORE deduplication.** Post-MarkDuplicates BAMs report meaningless library complexity (NRF ~ 1.0).
- **FRiP varies with peak caller stringency.** Tighter q-value -> fewer peaks -> higher FRiP per peak but may understate enrichment quality. Use ENCODE `-p 1e-2` calls for consistent FRiP reporting.
- **phantompeakqualtools is depth-sensitive.** Subsample to 15-25M reads (`samtools view -s`) for consistent NSC/RSC across replicates of different depth.
- **plotFingerprint JS distance has mark-specific thresholds.** TFs > 0.3; broad histones can be 0.05-0.15 and still be high-quality.
- **Replicate correlation is not transitive.** Three-replicate experiments need pairwise Spearman > 0.8 for ALL pairs, not just the average.
- **Custom blacklist beats ENCODE blacklist alone** for cell-type-specific hyper-ChIPable regions. ENCODE v2 catches repeats; cell-type-specific high-input regions catch transcription artifacts.
- **ChIPQC R package gives the full battery in one call** for 4-12 samples; phantompeakqualtools is the canonical source for individual NSC/RSC values.
- **For CUT&RUN/CUT&Tag, QC thresholds differ** (FRiP > 25%, spike-in % aligned 0.5-5%, no input control). See cut-and-run-tag skill.
- **Document everything.** Methods sections should include: antibody catalog + lot + validation method, fragment-size summary, FRiP, NSC, RSC, library complexity, replicate concordance metric, IDR/overlap threshold, blacklist version.

## Troubleshooting

### Low FRiP (< 1% TF / < 5% histone)

The single most common ChIP failure mode. Causes in order of frequency:
1. Antibody specificity / lot issue -> re-validate with KO/KD
2. IP efficiency low -> check IP-Western
3. Fragment size wrong -> check sonication, library prep
4. Sequencing depth too shallow -> typically need ≥ 20M unique mapped
5. Hyper-ChIPable artifact dominance -> build custom blacklist and recompute

Increasing sequencing depth on a failed-FRiP sample does not fix the underlying issue.

### Low NSC/RSC (< 1.05 NSC, < 0.8 RSC)

Cross-correlation reflects strand-shift enrichment. Causes:
1. Weak ChIP signal -> check FRiP
2. Wrong fragment length estimate -> check `_cc.pdf` plot for clear peak
3. Phantom peak dominates (high signal at read length) -> bad library
4. phantompeakqualtools R compatibility issue -> use kundajelab fork

For broad histone marks (H3K27me3, H3K9me3), NSC marginal at 1.05 with strong FRiP > 15% is acceptable.

### Poor replicate concordance (Spearman < 0.6)

Causes:
1. Sample swap -> verify with sequencing index, library prep batch records
2. Batch effect -> check date of IP, library prep, sequencer
3. Biological variability -> expected for some cell states (differentiation time course)
4. One bad replicate -> check per-rep FRiP / NSC; drop and repeat

Do NOT average across replicates to fix poor concordance.

### IDR returns 0 reproducible peaks

1. Sorted by wrong column -> use `-k8,8nr` (p-value descending)
2. Library size imbalance > 2× -> downsample to common depth
3. One replicate dominated -> check per-rep peak count
4. Truly different conditions labeled as replicates -> check metadata

### Hyper-ChIPable region detected

1. Build custom blacklist from top-1% input signal
2. Recompute FRiP and peak counts after blacklist filter
3. For specific high-stakes loci (rRNA, mtDNA), require motif + KO/KD validation before publishing

## Related Skills

- chip-seq/peak-calling - Use QC outputs to set peak-calling parameters; FRiP/NSC inform threshold choice
- chip-seq/cut-and-run-tag - Different QC battery (spike-in % aligned, FRiP > 25%)
- chip-seq/spike-in-normalization - Spike-in QC (E. coli or Drosophila read depth and consistency)
- chip-seq/differential-binding - Replicate concordance required before differential testing
- atac-seq/atac-qc - Parallel ATAC QC (no input control; TSS enrichment instead of NSC)
- alignment-files/duplicate-handling - MarkDuplicates is upstream; NRF computed before dedup filter
- alignment-files/bam-statistics - Library size, mapping rate, flagstat
