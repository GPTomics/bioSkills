# Nucleosome Positioning - Usage Guide

## Overview

Map nucleosome center positions, occupancy, and fuzziness from ATAC-seq fragment-size patterns. Covers V-plot interpretation, +1 nucleosome calling at TSSs, NucleoATAC vs ATACseqQC vs DANPOS3 vs scprinter selection, NRL estimation, and differential positioning between conditions. Provides the structural complement to peak-level (atac-peak-calling) and TF-level (footprinting) analyses.

## Prerequisites

```bash
# Dedicated NucleoATAC env (pinned because last release 2018)
conda create -n nucleoatac python=3.7 numpy=1.18 scipy=1.5 pysam
conda run -n nucleoatac pip install nucleoatac

conda install -c bioconda samtools
# DANPOS3 (python danpos.py) installs from github.com/sklasfeld/DANPOS3 (bioconda `danpos` is DANPOS2)
# scPrinter installs from source: git clone https://github.com/buenrostrolab/scPrinter && cd scPrinter && pip install ./
pip install pysam pyBigWig matplotlib
```

```r
BiocManager::install(c('ATACseqQC', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
                       'BSgenome.Hsapiens.UCSC.hg38', 'MotifDb', 'motifmatchr'))
```

Inputs: deduplicated, MAPQ-filtered, chrM-stripped paired-end BAM with >= 30M nuclear reads.

## Quick Start

Tell your AI agent what you want to do:
- "Run NucleoATAC on accessible regions to call per-base nucleosome occupancy"
- "Generate a V-plot at TSSs to verify nucleosome positioning is recoverable"
- "Estimate the nucleosome repeat length (NRL) from fragment-size autocorrelation"
- "Call +1 nucleosomes downstream of every protein-coding TSS"
- "Run DANPOS3 dpos with ATAC-tuned parameters for differential positioning"
- "Use scprinter for single-cell nucleosome positioning per cell type"

## Example Prompts

### Per-Region Nucleosome Calling
> "Run NucleoATAC on consensus peaks merged with bedtools (regions >= 500 bp). Use a dedicated conda env with Python 3.7 because NucleoATAC is pinned to old NumPy. Report nucpos.bed and occupancy bigWig."

### V-Plot Diagnostic
> "Generate a V-plot at TSSs (+/- 1 kb) to assess whether nucleosome positioning is recoverable. Classic V or W shape = positioning visible; flat horizontal band = library over-transposed and positioning is gone."

### +1 Nucleosome Calling
> "Define gene bodies starting at TSS+50 to TSS+1000 (where +1 nucleosome lives in metazoa); run NucleoATAC; for each gene, report the most-downstream nucleosome center as the +1 position."

### Differential Positioning
> "Compare two conditions with DANPOS3 dpos using `--width 145 --smooth_width 80` (ATAC-tuned). Filter to nucleosome shifts >= 30 bp at FDR < 0.05."

### NRL Estimation
> "Calculate the autocorrelation of cumulative Tn5 cut signal across a 2 kb window; the first non-zero peak position is the NRL. Compare against expected NRL for the cell type."

### Single-Cell Nucleosomes
> "Run scprinter on the scATAC fragment file to get per-cell nucleosome activity tracks; aggregate by cluster and plot V-plots per cluster."

## What the Agent Will Do

1. Verify input BAM is paired-end (mono-nuc analysis requires it)
2. Generate V-plot diagnostic at TSSs (or feature of choice)
3. If V-plot shows clear V/W -- proceed with calling. If flat band -- abort with diagnosis.
4. Define analysis regions (consensus peaks merged, or TSS-flanking, or whole genome)
5. Choose tool (NucleoATAC for ATAC-specific calls; DANPOS3 for differential; scprinter for single-cell)
6. Run nucleosome calling; output nucleosome positions, occupancy, fuzziness
7. For differential: DANPOS3 dpos with ATAC parameters, filter at >= 30 bp shift
8. Optionally: estimate NRL from autocorrelation; report per-region positioning summary
9. Annotate +1 / -1 / +2 / etc. nucleosomes relative to TSS

## Tool Selection Quick Reference

| Goal | Tool |
|------|------|
| Per-base occupancy track | NucleoATAC (caveat: unmaintained); scprinter alternative |
| V-plot at features | ATACseqQC vPlot |
| Differential positioning between conditions | DANPOS3 dpos (ATAC-tuned params) |
| +1 calling | NucleoATAC + downstream filter |
| Single-cell | scprinter |
| NRL estimation | Custom autocorrelation script |

## V-Plot Pattern Quick Reference

| Pattern | Meaning |
|---------|---------|
| Classic V at center | NFR with flanking nucleosomes; positioning recoverable |
| W (two V's) | Strong NFR + clear +1 / -1 |
| Inverted V (mountain) | Fragment fully enclosed; possibly nucleosome-bound TF |
| Flat horizontal band | No positioning info; library over-transposed |
| 10.4 bp helical sub-peaks | High-quality library with helical phasing visible |

## Fragment-Size Class Quick Reference (Buenrostro 2013)

| Class | Range | Use |
|-------|-------|-----|
| NFR | < 100 bp | Footprinting input; TF binding |
| Mono-nuc | 180-247 bp | Standard nucleosome positioning |
| Di-nuc | 315-473 bp | Phasing, NRL estimation |
| Tri-nuc | 558-615 bp | Heterochromatin, dense phasing |

## Tips

- V-plot is the most useful diagnostic before nucleosome calling. Always plot it first; flat-band patterns mean positioning is not recoverable and downstream calling is invalid.
- NucleoATAC is unmaintained since 2018 but remains the canonical ATAC nucleosome caller. Pin Python to 3.7 in a dedicated env if installing.
- Mono-nucleosome window is 180-247 bp (NOT 147 bp). The protected nucleosomal DNA is 147 bp, but the fragment includes flanking linker.
- NRL is species- and cell-type-specific (165 in yeast, ~196-200 in human somatic, 211 in cortical neurons). Default mono-nuc window may miss some libraries; verify with fragment-size distribution.
- +1 nucleosome is ~50-60 bp downstream of TSS in metazoa (~120 bp in yeast). It's the most-studied positioning feature; its absence in aggregate suggests TSS annotation is wrong.
- DANPOS3's MNase defaults are wrong for ATAC. Use `--width 145 --smooth_width 80`.
- scprinter is GPU-recommended; CPU runs are slow on large datasets.
- Pioneer TFs (FOXA1, GATA, OCT4) produce asymmetric V-plots because one DNA face is on the histone. This is biology, not artefact.
- Differential positioning calls (DANPOS3 dpos) require >= 30 bp shift at FDR < 0.05 to be biologically meaningful; smaller shifts are within fragment-size noise.
- Nucleosome fuzziness 20-50 bp is "well-positioned" in metazoan; > 100 bp is effectively unpositioned.

## Related Skills

- atac-seq/atac-qc - Fragment-size periodicity QC
- atac-seq/atac-peak-calling - HMM-based nucleosome-aware peak calling
- atac-seq/footprinting - Per-TF flanking nucleosome analysis
- atac-seq/single-cell-atac - scprinter for sc nucleosome positioning
- chip-seq/peak-annotation - Annotate nucleosome positions to genes
- alignment-files/bam-statistics - Insert-size statistics upstream
