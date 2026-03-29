---
name: bio-chipseq-peak-calling
description: ChIP-seq peak calling using MACS3 and HOMER findPeaks. Call narrow peaks for transcription factors or broad peaks for histone modifications. Supports single-caller and multi-caller consensus approaches, input control, fragment size modeling, and various output formats. Use when calling peaks from ChIP-seq alignments.
tool_type: cli
primary_tool: macs3
---

## Version Compatibility

Reference examples tested with: MACS2 2.2+, MACS3 3.0+, HOMER 4.11+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Peak Calling

**"Call peaks from my ChIP-seq data"** → Identify significantly enriched regions (narrow peaks for TFs, broad peaks for histone marks) by comparing IP signal to input control.
- CLI (MACS3): `macs3 callpeak -t chip.bam -c input.bam -f BAM -g hs -n sample`
- CLI (HOMER): `makeTagDirectory tags/ chip.bam` then `findPeaks tags/ -style factor -i input_tags/ -o peaks.txt`

## Choosing a Peak Caller

MACS3 and HOMER use fundamentally different statistical approaches. MACS3 builds a dynamic local Poisson model (taking the maximum of genome-wide, 1kb, 5kb, and 10kb background estimates). HOMER applies three independent sequential filters: control enrichment, local enrichment, and clonal signal complexity. Neither is universally superior — the choice depends on the analysis context.

| Context | Recommended | Why |
|---------|-------------|-----|
| ENCODE compliance / IDR workflows | MACS3 | ENCODE standard; IDR designed around MACS2/SPP output |
| Broad marks (H3K27me3, H3K36me3) | MACS3 `--broad` | Two-tier Poisson produces gapped peaks with internal structure |
| Paired-end data | MACS3 with `-f BAMPE` | Native fragment pileup without read extension |
| Integrated motif + annotation workflow | HOMER | Peak calling, motif discovery, and annotation share one tag directory |
| Super-enhancer identification | HOMER `-style super` | Native ROSE-like algorithm from H3K27ac data |
| High-confidence peak sets | Both + consensus | Multi-caller intersection removes tool-specific false positives |
| Repetitive genomes or PCR artifacts | HOMER | Clonal filter (`-C`) catches false positives over repeats |

When both tools are available and runtime permits, running both and taking the intersection produces a higher-confidence peak set than either caller alone. See Multi-Caller Consensus below.

MACS3 is the actively developed successor to MACS2. Commands are identical except the binary name. MACS2 is in maintenance mode.

## Decision Framework

Before running any peak calling command, assess the data and choose parameters accordingly. The three key decisions are peak mode, model building strategy, and genome size.

### 1. Peak Mode: Narrow vs Broad

The choice depends on the biology of the target, not a preference:

| Target type | Examples | Mode | Why |
|-------------|----------|------|-----|
| Transcription factors | CTCF, p53, GATA1 | Narrow (default) | TFs bind discrete motif sites, producing sharp peaks |
| Active promoter marks | H3K4me3, H3K27ac | Narrow (default) | These marks localize to narrow regulatory elements |
| Elongation/body marks | H3K36me3, H3K79me2 | `--broad` | Deposited across gene bodies, forming wide domains |
| Repressive marks | H3K27me3, H3K9me3 | `--broad` | Spread across large chromatin domains |

If uncertain about the target's enrichment pattern, check published ENCODE data for that mark or search current literature before proceeding.

**HOMER note:** The MACS3 modes above do not map directly to HOMER styles. For HOMER, use `-style factor` only for transcription factors. All histone marks — including narrow marks like H3K4me3 and H3K27ac — perform better with `-style histone` (benchmarked in Omnipeak, NAR 2025). See the HOMER section below.

### 2. Model Building vs --nomodel

MACS3 estimates fragment size by finding paired plus/minus strand peaks. This requires at least 100 such pairs within the `--mfold` enrichment range. Use `--nomodel` when model building is expected to fail:

| Condition | Use --nomodel? | Why |
|-----------|---------------|-----|
| Whole-genome, >1M treatment reads | No — let MACS3 model | Enough data for reliable fragment size estimation |
| Single chromosome or small region | Yes | Too few enriched regions for 100 paired peaks |
| Low read count (<500k treatment) | Yes | Sparse signal makes modeling unreliable |
| ATAC-seq or DNase-seq | Yes, with `--extsize 150 --shift -75` | Open chromatin has no directional shift to model |
| Paired-end data (`-f BAMPE`) | N/A | Fragment size comes from mate pairs, no modeling needed |

When using `--nomodel`, set `--extsize` based on the mark type (see table below). When model building fails unexpectedly, try widening `--mfold` (e.g., `--mfold 3 50`) before falling back to `--nomodel`.

### 3. Fragment Size and Extension

When using `--nomodel`, the `--extsize` value should come from data when possible, not a generic default:

1. **Best:** Estimate from strand cross-correlation using phantompeakqualtools (SPP). This is the ENCODE standard — the same analysis produces QC metrics (NSC/RSC). See chipseq-qc for cross-correlation details.
2. **Good:** Estimate using `macs3 predictd -i chip.bam -g hs --outdir .` and read the fragment size from stderr. Requires sufficient enriched regions (same constraint as model building).
3. **Fallback:** Use mark-type defaults from the table below.

Fragment sizes vary across experiments — ENCODE H3K4me3 datasets show cross-correlation estimates from 95-245bp within a single experiment. Using an estimated value is always preferable to a default.

| Mark type | Fallback --extsize | Rationale |
|-----------|-------------------|-----------|
| TFs (CTCF, p53) | 150-200 | Typical ChIP fragment size range |
| H3K4me3, H3K27ac | 147 | Nucleosome core particle is 147bp; these marks are nucleosome-proximal |
| H3K4me1 | 200 | Slightly broader enhancer marks |
| H3K36me3, H3K79me2 | 200-300 | Gene body marks, broader signal |
| H3K27me3, H3K9me3 | 200 | Extension less critical since `--broad` links subpeaks |

The ENCODE pipeline always uses `--nomodel --shift 0 --extsize {fraglen}` where fraglen comes from cross-correlation, even when model building would succeed. This ensures consistency between peak calling and signal track generation.

### 4. Input Format

| File type | Extension | -f flag | Notes |
|-----------|-----------|---------|-------|
| Aligned BAM (SE) | .bam | BAM | Most common |
| Aligned BAM (PE) | .bam | BAMPE | Uses actual fragment sizes from mate pairs |
| tagAlign / BED6 | .tagAlign.gz, .bed.gz | BED | ENCODE standard format; stores chr, start, end, name, score, strand |
| SAM | .sam | SAM | Rarely used directly |

TagAlign is a BED6 format widely used by ENCODE. MACS3 reads it with `-f BED` (not a special tagAlign flag).

## Basic Peak Calling

```bash
macs3 callpeak -t chip.bam -c input.bam -f BAM -g hs -n sample --outdir peaks/
```

## Without Input Control

```bash
macs3 callpeak -t chip.bam -f BAM -g hs -n sample --outdir peaks/
```

Calling without input control uses genomic background only. This is less accurate — always prefer matched input/IgG when available.

## Narrow Peaks (TF, H3K4me3, H3K27ac)

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    -n sample_narrow \
    --outdir peaks/ \
    -q 0.05
```

## Broad Peaks (H3K36me3, H3K27me3, H3K9me3)

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    -n sample_broad \
    --outdir peaks/ \
    --broad \
    --broad-cutoff 0.1
```

## Paired-End Data

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAMPE \
    -g hs \
    -n sample_pe \
    --outdir peaks/
```

BAMPE uses actual fragment sizes from mate pairs — no model building or `--extsize` needed.

## tagAlign / BED Input

```bash
macs3 callpeak \
    -t treatment.tagAlign.gz \
    -c control.tagAlign.gz \
    -f BED \
    -g hs \
    -n sample \
    --outdir peaks/
```

MACS3 accepts gzipped tagAlign directly. No decompression needed.

## Multiple Replicates

```bash
macs3 callpeak \
    -t rep1.bam rep2.bam rep3.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    -n pooled \
    --outdir peaks/
```

## Genome Size

Use built-in shortcuts for model organisms, or provide a numeric value for non-model organisms, single-chromosome data, or targeted captures:

| Genome | Flag | Effective Size |
|--------|------|----------------|
| Human (whole) | hs | 2.7e9 |
| Mouse (whole) | mm | 1.87e9 |
| C. elegans | ce | 9e7 |
| D. melanogaster | dm | 1.2e8 |
| Human chr21 only | 46700000 | ~46.7M mappable bases |
| Custom/targeted | numeric | Sum of mappable bases in target regions |

For single-chromosome or targeted data, always provide the numeric effective genome size rather than the whole-genome shortcut. MACS computes genome-wide background as lambda_BG = (control_reads x fragment_size) / genome_size. Using `hs` (2.9e9) on chr21-only data (~46.7M) deflates lambda_BG by ~62x, making the background floor artificially low and allowing false positives in regions where local lambda is also low.

## Fixed Fragment Size (--nomodel)

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    --nomodel \
    --extsize 150 \
    -n sample \
    --outdir peaks/
```

See the decision framework above for when to use `--nomodel` and how to choose `--extsize`.

## Rescuing Failed Model Building

If model building fails with "needs at least 100 paired peaks":

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    --mfold 3 50 \
    -n sample \
    --outdir peaks/
```

Widening `--mfold` from the default `[5, 50]` to `[3, 50]` lowers the enrichment threshold for selecting model-building regions, which can recover enough paired peaks. If this still fails, use `--nomodel` with an appropriate `--extsize`.

## Generate Signal Tracks

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    -n sample \
    --outdir peaks/ \
    -B \
    --SPMR

sort -k1,1 -k2,2n peaks/sample_treat_pileup.bdg > peaks/sample.sorted.bdg
bedGraphToBigWig peaks/sample.sorted.bdg chrom.sizes peaks/sample.bw
```

## Cutoff Analysis

```bash
macs3 callpeak \
    -t chip.bam \
    -c input.bam \
    -f BAM \
    -g hs \
    --cutoff-analysis \
    -n sample \
    --outdir peaks/
```

Generates a table of peak counts at various q-value cutoffs. Useful for choosing a threshold when expected peak counts are unclear.

## HOMER Peak Calling

HOMER peak calling requires two steps: creating a tag directory (which estimates fragment size via strand autocorrelation), then calling peaks against it.

### Tag Directories

```bash
makeTagDirectory chip_tags/ chip.bam
makeTagDirectory input_tags/ input.bam
```

For tagAlign/BED input, specify format explicitly:

```bash
makeTagDirectory chip_tags/ treatment.tagAlign.gz -format bed
makeTagDirectory input_tags/ control.tagAlign.gz -format bed
```

Tag directories store pre-processed alignment data and compute QC metrics (fragment length, strand enrichment) automatically.

### Transcription Factors

```bash
findPeaks chip_tags/ -style factor -i input_tags/ -gsize 2.7e9 -o peaks.txt
```

`-style factor` uses fixed-width peaks with all three filters active (control enrichment `-F 4`, local enrichment `-L 4`, clonal filtering `-C 2`). Peak width is auto-estimated from tag autocorrelation. Use for point-source TF binding (CTCF, p53, GATA1).

### Histone Marks (H3K4me3, H3K27ac, H3K27me3, H3K36me3)

```bash
findPeaks chip_tags/ -style histone -i input_tags/ -gsize 2.7e9 -o regions.txt
```

`-style histone` uses variable-width region stitching (500bp building blocks stitched within 1000bp) and disables local enrichment filtering (`-L 0`). Use for **all histone marks**, including narrow marks like H3K4me3 and H3K27ac. Benchmarking (Omnipeak, NAR 2025) shows histone mode outperforms factor mode for H3K4me3 because it captures the variable-width enrichment around modified nucleosomes rather than clipping to a fixed-width center. For broad marks (H3K27me3, H3K36me3), `-L 0` is essential since local enrichment filtering would eliminate the signal being sought.

### Converting HOMER Output to BED

```bash
pos2bed.pl peaks.txt > peaks.bed
```

HOMER peak files use 1-indexed inclusive coordinates; `pos2bed.pl` converts to 0-indexed BED format. For manual conversion:

```bash
awk 'BEGIN{OFS="\t"} !/^#/ && NF>=5 {print $2, $3, $4, $1, $8}' peaks.txt > peaks.bed
```

### HOMER Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| -style | factor | Peak mode: factor, histone, super, groseq, tss, dnase |
| -i | none | Control tag directory |
| -gsize | 2e9 | Effective genome size |
| -size | auto | Peak size (auto for factor; 500 for histone) |
| -F | 4.0 | Fold enrichment over control |
| -L | 4.0 | Fold enrichment over local region (0 in histone mode) |
| -C | 2.0 | Clonal signal filter (0 to disable) |
| -fdr | 0.001 | FDR threshold |
| -o | stdout | Output file (`auto` writes to tag directory) |

## Multi-Caller Consensus

When both MACS3 and HOMER are available, running both and intersecting results is the recommended approach for final peak sets. MACS3's local Poisson model and HOMER's three-filter approach have different error profiles, so their intersection removes tool-specific false positives. For quick exploration with a single tool, MACS3 alone is sufficient.

```bash
bedtools intersect -a macs3_peaks.narrowPeak -b homer_peaks.bed -wa -u > consensus_peaks.bed
```

For lenient matching (peaks within 500bp), use `bedtools window`:

```bash
bedtools window -a macs3_peaks.narrowPeak -b homer_peaks.bed -w 500 | cut -f1-10 | sort -k1,1 -k2,2n | uniq > consensus_peaks.bed
```

**Single-caller vs consensus:** Single-caller (MACS3 alone) is appropriate for exploratory analysis, quick surveys, or when only one tool is installed. Consensus is preferred for final peak sets used in downstream motif analysis, functional enrichment, cross-condition comparison, or publication.

**Consensus vs IDR vs naive overlap:** These address different questions:
- **IDR** measures replicate reproducibility (same caller, multiple replicates). ENCODE standard for TF ChIP-seq.
- **Naive overlap** identifies peaks found across replicates. ENCODE standard for histone marks — IDR is too conservative for the variable dynamic range of histone signal.
- **Multi-caller consensus** measures algorithmic agreement (multiple callers, same data). Widely used for higher-confidence peak sets.

For publication-quality analysis, consider both replicate thresholding (IDR or naive overlap) and multi-caller consensus.

## Sanity-Checking Results

After peak calling, verify that results are biologically plausible before proceeding to downstream analysis:

| Target | Typical peak count (whole genome) | Typical peak width |
|--------|-----------------------------------|--------------------|
| TFs (CTCF, p53) | 10,000-80,000 | 200-500 bp |
| H3K4me3 | 20,000-50,000 | 500-2,000 bp |
| H3K27ac | 30,000-80,000 | 500-3,000 bp |
| H3K27me3 | 5,000-30,000 broad domains | 10-100+ kb |
| H3K36me3 | 10,000-30,000 broad domains | 5-50 kb |

For subset data (e.g., single chromosome), scale expectations proportionally — chr21 is ~1.5% of the human genome, so expect ~1.5% of whole-genome peak counts.

Red flags that indicate parameter problems:
- **Zero peaks:** Check genome size flag (using `hs` on subset data?), input format (`-f`), and whether input control is properly matched
- **Orders of magnitude too many peaks:** Control sample may be swapped with treatment, or q-value cutoff is too permissive
- **Orders of magnitude too few peaks:** Genome size may be too large for the data, or model building failed silently (check stderr)

When results are unexpected, re-read the MACS3 stderr output — it reports tag counts, redundancy rates, fragment size estimates, and model building status. Search current MACS3 documentation or literature if the issue is not covered here.

## Output Files

| File | Description |
|------|-------------|
| *_peaks.narrowPeak | Peak coordinates (BED6+4) |
| *_peaks.broadPeak | Broad peak coordinates |
| *_summits.bed | Peak summit positions |
| *_model.r | R script for model visualization |
| *_treat_pileup.bdg | Treatment signal (with -B) |
| *_control_lambda.bdg | Control signal (with -B) |

## narrowPeak Format

```
chr1  100  200  peak_1  100  .  5.2  10.5  8.3  50
```
Columns: chr, start, end, name, score, strand, signalValue, pValue, qValue, peak

## Convert narrowPeak to BED

Many downstream tools and pipelines expect simple BED format (chr, start, end, name, score). Extract the first 5 columns from narrowPeak:

```bash
cut -f1-5 peaks.narrowPeak > peaks.bed
```

## Filter Peaks

```bash
awk '$9 > 2' peaks.narrowPeak > peaks.filtered.narrowPeak  # -log10(q) > 2 means q < 0.01
sort -k7,7nr peaks.narrowPeak > peaks.sorted.narrowPeak
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| -t | required | Treatment file(s) |
| -c | none | Control file(s) |
| -f | AUTO | Format: BAM, BAMPE, BED, SAM |
| -g | hs | Effective genome size (shortcut or numeric) |
| -n | NA | Output prefix |
| -q | 0.05 | Q-value cutoff |
| -p | none | P-value cutoff (overrides -q) |
| --broad | false | Broad peak calling |
| --broad-cutoff | 0.1 | Q-value cutoff for broad regions |
| --nomodel | false | Skip fragment size modeling |
| --extsize | 200 | Read extension size (with --nomodel) |
| --mfold | 5 50 | Enrichment range for model building |
| -B | false | Generate bedGraph |
| --SPMR | false | Signal per million reads |
| --nolambda | false | Use local background only (skip genome-wide lambda) |

## ENCODE Pipeline Reference

The ENCODE ChIP-seq pipeline uses these MACS2 parameters as its standard:

```bash
macs2 callpeak -t chip.tagAlign.gz -c input.tagAlign.gz \
    -f BED -g {gsize} -p 1e-2 \
    --nomodel --shift 0 --extsize {fraglen} \
    --keep-dup all -B --SPMR
```

Key differences from default MACS3 settings: uses `-p 0.01` (p-value, not q-value) because IDR or naive overlap thresholding is applied downstream to filter peaks across replicates; uses `--keep-dup all` because duplicates are handled upstream during BAM processing; fragment length comes from phantompeakqualtools cross-correlation, not MACS model building. For histone marks, ENCODE applies naive overlap across replicates rather than IDR.

## Related Skills

- peak-annotation - Annotate peaks to genes (ChIPseeker or HOMER annotatePeaks.pl)
- motif-analysis - De novo and known motif enrichment (HOMER findMotifsGenome.pl, MEME)
- differential-binding - Compare peaks between conditions
- chipseq-qc - QC metrics including FRiP, NSC/RSC, IDR
- super-enhancers - Super-enhancer identification with ROSE or HOMER -style super
- alignment-files/sam-bam-basics - BAM file preparation
- chipseq-visualization - Visualize peaks and signal tracks
- genome-intervals/interval-arithmetic - BED intersection, merging, and window operations
