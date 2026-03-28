# Peak Calling - Usage Guide

## Overview

Call peaks from ChIP-seq data using MACS3 and/or HOMER to identify transcription factor binding sites or histone modification regions. MACS3 uses a dynamic local Poisson model; HOMER applies three independent filters (control enrichment, local enrichment, clonal complexity). Both are widely used, and running both for a multi-caller consensus produces higher-confidence results than either alone.

## Prerequisites

```bash
# MACS3
conda install -c bioconda macs3
# or: pip install macs3

# HOMER (optional, for multi-caller consensus or motif integration)
conda install -c bioconda homer
```

## Quick Start

Tell your AI agent what you want to do:
- "Call peaks from my ChIP-seq BAM file with input control"
- "Identify H3K27ac enriched regions in my treatment vs input"
- "Call peaks from my tagAlign files for H3K4me3 on chr21"
- "Run both MACS3 and HOMER and give me a high-confidence consensus peak set"
- "Find transcription factor binding sites — I have a small dataset, may need --nomodel"

## Example Prompts

### Narrow Peaks (TFs, H3K4me3)
> "Call narrow peaks for my transcription factor ChIP-seq"

> "Run MACS3 on my H3K27ac ChIP-seq with input control"

### Broad Peaks (H3K27me3, H3K36me3)
> "Call broad peaks for my H3K27me3 ChIP-seq data"

> "Identify H3K9me3 domains using broad peak calling"

### Multi-Caller
> "Call peaks with both MACS3 and HOMER, then intersect for high-confidence peaks"

> "I want a robust peak set — use multiple callers"

### Subset or Small Data
> "Call H3K4me3 peaks from my chr21-only tagAlign files — about 400k reads"

> "I have a small pilot ChIP-seq dataset, call peaks with appropriate settings"

### ATAC-seq
> "Call peaks from my ATAC-seq data"

### Troubleshooting
> "My peak calling found no peaks, help me troubleshoot"

> "MACS3 model building failed — not enough paired peaks"

> "I'm getting too many peaks, how can I increase stringency?"

## What the Agent Will Do

1. Assess the data: identify input format (BAM, tagAlign/BED), sequencing layout (SE vs PE), target mark type, genome scope, and read count
2. Choose peak caller(s) based on available tools and analysis goals (MACS3, HOMER, or both)
3. Choose peak mode (narrow vs broad) based on the target mark's biology
4. Decide whether MACS3 model building is feasible or `--nomodel` is needed based on read count and genome scope
5. Set genome size appropriately — numeric value for subset data, shortcut for whole genomes
6. Run peak calling, review stderr/output for warnings
7. If both callers are available, consider intersecting results for a high-confidence consensus
8. Sanity-check peak counts against expected ranges for the mark type and data scope
9. Convert output to the requested format (e.g., narrowPeak to BED) if needed

## Troubleshooting

### Model building fails ("needs at least 100 paired peaks")

MACS3 estimates fragment size by finding enriched regions with reads on both strands. It needs at least 100 such regions within the `--mfold` enrichment range. Common causes of failure:

- **Single-chromosome or targeted data:** Too few enriched regions genome-wide. Use `--nomodel` with an appropriate `--extsize`.
- **Low read count:** Sparse signal means fewer detectable enriched regions. Use `--nomodel`.
- **Very narrow enrichment range:** Try widening `--mfold` from `[5, 50]` to `[3, 50]` before falling back to `--nomodel`.

### Zero peaks called

- **Wrong genome size:** Using `hs` (2.7e9) on single-chromosome data inflates the background model. Use the actual effective size (e.g., 46700000 for human chr21).
- **Wrong format flag:** Ensure `-f` matches the input: BAM for .bam, BAMPE for paired-end .bam, BED for .tagAlign/.bed files.
- **Swapped files:** Verify treatment and control are not reversed.
- **Too few reads:** Very low read counts may not produce statistically significant peaks at `-q 0.05`. Try `-q 0.1` for exploratory analysis, but interpret with caution.

### Too many peaks

- **Permissive threshold:** Tighten from `-q 0.05` to `-q 0.01` or `-q 0.001`.
- **Missing input control:** Without control, genomic biases (mappability, GC content) produce false positives. Add a matched input sample.
- **Swapped files:** Check that the enriched sample is `-t` and the control is `-c`.

### Unexpected peak widths

- **Narrow peaks for a broad mark:** Add `--broad` for H3K27me3, H3K9me3, H3K36me3.
- **Broad peaks for a TF:** Remove `--broad` — TFs produce sharp peaks by nature.

When troubleshooting does not resolve the issue, read the MACS3 stderr output carefully (it reports tag counts, redundancy, fragment size, and model status) and search current MACS3 documentation or literature for guidance.

## Tips

- Always use input/IgG control when available — peaks without control have higher false positive rates
- Use `--broad` (MACS3) or `-style histone` (HOMER) for histone marks that form wide domains
- Default q-value threshold is 0.05; use `-q 0.01` for higher stringency
- For subset data (single chromosome, targeted capture), always provide numeric genome size to both MACS3 (`-g`) and HOMER (`-gsize`)
- TagAlign files are BED6 format — use `-f BED` with MACS3 or `-format bed` with HOMER's makeTagDirectory
- MACS3 reads gzipped files directly — no need to decompress .tagAlign.gz or .bed.gz
- When both MACS3 and HOMER are installed, running both and intersecting peaks is a low-cost way to improve confidence
- Check peak counts against known ranges: H3K4me3 typically gives 20,000-50,000 peaks genome-wide; scale proportionally for subset data
- HOMER's tag directory computes fragment size and QC metrics automatically — check `tagInfo.txt` for diagnostics

## Related Skills

- peak-annotation - Annotate peaks to genes (ChIPseeker or HOMER annotatePeaks.pl)
- motif-analysis - De novo and known motif enrichment (HOMER findMotifsGenome.pl, MEME)
- differential-binding - Compare peaks between conditions
- chipseq-qc - QC metrics including FRiP, NSC/RSC, IDR
- super-enhancers - Super-enhancer identification with ROSE or HOMER -style super
- alignment-files/sam-bam-basics - BAM file preparation
- chipseq-visualization - Visualize peaks and signal tracks
- genome-intervals/interval-arithmetic - BED intersection, merging, and window operations
