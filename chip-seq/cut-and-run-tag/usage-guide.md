# CUT&RUN / CUT&Tag - Usage Guide

## Overview

Analyze CUT&RUN (Skene & Henikoff 2017) and CUT&Tag (Kaya-Okur et al 2019) chromatin profiling experiments. These methods use antibody-tethered MNase (CUT&RUN) or Tn5 (CUT&Tag) to profile protein-DNA binding with 10-100× lower background than traditional ChIP, 100-1000× lower cell input, and automatic E. coli spike-in carryover from bacterially-produced enzymes. Embeds SEACR vs MACS2 peak calling logic (per 2025 Bioinformatics btaf375 benchmark), pA vs pAG-Tn5 chimera selection, characteristic 25-75 bp CUT&Tag fragment signatures, IgG-only control logic, and `--keep-dup all` requirement for CUT&Tag.

## Prerequisites

```bash
# Aligner + alignment utilities
conda install -c bioconda bowtie2 samtools bedtools

# Peak callers
git clone https://github.com/FredHutch/SEACR.git    # SEACR_1.4+
conda install -c bioconda macs2 macs3

# Optional alternatives
conda install -c bioconda gopeaks
pip install lanceotron

# Adapter trimming (important for 25-75 bp CUT&Tag fragments)
conda install -c bioconda cutadapt

# QC and visualization
conda install -c bioconda deeptools
```

## Quick Start

Tell the agent what to do:
- "Align my CUT&Tag paired-end reads with bowtie2 using Henikoff lab parameters"
- "Run SEACR with stringent + norm + IgG control on H3K4me3 CUT&Tag bedgraphs"
- "Call CUT&RUN peaks with both SEACR and MACS2 and take the consensus (per 2025 benchmark)"
- "Verify my CUT&Tag fragment-size distribution shows the expected 25-75 bp Tn5 peak"
- "Compute E. coli spike-in fractions per sample and calculate scaling factors"
- "Why are my H3K27me3 CUT&RUN peaks called as narrow by MACS2 instead of broad?"
- "Compare CUT&Tag peak counts to a published ChIP-seq for the same TF; they shouldn't match exactly"

## Example Prompts

### Standard CUT&Tag pipeline
> "Run the Henikoff lab CUT&Tag pipeline: bowtie2 align with `--local --very-sensitive --no-mixed --no-discordant -I 10 -X 700`, sort, generate fragment bedgraph, run SEACR with `norm stringent` and IgG control. Use `--keep-dup all` if also running MACS2."

### CUT&RUN with H3K27me3
> "H3K27me3 CUT&RUN is broad. Use SEACR (better than MACS2 for broad CUT&RUN per the 2025 benchmark) with `norm stringent` mode and IgG control."

### Consensus peak set
> "Run MACS2 with `-f BAMPE --keep-dup all` AND SEACR with stringent + norm + IgG. Take the intersection for a publication-grade peak set."

### Spike-in normalization
> "Align reads to a combined hg38 + E. coli K12 genome. Report E. coli read counts per sample. Verify target samples have 0.5-2% E. coli alignment and IgG has 2-5%. Compute spike-in scaling factors as min(E_coli) / per_sample_E_coli."

### Fragment size diagnostic
> "Extract fragment-size histogram from my CUT&Tag BAM. Expected: sharp peak at 25-75 bp (Tn5 staggered cuts) + secondary at ~150-200 bp (mono-nucleosomal for H3K4me3 etc). If flat above 200 bp, suspect over-amplification or poor enzyme activity."

### Mouse antibody pitfall
> "My primary antibody is mouse anti-MYC. I'm using pA-Tn5 from the Henikoff protocol. The signal is weak; diagnose. (Answer: pA binds rabbit IgG preferentially; switch to pAG-Tn5 which binds both rabbit and mouse.)"

### Quality vs ChIP comparison
> "I have CUT&Tag H3K4me3 with 50k peaks and published ChIP-seq with 200k peaks. Is this a CUT&Tag failure? (Answer: usually not; CUT&Tag's higher specificity often gives fewer but cleaner peaks. Verify with motif enrichment and known target gene coverage.)"

## What the Agent Will Do

1. **Alignment**: bowtie2 with Henikoff parameters; combined hg38 + E. coli genome for spike-in counting
2. **Pre-peak-calling QC:**
   - Fragment-size distribution (expect 25-75 bp peak for CUT&Tag)
   - E. coli spike-in fraction (0.5-2% target, 2-5% IgG)
   - Library complexity (low PCR cycles, NRF may be lower than ChIP)
3. **Peak calling decision:**
   - Sharp marks (H3K4me3, TFs): MACS2 with `-f BAMPE --keep-dup all`
   - Broad marks (H3K27me3): SEACR with `norm stringent` + IgG
   - Publication-grade: MACS2 + SEACR consensus (per btaf375 2025)
4. **Generate bedgraph for SEACR** from paired-end fragments (NOT bigWig)
5. **Run SEACR** with appropriate mode:
   - `norm stringent` + IgG (recommended default)
   - Top-X% + `non stringent` if no IgG
   - `non` only if upstream spike-in scaling applied
6. **Spike-in normalization** if cross-condition comparison: scaling factor = min(E_coli) / per_sample_E_coli; apply via `bamCoverage --scaleFactor` for tracks; via DiffBind for differential
7. **Output**: peak BED, bedgraph, optional spike-in scaled bigWig, fragment-size diagnostic plot
8. **Document**: chimera variant (pA vs pAG vs pA-MNase), bead type, digitonin concentration, antibody catalog + lot, sequencing depth, alignment %, peak-call command exactly

## Tips

- **CUT&Tag duplicates contain biology**: always use `--keep-dup all` in MACS2 for CUT&Tag. CUT&RUN allows standard dedup (more PCR cycles typically).
- **`-f BAMPE` is mandatory for MACS2 on CUT&Tag/CUT&RUN PE data.** `-f BAM` modeling fails on 25-75 bp fragments.
- **SEACR mode matters:**
  - `norm` (recommended) scales target to IgG
  - `non` requires upstream spike-in scaling
  - `stringent` (recommended) takes top-half of signal distribution
  - `relaxed` includes the full distribution (use only for very sparse signal)
- **E. coli carryover is variable across enzyme production batches.** For high-stakes cross-condition claims, add deliberate Drosophila spike-in (ChIP-Rx style).
- **pAG-Tn5 binds both rabbit AND mouse IgG.** pA-Tn5 prefers rabbit. Choose chimera per primary antibody species.
- **Adapter readthrough is common for short fragments.** Trim aggressively with cutadapt `-e 0.1 -O 5 --minimum-length 25`. Consider 50 bp PE sequencing instead of 150 bp for CUT&Tag.
- **Digitonin concentration must be titrated per cell type.** 0.05% default fails for some cell lines.
- **Sequencing depth requirement is much lower**: 3-10M reads typical, vs 20-50M for ChIP. Sequencing more than ~10M reads often hits saturation.
- **FRiP > 0.25 is normal for CUT&Tag**; > 0.10 for TF CUT&Tag. Higher than ChIP because background is so low.
- **5,000-100,000 cells is sufficient input.** Standard ChIP needs 1-10M cells. This is a major advantage for clinical samples.

## Troubleshooting

### Very low E. coli alignment (< 0.1%)

Spike-in carryover lost. Either the enzyme was over-purified or the prep batch had unusual carryover. Cross-condition normalization unreliable; add deliberate Drosophila spike-in next time.

### Very high E. coli alignment (> 10%)

Target ChIP failed; sample is mostly E. coli. Re-do experiment with fresh enzyme and titrated antibody. Check antibody validation.

### Peak count much lower than ChIP for same antibody

Often expected; CUT&Tag's lower background reveals signal-vs-noise that ChIP's higher background obscures. Verify with:
1. Motif enrichment at peaks (should match known TF if validated antibody)
2. Coverage at known target genes
3. ENCODE / public CUT&Tag for same target if available

If genuinely failing, check digitonin concentration, antibody titration, and bead type.

### No clean 25-75 bp peak in CUT&Tag fragment size

Indicates poor Tn5 activity or over-tagmentation. Re-do with titrated Tn5; verify enzyme storage and freshness.

### MACS2 peak calling produces few peaks

1. `-f BAM` instead of `-f BAMPE` → switch
2. `--keep-dup auto` removed CUT&Tag biology → use `--keep-dup all`
3. `-q 0.05` too lenient for very-low-background CUT&Tag → tighten to `-q 0.01`
4. Library too shallow → 3-5M reads minimum

### Adapter content high in FASTQ

100-150 bp reads on 25-75 bp fragments read through adapters. Trim aggressively or use 50 bp PE.

### SEACR output peaks overlap genome regions outside chromatin

SEACR processes the input bedgraph as-is. Filter blacklist before generating bedgraph; ensure chromosome naming matches between target and IgG.

### IgG signal nearly as strong as target

Antibody specificity failure OR over-permeabilization OR too much pA-Tn5 added.

## Related Skills

- chip-seq/peak-calling - Traditional ChIP peak calling (MACS3 + IDR vs naive overlap)
- chip-seq/chipseq-qc - QC battery (different thresholds for CUT&RUN/Tag)
- chip-seq/spike-in-normalization - Deliberate Drosophila spike-in beyond E. coli carryover
- chip-seq/differential-binding - DiffBind / csaw differential
- chip-seq/peak-annotation - Annotate CUT&RUN/Tag peaks
- chip-seq/super-enhancers - SE calling on CUT&Tag H3K27ac
- read-qc/adapter-trimming - Aggressive adapter trimming for short fragments
- read-alignment/bowtie2-alignment - Bowtie2 with Henikoff parameters
- alignment-files/sam-bam-basics - BAM preparation
