# CLIP-seq QC - Usage Guide

## Overview

Comprehensive QC for CLIP libraries (eCLIP, iCLIP, iCLIP2, PAR-CLIP). Five gates: preprocessing retention -> alignment rate -> library complexity -> FRiP/IP enrichment -> replicate reproducibility (IDR). Each gate must pass for the library to be usable. ENCODE eCLIP thresholds: >= 1M unique fragments per rep, FRiP >= 0.005 (narrow RBPs), IDR rescue + self-consistency both < 2. SMInput is the canonical control - NOT IgG, NOT empty beads, NOT RNA-seq.

## Prerequisites

```bash
conda install -c bioconda preseq picard samtools bedtools deeptools idr multiqc rseqc fastp
pip install pysam
```

## Quick Start

Tell your AI agent:
- "Run library complexity with preseq lc_extrap; flag if predicted unique < 1M at 100M reads"
- "Compute FRiP against stringent peaks; ENCODE narrow-binding minimum 0.005"
- "Run IDR true-replicate and pseudo-replicate; report rescue and self-consistency ratios"
- "Read distribution metagene - HuR should be >50% 3' UTR"
- "Diagnose why my CLIP library has 95% PCR duplication"
- "Why is my FRiP so low - failed IP or wrong RBP class threshold?"
- "Antibody sanity check: do top peaks match expected biology?"

## Example Prompts

### Library Complexity

> "preseq lc_extrap on pre-dedup BAM; what unique-fragment count is predicted at 100M reads?"

> "preseq curve flattens early - over-amplified library; nothing to do but re-prep"

### FRiP

> "FRiP against stringent peak set; expected 0.01-0.10 for splicing factor"

> "FRiP < 0.005 - check antibody; not an atypical-binding RBP"

### IDR

> "IDR with --rank 5 (log2 FC) on signal-sorted replicate beds; threshold 0.05 true rep, 0.10 pseudoreplicate"

> "Rescue ratio > 2 - one replicate failed; identify and re-sequence"

### Diagnostics

> "Top peaks are dominated by GAPDH and ACTB - failed IP or expression bias?"

> "60% PCR duplication - this is normal CLIP, not low complexity"

> "Read distribution shows 50% intergenic - rRNA contamination or wrong TxDb?"

### Antibody Sanity

> "Verify antibody on KD lysate WB before claiming the CLIP worked"

> "Top motif AU-rich and GO terms generic - antibody non-specific"

## What the Agent Will Do

1. Inspect preprocessing logs (cutadapt retention >= 70%)
2. STAR Log.final.out (unique alignment >= 60% eCLIP, 70% iCLIP)
3. `preseq lc_extrap` on pre-dedup BAM (>= 1M unique predicted)
4. FRiP calculation against stringent peaks (>= 0.005 narrow-binding)
5. log2(IP/SMInput) global enrichment check
6. IDR true-replicate and pseudo-replicate; ENCODE rescue/self-consistency rules
7. RSeQC read_distribution + geneBody_coverage for metagene QC
8. Antibody sanity: top peaks match expected RBP biology
9. MultiQC for consolidated report

## Tips

- **Five gates in order.** Don't skip earlier ones.
- **40-70% duplication is normal CLIP.** UMI dedup recovers unique molecules.
- **preseq on pre-dedup BAM.** Post-dedup gives nonsense.
- **FRiP is RBP-class-specific.** Atypical-binding RBPs exempt.
- **IDR ranks by log2 FC or p-value, not raw score.** CLIPper score is tied.
- **SMInput is the canonical control.** Not IgG, not RNA-seq.
- **Antibody is the single most common failure point.** 30-50% commercial "IP-grade" antibodies fail.
- **Validate antibody on KD lysate WB.** Knockdown -> WB signal decrease is the test.
- **MultiQC unifies the report.** Standard practice for cross-lab work.
- **Failed at gate 3 (complexity) = no rescue.** Re-prep.

## Related Skills

- clip-seq/clip-preprocessing - Gate 1
- clip-seq/clip-alignment - Gate 2
- clip-seq/clip-peak-calling - Gates 4-5
- clip-seq/differential-clip - Differential QC
- read-qc/quality-reports - General QC
- read-qc/contamination-screening - Contamination
- chip-seq/chipseq-qc - ChIP-seq analogue
