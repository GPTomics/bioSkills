# Basecalling - Usage Guide

## Overview
Basecalling converts raw Oxford Nanopore signal (POD5, or legacy FAST5) into reads using Dorado, ONT's production basecaller. The model and accuracy tier chosen are not a neutral preprocessing choice - they are baked irreversibly into the output: methylation must be requested at basecall time or it cannot be recovered, and downstream polishing/variant models must match the basecaller model and version. This skill covers model selection, modified-base calling, duplex, demultiplexing, trimming, and HERRO read correction. Guppy is end-of-life; Dorado is used for all new work. PacBio HiFi reads are basecalled on-instrument and are out of scope here.

## Prerequisites
```bash
# Dorado (from ONT) - download the prebuilt binary
# https://github.com/nanoporetech/dorado

# POD5 tooling and post-basecall QC/filtering
pip install pod5
conda install -c bioconda chopper nanoplot samtools

# A CUDA GPU is strongly recommended; sup is GPU-bound. Models download on first use:
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0
```

## Quick Start
Tell your AI agent what you want to do:
- "Basecall my POD5 files with Dorado using the sup model"
- "Basecall and call 5mC CpG methylation, and remind me to keep the POD5"
- "Demultiplex my barcoded Nanopore run"
- "Correct my Nanopore reads with HERRO before assembly"

## Example Prompts

### Model selection and chemistry matching
> "I have R10.4.1 POD5 data and want to call variants downstream. Pick the right Dorado model, pin its version, and basecall at super-accuracy so the downstream Clair3 model matches."

### Methylation at basecall time
> "Basecall my POD5 with 5mC and 5hmC CpG methylation so I can run modkit later. Confirm the MM/ML tags are present and tell me whether I still need the POD5."

### Duplex
> "These reads were prepared for duplex. Run duplex basecalling, but tell me the realistic duplex yield and how to avoid double-counting the simplex parents."

### Demultiplexing
> "Basecall and demultiplex my SQK-NBD114-24 barcoded run. Avoid the trimming-before-demux trap so barcodes are not lost."

### Read correction for assembly
> "Correct my super-accuracy R10 reads with HERRO so I can build a phased diploid assembly."

## What the Agent Will Do
1. Identify flow-cell chemistry from the POD5 metadata (R10.4.1 vs legacy R9.4.1; DNA vs RNA004).
2. Select the accuracy tier (sup for analysis, hac for routine, fast for previews) and pin the model version.
3. Decide whether to request modified bases now (and warn that mods cannot be recovered later).
4. Run `dorado basecaller` (or `duplex`), producing an unaligned BAM with quality and any MM/ML tags.
5. For barcoded runs, basecall with `--no-trim` then `dorado demux`.
6. Optionally convert to FASTQ, filter with chopper, and QC with NanoPlot, or correct with `dorado correct` for assembly.

## Tips
- Use sup for anything you will analyze; hac for routine throughput; fast only for live/QC previews.
- Keep the POD5. It is the only way to re-call methylation, run duplex, or re-basecall with a newer model.
- Propagate the exact model string (e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`) to medaka and Clair3 so their models match.
- Re-basecall an entire cohort with one model version; mixing versions creates a batch effect.
- The per-read mean qscore is an overconfident posterior, not a calibrated accuracy - use it for relative filtering, not as ground truth. Real accuracy needs alignment to a reference (see long-read-qc).
- cDNA libraries cannot call RNA modifications; only direct RNA (RNA004) can.
- Convert FAST5 to POD5 first; basecalling FAST5 directly is slow.

## Related Skills

- long-read-qc - Read length/quality and run-health QC after basecalling
- nanopore-methylation - Pile up the MM/ML tags requested here into bedMethyl
- long-read-alignment - Map reads, carrying MM/ML tags through with `-y`
- medaka-polishing - Consensus model that must match this basecaller model+version
- clair3-variants - Variant model that must match this basecaller model+version
- genome-assembly/long-read-assembly - Assemble the (optionally HERRO-corrected) reads
- workflows/longread-sv-pipeline - End-to-end basecall -> align -> SV call

## Resources
- [Dorado GitHub](https://github.com/nanoporetech/dorado)
- [Dorado documentation](https://software-docs.nanoporetech.com/dorado/)
- [POD5 format](https://pod5-file-format.readthedocs.io/)
