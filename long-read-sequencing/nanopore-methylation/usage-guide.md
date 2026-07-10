# Nanopore Methylation Calling - Usage Guide

## Overview
This skill calls DNA base modifications (5mC, 5hmC, 6mA, 4mC) directly from Oxford Nanopore and PacBio HiFi long reads and piles them into per-site methylation with modkit (or pb-CpG-tools for PacBio). The fact that governs the whole pipeline is that methylation is a basecalling decision: modifications are inferred from raw signal at basecall time and stored as MM/ML SAM tags, so a BAM basecalled without a modification model can never yield methylation, and the tags silently die in a normal alignment workflow unless preserved. The skill covers the modBAM pipeline, the MM/ML tag spec, modkit's 18-column bedMethyl and its 10th-percentile auto-threshold, 5mC vs 5hmC handling relative to bisulfite, phased allele-specific methylation, and the PacBio path. RNA modifications are out of scope (see epitranscriptomics).

## Prerequisites
```bash
conda install -c bioconda modkit samtools minimap2 htslib
# PacBio path: pbmm2, pb-CpG-tools, jasmine
# Methylation must be requested at basecall time (see basecalling):
#   dorado basecaller sup,5mCG_5hmCG pod5/ > calls.bam
```

## Quick Start
Tell your AI agent what you want to do:
- "Call 5mC CpG methylation from my modBAM and make a bedMethyl"
- "Align my Dorado methylation BAM without losing the MM/ML tags"
- "Get per-haplotype allele-specific methylation"
- "Compare methylation between two samples"

## Example Prompts

### Tag check and pileup
> "Check whether my BAM actually has MM/ML methylation tags, then pile up 5mC at CpG sites with strands combined into a bedMethyl."

### Tag preservation through alignment
> "I have a Dorado modification BAM. Align it to the genome without dropping the MM/ML tags, and confirm they survived before I run modkit."

### Comparison to WGBS
> "I want to compare my Nanopore methylation to WGBS. Combine 5mC and 5hmC so the comparison is fair, and tell me why."

### Allele-specific methylation
> "Phase my sample, haplotag the BAM, and give me per-haplotype methylation at imprinted loci."

### PacBio
> "My data is PacBio HiFi with kinetics. Call 5mC and make a bedMethyl - and do not suggest primrose."

## What the Agent Will Do
1. Confirm the BAM carries MM/ML tags (and warn that if it does not, the only fix is re-basecalling from POD5).
2. Preserve the tags through alignment (dorado aligner, or `samtools fastq -T MM,ML | minimap2 -y -Y`) and verify they survived.
3. Pile up with modkit into an 18-column bedMethyl (auto-thresholded at the 10th percentile, not 0.5).
4. For WGBS comparison, combine 5mC and 5hmC; for 5hmC biology, keep them split.
5. For allele-specific methylation, phase and haplotag first, then pileup with `--partition-tag HP`.
6. For PacBio, run ccs --hifi-kinetics -> jasmine -> pb-CpG-tools (or modkit), and hand DMR statistics to methylation-analysis.

## Tips
- The first move is always to check the tags exist: `samtools view in.bam | head | grep MM:Z`. No tags means re-basecall, not re-analyze.
- The tags silently die unless you use `dorado aligner` or `samtools fastq -T MM,ML | minimap2 -y -Y`; always re-check after alignment.
- modkit auto-thresholds at the 10th percentile of the ML distribution, not 0.5; there is no `--min-coverage` flag on pileup - filter on Nvalid_cov afterward.
- WGBS conflates 5mC and 5hmC; combine them (`--combine-mods`/`--preset traditional`) to compare, keep them split to study 5hmC.
- The MM `?` modifier means no-call (Nnocall), not unmethylated; misreading it deflates methylation.
- For count-based DMR, pass Nmod and Nvalid_cov, never percent_modified.
- PacBio primrose is deprecated; use jasmine (or Revio on-instrument 5mC).

## Related Skills

- basecalling - Methylation must be requested at basecall time
- long-read-alignment - Carry MM/ML through alignment (`-y -Y`)
- haplotype-phasing - Phase + haplotag for allele-specific methylation
- methylation-analysis/dmr-detection - DMR statistics downstream
- methylation-analysis/methylkit-analysis - methylKit differential methylation
- epitranscriptomics/m6anet-analysis - Direct-RNA m6A (out of scope here)

## Resources
- [modkit](https://github.com/nanoporetech/modkit)
- [modkit docs](https://nanoporetech.github.io/modkit/)
- [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)
- [SAM tags spec (MM/ML)](https://samtools.github.io/hts-specs/SAMtags.pdf)
