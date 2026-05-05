# Alignment Amplicon Clipping

Trim PCR primers from aligned reads in amplicon-panel BAMs.

## Overview

Amplicon panels (SARS-CoV-2 ARTIC, hereditary cancer panels, ctDNA hot-spot panels, fusion panels, 16S rRNA) use designed PCR primers for enrichment. The primer footprint at read 5' ends does not represent biological sequence; without trimming, downstream variant calls suppress true variants under primers and falsely confirm reference at primer locations.

This skill covers post-alignment primer clipping with `samtools ampliconclip`, the required tag-repair steps (`fixmate`, `calmd`), and the relationship to duplicate marking (do not run markdup on amplicon BAMs without UMIs).

## Prerequisites

- samtools 1.11 or later (ampliconclip introduced in 1.11)
- Coordinate-sorted, indexed amplicon BAM
- Primer BED file with strand information in column 6 (for `--strand` mode)

```bash
conda install -c bioconda samtools
# Verify ampliconclip is available
samtools ampliconclip --help
```

## Quick Start

Tell your AI agent what you want to do:
- "Trim ARTIC SARS-CoV-2 primers from my amplicon BAM"
- "Soft-clip primers using my primer BED"
- "Why is markdup marking everything as duplicate on my amplicon panel?"
- "Repair the MD/NM tags after primer clipping"
- "Hard-clip primers for archival storage"

## Example Prompts

### SARS-CoV-2 / ARTIC

> "I have a SARS-CoV-2 ARTIC BAM aligned with bwa-mem. Soft-clip primers using artic_v3 primer BED, then prepare for consensus calling."

> "Compare samtools ampliconclip vs iVar trim for the ARTIC SARS-CoV-2 pipeline."

### Cancer Panels

> "ctDNA hot-spot panel with TruSeq UMIs. Trim primers, then run UMI consensus rather than markdup."

> "Hereditary cancer panel from a Roche AVENIO design. Why is my variant caller missing variants at primer-overlapping positions?"

### General Amplicon Workflows

> "My amplicon BAM ran through markdup and now flagstat shows 99% duplicates. What went wrong?"

> "After ampliconclip, my downstream mpileup is reporting wrong reference bases. What tag did clipping invalidate?"

## What the Agent Will Do

1. Inspect the BAM to confirm it is amplicon-style (read coordinates clustered at amplicon start/end positions; check `samtools view -H` for `@PG` indicating an amplicon panel kit).
2. Locate or generate the primer BED file matching the kit.
3. Run `samtools ampliconclip --both-ends --strand --soft-clip -b primers.bed`.
4. Re-fixmate and re-calmd to repair invalidated tags.
5. Confirm with a region inspection (`samtools view -h clipped.bam chr1:start-end`) that primer footprint has soft-clip CIGAR ops.

## Decision Points

- **Soft-clip vs hard-clip**: default to soft-clip (reversible). Hard-clip only when archiving and the trimmed bases will never be needed.
- **`--strand` on or off**: on for primer-strand-specific trimming (default for amplicon panels with single-stranded primers); off only if both strands carry primer-derived sequence.
- **`--both-ends` on or off**: on when reads can read through the entire amplicon and primers can appear at either 5' or 3'. Standard for short-amplicon panels (<300 bp).
- **Re-calmd**: required if downstream uses BAQ in `bcftools mpileup`, IGV mismatch visualization, or any tool that reads NM/MD.

## Tips

- Do NOT run `samtools markdup` on amplicon BAMs without UMIs -- amplicon reads share start coordinates by design and markdup will erase the dataset. See duplicate-handling for the assay-aware decision matrix.
- For UMI amplicon panels (Twist UMI, IDT xGen UMI, Roche AVENIO): use `fgbio GroupReadsByUmi` -> `CallMolecularConsensusReads` instead of markdup.
- For ARTIC SARS-CoV-2 specifically, modern pipelines often pair ampliconclip with `samtools consensus --config illumina --ambig` for IUPAC heterozygote handling.
- Strand bias at every amplicon end usually means `--strand` was forgotten.
- After clipping, mpileup with default BAQ may be wrong (MD tag stale); always re-calmd.
- For long-read amplicon (ONT viral, PacBio HiFi 16S): tool support varies; `samtools ampliconclip` works but verify CIGAR semantics with the long-read team's recommended pipeline.

## Related Skills

- duplicate-handling - Why amplicon BAMs need ampliconclip, not markdup
- alignment-filtering - Post-clip filtering for variant calling
- alignment-sorting - Re-sort after fixmate
- pileup-generation - mpileup flags for amplicon
- reference-operations - Consensus generation from amplicon BAMs
- read-qc/quality-reports - Pre-alignment QC
