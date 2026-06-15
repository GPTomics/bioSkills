# Primer Specificity - Usage Guide

## Overview

This skill checks whether a PCR primer PAIR amplifies only the intended target across a whole genome (and transcriptome), using pair-aware in-silico PCR. The central idea it enforces: plain BLAST is the wrong tool, because specificity is a property of a predicted amplicon (forward and reverse primer both anchored, convergent, and within range), not of one primer's similarity to the genome, and the 3' terminus that governs priming is exactly what a similarity search ignores. It covers why intron-spanning RT-qPCR is defeated by processed pseudogenes (which force a genome search, not just the transcriptome), how to read a Primer-BLAST report without over-trusting an empty section, and why in-silico checking reduces but never replaces empirical validation. Use it to confirm primers are specific, screen for off-target amplicons, and avoid paralog/pseudogene/repeat amplification; route design to primer-basics and dimer checks to primer-validation.

## Prerequisites

```bash
# offline 3'-anchor prefilter
pip install primer3-py
# pair-aware in-silico PCR (pick one): MFEprimer (local), UCSC isPcr (local), or NCBI Primer-BLAST (web)
conda install -c bioconda mfeprimer blast
# plus a genome (and, for RT-qPCR, transcriptome) FASTA to search against
```

- Input is a chosen primer PAIR (and the intended amplicon location), plus the correct sequence database.
- For RT-qPCR, search the GENOME (with pseudogenes and alt/unplaced contigs), not the transcriptome only; intron-spanning does not escape processed pseudogenes.
- The in-silico tools need their binaries and a database; the genome step is not runnable without them. The included Python example (3'-anchor logic) is offline-runnable.
- A passing in-silico check is necessary but not sufficient: empirical validation (gradient PCR, single band, melt curve, sequencing) still licenses the assay.

## Quick Start

Tell your AI agent what you want to do:
- "Check whether this primer pair amplifies anything besides my target in the human genome"
- "My RT-qPCR primers span an exon junction; do they still hit a processed pseudogene?"
- "Run in-silico PCR for this pair on hg38 and tell me every predicted amplicon"
- "Is there a common SNP under the 3' end of either primer?"

## Example Prompts

### Genome-wide off-target check
> "Run pair-aware in-silico PCR for this forward/reverse pair against the human genome and report all predicted amplicons, their sizes, and which is the intended one. Flag any off-target with a stable 3' anchor."

### RT-qPCR pseudogene trap
> "My primers span the exon 2-3 junction of this gene. Confirm they are cDNA-specific by searching the genome for processed pseudogenes that carry the junction, not just the transcriptome."

### Reading a Primer-BLAST report
> "Here is my Primer-BLAST result with an empty unintended-products section. Explain whether that proves the primers are unique, and how the default 3'-mismatch filter could be hiding an off-target."

### Allele dropout risk
> "Check the 3' ends of both primers against gnomAD common variants so I do not get allele dropout in some individuals."

## What the Agent Will Do

1. Confirm the assay type and pick the correct database (genome for genotyping; genome + transcriptome for RT-qPCR, including alt/unplaced contigs).
2. Run pair-aware in-silico PCR (MFEprimer, isPcr, or Primer-BLAST) on the pair and enumerate predicted amplicons.
3. Require exactly one intended amplicon of expected size and no qualifying off-target; weight 3'-anchored sites most.
4. Use the offline 3'-anchor prefilter (calc_end_stability) to explain why an overall-similar site may or may not actually prime.
5. Check primer 3' ends against common variants to flag allele-dropout risk.
6. State the remaining empirical validation (gradient PCR, single band, melt curve, sequencing) and route a redesign to primer-basics.

## Tips

- Use a pair-aware tool (MFEprimer / isPcr / Primer-BLAST), never plain BLAST, for the final specificity decision; reserve blastn-short (word_size 7, dust off) for a quick single-primer repeat scan.
- For RT-qPCR, search the genome; processed pseudogenes carry the exon junction and amplify like cDNA.
- Treat an empty Primer-BLAST unintended-products section as "none passed the filter," not "none exist"; loosen the mismatch settings to stress-test.
- Match the database to the assay and include alt/unplaced contigs; a "primary assembly only" search hides off-targets.
- Weight the 3' end: a site with high overall similarity but a disrupted 3' anchor will not prime, and one with internal mismatches but an intact 3' anchor will.
- In-silico checking reduces bad designs; confirm a quantitative assay empirically before trusting it.

## Related Skills

- primer-basics - Design (or redesign) primers when specificity fails
- primer-validation - Intramolecular dimers/hairpins of the chosen oligos
- qpcr-primers - qPCR assays where specificity and gDNA exclusion are mandatory
- read-alignment/bwa-alignment - Align candidate amplicons / reads to a genome
- database-access/blast-searches - Build/query BLAST databases for candidate finding
