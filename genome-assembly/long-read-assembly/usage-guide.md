# Long-Read Assembly (Noisy Reads) - Usage Guide

## Overview

This skill builds de novo genome assemblies from noisy long reads - Oxford Nanopore (R9.4.1, R10.4.1, Dorado simplex/duplex) and PacBio CLR. It covers the modern long-read assemblers (Flye, Canu, NextDenovo, Shasta, Raven, wtdbg2, miniasm) and the bacterial multi-assembler consensus tools (Trycycler, Autocycler). The two decisions it exists to get right are: matching the assembler input flag to the basecaller error regime (the wrong flag silently over-collapses repeats while raising N50), and treating the assembler as the halfway point because a raw long-read assembly is contiguous but low-QV and is not finished until it is polished and its QV is measured. PacBio HiFi reads are accurate, not noisy - use the hifi-assembly skill for those.

## Prerequisites

```bash
conda install -c bioconda flye canu nextdenovo shasta raven-assembler wtdbg miniasm minimap2 purge_dups porechop_abi trycycler
```

Genome-size/coverage estimation (GenomeScope2 + jellyfish/KMC) is best run on accurate short reads, not raw ONT, because errors inflate unique k-mers. Polishing (Medaka/Racon) and QV assessment (Merqury) are separate downstream tools - see the assembly-polishing and assembly-qc skills.

## Quick Start

Tell your AI agent what you want to do:
- "Assemble my ONT R10 reads with Flye"
- "I have Nanopore reads basecalled with Dorado SUP - which Flye flag do I use?"
- "Assemble a bacterial genome reliably from Nanopore reads"
- "My assembly is twice the expected size - what went wrong?"
- "Assemble a large plant genome from noisy long reads"

## Example Prompts

### Matching the flag to the basecaller era
> "I have Nanopore reads from an R10.4.1 flow cell basecalled with Dorado SUP. Assemble them with Flye using the correct input mode and explain why --nano-raw would be wrong here."

### Bacterial reliability
> "Assemble my bacterial isolate from ONT reads as reliably as possible. Use a multi-assembler consensus rather than trusting one assembler, and tell me which contigs are circular."

### Large eukaryote
> "Assemble a 3 Gb repetitive plant genome from noisy long reads. Pick the most memory-efficient contiguity-first assembler and downsample to the longest reads for the initial assembly."

### Diagnosing duplication
> "My Flye assembly of a diploid genome is 1.8x the expected genome size and BUSCO-Duplicated is high. Diagnose whether this is haplotig duplication and run purge_dups, checking the coverage cutoffs before trusting them."

### Pre-assembly cleanup
> "Before assembling my Nanopore reads, check for and remove mid-read adapter chimeras so I don't get structural mis-joins."

## What the Agent Will Do

1. Establish the platform and basecaller era (ONT R9 vs R10/Dorado vs PacBio CLR) - this is an assembly parameter, not optional metadata.
2. De-chimerize internal adapters with Porechop_ABI and confirm read length/quality and coverage (~30-60x) and read-N50 are adequate.
3. Select an assembler and the input mode matching the error regime (Flye `--nano-hq` for modern ONT; Trycycler/Autocycler for bacterial reliability; NextDenovo for large repetitive genomes).
4. Run the assembly, downsampling to the longest reads (`--asm-coverage`) on deep large-genome data.
5. Check for haplotig false-duplication (size vs expected, half-coverage depth peak, BUSCO-Duplicated) and run purge_dups if warranted.
6. State that the assembly is unpolished and low-QV, and hand off to polishing then QV assessment - the deliverable is a polished assembly with a measured QV, not a high-N50 FASTA.

## Tips

- The most dangerous mistake is telling the assembler the reads are noisier than they are: `--nano-raw` on R10/Dorado data over-collapses repeats and *raises* N50 as the assembly gets worse - nothing crashes.
- Record pore + kit + basecaller model before assembling; if the chemistry is unknown the flag cannot be chosen.
- A high N50 is not a finished assembly. Raw consensus is Q20-Q30; indels frameshift genes. Polish, then measure QV (Merqury Q40 ~ 1 error/10 kb).
- For bacteria, prefer a multi-assembler consensus (Trycycler/Autocycler) over agonizing about which single assembler - any one makes structural errors polishing cannot fix.
- Coverage and read-N50 are non-substitutable: 200x of 5 kb reads gives great consensus and a fragmented assembly. Spend effort on read length, not just gigabases.
- Filter conservatively - the longest reads are often the lowest quality, and aggressive length filters discard the reads that span hard repeats.
- wtdbg2 needs the separate `wtpoa-cns` consensus step; miniasm does no consensus at all and its output must be polished.
- PacBio CLR is legacy (superseded by HiFi). Do not recommend generating new CLR; treat existing CLR like noisy ONT and polish hard.
- Pre-2022 tutorials assume R9 and tend to default to `--nano-raw` plus an unnecessary short-read polish - check the chemistry the advice was written for before applying it.

## Related Skills

- hifi-assembly - PacBio HiFi (Q30+) phased haplotype-resolved assembly with hifiasm
- assembly-polishing - Polishes the low-QV contigs this skill produces
- assembly-qc - QUAST/BUSCO/Merqury QV assessment; the QV is the real deliverable
- long-read-sequencing/long-read-qc - Read QC and conservative filtering before assembly
- long-read-sequencing/basecalling - Basecaller era and model that set the assembler flag
- read-qc/contamination-screening - Screen reads for contamination before assembling
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
