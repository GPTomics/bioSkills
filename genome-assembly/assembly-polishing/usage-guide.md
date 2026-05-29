# Assembly Polishing - Usage Guide

## Overview

Polishing raises an assembly's *consensus accuracy* (per-base QV) by correcting residual substitutions and, above all, the homopolymer and low-complexity indels characteristic of noisy long reads. It does not change contiguity, join contigs, or fill N-gaps - those are scaffolding and gap-filling. The single most important judgment this skill encodes is *when not to polish*: HiFi and high-accuracy ONT assemblies start near-perfect, and mapping-based polishers routinely lower QV and collapse haplotypes on them. The correct workflow is read-type-matched (ONT self-polish with Racon then medaka; hybrid with Polypolish + pypolca; HiFi usually nothing), and every polish is validated with reference-free, held-out Merqury k-mers before and after - never the reads polished with, and never BUSCO.

## Prerequisites

```bash
# Long-read polishing
conda install -c bioconda racon medaka minimap2 samtools
# Short-read / hybrid polishing
conda install -c bioconda polypolish pypolca bwa pilon
# Measuring the gain (QV)
conda install -c bioconda merqury meryl
```

Notes: the medaka consensus model must match the basecaller/chemistry/version (`medaka tools list_models`); ONT increasingly steers toward `dorado polish`. Merqury needs a meryl k-mer DB built from an independent read set (k from `best_k.sh`). Pilon needs substantial heap (~1 GB per Mb of genome).

## Quick Start

Tell your AI agent what you want to do:
- "Should I polish my HiFi assembly?"
- "Polish my ONT Flye assembly with Racon then medaka"
- "Polish my hybrid bacterial assembly with Polypolish and pypolca"
- "Measure whether polishing actually improved my assembly QV"

## Example Prompts

### Deciding whether to polish
> "I have a hifiasm assembly from PacBio HiFi reads. Should I polish it with Pilon, and if not, how do I confirm it doesn't need polishing?"
> "My assembly is fragmented into hundreds of contigs - will polishing help, or is that a different problem?"

### ONT polishing
> "Polish my Flye ONT assembly: run two to four rounds of Racon, measure QV each round and stop at the plateau, then a single medaka pass with the model that matches my R10.4.1 sup basecalling."
> "My medaka output has a lower QV than the input - help me figure out if the model was wrong."

### Hybrid and short-read polishing
> "I have an ONT draft and 30x Illumina reads. Run the long-read polish then Polypolish plus pypolca, using the Bouras depth tiers."
> "Polish my heterozygous diploid assembly with short reads without collapsing the two haplotypes."

### Measuring the gain
> "Build a Merqury meryl database from my held-out Illumina reads and report QV before and after polishing so I know if it helped."

## What the Agent Will Do

1. Establish the read type and starting QV, and check whether the request is actually a polishing problem (vs scaffolding/gap-filling).
2. For HiFi or high-accuracy ONT, default to not polishing and measure Merqury QV first; only proceed if there is a real deficit.
3. Select a read-type-matched chain: Racon -> medaka (or dorado polish) for ONT; Polypolish + pypolca for hybrid/short reads; NextPolish2/DeepPolisher for HiFi deficits; Hapo-G for heterozygous genomes.
4. Map reads to the draft (bwa for short reads, minimap2 for long reads), confirming the medaka model matches the basecaller.
5. Run the polisher, iterating Racon 1-4 rounds and measuring QV per round, stopping at the plateau; run medaka exactly once.
6. Build a meryl k-mer DB from an independent / different-platform read set and report Merqury QV before and after; revert if QV dropped.

## Tips

- The default answer to "should I polish?" is increasingly "probably not, and definitely measure first." Reflexive Racon -> medaka -> Pilon on every assembly encodes 2018 practice.
- The medaka model mismatch is the #1 silent footgun: a wrong model degrades the consensus with no warning. Let medaka auto-detect from the FASTQ when possible; otherwise confirm the model in `medaka tools list_models`.
- "Changes made" is a risk signal, not a success signal - a polisher making thousands of edits to a HiFi assembly is alarming, not reassuring.
- Validate with held-out / different-platform k-mers; measuring with the polishing reads is circular and always looks good.
- For repeat-rich genomes prefer Polypolish (uses all alignments via `bwa mem -a`) over Pilon (best-placement mismaps).
- Run medaka exactly once after the Racon rounds; running it twice signals a misunderstanding of the chain.
- Keep the lanes separate: fragmentation is scaffolding, N-gaps are gap-filling, frameshifts/low QV/homopolymer indels are polishing.

## Related Skills

- long-read-assembly - Produces the contiguous-but-error-prone contigs this skill polishes
- short-read-assembly - Source of Illumina reads for hybrid/short-read polishing
- hifi-assembly - HiFi assemblies that usually should NOT be polished
- assembly-qc - Merqury QV before/after is the polishing stop signal
- read-alignment/bwa-alignment - Map short reads to the draft for Polypolish/Pilon
- long-read-sequencing/long-read-alignment - minimap2 mapping of long reads for Racon
- long-read-sequencing/medaka-polishing - medaka mechanics and model tables; this skill owns the strategy
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
