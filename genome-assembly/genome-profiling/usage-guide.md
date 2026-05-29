# Genome Profiling - Usage Guide

## Overview

Genome profiling is the reference-free, pre-assembly step that turns raw sequencing reads into a k-mer frequency spectrum and models it to estimate the genome's size, heterozygosity, repeat content, and ploidy. Those four numbers are not trivia: they set the expected assembly size (the NG50 denominator and the haplotig sanity check), the genome-size parameter for assemblers like Flye, the homozygous-coverage and purge settings for hifiasm, and the up-front decision about whether short reads can even produce the assembly being requested. This skill counts k-mers with KMC or Jellyfish, models the spectrum with GenomeScope2, and infers ploidy independently with Smudgeplot, while flagging the common traps that make the estimate wrong: too-low coverage, noisy ONT reads, contamination spikes, a too-small k, and confusing k-mer coverage with sequencing coverage.

## Prerequisites

Install the counters and modeling tools (Bioconda):

```bash
conda install -c bioconda kmc jellyfish genomescope2 smudgeplot merqury kat
```

- KMC or Jellyfish produces the k-mer histogram; only one is required.
- GenomeScope2 needs R; the Bioconda package pulls it in.
- Merqury provides `best_k.sh` (k selection) and meryl; `$MERQURY` should point at the install.
- Accurate reads are required for counting: Illumina or PacBio HiFi. Noisy ONT reads cannot be profiled reliably.
- A rough expected genome size (even an order of magnitude) helps choose k and validate the result.

## Quick Start

Tell your AI agent what you want to do:
- "Profile my genome from Illumina reads before I assemble"
- "Estimate genome size and heterozygosity from these reads"
- "Build a k-mer spectrum and run GenomeScope2"
- "Figure out the ploidy of this non-model organism"
- "Check whether my assembly is 2x too big because of haplotigs"

## Example Prompts

### Pre-assembly profiling
> "I have 60x Illumina reads from an unknown beetle. Count 21-mers with KMC, run GenomeScope2 as diploid, and tell me the estimated genome size, heterozygosity, and repeat content, plus whether short reads will struggle with the heterozygosity."

### Choosing k and counting
> "My expected genome size is about 12 Mb. Use Merqury best_k.sh to pick k, count k-mers with that k, and run GenomeScope2 so I can set the genome-size parameter for Flye."

### Ploidy inference
> "GenomeScope2 fits poorly as diploid for this plant. Run Smudgeplot on the het k-mer pairs and tell me whether it looks triploid or tetraploid, and whether that disagrees with GenomeScope."

### Diagnosing a bad spectrum
> "The GenomeScope size estimate is twice what I expected and the spectrum has an extra peak. Help me decide whether this is contamination, an organelle spike, or a too-small k."

### Sanity-checking an assembly
> "My HiFi assembly is 1.8x the GenomeScope estimate and BUSCO Duplicated is high. Confirm whether this is uncollapsed haplotigs and what to do before reporting the genome size."

## What the Agent Will Do

1. Confirm the reads are accurate (Illumina/HiFi) and adequate depth (homozygous peak >= ~15-20x); route noisy-ONT-only cases to assemble-then-estimate.
2. Choose k with Merqury `best_k.sh` from the expected genome size (k=21 for ~Gb genomes), using the same k for counting and modeling.
3. Count k-mers with KMC (or Jellyfish with `-C`), keeping singletons (`-ci1`) and capping the tail (`-cs10000`) to build the histogram.
4. Run GenomeScope2 with the matching `-k` and a starting ploidy `-p`, then read the model summary: haploid genome length, heterozygosity, and percent unique (repeat content).
5. Interpret the spectrum: single vs het (AB) peak, repeat tail, organelle/contamination spikes; distinguish k-mer coverage (lambda) from sequencing depth.
6. Run Smudgeplot to infer ploidy from het k-mer-pair coverage ratios and reconcile it with GenomeScope's `-p`.
7. Hand the estimate downstream: the NG50 denominator and haplotig sanity check (assembly-qc), Flye `-g`, hifiasm `--hom-cov`/purge level, and the short-vs-long-read decision.

## Tips

- Always use the same k for counting and for GenomeScope2's `-k`; a mismatch silently mis-scales the model.
- Count canonical k-mers: KMC does this by default, Jellyfish needs `-C`. WGS is unstranded, so forward and reverse k-mers must be merged.
- Read the homozygous peak position off the plot rather than assuming it equals sequencing depth; they differ by the `(L-k+1)/L` factor.
- A distinct peak at roughly half the homozygous coverage is the heterozygous (AB) signal; its area scales with heterozygosity.
- A spike at very high multiplicity is usually organelle or high-copy repeat DNA, not the nuclear genome; cap the histogram tail so it does not dominate the fit.
- If GenomeScope's diploid fit is poor, treat ploidy as unknown and let Smudgeplot decide before forcing a `-p`.
- Never report an assembly size larger than the profiled estimate as the genome size; a 1.5-2x excess is almost always uncollapsed haplotigs.

## Related Skills

- short-read-assembly - High het from the profile predicts short-read fragmentation and haplotig inflation
- hifi-assembly - The profiled size/het sets hifiasm --hom-cov and the purge level
- long-read-assembly - The profiled genome size feeds Flye -g and read-chemistry expectations
- assembly-qc - The estimate is the NG50 denominator and the 1.5-2x haplotig sanity check
- read-qc/quality-reports - QC and contamination-screen reads before profiling and assembling
- workflows/genome-assembly-pipeline - Profiling is the first step of the end-to-end assembly workflow
