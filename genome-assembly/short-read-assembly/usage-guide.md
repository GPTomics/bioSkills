# Short-Read Assembly - Usage Guide

## Overview

This skill assembles a genome de novo from Illumina short reads. SPAdes is the field default for bacterial, fungal, and small-eukaryotic isolates; MEGAHIT handles huge or memory-constrained datasets and metagenomes; Unicycler finishes bacterial genomes (short-only or hybrid); Platanus targets highly heterozygous diploids. The defining reality is that short reads cannot span repeats longer than the read or paired insert, so a short-read assembly is structurally fragmented at every long repeat - its contiguity ceiling is set by the genome's repeat structure and the input DNA quality, not by the sequencing depth or the assembler. The skill profiles the reads first (GenomeScope2 on a k-mer histogram) to know that ceiling before touching parameters, and reports NG50/auN/BUSCO rather than a vanity N50. For anything that needs to be finished or complete, it routes to long-read assembly.

## Prerequisites

```bash
conda install -c bioconda spades megahit unicycler abyss genomescope2 kmc
# SPAdes 4.0+ requires Python >=3.8
# GenomeScope2 also needs R; smudgeplot ships alongside it
```

## Quick Start

Tell your AI agent what you want to do:
- "Profile my reads with GenomeScope2 before I assemble, to estimate genome size and heterozygosity"
- "Assemble my bacterial isolate from paired-end Illumina reads with SPAdes"
- "Assemble a finished, circularized bacterial genome from Illumina plus Nanopore reads"
- "My eukaryote assembly is twice the expected size - tell me whether that is heterozygosity"
- "Assemble a huge metagenomic dataset with limited RAM"

## Example Prompts

### Pre-Assembly Profiling
> "Build a k=21 k-mer histogram with KMC and fit GenomeScope2 with ploidy 2, then tell me the estimated genome size, heterozygosity, and repeat content, and whether short reads can give me a usable assembly for this genome."

### Bacterial Isolate
> "Assemble my paired-end reads R1.fastq.gz and R2.fastq.gz with SPAdes in isolate mode at 16 threads, then report contig count, NG50 against the expected genome size, and whether the contig count looks like rRNA-operon fragmentation rather than a problem."

### Hybrid Finishing
> "Combine my Illumina reads with ONT long reads in Unicycler to get a finished, circularized chromosome plus plasmids, rotated to dnaA."

### Heterozygosity Diagnosis
> "My SPAdes assembly is 1.8x the GenomeScope estimate and BUSCO-Duplicated is high - tell me whether this is un-purged haplotigs and what to do about it."

### Memory-Constrained / Metagenome
> "Assemble this large short-read metagenome with MEGAHIT because I am RAM-limited, and explain why MEGAHIT instead of metaSPAdes here."

## What the Agent Will Do

1. Profile the reads first - build a k-mer histogram (KMC) and fit GenomeScope2 to estimate genome size, heterozygosity, repeat content, and ploidy, setting the contiguity/size ceiling before any parameter choice.
2. Confirm reads are QC'd - adapter/light-quality trimmed (not over-trimmed) and screened for contamination.
3. Select the assembler and mode from the data: SPAdes `--isolate` for a clean isolate, `--careful` only for small genomes, `--sc` for single-cell, MEGAHIT for huge/low-memory, Unicycler for finishing/hybrid, Platanus-allee for high heterozygosity.
4. Run the assembly with multi-k (auto-selected), avoiding hand-picked single k and the `--careful`-on-a-eukaryote anti-pattern, and using `--only-assembler` for heterozygous/pooled data to preserve rare alleles.
5. Report contig metrics (NG50, auN, contig count, total size vs the GenomeScope estimate) and flag heterozygosity inflation, GC dropout, or contamination.
6. Hand off to assembly-qc (BUSCO, Merqury QV) and, where a finished genome is needed, to long-read assembly.

## Tips

- Run GenomeScope2 before assembling - it indicates whether short reads can even produce the requested assembly, and gives the expected size to report against.
- More coverage past ~50-100x (isolate) or ~50-60x (eukaryote) does not lengthen contigs; it only adds error-correction and, above ~150x, amplifies GC and duplicate bias.
- Do not hand-pick a single k; let multi-k auto-select. If overriding, give a small-to-large list of odd values <= read length.
- `--careful` is small-genome isolate-only and incompatible with `--meta`/`--rna`; never `--careful` a eukaryote, and avoid the redundant `--isolate --careful` combination.
- Do not error-correct heterozygous, pooled, or metagenomic data (use `--only-assembler`) - k-mer-spectrum correction erases the rare alleles or strains that are the signal.
- An assembly ~1.5-2x the expected size with high BUSCO-Duplicated is un-purged haplotigs, not a large genome; purge or move to long reads.
- A clean bacterial isolate fragmenting into 30-100 contigs is normal (rRNA operons and IS elements); for a single finished contig, use Unicycler hybrid or long reads.
- Report NG50 + auN + BUSCO + a reference-free correctness check, never N50 alone - N50 rises when sequence is thrown away and inflates on mis-joins.
- Pin the assembler version in methods; SPAdes default behavior changed across 3.x to 4.0.

## Related Skills

- assembly-qc - NG50 + auN + BUSCO + Merqury QV; never report N50 alone
- assembly-polishing - Base-accuracy polishing belongs here, not in `--careful` reflexes
- metagenome-assembly - Community co-assembly and MAG recovery (metaSPAdes/MEGAHIT)
- long-read-assembly - The real fix for repeats short reads cannot span
- hifi-assembly - Phased, repeat-spanning assembly for heterozygous/complex eukaryotes
- read-qc/quality-reports - Garbage-in caps assembly quality; QC reads before assembling
- read-qc/adapter-trimming - Light adapter/quality trimming; over-trimming fragments the graph
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
