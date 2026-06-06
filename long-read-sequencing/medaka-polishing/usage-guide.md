# Medaka Polishing - Usage Guide

## Overview
medaka polishes Oxford Nanopore draft assemblies to higher consensus accuracy, produces haploid variant calls (VCF) for microbial, mitochondrial, and viral samples, and generates amplicon/viral consensus sequences. It is a basecaller-model-specific neural consensus net: the single most important rule is that the model must match the basecaller (pore, chemistry, speed, mode, version) or the consensus silently degrades. medaka runs directly on the assembler output as a single pass - the old "Racon then medaka" chain is obsolete for contemporary models. medaka is ONT-only; PacBio HiFi/CLR must never be fed to it. For ONT diploid/germline small-variant calling, use Clair3; for read-level/human polishing, dorado polish is the emerging successor.

## Prerequisites
```bash
conda install -c bioconda medaka minimap2 samtools seqkit
# Models download on first use; a GPU accelerates inference but is not required.
medaka tools list_models    # see installed/available models
```

## Quick Start
Tell your AI agent what you want to do:
- "Polish my ONT Flye assembly with medaka, auto-detecting the model"
- "Polish my bacterial isolate assembly with the methylation-aware model"
- "Make a haploid consensus VCF for my viral sample"
- "Check whether medaka actually improved my assembly"

## Example Prompts

### Assembly polishing
> "I have an ONT-only Flye assembly and the Dorado-basecalled reads. Polish it with medaka, let it auto-detect the model from the reads, and run a single pass. Do not run Racon first."

### Bacterial isolate
> "Polish my native bacterial isolate assembly with medaka's methylation-aware bacterial model so methylation-motif errors are corrected."

### Model matching
> "These reads were basecalled with an older Dorado model. Tell me whether the medaka model will match, and whether I should re-basecall before polishing."

### Haploid consensus
> "Generate a haploid consensus and VCF for my mitochondrial Nanopore reads against the reference."

### Honest measurement
> "I polished with medaka. Measure whether it helped using held-out k-mers, not the reads I polished with."

## What the Agent Will Do
1. Establish what basecaller model and version produced the reads, so the medaka model matches (preferring auto-detection from the BAM).
2. Confirm the reads are ONT (and refuse to ONT-polish PacBio HiFi/CLR).
3. Run `medaka_consensus` as a single pass directly on the assembler output (no Racon pre-step).
4. For native bacterial isolates, add `--bacteria` to use the methylation-aware model.
5. For haploid samples, run `medaka_variant` to produce a VCF (apply it to the reference with `bcftools consensus` for a consensus FASTA).
6. Recommend a reference-free Merqury QV before/after measurement on held-out data, not the polishing reads.

## Tips
- Model matching is the number-one footgun: a mismatched model degrades the consensus silently. Auto-detect; treat a stale model name as a reason to re-basecall.
- Run medaka directly on the Flye output. Do not pre-run Racon, and do not run medaka twice.
- Never feed PacBio HiFi/CLR to medaka - it is ONT-only. Route PacBio to genome-assembly/assembly-polishing.
- medaka's "N changes made" is a risk signal, not a success signal. Validate with held-out Merqury QV.
- For ONT diploid small-variant calling, use Clair3; medaka diploid calling was deprecated in v2.
- v2 renamed subcommands: `inference` (was `consensus`), `sequence` (was `stitch`), `vcf` (was `variant`). The `medaka_consensus` wrapper kept its name.
- Make the assembly structurally complete before polishing; a missing plasmid makes medaka corrupt the chromosome.

## Related Skills

- basecalling - The basecaller model+version medaka must match
- clair3-variants - ONT diploid/germline small-variant calling
- long-read-alignment - minimap2 map-ont alignment medaka wraps internally
- genome-assembly/assembly-polishing - Polishing strategy, HiFi doctrine, Merqury QV design
- genome-assembly/long-read-assembly - Produces the draft medaka polishes
- genome-assembly/assembly-qc - QV/BUSCO before-vs-after measurement

## Resources
- [medaka GitHub](https://github.com/nanoporetech/medaka)
- [ARTIC fieldbioinformatics](https://github.com/artic-network/fieldbioinformatics) (medaka amplicon consensus)
