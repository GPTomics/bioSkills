# Long-Read SV Pipeline - Usage Guide

## Overview

This workflow detects structural variants (deletions, insertions, inversions, duplications, translocations) from Oxford Nanopore or PacBio HiFi long-read data, and it exists for one physical reason: a single long read spans the SV and both flanks in one molecule, so it resolves the insertions and repeat-mediated events short reads cannot. It is an orchestration skill - it makes the stage-to-stage decisions (basecalling, alignment preset, caller, cohort merge, benchmark) and delegates the SV mechanism itself to the component skills it chains.

## Prerequisites

```bash
conda install -c bioconda minimap2 samtools sniffles cutesv pbsv nanoplot bcftools truvari
# ONT basecalling (Dorado) and assembly-based calling (dipcall) install separately - see their skills.
```

## Quick Start

Tell your AI agent what you want to do:
- "Detect structural variants from my Nanopore data"
- "Run the long-read SV pipeline on my PacBio HiFi reads"
- "Find insertions and deletions my short-read caller missed"
- "Joint-call SVs across my long-read cohort"
- "Benchmark my SV calls against GIAB HG002"

## Example Prompts

### Platform and preset
> "My reads are ONT R10.4 - which minimap2 preset and SV caller should I use?"

> "Should I switch from short-read to long-read SV calling for insertions?"

### SV calling
> "Call SVs from my aligned long reads with Sniffles2 and a tandem-repeat BED"

> "Use cuteSV with the right parameters for my PacBio HiFi data"

> "Call the highest-quality SV set from my phased diploid assembly"

### Cohort and benchmarking
> "Joint-genotype SVs across my three long-read samples"

> "Benchmark my SV calls against GIAB Tier 1 and CMRG with Truvari"

> "Report an SV F1 that is actually reproducible"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Reads | POD5/FAST5 (ONT) or FASTQ/BAM (HiFi) | ONT signal to basecall, or HiFi reads |
| Reference | FASTA | Reference genome (GRCh38 or T2T-CHM13) |
| Tandem-repeat BED | BED | reference-matched TR annotation for the caller (biggest FP lever) |
| Coverage | >=15x | higher resolves more and rarer SVs |
| Truth set (optional) | VCF | GIAB HG002 Tier 1 + CMRG for benchmarking |

## What the Agent Will Do

1. (ONT only) Basecall POD5/FAST5 with Dorado using a sup model, requesting methylation at basecall time if it will ever be needed.
2. QC read length and quality with NanoPlot and gate on read N50 and mean quality.
3. Align with minimap2 using the platform-matched preset (lr:hq / map-ont / map-hifi) and keep -Y so split reads retain breakpoint sequence.
4. Call SVs with Sniffles2 (or cuteSV/pbsv), supplying a reference-matched tandem-repeat BED.
5. Optionally call the highest-quality set from a phased diploid assembly with dipcall or PAV.
6. Build a joint-genotyped cohort with the two-step Sniffles2 .snf design rather than unioning discovery VCFs.
7. Filter and annotate, then benchmark with Truvari against Tier 1 and CMRG, reporting the full parameter set.

## Platform choice: ONT vs PacBio HiFi

| Feature | ONT (R10.4.1) | PacBio HiFi |
|---------|---------------|-------------|
| Accuracy | ~Q20+ simplex, higher duplex | ~Q30+ |
| Length | tens of kb, ultralong >Mb | ~15-25 kb |
| minimap2 preset | lr:hq (R10) or map-ont (R9) | map-hifi |
| SV caller | Sniffles2 / cuteSV | Sniffles2 / cuteSV / pbsv |
| Bonus | native methylation if requested at basecall | native phasing from read length |

## Sniffles2 vs cuteSV

| Feature | Sniffles2 | cuteSV |
|---------|-----------|--------|
| Speed | Moderate | Fast |
| Recall on noisy ONT | High | Highest |
| Multi-sample | Built-in two-step .snf | External merge |
| Defaults | sensible | NOT platform-appropriate; must tune + enable --genotype |
| Best for | general use, cohorts, mosaic | high-recall ONT |

## Tips

- 15-30x coverage is recommended; below ~10x callers drift toward false negatives.
- Longer reads detect larger SVs and resolve complex rearrangements better.
- Provide a reference-matched tandem-repeat BED to the caller - it is the single biggest false-positive lever in repeats.
- cuteSV defaults are not platform-appropriate; use the ONT, HiFi, or CLR parameter set and enable --genotype (off by default).
- Map with minimap2 >= 2.28 and keep -Y (soft-clipped supplementaries) so split reads retain breakpoint sequence; use lr:hq for accurate R10 ONT, map-ont only for older R9.
- Request methylation at basecall time - it cannot be recovered later and rides through alignment as a free per-haplotype channel.
- Build cohorts with the two-step .snf design (per-sample .snf, then combine); a union of discovery VCFs produces false 0/0 genotypes and wrong allele frequencies.
- An SV F1 is meaningless without its Truvari parameters; report refdist/pctsize/pctseq/sizemin and never disable pctseq to inflate insertion scores.
- Benchmark on GIAB HG002 Tier 1 AND CMRG - Tier 1 excludes the medically relevant repetitive genes.
- For tumor-normal somatic SVs use a paired caller (Severus/nanomonsv), not Sniffles --mosaic.

## Related Skills

- long-read-sequencing/basecalling - Dorado model choice and requesting methylation at basecall time
- long-read-sequencing/long-read-alignment - minimap2 preset selection and -Y soft-clipping
- long-read-sequencing/long-read-qc - read QC and chimera screening before alignment
- long-read-sequencing/structural-variants - caller tuning, tandem-repeat BED, Truvari (the long-read SV mechanism)
- long-read-sequencing/haplotype-phasing - haplotag the BAM for phased/somatic SVs
- variant-calling/structural-variant-calling - the SV signal model, VCF representation, force-genotyping, merging (also short-read SV)
- variant-calling/consensus-sequences - why symbolic SV alleles are not directly consensus-able
