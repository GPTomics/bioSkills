# Read-Backed Haplotype Phasing - Usage Guide

## Overview
This skill performs read-backed (physical) haplotype phasing of a single sample's Oxford Nanopore or PacBio long reads and haplotags the BAM for allele-resolved downstream analysis. Because a long read physically spans multiple heterozygous sites, haplotypes are reconstructed directly from the reads with no reference panel - so phasing works on rare and private variants, but phase blocks are capped by read length and heterozygosity and break wherever no read spans two adjacent hets. The central lesson is that phasing the VCF (GT pipe + PS) is not enough: every read-level downstream tool (allele-specific methylation, phased SVs, IGV coloring, read splitting) needs the per-read HP tag that only the separate haplotag step writes. This skill covers WhatsHap, LongPhase, and HiPhase, trio phasing as the gold standard, the phasing-quality metrics, and the diploid-assumption traps. Statistical/panel phasing for imputation is a separate skill (phasing-imputation/haplotype-phasing).

## Prerequisites
```bash
conda install -c bioconda whatshap longphase samtools htslib
# HiPhase (PacBio HiFi) and HapCUT2 (multi-tech) optionally
# Inputs: a het VCF (Clair3/DeepVariant) + the aligned long-read BAM + reference
```

## Quick Start
Tell your AI agent what you want to do:
- "Phase my Clair3 variants with WhatsHap and haplotag the BAM"
- "Phase my whole ONT genome fast and co-phase the SVs with LongPhase"
- "Trio-phase my family"
- "Report my phasing quality (block N50 and switch error)"

## Example Prompts

### Phase and haplotag
> "I have a Clair3 het VCF and my ONT BAM. Phase the variants with WhatsHap using realignment for indels, then haplotag the BAM and confirm the HP tags are present so I can run allele-specific methylation."

### Whole-genome speed with SV co-phasing
> "Phase my 30x ONT genome quickly with LongPhase, co-phasing the Sniffles SVs so I get long phase blocks, then haplotag the BAM."

### Trio
> "I sequenced a mother-father-child trio. Trio-phase the child with WhatsHap using the pedigree, since that is the gold standard."

### Quality assessment
> "Assess my phasing: report block N50, phased fraction, and the switch error rate against the GIAB trio-phased truth, and explain the flip decomposition."

### Allele-specific methylation handoff
> "Haplotag my modBAM and give me per-haplotype methylation at imprinted loci."

## What the Agent Will Do
1. Take the het VCF and the aligned BAM and choose a phaser (WhatsHap default; LongPhase for speed/SV co-phasing; HiPhase for HiFi; `--ped` for trios).
2. Phase the variants (`--reference --indels` on long reads), writing GT pipe + PS into the VCF.
3. Haplotag the BAM (HP/PS read tags) and verify the tags are present.
4. Report block N50 together with switch error (not N50 alone).
5. Hand the haplotagged BAM to downstream consumers (modkit ASM, Severus phased SVs, IGV).
6. Flag haploid/CNV/segdup regions where phasing is unreliable.

## Tips
- `phase` writes the VCF; `haplotag` writes the BAM - run both, and verify HP tags with `samtools view ... | grep HP:i:`.
- Always pass `--reference` on long reads (realignment mode) and `--indels` to phase indels.
- The trio flag is `--ped` (a PED file), not `--trio`; trio phasing is the gold standard when parents are sequenced.
- Report block N50 AND switch error together; N50 alone is gameable by over-joining blocks.
- Short blocks on a homozygosity-rich sample are biology, not tool failure - co-phase SVs (LongPhase) or use ultra-long reads to bridge sparse-het gaps.
- LongPhase uses bare `--ont`/`--pb` and subcommands `phase`/`haplotag`/`modcall`; `--max-coverage 15` in WhatsHap is a runtime cap, not a minimum depth.
- Phasing of chrX/Y/MT, CNV regions, and segdups is unreliable (the diploid assumption is false there).

## Related Skills

- clair3-variants - Produces the het VCF to phase
- long-read-alignment - Produces the BAM
- nanopore-methylation - Allele-specific methylation via `--partition-tag HP`
- structural-variants - Severus consumes a haplotagged BAM
- phasing-imputation/haplotype-phasing - Statistical/panel phasing for imputation
- genome-assembly/hifi-assembly - Phased de novo haplotype contigs

## Resources
- [WhatsHap docs](https://whatshap.readthedocs.io/)
- [LongPhase](https://github.com/twolinin/longphase)
- [HiPhase](https://github.com/PacificBiosciences/HiPhase)
