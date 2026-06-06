# Long-Read Alignment - Usage Guide

## Overview
Long-read alignment maps Oxford Nanopore and PacBio reads (or whole assemblies) to a reference with minimap2. The central decision is the `-x` preset, which is not a label but a bundle that rewrites k-mer/window and the entire scoring and chaining model - so the right preset depends on the reads' error rate (noisy R9 vs accurate Q20/R10/HiFi), not just the platform. This skill also covers the tags downstream callers silently depend on (`-Y` for SV split reads, `-y` for methylation MM/ML, `--MD`/`--cs` for variant tools), the multi-part-index MAPQ trap on large references, and when to swap minimap2 for Winnowmap (repeats/centromeres), VACmap/lra (complex SVs), or pbmm2 (PacBio-native).

## Prerequisites
```bash
conda install -c bioconda minimap2 samtools
# Optional specialized aligners:
conda install -c bioconda winnowmap pbmm2 lra
# VACmap: see https://github.com/micahvista/VACmap
```

## Quick Start
Tell your AI agent what you want to do:
- "Align my accurate R10 Nanopore reads with the lr:hq preset and give me a sorted BAM"
- "Map PacBio HiFi reads and keep the tags Sniffles needs"
- "Align my Dorado methylation BAM without losing the MM/ML tags"
- "Spliced-align my direct RNA reads"

## Example Prompts

### Error-rate-matched preset
> "I have ONT R10.4.1 reads basecalled with Dorado sup. Map them to GRCh38 with the correct preset, keep supplementary alignments soft-clipped for SV calling, add the MD tag, and give me a sorted+indexed BAM ready for Clair3 and Sniffles."

### Methylation passthrough
> "These are Dorado-basecalled reads with 5mC MM/ML tags. Align them to the genome without dropping the methylation tags so I can run modkit pileup."

### Repeats and centromeres
> "My reads are mismapping in centromeric and segmental-duplication regions of a T2T reference. Align them with something that handles repetitive sequence correctly."

### Spliced cDNA / Iso-Seq
> "These are PacBio Iso-Seq full-length cDNA reads. Spliced-align them to the genome so I can feed the BAM to isoform collapse."

### Assembly to reference
> "Align my assembly to the reference at ~1% divergence and call the variants from the alignment."

## What the Agent Will Do
1. Determine read type and error rate (R9 noisy vs R10/Q20/duplex accurate vs HiFi vs CLR vs RNA).
2. Choose the matching preset (map-ont / lr:hq / map-hifi / map-pb / splice / asm*).
3. Add the tags the downstream tool needs (`-Y`, `-y`, `--MD`, `--cs`) and a read group.
4. Pipe to `samtools sort` and index, or use pbmm2 for a one-call sorted+indexed PacBio BAM.
5. For repeat-heavy references or complex SVs, route to Winnowmap2 or VACmap/lra.
6. For very large references, set `-I`/`--split-prefix` to avoid the multi-part-index MAPQ problem.

## Tips
- Accurate ONT (Q20+/duplex/R10 sup) belongs on `lr:hq`, not `map-ont` - it is ~4x faster and at least as accurate.
- `map-pb` is PacBio CLR only; never use it for HiFi (use `map-hifi`).
- For SV calling, keep `-Y` (soft-clip supplementaries) so breakpoint sequence survives; `--secondary=no` is fine.
- Methylation needs both `samtools fastq -T MM,ML` and `minimap2 -y`, or the tags vanish with no error.
- A prebuilt `.mmi` bakes in k/w/H/I; build it with the same preset you will align with.
- Direct RNA needs `-uf` (stranded); without it minimap2 invents junctions on the wrong strand.
- Use minimap2 >= 2.28: it has `lr:hq`/`lr:hqae` and fixes the 2.27 `--MD` regression.

## Related Skills

- basecalling - Basecaller chemistry/error rate that picks the preset
- long-read-qc - Read QC before alignment; percent identity from the aligned BAM
- structural-variants - SV calling from the supplementary split-read signal
- clair3-variants - Small-variant calling on the aligned BAM
- nanopore-methylation - Pileup of MM/ML tags carried through with `-y`
- isoseq-analysis - Spliced alignment of full-length transcripts
- alignment-files/sam-bam-basics - Sort/index/inspect the BAM
- genome-assembly/long-read-assembly - De novo assembly as an alternative to mapping

## Resources
- [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html)
- [minimap2 cookbook](https://github.com/lh3/minimap2/blob/master/cookbook.md)
- [Winnowmap](https://github.com/marbl/Winnowmap)
- [pbmm2](https://github.com/PacificBiosciences/pbmm2)
