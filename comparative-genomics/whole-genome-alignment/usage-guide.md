# Whole Genome Alignment - Usage Guide

## Overview

Whole-genome alignment (WGA) is the substrate for comparative genomics at the base-pair level. Modern production WGA uses Progressive Cactus (Armstrong 2020 Nature 587:246) for reference-free multi-genome alignment and Minigraph-Cactus (Hickey 2024 Nat Biotech 42:663) for pangenome graph construction. LASTZ chain/net (UCSC pipeline) and minimap2 remain standards for pairwise reference-anchored alignment. The HAL (Hierarchical Alignment Format; Hickey 2013) is the standard multi-genome alignment container, supporting downstream tools like TOGA, halSynteny, and halLiftover.

The fundamental decision is reference-free (Cactus, all genomes equally weighted) vs reference-anchored (LASTZ chains/nets, one privileged reference). For multi-species comparative genomics, Cactus is the modern standard; for pairwise pipelines like UCSC genome browser tracks, LASTZ chain/net remains canonical.

## Prerequisites

```bash
# Progressive Cactus + Minigraph-Cactus
docker pull quay.io/comparative-genomics-toolkit/cactus:latest
# Or Python venv: pip install cactus

# HAL toolkit
git clone https://github.com/ComparativeGenomicsToolkit/hal

# LASTZ + UCSC kentUtils
conda install -c bioconda lastz ucsc-axt-chain ucsc-chain-net

# Pairwise tools
conda install -c bioconda minimap2 mummer4 anchorwave winnowmap

# Repeat masking
conda install -c bioconda repeatmasker repeatmodeler

# Toil for HPC / cloud
pip install toil[all]
```

## Quick Start

Tell the AI agent what to do:
- "Build a Progressive Cactus alignment of these 20 mammal genomes with the provided guide tree"
- "Run pairwise LASTZ + chain/net for query genome against the human reference for UCSC browser tracks"
- "Use Minigraph-Cactus to build a pangenome graph from these 10 haplotype-resolved assemblies"
- "Use minimap2 -x asm5 to align two closely related strains then run SyRI for structural variant calls"

## Example Prompts

### Multi-Genome Reference-Free WGA

> "Build a Progressive Cactus alignment of these 30 vertebrate genomes (FASTA files + guide tree provided). Pre-mask repeats with RepeatModeler2 + RepeatMasker softmask. Use `--branchScale 1.0` for vertebrates. Produce a HAL file and extract syntenic blocks via halSynteny for downstream TOGA annotation projection."

### Pangenome Graph Construction

> "Build a pangenome graph from 90 haplotype-resolved human assemblies using Minigraph-Cactus with GRCh38 as reference. Produce GFA, VCF (relative to GRCh38), and GBZ outputs. Verify the graph contains all input haplotypes and report the number of short variants and structural variants."

### Pairwise Synteny Comparison

> "Align two plant genomes (Arabidopsis thaliana and Arabidopsis lyrata) using minimap2 `-x asm5` followed by SyRI for structural variant calling. Plot the result with plotsr. Identify inversions and translocations between the two genomes."

## What the Agent Will Do

1. **Validate inputs**: BUSCO/Compleasm completeness, N50, FASTA validity, guide-tree compatibility
2. **Pre-mask repeats** with RepeatModeler2 species-specific library + RepeatMasker softmask
3. **Configure Cactus seqFile** with guide tree + per-species FASTA paths
4. **Run alignment** with appropriate tool:
   - Cactus for multi-genome reference-free
   - Minigraph-Cactus for pangenome graph
   - LASTZ chain/net for pairwise reference-anchored
   - minimap2 for fast pairwise close species
5. **Quality-check output**: halStats for HAL; chain count for LASTZ; alignment fraction
6. **Extract downstream data**: synteny blocks, MAF for phyloP, VCF from pangenome
7. **Document**: tool versions, branch scale, masking strategy, reference choice
8. **Caveats**: reference bias, assembly fragmentation, repeat-driven artifacts

## Tips

- Progressive Cactus is the modern standard for multi-species WGA; LASTZ chain/net for pairwise UCSC tracks
- Always provide a guide tree for Cactus; without it scaling is O(N^2) instead of O(N)
- Pre-mask repeats softmasked (lowercase) before WGA; hardmask (N's) breaks LASTZ scoring
- Branch scale: 1.0 vertebrates; 0.5-0.7 plants; 0.3-0.5 bacteria
- minimap2 sensitivity loss below ~70% identity; switch to LASTZ for distant comparisons
- Toil for Cactus on cluster; use `--batchSystem slurm` for HPC
- For pangenome graphs (intra-species or sister-species), Minigraph-Cactus is HPRC standard
- PGGB is reference-free; use for transparent pangenome comparison; max ~20 large genomes per run
- Pin HAL toolkit version; HAL file format is versioned
- Cactus jobStore checkpointing requires shared storage accessible to all nodes
- Verify guide-tree leaf names match seqFile genome names exactly
- For pangenome graph variant calling, use vg + Giraffe; for SV genotyping, PanGenie
- Centromeric / telomeric regions usually unalignable; use Winnowmap2 for explicit T2T work
- Reference-guided assemblies create circular reasoning if used for self-comparison
- Minimum N50 for reliable WGA: 1 Mb; chromosome-level preferred

## Related Skills

comparative-genomics/synteny-analysis - Synteny detection from WGA
comparative-genomics/comparative-annotation-projection - TOGA / CESAR uses Cactus HAL
comparative-genomics/pangenome-analysis - Minigraph-Cactus / PGGB / PGR-TK for pangenome graphs
alignment/multiple-alignment - Protein MSA from WGA orthologs
alignment/pairwise-alignment - LASTZ / minimap2 sequence-level alignment
alignment/structural-alignment - WGA-extracted protein orthologs
genome-assembly/assembly-qc - BUSCO checks before WGA
variant-calling/structural-variant-calling - WGA-based SV detection
