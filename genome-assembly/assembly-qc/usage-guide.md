# Assembly QC - Usage Guide

## Overview

Assembly QC judges whether a genome assembly is good enough to annotate, submit, or publish by measuring three orthogonal axes: contiguity (how fragmented), completeness (is everything there), and correctness (is what is there right). The central lesson is that no single number is quality - least of all N50, the most-gamed metric in genomics. The modern, reviewer-expected standard reports reference-free metrics: auN/NGx for contiguity, BUSCO or compleasm gene-space completeness plus Merqury k-mer completeness, and a Merqury QV for per-base accuracy, with reference-based QUAST reserved for the special case of a same-organism reference. The agent picks methods by whether a trusted reference exists, the read chemistry, the ploidy, and the downstream goal.

## Prerequisites

```bash
conda install -c bioconda quast busco compleasm merqury meryl minimap2 inspector merfin genomescope2
# k8 + calN50.js/paftools.js ship with minimap2; CRAQ installs from its GitHub repo
# BUSCO/compleasm lineage datasets and the GenomeScope2 model download separately
```

Large downloads to budget for: BUSCO/compleasm lineage datasets (use `--offline` + a local path on HPC), and a meryl k-mer database built from the raw reads (can be large for big genomes).

## Quick Start

Tell your AI agent what you want to do:
- "Assess the quality of my genome assembly across contiguity, completeness, and correctness"
- "Compute a reference-free QV for my assembly with Merqury"
- "Run BUSCO and compleasm on my assembly and tell me which to trust"
- "My BUSCO Duplicated is high - is that haplotigs or a real whole-genome duplication?"
- "Compute auN and NG50 instead of plain N50"

## Example Prompts

### Reference-free QC of a novel genome
> "I assembled a non-model insect from HiFi reads and have no trusted reference. Give me a reference-free QC: Merqury QV and k-mer completeness, compleasm completeness, and auN/NG50, and tell me whether it clears the EBP 6.C.Q40 bar."

> "Run Inspector on my Nanopore assembly to find structural errors without a reference, and explain why I should not use QUAST against the congener genome I have."

### Diagnosing a duplicated or fragmented assembly
> "My assembly is about 1.8x the GenomeScope2 size estimate and BUSCO Duplicated is 14%. Walk through whether this is uncollapsed haplotigs or a real WGD, and what to check before running purge_dups."

> "BUSCO Complete is 91% on my T2T-grade HiFi assembly, which seems low. Re-run completeness with compleasm and explain whether the gap is the assembly or BUSCO's predictor."

### QV and the circularity trap
> "I polished my assembly with Illumina reads. Build a Merqury k-mer database and compute QV, and make sure the QV is not measured on the same reads I polished with."

### Comparing assemblers
> "Compare my SPAdes, Flye, and hifiasm assemblies on contiguity (auN/NG50), completeness (compleasm), and accuracy (Merqury QV), and tell me which is best and why N50 alone would mislead me."

### Phased diploid QC
> "I have a trio-phased hifiasm assembly with parental reads. Compute switch and hamming error with Merqury hap-mers and interpret the hap-mer blob plot."

## What the Agent Will Do

1. Establish whether a trusted same-organism reference and a genome-size estimate exist, and choose reference-free methods by default.
2. Compute contiguity with auN/NGx (calN50/paftools) and QUAST, reporting contig and scaffold N50 separately.
3. Score gene-space completeness with BUSCO and/or compleasm at the deepest applicable clade lineage, recording the OrthoDB generation.
4. Build a meryl k-mer database from accurate reads at the best_k.sh-derived k and run Merqury for QV, k-mer completeness, and the spectra-cn copy-number plot.
5. For novel genomes, map raw reads back with Inspector or CRAQ for reference-free structural errors; for same-organism references, add QUAST structural metrics.
6. For phased assemblies, compute switch/hamming error and read the hap-mer blob plot.
7. Triangulate high BUSCO-Duplicated against size, spectra-cn, and asmgene to call haplotigs vs WGD, and recommend purge_dups or not.
8. Report all three axes with method and read-source, compare against the EBP/VGP bar appropriate to the organism, and flag any axis that is missing.

## Tips

- Distrust any assembly described by N50 alone; demand a QV, completeness, and auN/NGx.
- Never compute QV from the reads used to polish - it is circular and inflated. Build the k-mer DB from accurate, ideally independent reads (HiFi/Illumina, not noisy ONT).
- NG50, auNG, and the spectra-cn plot all need a genome-size estimate - profile the reads first (GenomeScope2).
- Use the deepest applicable BUSCO/compleasm lineage; a 99% on shallow eukaryota_odb10 (~255 genes) is not comparable to a deep clade set, and odb10 vs odb12 are not comparable across the version boundary.
- On high-quality HiFi/T2T-era genomes, prefer compleasm - BUSCO under-reports because its own predictor misses present genes.
- Read BUSCO C, F, and M together; high Fragmented behind a good Complete signals a contiguity/base-quality problem.
- QUAST is a comparison tool, not a quality tool, unless the reference is the same (or near-isogenic) organism.
- High BUSCO-Duplicated with size above the estimate is uncollapsed haplotigs until proven WGD; but watch the over-purge tightrope - purging below the estimate deletes real segmental duplications.
- Always report contig AND scaffold N50; scaffold >> contig means the contiguity is gap-Ns and every join is a misassembly hypothesis.

## Related Skills

- short-read-assembly - Short-read assemblies plateau at the repeat structure; QC shows it
- long-read-assembly - Produces the contiguous-but-error-prone contigs this QC evaluates
- hifi-assembly - Phased diploid output whose false duplication and switch/hamming error this QC checks
- assembly-polishing - Merqury QV plateau is the honest stop signal; never QV on the polishing reads
- scaffolding - Validate chromosome-scale joins before trusting scaffold NG50
- contamination-detection - MAG completeness/contamination (CheckM2/GUNC/MIMAG) is a separate problem
- genome-profiling - GenomeScope2 genome-size estimate that NG50/auNG/spectra-cn require
- genome-annotation/annotation-qc - Assembly-side completeness; purge haplotigs before annotating
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> QC
