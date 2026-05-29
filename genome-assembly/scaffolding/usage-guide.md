# Scaffolding - Usage Guide

## Overview

Scaffolding orders and orients assembled contigs into chromosome-scale scaffolds using long-range linking data, inserting runs of N as gap spacers between joins. It adds no new sequence - so a scaffold N50 can be far larger than the underlying contig N50, and the two must be reported separately. Hi-C / Omni-C proximity ligation is the dominant modern modality (YaHS, SALSA2, 3D-DNA); genetic/linkage maps (ALLMAPS), Bionano optical maps, and reference homology (RagTag) are alternatives keyed to the data available. The single most important point: automated scaffolding produces a *draft* chromosome structure. It becomes trustworthy only after manual curation of the contact map (the VGP/Darwin Tree of Life standard), where misjoins, inversions, and false duplications are fixed by eye. Reference-guided scaffolding additionally imposes the reference's karyotype, erasing real fusions/translocations - so it must never be used on a genome whose chromosome structure is a biological question.

## Prerequisites

```bash
# Hi-C scaffolders
conda install -c bioconda yahs salsa2

# Hi-C read mapping + dedup
conda install -c bioconda bwa chromap samtools

# Reference-guided + linkage maps
conda install -c bioconda ragtag jcvi   # jcvi provides ALLMAPS

# QC: telomeres, stats
conda install -c bioconda tidk seqkit assembly-stats

# 3D-DNA / Juicer (clone; juicer_tools is a jar)
git clone https://github.com/aidenlab/3d-dna.git
wget https://github.com/aidenlab/juicer/releases/download/v2.0/juicer_tools.2.0.jar

# Curation viewers (no conda; from GitHub releases)
# PretextMap / PretextView / PretextSnapshot (wtsi-hpag), Juicebox Assembly Tools (aidenlab)
```

## Quick Start

Tell your AI agent what you want to do:
- "Scaffold my draft contigs to chromosome level with Hi-C using YaHS"
- "Map my Hi-C reads correctly before scaffolding (each end separately, dedup)"
- "Prepare a .hic file and AGP so I can curate the scaffolds in Juicebox"
- "Scaffold against a closely related reference, but only for gene coordinates"
- "Order my scaffolds with a genetic linkage map using ALLMAPS"
- "Check telomeres and report contig N50 vs scaffold N50 separately"

## Example Prompts

### Hi-C / Omni-C scaffolding
> "I have HiFi contigs and an Arima Hi-C library. Map the Hi-C reads with each end aligned separately and PCR duplicates removed, run YaHS with the correct enzyme, and give me the AGP plus a .hic file for Juicebox curation."
> "My Hi-C is Omni-C (enzyme-free). Scaffold with SALSA2 using -e DNASE and -m yes so chimeric contigs get broken first, and tell me why -m yes matters."

### Reading and curating the contact map
> "Render the contact map with PretextMap and walk me through what a misjoin, an inversion bowtie, and a leaked-haplotig stripe look like so I can decide what to break in PretextView."
> "The automated scaffolds_final looks chromosome-scale but a downstream synteny analysis disagrees - explain whether I am shipping an uncurated draft and what curation would change."

### Reference-guided and linkage-map scaffolding
> "Scaffold my non-model genome against a related reference with RagTag - but first tell me whether this will erase any real chromosome fusions, because I care about karyotype evolution."
> "I have two genetic maps. Integrate them with ALLMAPS to order and orient my scaffolds and use it to validate my Hi-C chromosome assignment."

### QC and finishing
> "Report contig N50, scaffold N50, gap count, and total N bases separately, then run tidk to find which scaffolds have telomeres at both ends."
> "Gap-fill the scaffolds with my ONT reads, but only after curation, and flag any filled gap whose length is far from the estimate."

## What the Agent Will Do

1. Establish whether the job is haploid scaffolding or diploid phasing-then-scaffolding, and which linking modality is available.
2. For Hi-C: map each end separately with no mate rescue (`bwa mem -5SP` or the per-end pipeline, or chromap), MAPQ-filter, and mark/remove duplicates.
3. Break chimeric contigs before scaffolding (SALSA2 `-m yes` / YaHS contig error-correction on by default).
4. Run the scaffolder appropriate to the data (YaHS default; SALSA2 for graph-aware conservatism; 3D-DNA for interactive curation; RagTag/ALLMAPS/Bionano for their modalities).
5. Produce the AGP + FASTA and a contact map (.hic / Pretext), and prepare files for manual curation in Juicebox/PretextView.
6. Treat the automated output as a draft - inspect the contact-map diagonal for misjoins/inversions/false duplications and recommend curation before large-scale-order claims.
7. QC: report contig N50 vs scaffold N50 separately, gap/N counts, telomere recovery (tidk), and optionally gap-fill after curation.

## Tips

- The output file named `_scaffolds_final` is a first pass, not a finished genome - curate the contact map before trusting chromosome structure or claiming fusions/synteny.
- Hi-C reads are not a normal paired-end library; mapping them with proper-pair/insert-size logic corrupts the signal. Map each end separately, no mate rescue, then dedup.
- Set the enzyme correctly: `-e GATC` for DpnII/Arima, omit `-e` (YaHS) or `-e DNASE` (SALSA2) for Omni-C/enzyme-free. Wrong enzyme modeling degrades joins.
- Never reference-scaffold a genome whose karyotype matters - RagTag imposes the reference's structure and silently erases real rearrangements.
- A chimeric contig cannot be fixed by reordering; keep contig error-correction on so it is broken before scaffolding.
- Report contig N50 and scaffold N50 separately; a large ratio means contiguity is borne by N-gaps, not finished sequence.
- Gap-filling is a separate step after curation, and it can mis-fill at repeats - sanity-check filled lengths and report which gaps were closed.
- Hi-C coverage matters: roughly 50-100M valid pairs for a vertebrate-size genome; too little leaves a sparse, noisy map with weak joins.

## Related Skills

- long-read-assembly - Produces the contigs this skill orders into chromosomes
- hifi-assembly - hifiasm phases haplotypes with Hi-C before each is scaffolded
- assembly-polishing - Polish contigs before scaffolding
- assembly-qc - Contig-vs-scaffold N50, BUSCO, and Merqury QV on the scaffolded result
- hi-c-analysis/hic-data-io - Hi-C read/pairs handling feeding the scaffolder
- comparative-genomics/synteny-analysis - Validate scaffold order/orientation against a related genome
- workflows/genome-assembly-pipeline - End-to-end QC -> assemble -> polish -> scaffold -> curate -> QC
