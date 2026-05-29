# Annotation QC - Usage Guide

## Overview

Assess whether a genome annotation is good enough to publish or submit, and diagnose what is wrong when it is not, using BUSCO (conserved single-copy ortholog recovery), OMArk (proteome completeness, consistency, contamination), CheckM2 (prokaryotic completeness/contamination), and a gene-set sanity panel (gene count, mono-exonic fraction, protein-length distribution, mRNA:gene ratio, coding density). The orienting principle: BUSCO measures only the easy conserved core, so the diagnostic that matters is assembly-BUSCO vs proteome-BUSCO on the same assembly, gene count is a vanity metric, and a green BUSCO can sit on top of a proteome full of chimeric, fragmented, or contaminant models.

## Prerequisites

```bash
# Completeness and consistency
conda install -c bioconda busco omark compleasm

# Prokaryotic completeness/contamination
conda install -c bioconda checkm2

# Python sanity panel
pip install gffutils matplotlib pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Is my annotation good enough to publish - run BUSCO on the proteome and the genome and compare"
- "My vertebrate annotation has 45,000 genes - is that haplotigs, unmasked repeats, or split models?"
- "Run CheckM2 on my MAG before I trust the annotation numbers"
- "Check whether my proteome has contaminant or chimeric models with OMArk"

## Example Prompts

### Completeness

> "Run BUSCO in proteins mode on my delivered proteome and in genome mode on the assembly with the deep clade dataset, and tell me which one is the limiting factor."

> "My BUSCO is 98% but the annotation looks rough - run OMArk and the gene-set sanity panel to find the over-prediction."

### Diagnosing a Suspect Annotation

> "BUSCO-Duplicated is 12% - is this real WGD or uncollapsed haplotigs? Check synteny and ploidy."

> "Compute the mono-exonic fraction and protein-length distribution - is this annotation full of TE ORFs?"

### Prokaryotic and Transferred

> "Gate my bacterial annotation on CheckM2 completeness and contamination, then check coding density."

> "Run BUSCO on my lifted annotation vs the reference to see how many conserved genes the transfer lost."

## What the Agent Will Do

1. For eukaryotes, run BUSCO `-m proteins` AND `-m genome` on the same assembly and compare (the diagnostic fork)
2. Use the deepest applicable clade dataset and record it
3. Run OMArk for proteome consistency and contamination; CheckM2 for prokaryotic completeness/contamination
4. Compute the gene-set sanity panel (gene count vs relative, mono-exonic fraction, protein length, mRNA:gene)
5. Interpret BUSCO-Duplicated against synteny and clade ploidy (haplotigs vs WGD)
6. Conclude whether the limit is the assembly or the predictor, and route accordingly

## Tips

- **Genome-vs-proteome BUSCO is the diagnostic** - large gap (assembly high, proteome low) = predictor missed present genes (fix evidence/masking); both low = fix the assembly. Almost nobody runs this fork.
- **Report `-m proteins` on the delivered proteome** - genome-mode BUSCO does its own gene-finding and is not your annotation's score.
- **Gene count is a vanity metric** - read it against the nearest relative and ploidy; the signal is in mono-exonic fraction, protein length, mean exons/gene, and mRNA:gene ratio.
- **High Duplicated has three causes** - haplotigs (purge_dups before annotating), real WGD (keep, confirm via synteny/Ks), or split models. >5-8% with no known WGD means purge first.
- **Use the deep clade dataset** - `eukaryota_odb10` (~255 genes) is trivially easy; a 99% there is not a 99% on a 5,500-gene clade set.
- **BUSCO alone is not enough** - add OMArk (consistency + contamination); BUSCO's single-copy lens is blind to over-prediction and chimeras.
- **Prokaryotes: CheckM2 first** - contamination >5% or completeness <90% makes the annotation numbers uninterpretable.

## Related Skills

- genome-annotation/eukaryotic-gene-prediction - The annotation this QC evaluates
- genome-annotation/prokaryotic-annotation - CheckM2 gate and coding-density sanity
- genome-annotation/annotation-transfer - BUSCO on a lifted set quantifies lost conserved genes
- genome-assembly/assembly-qc - Assembly-side completeness; purge haplotigs before annotating
