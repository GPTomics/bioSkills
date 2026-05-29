# Annotation Transfer - Usage Guide

## Overview

Transfer gene annotations between genome assemblies by matching the method to the divergence: coordinate liftover (UCSC liftOver, CrossMap) for same-species version updates and intervals, Liftoff for same/close-species gene models, miniprot for protein-level cross-species transfer, and TOGA/GeMoMa/CAT for distant clades. The orienting principle: a lift is geometry, not biology - a successful coordinate placement says nothing about whether the gene is intact, transfer can only reproduce what the reference annotated (reference bias), and the unmapped file is the most information-rich output. Validate every transfer against the target before trusting it.

## Prerequisites

```bash
# Feature/sequence projection
pip install liftoff lifton
conda install -c bioconda miniprot

# Coordinate liftover (same-species intervals / build updates)
conda install -c bioconda ucsc-liftover crossmap

# Python utilities
pip install gffutils biopython pandas
```

## Quick Start

Tell your AI agent what you want to do:
- "Transfer the reference annotation to my new same-species assembly and check ORF integrity"
- "Lift over my GRCh38 variant coordinates to T2T-CHM13"
- "Map proteins from a related species to my genome (cross-species)"
- "Did the transfer silently lose any gene families - read the unmapped file"

## Example Prompts

### Same-Species Transfer

> "Use Liftoff to transfer the GENCODE annotation to my new human assembly, then run the ORF-integrity check and classify the unmapped file."

> "Lift over my GRCh37 ClinVar coordinates to GRCh38 with the correct chain - watch the inverted regions and indels."

### Cross-Species Transfer

> "Align proteins from a related species to my genome with miniprot and flag any frameshift/stop tags as pseudogenized."

> "These species are tens of millions of years apart - should I use TOGA instead of coordinate liftover?"

### Quality Assessment

> "Compare the transferred annotation to the reference - transfer rate, identity distribution, and BUSCO on the lifted set."

> "How many transferred gene models have intact open reading frames vs internal stops?"

## What the Agent Will Do

1. Choose the paradigm by divergence (liftOver/Liftoff same-species; +miniprot for genus; TOGA/GeMoMa for family; de novo beyond)
2. Confirm the chain matches the exact build/patch pair and record source->target build
3. Run the transfer (target before reference for Liftoff)
4. Read and **classify** the unmapped file (deletion/split/duplicated), not just count it
5. Validate intactness: `-polish` + ORF check, identity distribution, BUSCO on the lifted set
6. Compare to a de novo annotation to expose reference bias; recommend de novo for regions that failed

## Tips

- **A successful lift is not biological confirmation** - coordinate placement can land in a pseudogene; check the ORF, start/stop, and splice sites before trusting a lifted gene.
- **Read the unmapped file** - liftOver exits 0 and prints a clean shorter GFF; entire gene families can vanish with no warning. The unmapped categories are distinct diagnoses.
- **Transfer is reference-biased** - it cannot discover target-specific genes; always pair with de novo + RNA-seq evidence (the target's novelty is exactly what transfer misses).
- **Match paradigm to divergence** - same species -> liftOver (intervals) / Liftoff (models); genus -> +miniprot; family/order -> TOGA/GeMoMa; beyond -> de novo. Cross-species coordinate liftover is the wrong paradigm, not a tuning problem.
- **Liftoff is target-then-reference** - swapping the positional args is a silent error.
- **Don't loosen `-a`/`-s` to rescue features** - a 35%-coverage hit is usually a paralog/pseudogene; default-threshold failure is often real signal.
- **Build identity is non-negotiable** - record source/target build; hg19 chrM is not rCRS; the chain must match the exact assembly pair. Same-species mapping <99% signals the wrong chain.

## Related Skills

- genome-annotation/eukaryotic-gene-prediction - De novo prediction; captures target-specific genes
- genome-annotation/annotation-qc - BUSCO on the lifted set and structure-integrity checks
- comparative-genomics/ortholog-inference - Orthology relationships across species
- comparative-genomics/synteny-analysis - Synteny context for transfer
- genome-intervals/gtf-gff-handling - Parse transferred annotation files
