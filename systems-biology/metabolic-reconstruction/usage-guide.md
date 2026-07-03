# Metabolic Reconstruction - Usage Guide

## Overview

Automated reconstruction turns an annotated genome into a draft genome-scale metabolic model. Two tools dominate microbial work and they are not interchangeable: CarveMe carves a draft top-down out of a curated BiGG universal model (fast, simulation-ready, biased toward well-studied organisms), while gapseq builds bottom-up from pathway and transporter evidence (slower, more transparent, better on non-model clades and fermentation phenotypes, ModelSEED namespace). The result is a hypothesis, not a finished model: both tools gap-fill specifically to force growth on a chosen medium, so a draft that grows was made to grow and proves nothing biological. Gap-filled reactions are the least-evidenced part of the model, the gap-fill medium determines what gets added, and reconstruction is the first step of curation, never a substitute for it.

## Prerequisites

```bash
pip install carveme          # CarveMe: also needs DIAMOND and an LP solver (CPLEX/Gurobi academic; SCIP fallback)
git clone https://github.com/jotech/gapseq   # gapseq: cloned, not pip; needs R + bedtools/hmmer etc.
pip install cobra            # to load and inspect the draft
```

Inputs: an annotated protein FASTA (`.faa`) for CarveMe (or nucleotide with `--dna`); a genome FASTA for gapseq. Raw/GenBank genomes are not accepted by CarveMe. For downloading inputs, see database-access/ncbi-datasets-cli and genome-annotation/prokaryotic-annotation.

## Quick Start

Tell your AI agent:
- "Build a metabolic model from this bacterial protein FASTA with CarveMe"
- "Reconstruct a model for this non-model soil bacterium with gapseq"
- "Gap-fill the draft so it grows on M9 minimal medium, and flag what was added"
- "Which tool should I use - CarveMe or gapseq - for my organism?"
- "Load my draft and tell me what to curate first"

## Example Prompts

### Tool Choice
> "I have a well-studied Gram-positive bacterium and need drafts for 50 strains quickly - recommend CarveMe or gapseq and justify it, then give the command with the right gram universe."

### Gap-Filling
> "Reconstruct with CarveMe and gap-fill for M9 glucose, then list which reactions were added by gap-filling so I can treat them as low-confidence."

### Non-Model Organism
> "Build a gapseq model for this environmental isolate, and explain why its ModelSEED IDs will complicate merging it with a CarveMe model later."

### Sanity Check
> "Load my draft SBML, confirm it grows on the gap-fill medium, and inventory orphan reactions and namespace so I know what curation it needs."

## What the Agent Will Do

1. Confirm the input format (protein FASTA for CarveMe; genome for gapseq) and pick the tool from the decision table.
2. Run reconstruction with the correct gram/universe setting and gap-fill to the intended medium.
3. Load the draft in COBRApy, report network size, gene coverage, exchanges, and orphan reactions.
4. Confirm growth on the gap-fill medium (true by construction) and flag gap-filled reactions as low-confidence.
5. Identify the namespace (BiGG vs ModelSEED) and note that cross-tool merge/comparison needs MetaNetX reconciliation.
6. Hand the draft to model-curation.

## Tips

- CarveMe takes a protein FASTA by default; use `--dna` for nucleotide. It rejects raw/GenBank genomes.
- Gram type and universe are values of `-u/--universe` (`grampos`, `gramneg`, `bacteria`, `archaea`, `cyanobacteria`), not `--grampos`/`--gramneg` flags.
- Always gap-fill to the medium you actually care about; the medium determines which reactions get added.
- Treat a draft that grows with suspicion, not satisfaction - it was gap-filled to grow. Flag gap-filled reactions low-confidence.
- gapseq's `find-transport` is its own subcommand, not `find -t`; its transporter file is `-Transporter.tbl` (singular).
- CarveMe outputs BiGG IDs, gapseq/ModelSEED outputs `seed.*` IDs; reconcile through MetaNetX/MNXref before merging or comparing models.
- Typical bacterial draft is ~1000-2500 reactions; a much smaller network usually means poor annotation or the wrong input file.
- No single tool dominates - different tools give different models from the same genome; consider consensus/ensemble reconstruction for important organisms.
- Reconstruction is the start of curation. Send the draft to model-curation before trusting any prediction.

## Related Skills

- systems-biology/model-curation - The required next step: curate, gap-fill deliberately, validate
- systems-biology/flux-balance-analysis - Predict growth/flux once the model is trustworthy
- systems-biology/community-metabolic-modeling - Combine reconstructions into a community model
- genome-annotation/prokaryotic-annotation - Produce the annotated protein FASTA these tools consume
- database-access/ncbi-datasets-cli - Fetch genome/proteome inputs
