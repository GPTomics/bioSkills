# AlphaFold Predictions - Usage Guide

## Overview

This skill retrieves precomputed models from the AlphaFold Protein Structure Database (AFDB) by UniProt accession and, more importantly, reads their confidence files correctly. An AFDB entry is a single AlphaFold2 prediction of one UniProt sequence modeled as an isolated chain - it has no ligands, ions, cofactors, PTMs, quaternary assembly, or alternative conformations, so it answers "what fold does this sequence adopt on its own?" not "what is the functional state in its complex?" The two confidence signals answer different questions: pLDDT is per-residue local confidence (stored in the B-factor column, opposite polarity to thermal motion), while PAE is per-residue-pair confidence that governs inter-domain placement. Low-pLDDT stretches are usually intrinsically disordered regions, not errors. The skill guides fetching via REST metadata (so version suffixes are discovered, not hard-coded), banding pLDDT, segmenting domains by PAE, and deciding when to trust AFDB versus run a fresh prediction.

## Prerequisites

```bash
pip install biopython numpy requests matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Download the AlphaFold model for UniProt P04637 and report mean pLDDT"
- "Which regions of this AlphaFold prediction are low confidence, and are they disordered?"
- "Read the PAE matrix and tell me if the two domains are confidently placed relative to each other"
- "Is this AlphaFold model suitable for docking, or should I run my own prediction?"

## Example Prompts

### Fetching Models
> "Retrieve the AlphaFold model for UniProt ID P04637 using the current database version"

> "This is a 3200-residue human protein - get all AlphaFold fragments for it"

> "Fetch the AlphaFold entry and its PAE file for accession Q9Y6K9"

### Reading pLDDT Correctly
> "Report the per-residue pLDDT for this model and classify each residue into confidence bands"

> "There is a long low-pLDDT stretch here - is it a disordered region or a modeling error?"

> "The mean pLDDT is 88 - can I trust the distance between the N- and C-terminal domains?"

### Reading PAE
> "Load the PAE matrix and identify the confidently predicted domains"

> "Are the relative orientations of these two domains reliable, or should I treat them as independent rigid bodies?"

> "Segment this multidomain model into blocks I can place independently for molecular replacement"

### Deciding Downstream Use
> "Is this AFDB monomer usable as a molecular-replacement search model?"

> "I need the ligand-bound state of this enzyme - can AFDB give me that?"

> "Should I use the AFDB entry or run AlphaFold3 for this protein-protein complex?"

## What the Agent Will Do

1. Query the AFDB REST metadata endpoint for the accession to discover current download URLs (the version suffix drifts, so it is never assembled by hand).
2. Download the coordinate file (mmCIF) and the PAE JSON, iterating fragments for long proteins.
3. Extract per-residue pLDDT from the B-factor column and band it (90 / 70 / 50 cutoffs), flagging low-pLDDT stretches as candidate intrinsically disordered regions rather than errors.
4. Load the compact PAE matrix and segment confident domains, judging inter-domain placement that mean pLDDT cannot certify.
5. Advise on downstream suitability (Foldseek search, molecular replacement, docking) and whether the static AFDB monomer answers the question or a fresh prediction is required.

## Database Access

The prediction metadata endpoint returns the current file URLs, so code should read them rather than hard-code a version token:

```
GET https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ACCESSION}
```

The response is a JSON list; each entry carries `cifUrl`, `pdbUrl`, `bcifUrl`, `paeDocUrl` (the PAE JSON), and `paeImageUrl`. Files follow the pattern `AF-{accession}-F{fragment}-model_v{N}.cif` where the version token is v6 as of 2025 and F1 is the only fragment for most proteins.

## Confidence Score Interpretation

| pLDDT | Band | Meaning | AFDB color |
|-------|------|---------|------------|
| 90-100 | Very high | Backbone and well-oriented side chains | Blue |
| 70-90 | Confident | Backbone generally correct | Cyan |
| 50-70 | Low | Backbone uncertain, cautionary zone | Yellow |
| < 50 | Very low | Placeholder ribbon, frequently an IDR | Orange |

## Tips

- pLDDT is in the B-factor column but is confidence, not motion - never color "by B-factor" to infer flexibility, and convert to a pseudo-B before any crystallographic use.
- A long low-pLDDT stretch is usually a genuine intrinsically disordered region; keep it and annotate it, and cross-check a sequence-based disorder predictor rather than trimming it as junk.
- High mean pLDDT does not certify inter-domain geometry - read PAE before measuring any inter-domain distance or seeding a rigid-body docking.
- PAE is asymmetric; average both off-diagonal blocks when scoring the relative placement of two domains.
- AFDB is a monomer apo model with no ligands, cofactors, PTMs, or assembly - do not read mechanism or catalytic geometry from it; get a holo experimental structure.
- Long non-human proteins over 2700 residues are often absent, and long human proteins are fragmented into F1, F2 with untrustworthy relative placement across fragments.
- For remote-homology detection, feed the confident core into Foldseek 3Di search; the low-pLDDT spaghetti degrades hits.
- Reach for a fresh prediction (modern-structure-prediction) the moment a complex, ligand, mutant, custom MSA depth, or a specific conformational state is needed.

## Related Skills

- structural-biology/modern-structure-prediction - run a new prediction for complexes, ligands, mutants, or a specific state
- structural-biology/structure-io - parse and convert the downloaded PDB/mmCIF
- structural-biology/geometric-analysis - RMSD and superposition against an experimental structure
- structural-biology/structure-modification - trim low-pLDDT regions or overload the B-factor column for coloring
- alignment/structural-alignment - Foldseek 3Di search over AFDB
- database-access/uniprot-access - resolve names and sequences to the UniProt accession AFDB is keyed on
