# Structure Navigation - Usage Guide

## Overview

This skill covers walking the Bio.PDB SMCRA hierarchy (Structure-Model-Chain-Residue-Atom) while surfacing the heterogeneity it hides by default. The tree is convenient but lossy: a DisorderedAtom silently forwards to one conformer, a residue key is a 3-tuple whose insertion code and hetero flag are routinely dropped, the sequence PPBuilder returns is the observed (ATOM) sequence with missing-density gaps concatenated away rather than the SEQRES or UniProt canonical sequence, and NMR files carry many models that naive iteration conflates. The skill turns each of those defaults into an explicit decision so the coordinates, sequence, and residue selection are statements about the molecule rather than artifacts of the parser.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "List all chains and how many residues each has"
- "Extract the observed sequence from chain A and tell me where the density gaps are"
- "Show every alternate conformation at residue 42 before I measure a distance"
- "Find all ligands but keep the catalytic metal"

## Example Prompts

### Accessing Structure Parts
> "Get residue 100A from chain H using the full residue id tuple"

> "What are the coordinates of every altloc of the CB atom in residue 50?"

> "Read this mmCIF with auth numbering so the residue numbers match the paper"

### Iterating and Counting
> "Count amino acids, ligands, and waters per chain"

> "List all C-alpha coordinates for model 0 only"

### Sequences
> "Extract the observed sequence and compare its length to SEQRES"

> "Build the one-letter sequence but map selenomethionine to M"

### Handling Heterogeneity
> "Enumerate all disordered atoms in chain A and their occupancies"

> "This is an NMR ensemble - compute per-model CA spread instead of using model 1"

## What the Agent Will Do

1. Parse the structure and establish how many models exist and why (NMR ensemble vs single asymmetric unit).
2. Navigate to the requested level using the full residue id tuple, keeping insertion codes and hetero flags intact.
3. Detect disordered atoms and residues and enumerate conformers rather than silently using the highest-occupancy child.
4. Extract sequence from the correct source (observed ATOM, declared SEQRES, or canonical UniProt) for the question asked, and flag missing-density gaps.
5. Filter waters, ligands, and metals by hetflag with an explicit policy, never a blanket hetero strip.
6. Return organized results with the numbering scheme (auth vs label) stated.

## SMCRA Hierarchy

```
Structure (whole PDB/mmCIF entry)
    Model (NMR conformer, or the single asymmetric unit)
        Chain (polypeptide, nucleic acid, or hetero groups)
            Residue (amino acid, nucleotide, ligand, water)
                Atom (individual atom; may be a DisorderedAtom)
```

## Tips

- A DisorderedAtom forwards uncaught calls to its HIGHEST-OCCUPANCY child, not literally altloc 'A' - enumerate with `disordered_get_list()` before any distance or RMSD.
- The residue id is `(hetflag, resseq, icode)`; index with the full tuple and never drop the insertion code (antibody 100A/100B) or the hetero flag.
- PPBuilder returns the observed sequence with density gaps concatenated away - use SEQRES (`pdb-seqres`) for true length and map to UniProt by residue number or SIFTS, not by string position.
- Water hetflag is `'W'`, not `'H_'`; a `startswith('H_')` water strip misses water and a blanket `!= ' '` strip deletes catalytic metals and mid-chain selenomethionine.
- CaPPBuilder bridges CA atoms within ~4.3 Angstroms so it handles broken backbones, but it can mis-join genuine gaps.
- NMR files have many models - select a representative or report per-model spread; never average coordinates across the ensemble.
- MMCIFParser defaults to auth numbering (matches the paper); label numbering is a gapless 1..N scheme where the same number is a different residue.

## Related Skills

- structure-io - Parse and write PDB/mmCIF; auth vs label numbering at the I/O layer
- geometric-analysis - Distances, angles, RMSD, SASA once heterogeneity is resolved
- structure-modification - Strip waters/hetero, edit coordinates and B-factors safely
- interface-analysis - Requires the biological assembly, not the deposited asymmetric unit
- sequence-manipulation/seq-objects - Work with the extracted Seq objects
- alignment/msa-parsing - Map SEQRES/ATOM sequences onto alignment columns
- database-access/uniprot-access - Fetch the canonical sequence for SIFTS-based mapping
