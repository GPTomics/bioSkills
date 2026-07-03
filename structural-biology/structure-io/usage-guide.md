# Structure I/O - Usage Guide

## Overview

This skill reads, writes, downloads, and converts macromolecular structure files with Biopython Bio.PDB, treating format choice and numbering scheme as decisions rather than defaults. mmCIF (PDBx) is the canonical modern format; the legacy fixed-column PDB format is frozen and cannot represent large assemblies (its hard limits are ~62 chains and 99,999 atoms), so big complexes exist only as mmCIF. Bio.PDB parses permissively and silently drops metadata, anisotropic B-factors, and alternate conformers, so a successful parse is not proof of data integrity. The skill covers auth_* vs label_* numbering, reading metadata through MMCIF2Dict, downloading from RCSB, and obtaining the biological assembly separately from the deposited asymmetric unit.

## Prerequisites

```bash
pip install biopython
# Optional, for assembly generation / large-structure fidelity:
pip install gemmi
```

## Quick Start

Tell the AI agent what the task is:
- "Download structure 4HHB from RCSB as mmCIF"
- "Parse this mmCIF and list the chains and residue count"
- "Why do the residue numbers not match the paper?"
- "Get the biological assembly for this entry, not the asymmetric unit"
- "Extract the resolution and R-free from this mmCIF"

## Example Prompts

### Reading Structures
> "Parse 1crn.pdb and report how many models, chains, and atoms it has"

> "Load this mmCIF with author numbering so the residue ids match the publication"

> "This structure downloaded as a .bcif file -- read it into Bio.PDB"

### Downloading Structures
> "Download PDB 4HHB as mmCIF into the current directory"

> "Fetch the biological assembly for 1ABC, since the deposited file is only the asymmetric unit"

> "Download these ten entries as mmCIF and tell me which failed"

### Converting and Writing
> "Save only chain A, protein residues, to a new PDB file"

> "Convert this large mmCIF to PDB and warn me if chains or serials would overflow"

> "Write this structure back out as mmCIF preserving the multi-character chain ids"

### Extracting Metadata
> "What is the resolution and experimental method of this structure?"

> "Bio.PDB says the resolution is None -- get it from the mmCIF dictionary instead"

> "Report R-work and R-free for this entry"

## What the Agent Will Do

1. Identify the format and pick the matching parser (MMCIFParser as the default, PDBParser for legacy fixed-column files, BinaryCIFParser for `.bcif`).
2. Choose the numbering scheme deliberately, defaulting to auth numbering so residue ids match the literature.
3. Download from RCSB with PDBList (noting the two-char divided subdirectory tree) or fetch a specific assembly file directly from files.rcsb.org.
4. Reach for MMCIF2Dict when metadata the object model omits (resolution, method, R-free, assembly operators) is requested.
5. Establish the biological assembly before any interface or oligomeric-state question, since the deposited coordinates are usually the asymmetric unit.
6. Write output with MMCIFIO for fidelity, or PDBIO with a Select subclass for filtered subsets, warning when a large mmCIF cannot fit the PDB format.

## Tips

- Default to mmCIF; only produce legacy PDB for a tool that demands it and only when the structure fits the format's limits.
- Converting a large mmCIF down to PDB silently renames multi-character chains and overflows atom serials -- stay in mmCIF or use gemmi's hybrid-36 writer.
- `MMCIFParser` defaults to `auth_residues=True`/`auth_chains=True`; flipping to label numbering makes residue 100 a different residue, which is the classic "selection points at the wrong residue" bug.
- Bio.PDB's `structure.header` is thin; read resolution, method, and R-free from `MMCIF2Dict` (values are lists of strings, index `[0]`).
- MMTF is dead: RCSB stopped serving it on 2 July 2024. Use BinaryCIF for compact binary; do not build a live MMTF download path.
- The deposited coordinates for an X-ray entry are usually the asymmetric unit, not the biological assembly; "one chain in the file" is never evidence of a monomer.
- Bio.PDB cannot apply assembly operators -- download the `-assembly1.cif.gz` file from RCSB or use gemmi's `transform_to_assembly`.
- Parse unfamiliar files with `QUIET=False` and inspect the `PDBConstructionWarning`s; a permissive parse can hide a discontinuous chain or a numbering gap.

## Related Skills

- structure-navigation - Walk the SMCRA tree, handle altlocs and insertion codes, extract observed vs SEQRES sequence
- structure-modification - Transform coordinates, strip waters/hetero safely, edit B-factors before writing
- geometric-analysis - Measure distances, angles, SASA, and superimpose once the correct assembly is loaded
- interface-analysis - Analyze interfaces that only exist in the biological assembly, not the asymmetric unit
- structure-validation - Read resolution/R-free/clashscore to judge whether the loaded model is trustworthy
- alignment/structural-alignment - Superpose sequence-different structures that Bio.PDB's ordered correspondence cannot handle
- database-access/uniprot-access - Map structure residues back to a UniProt reference sequence
