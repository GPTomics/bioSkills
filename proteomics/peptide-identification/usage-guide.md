# Peptide Identification - Usage Guide

## Overview
Match MS/MS spectra to peptide sequences by database search (or spectral-library search), then control false discovery rate with target-decoy competition. The deliverable an agent should act on is a q-value (list-level error) or PEP (per-ID error), never a raw engine score, and PSM-level FDR is a separate problem from protein-level FDR.

## Prerequisites
```bash
pip install pyopenms pandas numpy
# CLI engines: comet, msfragger, sage, MSGFPlus; rescoring: percolator, mokapot
# R alternative: BiocManager::install(c("mzID", "mzR", "PSMatch"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run a concatenated target-decoy database search on my mzML against UniProt human"
- "Configure trypsin with 2 missed cleavages, 10 ppm precursor and 0.02 Da fragment tolerance"
- "Annotate decoys and filter PSMs to 1% FDR with a proper q-value"
- "Rescore my Sage results with mokapot to gain IDs at the same FDR"
- "Explain why my PEP <= 0.01 cutoff kept so many fewer peptides than q <= 0.01"

## Example Prompts

### Database Search Setup
> "Configure a database search with trypsin, 2 missed cleavages, carbamidomethyl C fixed and oxidation M variable"

> "Set up an MSFragger open search from -150 to +500 Da to discover unknown modifications"

### Running Searches
> "Run a peptide search against a concatenated target-decoy human FASTA with pyOpenMS"

> "Pick a search engine for varied fragmentation across instruments and justify the choice"

### FDR and Rescoring
> "Annotate target/decoy from the DECOY_ prefix and filter to 1% peptide FDR"

> "Compute q-values from this PSM table assuming a concatenated competition search"

> "Rescore the search with Percolator and report the ID gain at q <= 0.01"

### Results Processing
> "Read the mzIdentML output in R and filter to 1% FDR"

> "Explain the difference between PEP and q-value for these PSMs"

## What the Agent Will Do
1. Confirm the search mode (concatenated vs separate target-decoy) so the correct FDR formula is used.
2. Configure search parameters (enzyme, missed cleavages, precursor/fragment tolerances, fixed/variable mods).
3. Run the database or spectral-library search and produce scored PSMs.
4. Annotate targets vs decoys, then estimate q-values by target-decoy competition.
5. Optionally rescore with Percolator/mokapot to boost IDs at the same FDR.
6. Filter to 1% peptide FDR on q-value, remove decoys, and report identification statistics.
7. Hand protein grouping and protein-level FDR to protein-inference.

## Tips
- Act on a q-value for list cutoffs and a PEP only for per-ID decisions; "PEP <= 0.01" is far stricter than "q <= 0.01".
- Never threshold on a raw score (XCorr, hyperscore, Andromeda, SpecEValue) or compare scores across engines.
- Concatenated competition uses FDR = decoys/targets; separate searches use the mix-max / 2x form -- do not mix them.
- Generate decoys at the protein level then digest, so decoy peptides obey the same enzyme rules.
- Below ~hundreds of PSMs the decoy FDR is unreliable; a single-protein pulldown is essentially uninformative.
- Open/mass-tolerant search is for discovery; follow with a closed search on the discovered mods before reporting clean FDR.
- PSM FDR at 1% does not give 1% protein FDR; estimate protein-level (picked) FDR separately.
- Do not impose a blanket "2 unique peptides" rule; let FDR and score decide (the two-peptide rule increases protein FDR).

## Related Skills
- protein-inference - Group peptides to protein groups and control protein-level (picked) FDR
- ptm-analysis - Open/variable-mod search follow-up and per-site PTM localization
- dia-analysis - DIA peptide-centric extraction and scoring; entrapment FDR validation
- quantification - FDR-filtered IDs feed label-free/TMT intensity quantification
- spectral-libraries - Empirical and predicted spectral-library search as an ID alternative
- data-import - Load mzML/raw MS data before identification
- database-access/uniprot-access - Build the target FASTA (canonical vs isoform, contaminants)
