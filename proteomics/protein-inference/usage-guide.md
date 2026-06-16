# Protein Inference - Usage Guide

## Overview
Protein inference decides which proteins are present from identified peptides. It is fundamentally underdetermined: bottom-up MS sees peptides, many peptides map to multiple proteins, so the protein set is a chosen explanation under an assumption (parsimony or a probability model), not a measurement. The honest reporting unit is a protein GROUP -- proteins indistinguishable by the observed peptides -- with a designated leading protein, not a flat list. This skill covers grouping, inference-method choice, and protein/group-level FDR done correctly (picked FDR, no two-peptide rule).

## Prerequisites
```bash
pip install pyopenms pandas
# CLI alternatives: ProteinProphet (TPP), Epifany (OpenMS), Philosopher (FragPipe)
```

## Quick Start
Tell your AI agent what you want to do:
- "Group these proteins by shared peptide evidence and give each group a leading protein"
- "Apply parsimony to find the minimal protein set explaining all my peptides"
- "Run EPIFANY for probabilistic inference with protein-group FDR"
- "Estimate protein-group FDR with picked FDR, not the PSM formula"
- "Quantify on unique peptides only to avoid razor-peptide artifacts"

## Example Prompts

### Protein Grouping
> "Parse MaxQuant proteinGroups.txt and explain the protein group structure and leading proteins"

> "Identify protein groups whose proteins share all peptide evidence (indistinguishable)"

> "List proteins with unique peptides versus those identified only by shared peptides"

### Inference Methods
> "Apply the parsimony principle to report the minimal protein set explaining all peptides"

> "Run EPIFANY for Bayesian protein inference with protein-group FDR control"

> "Use ProteinProphet to apportion shared peptides across candidate proteins"

### FDR Control
> "Estimate protein-group FDR with picked-group FDR and report groups at 1%"

> "Explain why my 1% PSM FDR left a much higher protein FDR on this deep dataset"

> "Flag single-peptide identifications and score them individually instead of dropping them"

### Reporting
> "Build a protein-group table with the leading accession, unique-peptide count, and group members"

> "Summarize how many protein groups pass at 1% picked FDR"

## What the Agent Will Do
1. Load the FDR-filtered peptide-to-protein evidence (idXML, proteinGroups.txt, or a PSM table)
2. Build protein groups, collapsing indistinguishable proteins and dropping subsumable ones
3. Apply the chosen inference method (parsimony or probabilistic) and assign a leading protein per group
4. Estimate protein/protein-group FDR with picked FDR, not the reused PSM formula
5. Report groups at the FDR cutoff, judging single-peptide IDs by score rather than a two-peptide rule
6. Hand the inferred groups to quantification, noting the razor-vs-unique choice

## Tips
- Report protein groups with a leading protein, never a flat protein list -- a flat list double-counts indistinguishable proteins and breaks protein-level FDR.
- Protein FDR is a separate estimation problem from peptide FDR and inflates on large data; use picked-group FDR (Savitski 2015; The & Kall 2016).
- Do not impose a two-peptide rule -- it raises protein FDR and discards real single-peptide proteins (Gupta & Pevzner 2009).
- "Unique" is database-relative; a peptide unique against SwissProt may be shared once isoforms and TrEMBL are added. Fix and document the database.
- For sensitive differential abundance, quantify on unique peptides only; razor assignment can flip between conditions and fake DE.
- A protein group is not a proteoform -- do not claim a specific isoform from shared-peptide evidence.

## Related Skills
- peptide-identification - Produces the FDR-filtered peptide list that feeds inference and shares the target-decoy machinery
- quantification - Consumes inferred groups; razor-vs-unique peptide choice lives here
- data-import - Loads idXML/mzML identification files
- database-access/uniprot-access - Canonical-vs-isoform databases drive uniqueness and the leading-protein convention
