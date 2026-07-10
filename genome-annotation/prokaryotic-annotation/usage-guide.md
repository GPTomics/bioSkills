# Prokaryotic Annotation - Usage Guide

## Overview

Annotate bacterial and archaeal genomes (isolates, MAGs, plasmids) with structural and functional gene predictions using Bakta (default) or Prokka (legacy), with NCBI PGAP and DFAST as the submission-grade paths. Produces GFF3, GenBank/EMBL, protein FASTA, and a feature table ready for downstream analysis or INSDC submission. The hard part is not finding genes (callers recover >95%) but assigning trustworthy function, picking the right translation table, and not comparing gene counts across mismatched tools or database versions.

## Prerequisites

```bash
# Bakta (default) + database (record the version)
conda install -c bioconda bakta
bakta_db download --output /path/to/bakta_db --type full   # ~30 GB; --type light for triage/CI

# Prokka (legacy)
conda install -c bioconda prokka

# QC gate (run before trusting annotation numbers)
conda install -c bioconda checkm2

# Python parsing
pip install gffutils biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Annotate my bacterial genome assembly with Bakta and report coding density"
- "Annotate my Mycoplasma genome with the correct translation table"
- "Re-annotate this collection of Klebsiella genomes uniformly for pangenomics"
- "Prepare my genome annotation for NCBI submission"

## Example Prompts

### Basic Annotation

> "Annotate this bacterial assembly with Bakta using the full database, gram-negative, locus tag prefix SAUR, and give me the GFF3 plus a coding-density and hypothetical-protein QC."

> "Run Bakta on my Staphylococcus genome and flag whether the coding density and gene count are in the expected range."

### Organism-Specific Cases

> "Annotate my Mycoplasma genitalium assembly - it uses genetic code 4 (UGA=Trp), so make sure the translation table is set correctly."

> "Annotate this archaeon with Bakta and confirm the archaeal tRNA model was used."

### Collections and Submission

> "Re-annotate all of these E. coli assemblies from FASTA with one Bakta version and DB so the pangenome accessory genome isn't inflated by annotation drift."

> "Prepare my Bakta annotation for NCBI GenBank submission with INSDC-compliant locus tags."

## What the Agent Will Do

1. Gate on assembly quality first (CheckM2 completeness/contamination); flag if contamination >5% or completeness <90%
2. Determine the correct translation table from taxonomy (GTDB-Tk if available), not a guess
3. Run Bakta (or Prokka/PGAP/DFAST per the target) with organism-specific parameters and a recorded DB version
4. Generate GFF3, GenBank, protein FASTA, and summary outputs
5. Report annotation statistics (gene count, coding density, tRNA/rRNA counts, hypothetical fraction)
6. Flag failure signatures (low coding density -> wrong table/fragmentation; pseudogene spike -> homopolymer/phase-variation/decay; near-0% hypothetical -> over-confident transfer)
7. Suggest functional-annotation follow-up (eggNOG-mapper/InterProScan) for hypotheticals

## Tips

- **Database version determines the result** - Record the Bakta DB version in methods; annotation content tracks the DB, not just the binary. Two runs months apart can differ purely from DB updates.
- **Never compare gene counts across tools/versions** - For any collection, re-annotate from FASTA with one pipeline + one pinned DB, then cluster with an error-aware tool (Panaroo). Merging published annotations inflates the accessory genome ~10x.
- **Translation table from taxonomy** - Default 11; table 4 for Mycoplasma/Mollicutes (UGA=Trp), table 25 for Gracilibacteria/SR1 (UGA=Gly). A wrong table is silent and looks like fragmentation.
- **"Hypothetical" is honest** - 25-50% hypothetical is healthy for a non-model organism; near 0% means over-confident transfer. Open Bakta's `.inference.tsv` to see why a product was assigned, and prefer the InterPro/Pfam domain architecture over a loose free-text product name.
- **Pseudogene spikes are ambiguous** - Could be reductive biology (symbionts), phase variation (homopolymer/SSR tracts - do NOT polish away), or ONT homopolymer indel error. Disambiguate with short-read/HiFi corroboration.
- **RefSeq won't match your GFF** - RefSeq re-annotates submissions with PGAP; report the stable WP_ protein accession alongside the locus tag.
- **macOS signal peptides** - DeepSig was dropped from Bakta's default mac conda env; signal-peptide calls can be silently absent. Run on Linux if they matter.
- **Complete vs draft** - `--complete` only for finished/circularized replicons; using it on draft contigs mislabels everything.

## Related Skills

- genome-annotation/functional-annotation - Add GO/KEGG/Pfam to hypothetical proteins
- genome-annotation/ncrna-annotation - Detailed ncRNA identification with Infernal/Rfam
- genome-annotation/annotation-qc - CheckM2 gate and gene-set sanity metrics
- genome-assembly/assembly-qc - Check assembly quality before annotation
- comparative-genomics/pangenome-analysis - Uniform re-annotation before pangenome clustering
- genome-intervals/gtf-gff-handling - Parse GFF3 output files
