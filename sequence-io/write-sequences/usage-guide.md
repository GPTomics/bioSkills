# Write Sequences - Usage Guide

## Overview

This skill enables AI agents to help you write biological sequence data to files (FASTA, FASTQ, GenBank, EMBL) using Biopython. The recurring theme is that what survives a write depends on which SeqRecord fields are populated before the call, since each format reads only the fields it understands.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Create a FASTA file with these sequences"
- "Save the modified sequences to a new file"
- "Convert these records to GenBank format with organism annotation"
- "Write FASTQ output with quality scores"
- "Append this new sequence to my existing FASTA file"

## Example Prompts

### Basic Writing
> "Create a FASTA file with sequences for gene1, gene2, and gene3"

### Modifying and Saving
> "Read input.fasta, uppercase all sequences, and save to output.fasta"

### Format-Specific
> "Write these reads to FASTQ with a constant quality of Q40"

> "Write these sequences to GenBank format with molecule_type and organism set"

### Appending
> "Add this new sequence to my existing FASTA file without overwriting it"

## What the Agent Will Do

1. Import Bio.SeqIO and build SeqRecord objects with the required fields
2. Set format-specific fields: `letter_annotations['phred_quality']` for FASTQ, `annotations['molecule_type']` for GenBank/EMBL
3. Make the FASTA description start with the id so the header is not mangled
4. Write records with `SeqIO.write()` (filename or open handle, 'a' mode to append)
5. Confirm the record count returned matches what was expected

## Tips

- FASTA writes `>` + `record.description`, not `record.id`. If you set both, make the description start with the id plus a space, or the id is dropped from the header.
- The FASTA writer wraps sequence lines at 60 characters by default.
- FASTQ requires a quality score per base. `letter_annotations` is length-locked to the sequence, so set the sequence first, then a quality list of matching length.
- GenBank and EMBL require `annotations['molecule_type']` (the alphabet that carried this was removed in BioPython 1.78); writing without it raises.
- When both Phred and Solexa qualities are present, the writer uses Phred. Writing 'fastq-solexa' from a Phred-only record is lossy and warns at high quality; prefer plain 'fastq' for modern data.
- Use a file handle in 'a' mode to append records to an existing file.

## Related Skills

- read-sequences - Read sequences before modifying and writing
- format-conversion - Direct format conversion without intermediate processing
- filter-sequences - Filter sequences before writing a subset
- fastq-quality - Phred/Solexa encodings and quality-score handling
- sequence-manipulation/seq-objects - Create SeqRecord objects to write
- alignment-files/sam-bam-basics - For SAM/BAM output, use samtools/pysam
