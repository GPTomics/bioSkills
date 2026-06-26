# Format Conversion - Usage Guide

## Overview

This skill enables AI agents to help you convert between biological sequence file formats using Biopython, and to re-encode FASTQ quality offsets. The central caution is that conversion is lossy whenever the target format cannot represent the source's information, and Biopython performs the loss silently.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Convert this GenBank file to FASTA"
- "Convert all .gb files in this folder to FASTA"
- "Change my FASTQ to FASTA format"
- "Convert FASTA to GenBank with DNA molecule type"
- "Re-encode this old Illumina 1.3 FASTQ to Phred+33"
- "Will converting this GenBank to FASTA lose anything?"

## Example Prompts

### Simple Conversion
> "Convert sequence.gb to FASTA format"

### Batch Conversion
> "Convert all GenBank files in data/ to FASTA"

### Lossy-Conversion Check
> "I need a FASTA from this GenBank but I want to keep the gene features somewhere too"

### FASTQ Encoding
> "This FASTQ is Phred+64 from an old Illumina run, convert it to standard Phred+33"

### With Modifications
> "Convert to FASTA and uppercase all sequences"

## What the Agent Will Do
1. Parse sequences from the input format
2. Decide whether the target format can hold the source's information, and flag any silent loss
3. Use SeqIO.convert() for a streaming one-shot, or parse/write when records need modification
4. Add molecule_type or quality scores when converting up to a richer format
5. Write sequences in the target format and report the record count

## Tips

- Conversion is lossy whenever the target cannot represent the source: GenBank/EMBL to FASTA silently drops all features, qualifiers, annotations, and dbxrefs.
- Use SeqIO.convert() for plain conversions; it streams record-by-record and is more memory-efficient than parse + write. Use parse/write only when records change.
- FASTA to GenBank requires adding molecule_type; FASTA to FASTQ requires adding per-base quality scores (which are then fabricated, not measured).
- Never guess a FASTQ quality encoding: the wrong variant can silently shift every score by 31. Confirm Phred+33 vs Phred+64 vs Solexa from the pipeline before re-encoding.
- Phred to Solexa is lossy at low quality and warns when scores exceed the Solexa range; only re-encode to Solexa for legacy tools.
- Use AlignIO.convert() (not SeqIO) to convert alignments and keep gaps and columns.

## Related Skills

- read-sequences - Parse sequences and choose parse vs index for the input
- write-sequences - Write converted sequences with modifications
- fastq-quality - Phred/Solexa/Illumina encoding details and quality handling
- batch-processing - Convert many files across a directory
- compressed-files - Handle gzip/BGZF input and output during conversion
- alignment-files/sam-bam-basics - For SAM/BAM/CRAM conversion, use samtools view
