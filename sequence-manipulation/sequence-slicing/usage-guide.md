# Sequence Slicing - Usage Guide

## Overview

This skill enables AI agents to help you extract, slice, and concatenate biological sequences with Biopython. It covers bare-Seq slicing, carrying a sub-region of an annotated record (with quality scores and features) into a new record, coordinate-system conversion, and joining sequences. The central judgment it encodes is what metadata survives a SeqRecord slice and how to avoid the 0-based vs 1-based off-by-one.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Extract positions 100-200 from this sequence"
- "Pull out this CDS but keep its per-base quality scores"
- "The gene is at 1234-5678 in the GFF; extract it with the right coordinates"
- "Splice these exon regions into one transcript"
- "Join these sequences with an NNN linker"
- "Split this sequence into 100 bp chunks"

## Example Prompts

### Basic Extraction
> "Extract nucleotides 100 to 500 from my sequence and return a Seq."

### Coordinates From a File
> "My GFF says the feature spans 1234 to 5678 (1-based inclusive). Slice the matching subsequence without an off-by-one."

> "This interval came from a BED file (0-based). Extract it correctly."

### Keep Annotations and Qualities
> "Slice positions 100-400 out of this FASTQ record but keep the per-base quality scores."

> "Extract this region of a GenBank record and carry the organism and taxonomy into the new record."

### Splicing and Features
> "Splice these three exon coordinate pairs into a single mRNA sequence."

> "Extract every CDS feature from this GenBank file, strand-aware."

### Splitting and Joining
> "Split this coding sequence into codons."

> "Concatenate these sequences with NNNNNN between each."

> "Give me 50 bp on each side of position 1000."

## What the Agent Will Do

1. Decide whether the request needs a bare `Seq` slice or a `SeqRecord` slice (qualities/features to keep).
2. Convert coordinates to 0-based half-open, subtracting 1 from the start for 1-based GFF/GenBank input.
3. Apply Python slicing, or `feature.extract()` for strand-aware/joined locations.
4. Copy the `annotations` dict explicitly when carrying a sub-region into a new record.
5. Return the extracted, spliced, or concatenated sequence or record.

## Coordinate System Notes

Python and Biopython slice 0-based and half-open: `seq[0]` is the first base, `seq[0:3]` returns three bases (0, 1, 2), `seq[100:200]` returns 100 bases. File formats differ:

- GenBank / EMBL / GFF / GTF / VCF feature coordinates are 1-based and inclusive. Convert by subtracting 1 from the START only: `seq[start-1:end]`. The exclusive Python end already cancels the inclusive file end, so the end needs no change.
- BED coordinates are already 0-based half-open and slice directly.
- `Bio.SeqFeature` locations store 0-based start and Python-style end, so `feature.extract(record.seq)` slices correctly and handles strand and joined exons.

## What a SeqRecord Slice Keeps and Drops

`record[start:end]` preserves `id`, `name`, `description`, and `molecule_type`, and auto-slices per-letter data (PHRED quality in `letter_annotations`) to match. It keeps only features fully contained in the range, with recalculated locations. It SILENTLY drops the `annotations` dict (organism, taxonomy, references), `dbxrefs`, and any feature straddling the boundary. A stride such as `record[::2]` drops features entirely. Copy what you need: `sub.annotations = record.annotations.copy()`.

## Tips

- Always clarify whether coordinates are 0-based or 1-based before slicing; the off-by-one is silent.
- Slice the `SeqRecord` (not `record.seq`) when per-base qualities or contained features must ride along.
- Copy `annotations` after a SeqRecord slice; slicing always drops it without warning.
- `seq[::-1]` reverses but does NOT complement; use `reverse_complement()` for the opposite strand.
- Use `feature.extract(record.seq)` for any strand-aware or joined/compound location instead of a manual slice.
- `sum()` over slices needs a `Seq('')` start value; the default `0` cannot be added to a `Seq`.
- Guard `str(record.seq)` with `len()` first; undefined content raises `UndefinedSequenceError`.

## Related Skills

- seq-objects - Create Seq/SeqRecord objects and handle undefined sequence content
- reverse-complement - Reverse-complement an extracted region
- transcription-translation - Translate an extracted CDS or spliced transcript
- sequence-io/read-sequences - Parse GenBank/FASTQ records to slice
- genome-intervals/gtf-gff-handling - Read 1-based GFF/GTF feature coordinates before slicing
- alignment-files/sam-bam-basics - Extract sequences from BAM regions with samtools
