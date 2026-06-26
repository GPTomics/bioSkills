# Seq Objects - Usage Guide

## Overview

This skill enables AI agents to help create and manipulate biological sequence objects in Biopython. It covers the core Seq, MutableSeq, and SeqRecord classes, the breaking changes from the 1.78 alphabet removal, and the silent-corruption traps that follow from losing alphabet validation.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Create a Seq object from this DNA string and show its length"
- "I need a mutable sequence I can edit at position 5 in place"
- "Build a SeqRecord with id and description ready to write to FASTA"
- "Set the molecule_type so this record writes to GenBank without erroring"
- "Why does my old `from Bio.Alphabet import IUPAC` code raise ImportError?"

## Example Prompts

### Basic Sequence Creation
> "Create a Seq object from 'ATGCGATCGATCG', show its length, and count the G bases"

### Mutable Sequences
> "I need to replace position 5 with a G in this sequence and append three Ts"

### SeqRecord Creation
> "Create a SeqRecord with id 'gene1' and description 'Example gene' that I can write to a file"

### Annotated Records
> "Build a SeqRecord with organism 'E. coli' and molecule_type set so it writes to GenBank"

### Migration and Debugging
> "My script does `Seq('ACGT', IUPAC.unambiguous_dna)` and crashes on a new Biopython - fix it"
> "A record from a lazy parser raises UndefinedSequenceError when I print it - how do I handle that?"

### Batch Processing
> "Create SeqRecords from this list of sequences with sequential ids"

## What the Agent Will Do

1. Import the appropriate classes (Seq, MutableSeq, SeqRecord)
2. Choose the object type from the task (immutable read, in-place edit, or annotated record)
3. Create the object and set any required attributes (id, description, molecule_type)
4. Perform the requested operations, preserving metadata where methods drop it
5. Return or display the result

## When to Use Each Object

- **Seq**: Most common. Immutable, string-like. Reading, analyzing, translating, slicing.
- **MutableSeq**: Editing individual positions in place. Cast back to Seq before writing.
- **SeqRecord**: When metadata (id, description, features, annotations) or file output is needed.

## Tips

- Writing to a file needs a SeqRecord, not a bare Seq.
- For GenBank/EMBL output, set `record.annotations['molecule_type']`; slicing a record drops it, so restore it on the slice.
- Seq is immutable: `seq[0] = 'A'` raises TypeError. Use MutableSeq, or rebuild by slicing/concatenation.
- Since 1.78 nothing validates the alphabet - reverse-complementing or transcribing a protein returns silent garbage. Track molecule type yourself.
- Since 1.79 Seq is bytes-backed, not a str subclass: `isinstance(seq, str)` is always False. Use `str(seq)` for text and test `isinstance(x, (Seq, MutableSeq))`.
- `Seq(None, length=n)` from lazy parsers has a length but no content; reading it raises UndefinedSequenceError (a ValueError subclass).
- Copy records with `deepcopy()` - plain assignment shares the features and annotations.

## Related Skills

- sequence-io/read-sequences - Parse files to get SeqRecord objects
- sequence-io/write-sequences - Write SeqRecord objects to files
- transcription-translation - Transform Seq objects (DNA to protein)
- reverse-complement - Get reverse complement of Seq
- sequence-slicing - Slice and extract from Seq/SeqRecord
- database-access/entrez-fetch - Fetch sequences from NCBI as SeqRecords
