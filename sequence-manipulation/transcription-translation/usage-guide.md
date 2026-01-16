# Transcription and Translation - Usage Guide

## Overview

This skill enables AI agents to help you convert between DNA, RNA, and protein sequences using Biopython. It covers the central dogma operations: transcription, back-transcription, and translation with support for alternative codon tables.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Transcribe this DNA sequence to RNA"
- "Translate this coding sequence to protein"
- "Find all open reading frames in this sequence"
- "Translate using the mitochondrial codon table"

## Example Prompts

### Basic Transcription
> "Convert this DNA to RNA: ATGCGATCGATCG"

### Basic Translation
> "Translate this DNA sequence to protein"

### Stop Codon Handling
> "Translate this sequence but stop at the first stop codon"

### Alternative Genetic Codes
> "This is from E. coli, translate using the bacterial codon table"

### Mitochondrial Sequences
> "Translate this human mitochondrial sequence"

### ORF Finding
> "Find all open reading frames longer than 100 amino acids"

### Six-Frame Translation
> "Show me the translation in all six reading frames"

## What the Agent Will Do

1. Import Bio.Seq
2. Create a Seq object from your sequence
3. Apply the appropriate method (transcribe, translate, etc.)
4. Handle codon table selection if needed
5. Return the converted sequence(s)

## Codon Tables

Biopython includes all NCBI codon tables. Common ones:

- **1 (Standard)**: Most nuclear genes
- **2 (Vertebrate Mitochondrial)**: Human, mouse, etc. mitochondria
- **11 (Bacterial)**: E. coli, plastids, most prokaryotes
- **4 (Mold Mitochondrial)**: Fungi, some protists
- **5 (Invertebrate Mitochondrial)**: Insects, worms

## Tips

- Translation works on both DNA and RNA sequences
- Use `to_stop=True` to get clean protein without the `*` character
- Use `cds=True` for validated coding sequences (will error if invalid)
- For finding ORFs, search all six reading frames
- Mention the organism if using non-standard codon tables
