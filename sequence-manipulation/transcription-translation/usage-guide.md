# Transcription and Translation - Usage Guide

## Overview

This skill enables AI agents to help you convert between DNA, RNA, and protein sequences using Biopython. It covers transcription, back-transcription, and translation, with emphasis on the decisions that silently corrupt results: choosing the correct NCBI genetic code, validating coding sequences, and handling partial codons, gaps, and recoded stops.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Transcribe this DNA sequence to RNA"
- "Translate this coding sequence to protein"
- "Translate this human mitochondrial gene with the right genetic code"
- "Validate that this is a complete CDS and translate it"
- "Find all open reading frames longer than 100 amino acids"
- "Show me the translation in all six reading frames"

## Example Prompts

### Basic Transcription
> "Convert this DNA to RNA: ATGCGATCGATCG"

### Template-Strand Transcription
> "This is the template strand, give me the mRNA it produces"

### Basic Translation
> "Translate this DNA sequence to protein"

### Stop Codon Handling
> "Translate this sequence but stop at the first stop codon"

### Alternative Genetic Codes
> "This is from E. coli, translate it with the bacterial code"

### Mitochondrial Sequences
> "Translate this human mitochondrial sequence with the correct table"

### CDS Validation
> "Check that this is a valid complete CDS, then translate it"

### ORF Finding
> "Find all open reading frames longer than 100 amino acids"

### Six-Frame Translation
> "Show me the translation in all six reading frames"

## What the Agent Will Do

1. Confirm the molecule type and strand (no type checks exist since Biopython 1.78)
2. Create a Seq object from your sequence
3. Select the correct NCBI codon table for the organism if it is not standard
4. Apply transcribe, back_transcribe, or translate with the right options
5. Use `cds=True` to validate complete coding sequences and surface defects loudly
6. Return the converted sequence(s)

## Codon Tables

Biopython includes every NCBI genetic code. The choice changes the protein, and a valid-but-wrong table produces a plausible wrong result with no error:

- **1 (Standard)**: Most nuclear genes
- **2 (Vertebrate Mitochondrial)**: Human/vertebrate mtDNA; AGA/AGG are stops, AUA=Met, UGA=Trp
- **3 (Yeast Mitochondrial)**: CTG=Thr lives here (not table 12)
- **4 (Mold/Protozoan Mito + Mycoplasma)**: UGA=Trp
- **5 (Invertebrate Mitochondrial)**: AGA/AGG=Ser
- **6 (Ciliate Nuclear)**: UAA/UAG=Gln
- **11 (Bacterial/Archaeal/Plastid)**: expanded start codons
- **12 (Alternative Yeast Nuclear)**: CUG=Ser (not Thr)

## Tips

- Translation works directly on both DNA and RNA, so explicit transcription is rarely needed first.
- A wrong codon table is the most common silent bug; pick the organism's table deliberately.
- Use `cds=True` for complete coding sequences; it raises on a bad start, length, or stop instead of returning a wrong protein silently.
- An alternative start (GTG/TTG/ATT) under `cds=True` is translated as M, which is biologically correct.
- A length not divisible by 3 only warns and drops trailing bases; `cds=True` turns that into a loud error.
- `transcribe()` and `back_transcribe()` only swap T and U; they do not splice, cap, or add poly-A. For real transcription from the template strand, reverse-complement first.
- Selenocysteine (UGA) and pyrrolysine (UAG) are not recoded by Biopython; they appear as `*` or silently truncate with `to_stop=True`.

## Related Skills

- seq-objects - Create and inspect Seq objects before translation
- reverse-complement - Strand handling for six-frame translation and template-strand transcription
- codon-usage - Analyze codon bias and adaptation in coding sequences
- sequence-io/read-sequences - Parse GenBank/FASTA records and CDS features for translation
