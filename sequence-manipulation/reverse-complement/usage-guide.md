# Reverse Complement - Usage Guide

## Overview

This skill enables AI agents to help you generate complementary and reverse complementary sequences using Biopython. It covers the everyday operation (opposite strand, 5'->3'), the RNA-specific methods, IUPAC ambiguity codes, gapped alignment columns, in-place mutation, and the silent-corruption traps (protein sequences, double reverse-complementing minus-strand features, confusing complement with transcription).

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Get the reverse complement of this sequence"
- "Show me both strands of this DNA"
- "Is this a palindromic sequence?"
- "Design a reverse primer for this region"
- "Get the coding sequence for this minus-strand gene"
- "Reverse complement every sequence in this FASTA file"

## Example Prompts

### Basic Reverse Complement
> "What is the reverse complement of ATGCGATCGATCG?"

### Visualize Double-Stranded DNA
> "Show me this sequence as double-stranded DNA with both strands"

### Palindrome Check
> "Is GAATTC a palindromic sequence? What about ATGCAT?"

### RNA Strand
> "Give me the reverse complement of this RNA but keep it as U, not T"

### Primer Design
> "I need to amplify positions 100-500 of this sequence. What would the reverse primer look like?"

### Minus-Strand Feature
> "This CDS is annotated on the minus strand. Extract its coding sequence in the right orientation."

### Strand Conversion
> "Convert this template strand sequence to the coding strand"

## What the Agent Will Do

1. Wrap your input in a Bio.Seq.Seq object (or read records with SeqIO)
2. Confirm the molecule type so it does not silently reverse-complement a protein
3. Pick the method: reverse_complement() for DNA, reverse_complement_rna() to keep U, complement() for same-direction base pairs
4. For an annotated feature, use feature.extract() which already handles minus-strand orientation
5. Return the result in the conventional 5'->3' orientation

## Understanding the Methods

- **reverse_complement()**: What you usually want. Gives the opposite strand in 5'->3' direction. Runs in DNA mode, so any U is treated as T and emitted as T.
- **reverse_complement_rna()** / **complement_rna()**: Use these to keep RNA output (U, not T).
- **complement()**: Less common. Gives base pairs without reversing (3'->5' of the opposite strand).
- **transcribe()**: NOT a complement. It only swaps T->U on the same strand. To transcribe from a template strand: `template.reverse_complement().transcribe()`.

## Strand Terminology

```
Coding strand:     5'-ATGCGATCG-3'  (matches mRNA, except T instead of U)
Template strand:   3'-TACGCTAGC-5'  (used for transcription)
                      |||||||||
```

`reverse_complement()` of the coding strand = the template strand (written 5'->3').

## Tips

- For primer design, the reverse primer is the reverse complement of your target region's 3' end.
- Restriction enzyme sites are often palindromic - they equal their own reverse complement (S, W, and N bases are themselves self-complementary).
- When searching for motifs, search both the sequence and its reverse complement.
- Biopython handles all 15 IUPAC ambiguity codes correctly and case-insensitively. Never hand-roll the complement table - B/V and D/H are easy to swap and the error is silent.
- complement() and reverse_complement() do NOT validate the alphabet. Gaps (`-`) pass through unchanged (good for aligned sequences), but a protein passed in by mistake produces silent garbage with no warning - guard on molecule_type.
- After feature.extract() on a minus-strand feature, do NOT reverse-complement again; extract() already did it.
- inplace=True only works on a MutableSeq; on an immutable Seq it raises TypeError.

## Related Skills

- seq-objects - Create and mutate Seq/MutableSeq objects to complement
- transcription-translation - transcribe() vs complement(); six-frame translation uses the reverse complement
- motif-search - Search both strands by reverse-complementing the query or sequence
- sequence-io/read-sequences - Parse FASTA/GenBank records before reverse-complementing
- primer-design/primer-basics - Reverse primers are the reverse complement of the target 3' end
- restriction-analysis/restriction-sites - Restriction sites are often palindromic (self-complementary)
