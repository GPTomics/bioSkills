# Reverse Complement - Usage Guide

## Overview

This skill enables AI agents to help you generate complementary and reverse complementary sequences using Biopython. Essential for primer design, working with opposite DNA strands, and identifying palindromic sequences like restriction sites.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Get the reverse complement of this sequence"
- "Show me both strands of this DNA"
- "Is this a palindromic sequence?"
- "Design primers for this region"

## Example Prompts

### Basic Reverse Complement
> "What is the reverse complement of ATGCGATCGATCG?"

### Visualize Double-Stranded DNA
> "Show me this sequence as double-stranded DNA with both strands"

### Palindrome Check
> "Is GAATTC a palindromic sequence? What about ATGCAT?"

### Primer Design
> "I need to amplify positions 100-500 of this sequence. What would the forward and reverse primers look like?"

### Strand Conversion
> "Convert this template strand sequence to the coding strand"

## What the Agent Will Do

1. Import Bio.Seq
2. Create a Seq object from your sequence
3. Apply the appropriate method (reverse_complement, complement)
4. Return the result in the proper orientation

## Understanding the Methods

- **reverse_complement()**: What you usually want. Gives the opposite strand in 5' to 3' direction.
- **complement()**: Less common. Gives base pairs without reversing (3' to 5' of opposite strand).

## Strand Terminology

```
Coding strand:     5'-ATGCGATCG-3'  (matches mRNA, except T instead of U)
Template strand:   3'-TACGCTAGC-5'  (used for transcription)
                      |||||||||
```

`reverse_complement()` of coding strand = template strand (but written 5' to 3')

## Tips

- For primer design, the reverse primer should be the reverse complement of your target region's 3' end
- Restriction enzyme sites are often palindromic - they equal their own reverse complement
- When searching for motifs, search both the sequence and its reverse complement
- Biopython correctly handles IUPAC ambiguity codes (R, Y, N, etc.)
