# sequence-manipulation

## Overview

Working with sequence data programmatically using Biopython's Bio.Seq and Bio.SeqUtils modules. Handles transcription, translation, reverse complement, motif finding, and sequence property calculations.

**Tool type:** python | **Primary tools:** Bio.Seq, Bio.SeqUtils

## Skills

| Skill | Description |
|-------|-------------|
| seq-objects | Create and modify Seq, MutableSeq, and SeqRecord objects |
| transcription-translation | DNA to RNA to protein, codon tables, ORF finding |
| reverse-complement | Reverse complement, complement, palindrome detection |
| sequence-slicing | Slice, extract, concatenate, and manipulate sequences |
| motif-search | Find patterns with regex, PWM/PSSM, parse JASPAR/MEME files |
| sequence-properties | GC content, GC skew, molecular weight, melting temperature |
| codon-usage | Codon Adaptation Index, RSCU, codon optimization |

## Example Prompts

- "Create a Seq object from this DNA string"
- "Transcribe this DNA sequence to RNA"
- "Translate this coding sequence to protein"
- "Translate using the mitochondrial codon table"
- "Find all open reading frames in this sequence"
- "Get the reverse complement of this sequence"
- "Extract positions 100-200 from this sequence"
- "Join these sequences together with a linker"
- "Find all occurrences of GAATTC in my sequence"
- "Parse this JASPAR motif file and search my sequence"
- "Calculate the GC content of each sequence in my FASTA"
- "Plot GC skew along this sequence to find the origin"
- "Calculate the molecular weight of this protein"
- "Analyze this protein: pI, stability, hydropathy"
- "What is the codon usage bias in this gene?"
- "Calculate the CAI for E. coli expression"
- "Optimize this gene for yeast expression"

## Requirements

```bash
pip install biopython
```

## Related Skills

- **sequence-io** - Read sequences from files before manipulation
- **restriction-analysis** - Restriction enzyme analysis using Bio.Restriction
- **alignment** - Align sequences for comparison
- **database-access** - Fetch sequences from NCBI for analysis
