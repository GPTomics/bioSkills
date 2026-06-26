# sequence-manipulation

## Overview

Working with sequence data programmatically using Biopython's Bio.Seq and Bio.SeqUtils modules. Handles transcription, translation, reverse complement, motif finding, and sequence property calculations.

**Tool type:** python | **Primary tools:** Bio.Seq, Bio.SeqUtils

## Skills

| Skill | Description |
|-------|-------------|
| seq-objects | Create and modify Seq, MutableSeq, and SeqRecord objects; the post-1.78 no-alphabet model and molecule_type |
| transcription-translation | DNA to RNA to protein with NCBI codon-table selection, cds=True validation, and ORF finding |
| reverse-complement | Reverse complement and complement DNA/RNA, including IUPAC ambiguity codes and minus-strand features |
| sequence-slicing | Slice, extract, and concatenate sequences and annotated records, preserving what survives a SeqRecord slice |
| motif-search | Find degenerate IUPAC patterns and TF binding sites with regex and PWM/PSSM scoring; parse JASPAR/MEME files |
| sequence-properties | GC content/skew, molecular weight, and nearest-neighbor melting temperature; protein biophysical properties |
| codon-usage | Codon Adaptation Index (CAI), RSCU, Nc, and max-CAI optimization with its expression tradeoffs |

## Example Prompts

- "Create a Seq object from this DNA string"
- "Transcribe this DNA sequence to RNA"
- "Translate this coding sequence, validating it as a complete CDS"
- "Translate this mitochondrial gene with the correct codon table"
- "Find all open reading frames in this sequence"
- "Get the reverse complement of this sequence"
- "Extract this minus-strand CDS in the correct orientation"
- "Extract positions 100-200 from this record, keeping its quality scores"
- "Find all overlapping occurrences of GCGC in my sequence"
- "Build a PWM from these sites and scan both strands above a 1% false-positive threshold"
- "Calculate the GC content of each sequence, excluding ambiguous bases"
- "Find the replication origin from cumulative GC skew"
- "Compute the nearest-neighbor Tm for this primer with a Mg correction"
- "Analyze this protein: pI, stability, hydropathy"
- "Build a CAI index from highly expressed host genes and score my gene"
- "Codon-optimize this gene for the host, then flag the tradeoffs to screen"

## Requirements

```bash
pip install biopython
```

## Related Skills

- **sequence-io** - Read sequences from files before manipulation
- **restriction-analysis** - Restriction enzyme analysis using Bio.Restriction
- **alignment** - Align sequences for comparison
- **database-access** - Fetch sequences from NCBI for analysis
