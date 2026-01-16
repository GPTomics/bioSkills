# Alignment

Sequence alignment using Biopython's Bio.Align and Bio.AlignIO modules for pairwise and multiple sequence alignments.

## Overview

This category covers sequence-to-sequence alignment: pairwise alignment algorithms, reading/writing alignment files, parsing multiple sequence alignments, and calculating alignment statistics. Distinct from alignment-files which handles read-to-reference alignments (SAM/BAM).

**Tool type:** `python`
**Primary tools:** Bio.Align, Bio.AlignIO

## Skills

| Skill | Description |
|-------|-------------|
| [pairwise-alignment](pairwise-alignment/) | Global/local alignment using PairwiseAligner (Needleman-Wunsch, Smith-Waterman) |
| [alignment-io](alignment-io/) | Read, write, convert MSA files (Clustal, PHYLIP, Stockholm, FASTA) |
| [msa-parsing](msa-parsing/) | Parse and analyze MSA content: gaps, conservation, filtering, consensus |
| [alignment-statistics](alignment-statistics/) | Calculate identity, conservation scores, entropy, substitution patterns |

## Workflow

```
Sequences (FASTA, GenBank)
    |
    +---> [pairwise-alignment] - Compare two sequences
    |         |
    |         v
    |     Score, aligned sequences
    |
    v
[External MSA tool] - ClustalW, MUSCLE, MAFFT
    |
    v
MSA file (.aln, .phy, .sto)
    |
    v
[alignment-io] - Read, convert formats
    |
    +---> [msa-parsing] - Analyze gaps, filter sequences, consensus
    |         |
    |         v
    |     Cleaned alignment
    |
    v
[alignment-statistics] - Identity matrix, conservation, entropy
    |
    v
[phylogenetics] - Tree building, distance matrices
```

## Supported Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| `clustal` | .aln | ClustalW/Omega output |
| `fasta` | .fasta, .fa | Aligned FASTA |
| `phylip` | .phy | PHYLIP (interleaved) |
| `phylip-relaxed` | .phy | PHYLIP with long names |
| `stockholm` | .sto, .stk | Pfam/Rfam annotated |
| `nexus` | .nex | NEXUS for MrBayes/PAUP |
| `maf` | .maf | Multiple Alignment Format |

## Example Prompts

### Pairwise Alignment
- "Align these two DNA sequences and show the result"
- "Compare this protein to the reference using BLOSUM62"
- "Find the best matching region between these sequences"
- "What is the alignment score between seq1 and seq2?"

### Alignment I/O
- "Read this Clustal alignment and show sequence IDs"
- "Convert my PHYLIP alignment to FASTA format"
- "Extract columns 100-200 from the alignment"
- "Save this alignment as Stockholm format"

### MSA Analysis
- "Find conserved positions in this alignment"
- "Remove columns with more than 50% gaps"
- "Generate a consensus sequence"
- "Filter out sequences with too many gaps"

### Statistics
- "Calculate pairwise identity matrix"
- "Show conservation score at each position"
- "What is the transition/transversion ratio?"
- "Calculate Shannon entropy for each column"

## Requirements

```bash
pip install biopython numpy
```

## Notes

- **Deprecated `AlignInfo.SummaryInfo` avoided** - These skills use custom implementations for consensus, PSSM, and information content instead of the deprecated `Bio.Align.AlignInfo.SummaryInfo` class
- **Modern API** - Uses `Alignment.counts()` and `Alignment.substitutions` for pairwise alignment statistics
- **Both I/O modules documented** - `Bio.AlignIO` (legacy) and `Bio.Align` (modern) read/write functions are covered

## Alignment vs Alignment-Files

| This Category (alignment) | alignment-files |
|---------------------------|-----------------|
| Sequence-to-sequence | Read-to-reference |
| Pairwise/MSA | NGS alignments |
| Bio.Align, Bio.AlignIO | samtools, pysam |
| Clustal, PHYLIP, Stockholm | SAM, BAM, CRAM |
| Phylogenetics, conservation | Variant calling, coverage |

## Related Skills

- **alignment-files** - Process SAM/BAM/CRAM read alignments
- **phylogenetics** - Build phylogenetic trees from MSAs
- **sequence-io** - Read input sequences for alignment
- **sequence-manipulation** - Work with individual sequences

## References

- [Bio.Align documentation](https://biopython.org/docs/latest/api/Bio.Align.html)
- [Bio.AlignIO documentation](https://biopython.org/docs/latest/api/Bio.AlignIO.html)
