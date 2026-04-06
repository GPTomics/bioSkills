# alignment

## Overview

Pairwise and multiple sequence alignment: running MSA tools (MAFFT, MUSCLE5, ClustalOmega, T-Coffee), BioPython pairwise alignment, alignment I/O, and post-alignment analysis. Distinct from alignment-files which handles read-to-reference alignments (SAM/BAM).

**Tool type:** mixed | **Primary tools:** Bio.Align, Bio.AlignIO, MAFFT, MUSCLE5, ClustalOmega

## Skills

| Skill | Description |
|-------|-------------|
| multiple-alignment | Run MSA tools (MAFFT, MUSCLE5, ClustalOmega, T-Coffee) with algorithm selection guidance |
| pairwise-alignment | Global/local alignment using PairwiseAligner (Needleman-Wunsch, Smith-Waterman) |
| alignment-io | Read, write, convert MSA files (Clustal, PHYLIP, Stockholm, FASTA) |
| msa-parsing | Parse and analyze MSA content: gaps, conservation, filtering, consensus |
| msa-statistics | Calculate identity, conservation scores, entropy, substitution patterns |

## Example Prompts

- "Align these 50 protein sequences with the most accurate MSA method"
- "I have 5000 sequences, what MSA tool and settings should I use?"
- "Prepare a codon alignment for dN/dS analysis with codeml"
- "These sequences share about 30% identity, can I trust the alignment?"
- "Align these two DNA sequences and show the result"
- "Compare this protein to the reference using BLOSUM62"
- "Find the best matching region between these sequences"
- "Read this Clustal alignment and show sequence IDs"
- "Convert my PHYLIP alignment to FASTA format for RAxML"
- "Find conserved positions in this alignment"
- "Remove columns with more than 50% gaps"
- "Generate a consensus sequence"
- "Calculate pairwise identity matrix"
- "Show conservation score at each position"
- "Calculate Shannon entropy for each column"
- "Quantify alignment uncertainty before phylogenetic analysis"

## Requirements

```bash
pip install biopython numpy
conda install -c bioconda mafft muscle clustalo t-coffee pal2nal
```

## Related Skills

- **alignment-files** - Process SAM/BAM/CRAM read alignments
- **phylogenetics** - Build phylogenetic trees from MSAs
- **sequence-io** - Read input sequences for alignment
- **sequence-manipulation** - Work with individual sequences
