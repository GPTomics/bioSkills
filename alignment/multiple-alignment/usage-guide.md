# Multiple Sequence Alignment - Usage Guide

## Overview

This skill performs multiple sequence alignment (MSA) of three or more homologous sequences using command-line tools (MAFFT, MUSCLE5, ClustalOmega, T-Coffee). It guides tool selection based on dataset size, sequence divergence, and whether the alignment is destined for phylogenetics, selection analysis, or conservation studies.

## Prerequisites

```bash
# Install MSA tools (conda recommended)
conda install -c bioconda mafft muscle clustalo t-coffee pal2nal

# Or individually
conda install -c bioconda mafft        # Most versatile, required
conda install -c bioconda "muscle>=5"  # Highest accuracy
conda install -c bioconda clustalo     # Large datasets
conda install -c bioconda pal2nal      # Codon-aware alignment

pip install biopython  # For reading/writing alignment files
```

## Quick Start

Tell your AI agent what you want to do:
- "Align these protein sequences with the most accurate method available"
- "I have 5000 sequences to align, use something fast but reasonable"
- "Prepare a codon alignment for dN/dS analysis with codeml"

## Example Prompts

### Tool Selection
> "Align these 50 protein sequences with maximum accuracy"

> "I need to align 10,000 16S rRNA sequences, what tool and settings should I use?"

> "These sequences share about 30% identity, can I trust the alignment?"

### MAFFT
> "Run MAFFT L-INS-i on my protein FASTA file with 8 threads"

> "My sequences have conserved domains separated by variable linker regions, which MAFFT algorithm handles that?"

> "Add these 10 new sequences to my existing alignment without re-aligning everything"

### MUSCLE5
> "Align my sequences with MUSCLE5 and generate confidence scores for each column"

> "Use MUSCLE5 ensemble mode to quantify alignment uncertainty before phylogenetic analysis"

### Codon-Aware Alignment
> "I need to prepare a codon alignment for PAML. Align proteins first then thread the DNA"

> "My sequences may contain frameshifts, use MACSE for codon-aware alignment"

### Quality Assessment
> "Check if this alignment is reliable enough for phylogenetic inference"

> "Remove sequences that are making the alignment worse"

## What the Agent Will Do

1. Assess dataset characteristics (number of sequences, expected divergence, sequence type)
2. Select appropriate tool and algorithm (MAFFT mode, MUSCLE5, ClustalOmega)
3. Configure parameters (threads, gap penalties, output format)
4. Run the alignment
5. Validate output (check for excessive gaps, non-homologous sequences, alignment quality)
6. For codon alignments: align protein sequences first, then thread DNA with PAL2NAL

## Tool Selection Guide

| Scenario | Recommended Tool |
|----------|-----------------|
| <200 sequences, accuracy priority | MAFFT L-INS-i |
| <200 sequences, conserved motifs in variable regions | MAFFT E-INS-i |
| 200-10,000 sequences | MAFFT FFT-NS-2 or MUSCLE5 -super5 |
| >10,000 sequences | ClustalOmega or MAFFT --retree 2 |
| Need alignment confidence scores | MUSCLE5 ensemble mode |
| Structural information available | T-Coffee Expresso mode |
| Coding sequences for selection analysis | Protein MSA + PAL2NAL |
| Sequences with frameshifts | MACSE |

## Tips

- MAFFT L-INS-i is the most accurate general-purpose MSA method for datasets under 200 sequences; default to it unless dataset size prohibits
- Always align proteins rather than DNA when comparing sequences below ~70% nucleotide identity since protein alignment is 3x more sensitive for detecting distant homology
- For selection analysis (dN/dS), never align coding DNA directly; align proteins first and back-translate with PAL2NAL to preserve reading frame
- Guide tree effects are real: the same sequences aligned with different guide trees produce different alignments. For critical analyses, use MUSCLE5 ensemble or GUIDANCE2 to quantify uncertainty
- Sequences below ~20% protein identity are in the twilight zone where sequence alignment becomes unreliable; consider structural alignment (TM-align) or profile-profile methods (HHpred) instead
- After alignment, always check gap distribution before phylogenetic analysis. Columns with >50% gaps may indicate misalignment or inclusion of non-homologous sequences

## Related Skills

- pairwise-alignment - Compare two sequences using PairwiseAligner
- alignment-io - Read/write MSA files in various formats
- msa-parsing - Parse, filter, trim, and assess alignment quality
- msa-statistics - Calculate identity, conservation, entropy metrics
- phylogenetics/tree-building - Build phylogenetic trees from MSAs
