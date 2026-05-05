# alignment

## Overview

Sequence-to-sequence alignment: pairwise (Bio.Align, parasail, edlib, WFA), multiple sequence alignment (MAFFT, MUSCLE5, ClustalOmega, T-Coffee, BAli-Phy joint co-estimation), structural alignment (Foldseek, Foldseek-Multimer, TM-align, US-align, DALI, Foldmason), alignment I/O, post-MSA filtering / trimming / column-mapping (ClipKIT, trimAl, BMGE, PhyIN), and alignment analysis (conservation, identity, coevolution).

**Tool type:** mixed | **Primary tools:** Bio.Align, MAFFT, MUSCLE5, Foldseek, ClipKIT

## Skills

| Skill | Description |
|-------|-------------|
| pairwise-alignment | Global/local alignment using PairwiseAligner; library selection (parasail/edlib/WFA); Karlin-Altschul significance |
| multiple-alignment | Run MSA tools (MAFFT, MUSCLE5, ClustalOmega, T-Coffee) with algorithm-selection guidance and codon-aware modes |
| structural-alignment | Backbone-aware alignment with Foldseek 3Di, Foldseek-Multimer (complex search), TM-align, US-align, DALI Z-score, Foldmason; the dark-proteome alternative when sequence MSA fails |
| alignment-trimming | Post-MSA column trimming with ClipKIT, trimAl, BMGE, Divvier, HMMcleaner, PhyIN; goal-driven mode selection |
| msa-parsing | Parse and analyze MSA content: gaps, conservation, Henikoff weighting, Neff, MI-APC coevolution |
| msa-statistics | Calculate identity, Capra-Singh JSD conservation, KL information content, distance corrections |
| alignment-io | Read, write, convert MSA files (Clustal, PHYLIP, Stockholm, A2M/A3M, MAF) and stream Pfam-scale databases |

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
- "These proteins are 12% identical; align via predicted structures instead of sequence"
- "Search this AlphaFold model against AlphaFoldDB to find structural homologs"
- "Search this antibody-antigen complex against AFDB-Multimer with Foldseek-Multimer"
- "Build a structural MSA from these PDB structures"
- "Rank twilight-fold hits with DALI Z-score instead of TM-score"
- "Trim this MSA for IQ-TREE input using ClipKIT kpic-smart-gap"
- "Compare ClipKIT and trimAl on the same alignment and tell me which to use"
- "Apply PhyIN as a second-pass trimmer to flag phylogenetically incompatible columns"
- "Detect coevolving residue pairs with mutual-information APC correction"
- "Compute Henikoff sequence weights and effective sequence number (Neff)"
- "Co-estimate alignment and tree with BAli-Phy for these 50 sequences"

## Requirements

```bash
pip install biopython numpy pandas pyhmmer

conda install -c bioconda mafft muscle clustalo t-coffee pal2nal
conda install -c bioconda clipkit trimal
conda install -c bioconda foldseek tmalign usalign foldmason

# Foldseek-Multimer is included in foldseek 9.0+ (sub-command `easy-multimersearch` / `easy-multimercluster`); verify with `foldseek --help` if installed version is older.
# DALI server / DaliLite v5: http://ekhidna2.biocenter.helsinki.fi/dali/
# PhyIN (alignment-trimming, second-pass): https://github.com/wmaddisn/PhyIN
# BAli-Phy v3 (joint MSA+tree co-estimation): http://bali-phy.org/

# Optional: codon-aware tools
# MACSE: https://bioweb.supagro.inra.fr/macse/
# PRANK: https://ariloytynoja.github.io/prank-msa/
# GUIDANCE2: http://guidance.tau.ac.il/

# Optional: high-performance pairwise / database search
pip install parasail edlib pywfa mappy
conda install -c bioconda mmseqs2  # MMseqs2 (CPU); MMseqs2-GPU requires NVIDIA GPU + separate build

# Optional: pLM aligners (for dark proteome). fair-esm is archived (use ESM-2 weights only); for ESM3+, use `pip install esm` from EvolutionaryScale.
pip install fair-esm
# vcMSA:  https://github.com/clairemcwhite/vcmsa
# DEDAL:  https://github.com/google-research/google-research/tree/master/dedal
# TM-Vec: https://github.com/valentynbez/tmvec (active fork; original at tymor22/tm-vec is in limited maintenance)
```

## Related Skills

- **alignment-files** - Process SAM/BAM/CRAM read-to-reference alignments (post-mapping). For aligning reads to a reference in the first place, see **read-alignment**.
- **read-alignment** - Mapping short and long sequencing reads to a reference genome. Distinct from this category which handles sequence-to-sequence alignment.
- **phylogenetics** - Build phylogenetic trees from MSAs (uses trimmed output from alignment-trimming)
- **sequence-io** - Read input sequences for alignment
- **sequence-manipulation** - Work with individual sequences
- **rna-structure** - Secondary-structure-aware alignment (R-Coffee, Infernal cmalign)
