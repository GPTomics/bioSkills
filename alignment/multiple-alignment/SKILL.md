---
name: bio-alignment-multiple
description: Perform multiple sequence alignment using MAFFT, MUSCLE5, ClustalOmega, or T-Coffee. Guides tool and algorithm selection based on dataset size, sequence divergence, and downstream application. Use when aligning three or more homologous sequences for phylogenetics, conservation analysis, or evolutionary studies.
tool_type: mixed
primary_tool: MAFFT
---

## Version Compatibility

Reference examples tested with: MAFFT 7.520+, MUSCLE 5.1+, ClustalOmega 1.2.4+, T-Coffee 13+, PAL2NAL 14+, BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `mafft --version`, `muscle -version`, `clustalo --version`
- Python: `pip show biopython` then `help(module.function)` to check signatures

If code throws errors, introspect the installed tool and adapt the example to match the actual CLI flags rather than retrying.

# Multiple Sequence Alignment

**"Align multiple sequences"** -> Compute an optimal alignment of three or more homologous sequences using progressive, iterative, or consistency-based methods.
- CLI: `mafft` (most versatile), `muscle` (highest accuracy), `clustalo` (scales well), `t_coffee` (consistency-based)
- Python: `subprocess.run()` wrapping CLI tools; BioPython `Bio.Align.Applications` is deprecated

## Tool Selection

Choosing the right MSA tool is a consequential decision. Each tool has distinct strengths; the optimal choice depends on dataset size, divergence level, and whether accuracy or speed is prioritized.

| Tool | Best For | Max Sequences | Accuracy | Speed |
|------|----------|---------------|----------|-------|
| MAFFT L-INS-i | Highest accuracy, <200 seqs | ~200 | Highest | Slow |
| MAFFT FFT-NS-2 | Large datasets, good balance | ~50,000 | Good | Fast |
| MAFFT E-INS-i | Sequences with long unalignable internal regions | ~200 | High | Slow |
| MUSCLE5 | Benchmarked highest accuracy, uncertainty estimation | ~10,000 | Highest | Medium |
| ClustalOmega | Very large datasets, HMM-based profiles | ~100,000+ | Good | Fast |
| T-Coffee | Small datasets needing maximum accuracy, structural information | ~200 | Highest | Slowest |

**Default recommendation**: MAFFT L-INS-i for <200 sequences; MAFFT FFT-NS-2 or ClustalOmega for >1000 sequences; MUSCLE5 when alignment confidence estimates are needed.

## Critical Concepts

### Progressive Alignment Limitation

All major MSA tools use progressive alignment as their foundation: compute pairwise distances, build a guide tree, align sequences following the tree. The fundamental limitation is **once a gap is inserted, it is never removed** in the progressive phase. Early errors propagate. This is why iterative refinement (MAFFT -i modes, MUSCLE5) and consistency scoring (T-Coffee) exist; they mitigate but do not eliminate guide tree dependency.

### Guide Tree Dependency

The guide tree strongly influences the final alignment. Different guide tree topologies produce different alignments, and many gaps in progressive alignments reflect guide tree structure rather than true evolutionary events. Mitigations: MAFFT's iterative methods recompute guide trees; MUSCLE5 perturbs guide trees across its ensemble; GUIDANCE2 bootstraps over guide trees to quantify uncertainty.

### Sequence Divergence Thresholds

| Protein Identity | Signal Level | Recommendation |
|-----------------|--------------|----------------|
| >40% | Strong | Any MSA tool produces reliable alignment |
| 25-40% | Moderate (twilight zone begins) | Use iterative methods (L-INS-i, MUSCLE5); validate with GUIDANCE2 |
| 20-25% | Weak | Profile-profile methods (HHpred); consider structural alignment |
| <20% | Noise dominates signal | Sequence MSA is unreliable; use structural alignment (TM-align, DALI) or AlphaFold-based methods |

## Running MAFFT

### Algorithm Selection

MAFFT offers multiple algorithms with explicit accuracy/speed tradeoffs. Selecting the right mode is critical; the difference between L-INS-i and FFT-NS-1 can be the difference between a correct and incorrect downstream phylogeny.

| Algorithm | Flag | Strategy | Best For |
|-----------|------|----------|----------|
| FFT-NS-1 | `--retree 1` | Progressive only | Quick look, >10,000 seqs |
| FFT-NS-2 | `--retree 2` | Progressive + guide tree rebuild | Default balance, 200-10,000 seqs |
| FFT-NS-i | `--maxiterate 1000` | Iterative refinement | Moderate improvement, 200-2,000 seqs |
| G-INS-i | `--globalpair --maxiterate 1000` | Global pairwise + iterative | Sequences alignable over full length, <200 |
| L-INS-i | `--localpair --maxiterate 1000` | Local pairwise + iterative | Single alignable domain amid divergent flanks, <200 |
| E-INS-i | `--genafpair --maxiterate 1000` | Local with generalized affine gaps | Multiple conserved motifs separated by unalignable regions, <200 |
| Auto | `--auto` | Auto-selects based on dataset size | When unsure |

**Decision guide**: If sequences share a single conserved domain (most common case), use L-INS-i. If sequences are globally similar (e.g., ortholog set of similar length), use G-INS-i. If sequences have multiple conserved blocks separated by highly variable linker regions (e.g., multi-domain proteins with variable interdomain regions), use E-INS-i.

### Basic Usage

**Goal:** Run MAFFT on a FASTA file with appropriate algorithm selection.

**Approach:** Invoke MAFFT via command line or subprocess, selecting the algorithm based on dataset characteristics.

```bash
# Highest accuracy for <200 sequences (local pairwise iterative)
mafft --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Good balance for medium datasets
mafft --retree 2 input.fasta > aligned.fasta

# Auto-select algorithm based on dataset size
mafft --auto input.fasta > aligned.fasta

# Protein alignment with specific matrix (default BLOSUM62)
mafft --amino --localpair --maxiterate 1000 input.fasta > aligned.fasta

# DNA alignment (auto-detected, but can be explicit)
mafft --nuc --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Adjust gap penalties (op=gap open, ep=gap extension)
mafft --op 1.53 --ep 0.123 --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Multithreaded
mafft --thread 8 --localpair --maxiterate 1000 input.fasta > aligned.fasta
```

```python
import subprocess

def run_mafft(input_fasta, output_fasta, algorithm='linsi', threads=4):
    algo_flags = {
        'linsi': ['--localpair', '--maxiterate', '1000'],
        'ginsi': ['--globalpair', '--maxiterate', '1000'],
        'einsi': ['--genafpair', '--maxiterate', '1000'],
        'fftns2': ['--retree', '2'],
        'auto': ['--auto'],
    }
    cmd = ['mafft', '--thread', str(threads)] + algo_flags[algorithm] + [input_fasta]
    with open(output_fasta, 'w') as out:
        subprocess.run(cmd, stdout=out, check=True)

run_mafft('sequences.fasta', 'aligned.fasta', algorithm='linsi')
```

### Adding Sequences to an Existing Alignment

**Goal:** Add new sequences to an existing MSA without realigning the entire dataset.

**Approach:** Use MAFFT's `--add` flag, which aligns new sequences against the existing alignment profile.

```bash
# Add new sequences to existing alignment (keeps existing alignment fixed)
mafft --add new_seqs.fasta --keeplength existing_alignment.fasta > updated.fasta

# Add and allow existing alignment to adjust
mafft --add new_seqs.fasta existing_alignment.fasta > updated.fasta
```

## Running MUSCLE5

MUSCLE5 is a complete rewrite of MUSCLE3. It achieves the highest accuracy on standard benchmarks (59% correct columns on Balifam-10000 vs ClustalOmega 52%, MAFFT 47%) through an ensemble approach.

### Basic Usage

```bash
# Standard alignment
muscle -align input.fasta -output aligned.fasta

# Super5 algorithm for large datasets (divide-and-conquer)
muscle -super5 input.fasta -output aligned.fasta

# Specify number of threads
muscle -align input.fasta -output aligned.fasta -threads 8
```

### Ensemble Mode for Alignment Confidence

**Goal:** Quantify alignment uncertainty by generating multiple perturbed alignments and measuring consistency.

**Approach:** MUSCLE5's ensemble generates alignments with perturbed HMM parameters and guide trees. Column confidence equals the fraction of ensemble members supporting that column arrangement.

```bash
# Generate ensemble of alignments (default 4 replicates)
muscle -align input.fasta -output aligned.efa -perturb 4

# Disperse mode: output each replicate as separate file
muscle -align input.fasta -output rep -perturb 4 -disperse
```

This ensemble output enables downstream confidence assessment. Columns consistently aligned across replicates are reliable; columns that vary are uncertain. This directly addresses the problem of alignment uncertainty propagating to phylogenetic inference.

## Running ClustalOmega

ClustalOmega uses HMM-HMM alignment (via HHalign) for profile comparison, achieving better accuracy than ClustalW's profile-profile approach. The mBed algorithm enables O(N log N) guide tree construction, scaling to hundreds of thousands of sequences.

```bash
# Basic alignment
clustalo -i input.fasta -o aligned.fasta --auto

# Force overwrite output
clustalo -i input.fasta -o aligned.fasta --force

# Specify output format
clustalo -i input.fasta -o aligned.phy --outfmt=phylip

# Use more iterations for better accuracy
clustalo -i input.fasta -o aligned.fasta --iter=5

# Profile-profile alignment (align two existing MSAs)
clustalo --p1 profile1.fasta --p2 profile2.fasta -o merged.fasta

# Add sequences to existing alignment
clustalo -i new_seqs.fasta --profile1 existing.fasta -o updated.fasta

# Multithreaded
clustalo -i input.fasta -o aligned.fasta --threads=8
```

## Running T-Coffee

T-Coffee's consistency-based approach generates a library of pairwise alignments from multiple methods, then evaluates how consistently residue pairs align across all comparisons. Slower but integrates diverse information sources.

```bash
# Standard alignment (generates library from multiple pairwise methods)
t_coffee input.fasta -output fasta_aln -outfile aligned.fasta

# With structural information (dramatic accuracy improvement when structures available)
t_coffee input.fasta -template_file template.pdb -mode expresso

# Meta-mode: combine multiple MSA tools' outputs
t_coffee input.fasta -mode mcoffee

# RNA alignment with secondary structure
t_coffee input.fasta -mode rcoffee
```

**When to use T-Coffee**: Small datasets (<50 sequences) where maximum accuracy matters, especially when structural information (PDB templates) is available. Expresso/3D-Coffee mode provides ~15-49% improvement in correct columns when structures exist.

## Codon-Aware Alignment

### When Codon Alignment Is Required

Coding sequences destined for selection analysis (dN/dS with PAML/codeml, HyPhy) **must** be aligned respecting codon boundaries. Standard nucleotide MSA tools do not preserve reading frames and will produce incorrect dN/dS estimates.

### PAL2NAL: Protein-Guided Codon Alignment

**Goal:** Thread a nucleotide coding sequence alignment onto a protein alignment to preserve reading frame.

**Approach:** Align protein sequences first (higher sensitivity), then use PAL2NAL to map the protein alignment back to codons. This is the standard approach for preparing codeml input.

```bash
# Step 1: Align at protein level (more sensitive for divergent sequences)
mafft --localpair --maxiterate 1000 proteins.fasta > proteins_aligned.fasta

# Step 2: Thread DNA onto protein alignment
pal2nal.pl proteins_aligned.fasta codons.fasta -output fasta > codons_aligned.fasta

# Output as PAML format for codeml
pal2nal.pl proteins_aligned.fasta codons.fasta -output paml > codons_aligned.phy
```

### MACSE: Direct Codon-Aware Alignment

**Goal:** Align coding sequences directly at the nucleotide level while respecting codon structure.

**Approach:** MACSE scores based on amino acid translation while aligning DNA. Unlike PAL2NAL, it handles frameshifts and internal stop codons, making it suitable for pseudogenes and error-prone sequences.

```bash
# Align coding sequences (handles frameshifts)
java -jar macse_v2.jar -prog alignSequences -seq coding_seqs.fasta -out_NT aligned_nt.fasta -out_AA aligned_aa.fasta

# Add sequences to existing alignment
java -jar macse_v2.jar -prog enrichAlignment -align existing.fasta -seq new_seqs.fasta -out_NT updated.fasta
```

**PAL2NAL vs MACSE**: Use PAL2NAL when sequences are clean (no frameshifts, no premature stops) since it is simpler and faster. Use MACSE when working with pseudogenes, draft assemblies, or sequences that may contain sequencing errors causing frameshifts.

## Post-Alignment Validation Checklist

Before proceeding to downstream analysis, verify alignment quality:

1. **Visual inspection**: Scan for columns of mostly gaps with scattered residues (hallmark of misalignment)
2. **Gap distribution**: High gap fraction (>50% of columns with gaps) suggests problematic regions or inclusion of non-homologous sequences
3. **Sequence identity**: If average pairwise identity is <25% for proteins, alignment reliability is questionable
4. **Outlier sequences**: Sequences with excessive gaps relative to others may be non-homologous or fragments; consider removing and re-aligning
5. **Conservation pattern**: Functional domains should show clear conservation; absence of expected conserved motifs suggests alignment error or non-homology
6. **Run GUIDANCE2 or MUSCLE5 ensemble**: Quantify alignment confidence per column before phylogenetic inference

## When NOT to Run MSA

- **Non-homologous sequences**: MSA tools always produce an alignment, even for unrelated sequences; verify homology first (e.g., BLAST E-value < 1e-5)
- **Sequences below the twilight zone**: Below ~20% protein identity, sequence signal is lost in noise; structural alignment is needed
- **Different domain architectures**: Globally aligning multi-domain proteins with different domain orders produces meaningless results; align individual domains separately
- **Very different lengths without shared homology**: Aligning a 50-residue fragment against 1000-residue proteins globally forces biologically meaningless gaps; use local alignment or fragment-aware modes (E-INS-i)
- **Highly repetitive sequences**: Tandem repeats cause alignment ambiguity; specialized tools (e.g., TRUST for tandem repeats) may be needed

## Quick Reference

| Task | Command |
|------|---------|
| Best accuracy (<200 seqs) | `mafft --localpair --maxiterate 1000 in.fa > out.fa` |
| Large dataset | `mafft --retree 2 in.fa > out.fa` or `clustalo -i in.fa -o out.fa` |
| Uncertainty estimation | `muscle -align in.fa -output out.efa -perturb 4` |
| Codon-aware | Align protein first, then `pal2nal.pl prot.fa cds.fa -output fasta` |
| Add to existing MSA | `mafft --add new.fa existing.fa > updated.fa` |
| Profile merge | `clustalo --p1 msa1.fa --p2 msa2.fa -o merged.fa` |
| With structure info | `t_coffee in.fa -mode expresso` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| MAFFT runs out of memory | L-INS-i on too many sequences | Switch to FFT-NS-2 or auto mode |
| Alignment has many all-gap columns | Non-homologous sequences included | Filter input by BLAST hits first |
| ClustalOmega crashes on large input | Memory limits | Increase `--MAC` RAM or use MAFFT |
| PAL2NAL "length mismatch" | Protein/DNA sequences not from same genes | Verify correspondence with sequence IDs |
| MUSCLE5 slow on large datasets | Default algorithm not designed for >10K seqs | Use `-super5` mode |
| Poor alignment despite high identity | Wrong sequence type detection | Explicitly specify `--amino` or `--nuc` |

## Related Skills

- pairwise-alignment - Compare two sequences using PairwiseAligner
- alignment-io - Read/write MSA files in various formats
- msa-parsing - Parse, filter, trim, and assess MSA quality
- msa-statistics - Calculate identity, conservation, entropy metrics
- phylogenetics/tree-building - Build phylogenetic trees from MSAs
- sequence-io/read-sequences - Read input sequences for alignment
