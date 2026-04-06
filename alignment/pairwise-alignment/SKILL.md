---
name: bio-alignment-pairwise
description: Perform pairwise sequence alignment using Biopython Bio.Align.PairwiseAligner. Use when comparing two sequences, finding optimal alignments, scoring similarity, and identifying local or global matches between DNA, RNA, or protein sequences.
tool_type: python
primary_tool: Bio.Align
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pairwise Sequence Alignment

**"Align two sequences"** → Compute an optimal alignment between a pair of sequences using dynamic programming.
- Python: `PairwiseAligner()` (BioPython Bio.Align)
- CLI: `needle` (global) or `water` (local) from EMBOSS
- R: `pairwiseAlignment()` (Biostrings)

Align two sequences using dynamic programming algorithms (Needleman-Wunsch for global, Smith-Waterman for local).

## Required Import

**Goal:** Load modules needed for pairwise alignment operations.

**Approach:** Import the PairwiseAligner class along with sequence and I/O utilities from Biopython.

```python
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio import SeqIO
```

## Core Concepts

| Mode | Algorithm | Use Case |
|------|-----------|----------|
| `global` | Needleman-Wunsch | Full-length alignment, similar-length sequences |
| `local` | Smith-Waterman | Find best matching regions, different-length sequences |
| `global` + free end gaps | Semi-global | Overlap detection, fragment-to-reference alignment |

### Choosing the Right Mode

- **Global**: Both sequences are expected to be homologous over their full length (e.g., two orthologs of similar size). Forces end-to-end alignment.
- **Local**: Conserved domains or motifs within otherwise dissimilar sequences. BLAST uses local alignment internally. Preferred when protein termini are highly divergent (termini accumulate mutations faster than core regions).
- **Semi-global**: One sequence is a fragment or subsequence of the other (e.g., primer to template, read to reference, detecting overlap between shotgun reads). Free end gaps prevent penalizing unaligned flanking regions.

**Common mistake**: Global alignment of sequences with very different lengths forces biologically meaningless terminal gaps. If sequences differ substantially in length, use local or semi-global instead.

### DNA vs Protein Alignment

| Scenario | Align As | Rationale |
|----------|----------|-----------|
| Nucleotide identity >70% | DNA | Sufficient signal at nucleotide level |
| Nucleotide identity <70% | Protein | Codon degeneracy masks signal at DNA level; protein alignment is ~3x more sensitive |
| Noncoding sequences (UTRs, intergenic) | DNA | No protein translation possible |
| Coding sequences for dN/dS analysis | Protein first, then back-translate codons (PAL2NAL) | Preserves reading frame for selection analysis |

When in doubt, align at the protein level. It captures functional constraint better because 20 amino acids provide richer signal than 4 nucleotides.

## Creating an Aligner

**Goal:** Configure a PairwiseAligner with appropriate scoring for the sequence type.

**Approach:** Instantiate PairwiseAligner with mode, scoring parameters, or a substitution matrix depending on DNA vs protein input.

```python
# Basic aligner with defaults
aligner = PairwiseAligner()

# Configure mode and scoring
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)

# For protein alignment with substitution matrix
from Bio.Align import substitution_matrices
aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'))
```

## Performing Alignments

**"Align two sequences"** → Compute optimal alignment(s) between a pair of sequences, returning alignment objects or a score.

**Goal:** Align two sequences and retrieve the optimal alignment(s) or score.

**Approach:** Call `aligner.align()` for full alignment objects or `aligner.score()` for score-only (faster for large sequences).

```python
seq1 = Seq('ACCGGTAACGTAG')
seq2 = Seq('ACCGTTAACGAAG')

# Get all optimal alignments
alignments = aligner.align(seq1, seq2)
print(f'Found {len(alignments)} optimal alignments')
print(alignments[0])  # Print first alignment

# Get score only (faster for large sequences)
score = aligner.score(seq1, seq2)
```

## Alignment Output Format

```
target            0 ACCGGTAACGTAG 13
                  0 |||||.||||.|| 13
query             0 ACCGTTAACGAAG 13
```

## Accessing Alignment Data

**Goal:** Extract alignment properties including score, shape, aligned sequences, and coordinate mappings.

**Approach:** Access alignment object attributes and indexing to retrieve per-sequence aligned strings and coordinate arrays.

```python
alignment = alignments[0]

# Basic properties
print(alignment.score)                    # Alignment score
print(alignment.shape)                    # (num_seqs, alignment_length)
print(len(alignment))                     # Alignment length

# Get aligned sequences with gaps
target_aligned = alignment[0, :]          # First sequence (target) with gaps
query_aligned = alignment[1, :]           # Second sequence (query) with gaps

# Get coordinate mapping
print(alignment.aligned)                  # Array of aligned segment coordinates
print(alignment.coordinates)              # Full coordinate array
```

## Alignment Counts (Identities, Mismatches, Gaps)

**Goal:** Quantify identities, mismatches, and gaps in an alignment to calculate percent identity.

**Approach:** Use the `.counts()` method on the alignment object and derive percent identity from identity and mismatch totals.

```python
alignment = alignments[0]
counts = alignment.counts()

print(f'Identities: {counts.identities}')
print(f'Mismatches: {counts.mismatches}')
print(f'Gaps: {counts.gaps}')

# Calculate percent identity
total_aligned = counts.identities + counts.mismatches
percent_identity = counts.identities / total_aligned * 100
print(f'Percent identity: {percent_identity:.1f}%')
```

## Common Scoring Configurations

### DNA/RNA Alignment
```python
aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
```

### Protein Alignment
```python
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load('BLOSUM62')
aligner = PairwiseAligner(mode='global', substitution_matrix=blosum62, open_gap_score=-11, extend_gap_score=-1)
```

### Local Alignment (Find Best Region)
```python
aligner = PairwiseAligner(mode='local', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
```

### Semiglobal (Overlap/Fragment Alignment)
```python
# Free end gaps on query -- for aligning a fragment against a full-length reference
# or detecting overlap between reads
aligner = PairwiseAligner(mode='global')
aligner.query_left_open_gap_score = 0
aligner.query_left_extend_gap_score = 0
aligner.query_right_open_gap_score = 0
aligner.query_right_extend_gap_score = 0

# Free end gaps on BOTH sequences -- for overlap detection between two reads
aligner = PairwiseAligner(mode='global')
aligner.end_gap_score = 0.0
```

## Substitution Matrix Selection

**Goal:** Select the appropriate substitution matrix based on expected sequence divergence.

**Approach:** Match matrix to divergence level. BLOSUM and PAM number in **opposite directions**: higher BLOSUM = closer sequences; higher PAM = more distant sequences.

| Divergence Level | BLOSUM | PAM | When To Use |
|-----------------|--------|-----|-------------|
| Very close (<20% divergence) | BLOSUM80, BLOSUM90 | PAM30 | Recently duplicated genes, strain comparison |
| Moderate | BLOSUM62 (default) | PAM120 | General-purpose, most analyses |
| Distant (>50% divergence) | BLOSUM45, BLOSUM50 | PAM250 | Remote homology detection |

**BLOSUM62 is the universal default** (used by BLAST, most alignment tools). When in doubt, use BLOSUM62. Switch to BLOSUM80 for very similar proteins or BLOSUM45 for distant homologs.

**DNA matrices**: `NUC.4.4` (match=+5, mismatch=-4) handles IUPAC ambiguity codes. `HOXD70` is tuned for human-mouse whole-genome alignment from noncoding regions.

```python
from Bio.Align import substitution_matrices
print(substitution_matrices.load())  # List all 30 available matrices

blosum62 = substitution_matrices.load('BLOSUM62')  # General protein (default)
blosum80 = substitution_matrices.load('BLOSUM80')  # Close homologs
blosum45 = substitution_matrices.load('BLOSUM45')  # Distant homologs
nuc44 = substitution_matrices.load('NUC.4.4')      # DNA with IUPAC support
```

### Affine Gap Penalties: Biological Rationale

Gap penalties control how gaps (insertions/deletions) are scored. The **affine model** (`penalty = open + extend * (L-1)`) is almost always preferred over linear because it reflects indel biology: a DNA break introduces the first gap (costly), but extending an existing gap is mechanistically easier (less costly). This models the observation that indels in real sequences tend to occur as single contiguous events.

Typical values with BLOSUM62: gap open = -11, gap extend = -1 (BLASTP defaults). Setting gap open equal to gap extend (linear model) over-penalizes long indels and under-penalizes scattered single-residue gaps, producing biologically unrealistic alignments.

## Working with SeqRecord Objects

**Goal:** Align sequences loaded from FASTA files rather than hardcoded strings.

**Approach:** Parse SeqRecord objects from a FASTA file and pass their `.seq` attributes to the aligner.

```python
from Bio import SeqIO

records = list(SeqIO.parse('sequences.fasta', 'fasta'))
seq1, seq2 = records[0].seq, records[1].seq

aligner = PairwiseAligner(mode='global', match_score=1, mismatch_score=-1)
alignments = aligner.align(seq1, seq2)
```

## Iterating Over Multiple Alignments

```python
# Limit number of alignments returned (memory efficient)
aligner.max_alignments = 100

for i, alignment in enumerate(alignments):
    print(f'Alignment {i+1}: score={alignment.score}')
    if i >= 4:
        break
```

## Substitution Matrix from Alignment

**Goal:** Extract observed substitution frequencies from a completed alignment.

**Approach:** Access the `.substitutions` property to get a matrix of observed base/residue substitution counts.

```python
alignment = alignments[0]
substitutions = alignment.substitutions

# View as array (rows=target, cols=query)
print(substitutions)

# Access specific substitution counts
# substitutions['A', 'T'] gives count of A aligned to T
```

## Export Alignment to Different Formats

**Goal:** Convert an alignment to standard bioinformatics file formats for downstream tools.

**Approach:** Use Python's `format()` function with format specifiers (fasta, clustal, psl, sam) on the alignment object.

```python
alignment = alignments[0]

# Various output formats
print(format(alignment, 'fasta'))     # FASTA format
print(format(alignment, 'clustal'))   # Clustal format
print(format(alignment, 'psl'))       # PSL format (BLAT)
print(format(alignment, 'sam'))       # SAM format
```

## Quick Reference: Scoring Parameters

| Parameter | Description | Typical DNA | Typical Protein |
|-----------|-------------|-------------|-----------------|
| `match_score` | Score for identical bases | 1-2 | Use matrix |
| `mismatch_score` | Penalty for mismatches | -1 to -3 | Use matrix |
| `open_gap_score` | Cost to start a gap | -5 to -15 | -10 to -12 |
| `extend_gap_score` | Cost per gap extension | -0.5 to -2 | -0.5 to -1 |
| `substitution_matrix` | Scoring matrix | N/A | BLOSUM62 |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `OverflowError` | Too many optimal alignments | Set `aligner.max_alignments` |
| Low scores | Wrong scoring scheme | Use substitution matrix for proteins |
| No alignments in local mode | Scores all negative | Ensure `match_score` > 0 |

## Percent Identity: Definitions Matter

There are four common ways to calculate percent identity from the same alignment, producing different values:

| Method | Denominator | Best For |
|--------|-------------|----------|
| PID1 | Aligned positions + internal gaps | Gap-aware, conservative |
| PID2 | Aligned residue pairs (excluding gaps) | Always highest value |
| PID3 | Shorter sequence length | Length-normalized |
| PID4 | Mean sequence length | Best correlation with structural similarity |

**Practical impact**: Up to 11.5% difference between methods on a single alignment. Combined with different alignment algorithms, variation reaches 22%. Always report which method was used. The `counts()` method above uses aligned non-gap positions (similar to PID2).

## Interpreting Alignment Significance

Raw alignment scores are not directly interpretable across different scoring schemes. For database searches (BLAST), two normalized measures are used:

- **Bit score**: Normalized, database-size-independent. Scores >50 are generally reliable; <50 are suspect.
- **E-value**: Expected number of alignments with this score or better by chance. Database-size dependent. E < 1e-5 is a common significance threshold; E < 1e-50 is near-certain homology. E-values are NOT p-values (they can exceed 1).

**Compositional bias**: Low-complexity regions (poly-Q, proline-rich) and transmembrane segments inflate alignment scores. Use SEG (protein) or DUST (DNA) filtering to mask these regions before significance assessment.

## When Alignment Is NOT Appropriate

- **Twilight zone** (20-35% protein identity): Alignment reliability drops sharply. Below ~25% identity, >90% of detected sequence pairs are NOT structurally similar. Alignment length matters: 30% over 200 residues is more meaningful than 30% over 50.
- **Below 20% identity**: Sequence signal is lost in noise. Use profile-profile methods (HHpred), protein language model embeddings (ESM-2), or structural alignment (TM-align) instead.
- **Non-homologous sequences**: Alignment algorithms always produce an alignment, even for unrelated sequences. Statistical significance (E-value/bit score) is essential to distinguish true from false homology.
- **Highly repetitive sequences**: Tandem repeats cause alignment ambiguity and can produce artifactually high scores.

## Related Skills

- multiple-alignment - Align three or more sequences with MAFFT, MUSCLE5, ClustalOmega
- alignment-io - Save alignments to files in various formats
- msa-parsing - Work with multiple sequence alignments
- msa-statistics - Calculate identity, similarity metrics
- sequence-manipulation/motif-search - Pattern matching in sequences
