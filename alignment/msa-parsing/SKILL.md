---
name: bio-alignment-msa-parsing
description: Parse and analyze multiple sequence alignments using Biopython. Extract sequences, identify conserved regions, analyze gaps, work with annotations, and manipulate alignment data for downstream analysis. Use when parsing or manipulating multiple sequence alignments.
tool_type: python
primary_tool: Bio.AlignIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# MSA Parsing and Analysis

Parse multiple sequence alignments to extract information, analyze content, and prepare for downstream analysis.

## Required Import

**Goal:** Load modules for parsing, analyzing, and manipulating multiple sequence alignments.

**Approach:** Import AlignIO for reading, Counter for column analysis, and alignment classes for constructing modified alignments.

```python
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter
import numpy as np
import pandas as pd
```

Optional for streaming and Easel-based weighting:
```python
import pyhmmer
```

## Loading Alignments

**Goal:** Read an MSA file and inspect its dimensions.

**Approach:** Use `AlignIO.read()` specifying the file and format.

```python
from Bio import AlignIO

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'{len(alignment)} sequences, {alignment.get_alignment_length()} columns')
```

## Extracting Sequence Information

### Get All Sequence IDs
```python
seq_ids = [record.id for record in alignment]
```

### Get Sequences as Strings
```python
sequences = [str(record.seq) for record in alignment]
```

### Get Sequence by ID
```python
def get_sequence_by_id(alignment, seq_id):
    for record in alignment:
        if record.id == seq_id:
            return record
    return None

target = get_sequence_by_id(alignment, 'species_A')
```

### Access Descriptions and Annotations
```python
for record in alignment:
    print(f'ID: {record.id}')
    print(f'Description: {record.description}')
    print(f'Annotations: {record.annotations}')
```

## Column-wise Analysis

**Goal:** Analyze alignment content column by column to assess composition, conservation, and variability.

**Approach:** Use column indexing (`alignment[:, idx]`) and Counter to examine character frequencies at each position.

### Get Single Column
```python
column_5 = alignment[:, 5]  # Returns string of characters at position 5
print(column_5)  # e.g., 'AAAGA'
```

**API note:** `Bio.AlignIO` returns `MultipleSeqAlignment` objects whose `[:, idx]` returns a plain `str`; `[:, start:end]` returns another `MultipleSeqAlignment`. The newer `Bio.Align.Alignment` (from `Align.read` / `Align.parse`) uses numpy-backed slicing -- verify with `type(alignment[:, 0])` before assuming string methods work. For numpy-array access to the full alignment, use `np.array(alignment)`.

### Iterate and Count Columns
```python
for col_idx in range(alignment.get_alignment_length()):
    column = alignment[:, col_idx]
    counts = Counter(column)
```

### Find Conserved Positions
```python
def find_conserved_positions(alignment, threshold=1.0):
    conserved = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        counts = Counter(column)
        most_common_char, most_common_count = counts.most_common(1)[0]
        if most_common_char != '-':
            conservation = most_common_count / len(alignment)
            if conservation >= threshold:
                conserved.append((col_idx, most_common_char))
    return conserved

fully_conserved = find_conserved_positions(alignment, threshold=1.0)
mostly_conserved = find_conserved_positions(alignment, threshold=0.8)
```

## Gap Analysis

**Goal:** Quantify gap distribution across sequences and columns to identify problematic regions or sequences.

**Approach:** Count gap characters per sequence and per column, then identify positions exceeding a gap fraction threshold.

### Count Gaps Per Sequence
```python
gap_counts = [(record.id, str(record.seq).count('-')) for record in alignment]
for seq_id, gaps in gap_counts:
    print(f'{seq_id}: {gaps} gaps')
```

### Count Gaps Per Column
```python
def gaps_per_column(alignment):
    return [alignment[:, i].count('-') for i in range(alignment.get_alignment_length())]

gap_profile = gaps_per_column(alignment)
```

### Find Gappy Columns
```python
def find_gappy_columns(alignment, threshold=0.5):
    gappy = []
    num_seqs = len(alignment)
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        gap_fraction = column.count('-') / num_seqs
        if gap_fraction >= threshold:
            gappy.append(col_idx)
    return gappy

columns_to_remove = find_gappy_columns(alignment, threshold=0.5)
```

### Remove Gappy Columns
```python
def remove_gappy_columns(alignment, threshold=0.5):
    num_seqs = len(alignment)
    keep_columns = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        gap_fraction = column.count('-') / num_seqs
        if gap_fraction < threshold:
            keep_columns.append(col_idx)

    new_records = []
    for record in alignment:
        new_seq = ''.join(str(record.seq)[i] for i in keep_columns)
        new_records.append(SeqRecord(Seq(new_seq), id=record.id, description=record.description))
    return MultipleSeqAlignment(new_records)

cleaned = remove_gappy_columns(alignment, threshold=0.5)
```

## Alignment Trimming

Trimming controversy and tool selection (ClipKIT, trimAl, BMGE, Divvier, HMMcleaner, Noisy) is the subject of a dedicated skill. Use this short decision matrix for routing:

| Goal | First-line tool |
|------|-----------------|
| Phylogenetic-tree input | ClipKIT `kpic-smart-gap` (Steenwyk et al 2020 PLOS Bio) |
| HMM profile building | trimAl `-gappyout` (Capella-Gutierrez et al 2009 Bioinf) |
| Selection / dN/dS input | Avoid aggressive trimming; use TCS / GUIDANCE2 column masking |
| Deep prokaryotic phylogenomics | BMGE (Criscuolo & Gribaldo 2010 BMC Evol Biol) |
| Preserve column-mapping for residue-level analysis | trimAl `-colnumbering` |

See alignment/alignment-trimming for full mode comparisons, decision trees, and runnable examples.

## Gap Handling for Phylogenetics

How gaps are treated in downstream phylogenetic analysis significantly affects tree topology:

| Treatment | Method | Tradeoff |
|-----------|--------|----------|
| Missing data (default) | Gaps = unknown character | Most common; can be statistically inconsistent under ML |
| Fifth state | Gap = 5th nucleotide | Biologically problematic (gaps of different lengths treated equally) |
| Simple indel coding | Each unique indel coded as binary character | Most biologically realistic; adds phylogenetic signal |

For slow- to mid-rate datasets where indels are phylogenetically informative, prefer SIC indel coding or fifth-state treatment; for rapidly-evolving datasets (intra-species, ITS regions, retroelement-rich plant genomes), default to missing-data treatment because gap homology is unreliable. Run a sensitivity analysis comparing both treatments before drawing topological conclusions.

## Identifying Unreliable Alignment Regions

Columns exhibiting **both** high gap fraction AND low conservation are the strongest indicators of alignment uncertainty. These often reflect guide tree artifacts rather than true evolutionary events. Before phylogenetic analysis:

1. Flag columns with gap fraction >50%, which may be alignment artifacts
2. Check if gappy regions coincide with insertions in a single divergent sequence (remove that sequence and re-align)
3. For critical analyses, run GUIDANCE2 or MUSCLE5 ensemble to get per-column confidence scores; mask columns below the reliability threshold (default: 0.93 for GUIDANCE2)

## Consensus Sequence

**"Get consensus sequence"** → Derive a single representative sequence from an MSA based on majority-rule voting at each column.

**Goal:** Generate a consensus sequence from the alignment using a frequency threshold.

**Approach:** At each column, select the most common non-gap character if it exceeds the threshold; otherwise mark as ambiguous.

### Simple Majority Consensus
```python
def consensus_sequence(alignment, threshold=0.5, gap_char='-', ambiguous='N'):
    consensus = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        counts = Counter(column)
        most_common_char, most_common_count = counts.most_common(1)[0]
        if most_common_char == gap_char:
            counts.pop(gap_char, None)
            if counts:
                most_common_char, most_common_count = counts.most_common(1)[0]
            else:
                most_common_char = gap_char

        if most_common_count / len(alignment) >= threshold:
            consensus.append(most_common_char)
        else:
            consensus.append(ambiguous)
    return ''.join(consensus)

consensus = consensus_sequence(alignment, threshold=0.5)
```

### Note on Bio.Align.AlignInfo
The `AlignInfo.SummaryInfo` class is deprecated in recent Biopython versions. The custom `consensus_sequence()` function above is the recommended approach. When deprecation warnings appear from `AlignInfo`, callers should switch to the custom implementation.

## Extracting Regions

### Slice by Column Range
```python
region = alignment[:, 100:200]  # Columns 100-199
```

### Slice by Sequence Range
```python
subset = alignment[0:10]  # First 10 sequences
```

### Extract Ungapped Regions from Reference
```python
def extract_ungapped_regions(alignment, ref_idx=0):
    ref_seq = str(alignment[ref_idx].seq)
    ungapped_cols = [i for i, char in enumerate(ref_seq) if char != '-']

    new_records = []
    for record in alignment:
        new_seq = ''.join(str(record.seq)[i] for i in ungapped_cols)
        new_records.append(SeqRecord(Seq(new_seq), id=record.id, description=record.description))
    return MultipleSeqAlignment(new_records)

ungapped = extract_ungapped_regions(alignment, ref_idx=0)
```

## Sequence Filtering

**Goal:** Subset an alignment to retain only sequences matching specific criteria (ID pattern, gap content, uniqueness).

**Approach:** Iterate over alignment records, apply filter conditions, and reconstruct a new MultipleSeqAlignment from matching records.

```python
import re

def filter_by_id(alignment, pattern):
    regex = re.compile(pattern)
    return MultipleSeqAlignment([r for r in alignment if regex.search(r.id)])

def filter_by_gap_content(alignment, max_gap_fraction=0.1):
    return MultipleSeqAlignment(
        [r for r in alignment if str(r.seq).count('-') / len(r.seq) <= max_gap_fraction]
    )

def remove_duplicates(alignment):
    seen = set()
    return MultipleSeqAlignment([r for r in alignment if not (str(r.seq) in seen or seen.add(str(r.seq)))])
```

## Working with Annotations

Stockholm-derived alignments expose secondary-structure markup, per-column GC/GR annotations, and per-sequence metadata via `record.annotations`, `record.letter_annotations`, and `alignment.column_annotations`. See usage-guide.md "Working with Annotations" for the full access pattern.

## Position Mapping

**Goal:** Convert between alignment column coordinates and ungapped sequence coordinates.

**Approach:** For one-off lookups, walk the sequence tracking gap characters. For repeated queries on the same sequence, vectorize with `numpy.cumsum` over a gap-mask -- O(L) preprocessing, O(1) lookups.

### Vectorized Coordinate Mapping (Recommended)

```python
import numpy as np

def coordinate_map(record):
    chars = np.frombuffer(str(record.seq).encode('ascii'), dtype=np.uint8)
    is_residue = chars != ord('-')
    seq_to_aln = np.flatnonzero(is_residue)
    aln_to_seq = np.where(is_residue, np.cumsum(is_residue) - 1, -1)
    return seq_to_aln, aln_to_seq

seq_to_aln, aln_to_seq = coordinate_map(alignment[0])
column_for_residue_42 = seq_to_aln[42]
residue_at_column_100 = aln_to_seq[100]
```

`aln_to_seq[i] == -1` indicates a gap at alignment column `i`. This pattern handles 1 M-site genomic alignments in milliseconds compared to the loop-based version's seconds.

For single lookups without numpy, walk the sequence tracking a counter (`seq_pos += 1` for each non-`-` char). The vectorized form above subsumes both directions.

### Mapping Alignment Columns to PDB Residues

A column-to-PDB mapping requires THREE coordinate systems: alignment column -> SEQRES residue (ungapped FASTA) -> ATOM residue (resolved structure). The SEQRES-to-ATOM map is non-trivial because PDB structures have unmodelled loops, N-terminal tags, engineered mutations, and seleno-substitutions. Conservation scores mapped via the bare alignment-to-SEQRES path will be off-by-many residues whenever the structure has missing density. For SEQRES/ATOM extraction and the authoritative `_pdbx_poly_seq_scheme` mapping, see `structural-biology/structure-navigation`.

## Sequence Weighting and Neff

Compute sequence weights before any column-wise statistic on phylogenetically structured datasets (lots of closely-related sequences plus a few outliers); without weighting, every per-column metric is biased toward the over-represented clades.

### Henikoff Sequence Weights

**Goal:** Give each sequence a weight that reflects its non-redundant contribution.

**Approach:** Each column `c` contributes `1 / (k_c * n_{s,c})` to sequence `s`, where `k_c` is the number of distinct residues at column `c` and `n_{s,c}` is the count of `s`'s residue at that column (Henikoff & Henikoff 1994 JMB).

```python
def henikoff_weights(alignment):
    seq_array = np.array([list(str(r.seq)) for r in alignment])
    weights = np.zeros(len(alignment))
    for col_idx in range(seq_array.shape[1]):
        residues, inverse, counts = np.unique(seq_array[:, col_idx], return_inverse=True, return_counts=True)
        if '-' in residues:
            continue
        weights += 1.0 / (len(residues) * counts[inverse])
    return weights / weights.sum()
```

Full implementation: `examples/henikoff_weights.py`.

**Edge case:** The implementation skips columns containing any gap, so sequences whose non-gap residues fall ONLY in gap-rich columns receive weight zero. This typically affects fragmentary or terminal-truncated sequences. HMMER's `pb` weighting handles this by including gaps as a residue type. For Pfam-style alignments where many columns are gappy, switch to `pyhmmer.easel.MSAFile` + `compute_weights(method='pb')` rather than this implementation, or augment the pure-numpy version with an `ignore_gaps=False` branch.

### Effective Sequence Number (Neff)

**Goal:** Estimate effective non-redundant sequence count after similarity-based clustering.

**Approach:** Cluster sequences at an identity threshold (HMMER and AlphaFold use 0.62 for protein, 0.80 for nucleotide), assign each sequence weight `1 / cluster_size`, and sum. `Neff/L > 0.5` is the rule-of-thumb threshold for direct-coupling-analysis contact prediction. Full implementation: `examples/neff.py`.

The thresholds reflect convention: 0.62 traces back to BLOSUM62 derivation; 0.80 is HHsuite's default for cd-hit-style clustering of the BFD database.

**Neff is estimator-dependent.** Different tools report different "Neff" values for the same MSA, often differing by 2-3x:

| Tool | Method | Typical relationship |
|------|--------|---------------------|
| HMMER, pyhmmer (`method='pb'`) | Henikoff position-based | Baseline |
| HMMER (`method='gsc'`) | Gerstein-Sonnhammer-Chothia | ~0.5-1.5x of `pb` |
| HMMER (`method='blosum'`) | BLOSUM62-similarity clusters | ~0.7-1.2x of `pb` |
| AlphaFold2 / ColabFold | MMseqs2 cluster count at id=0.62 | ~1.5-3x of `pb` (counts unique sequences, not weights) |
| EVcouplings | HHfilter clusters at id=0.80 | ~0.5-1x of `pb` |
| HHsuite `hhmake` | inverse-purity count | Closer to `pb` |

Always report which estimator was used. Threshold rules ("Neff/L > 0.5 for DCA") are calibrated against a specific estimator -- Hopf et al 2017 use HHfilter-80 cluster count; Marks et al 2011 use position-based weighting. AlphaFold2's MSA-depth gating uses unweighted cluster count at 62% identity.

## Coevolution: Mutual Information with APC

**Goal:** Identify columns whose residue identities co-vary, indicating direct or indirect physical/functional coupling.

**Approach:** Compute pairwise mutual information across column pairs, then subtract the average-product correction (APC) from Dunn, Wahl & Gloor (2008 Bioinf) to remove per-column-entropy and phylogenetic background. This is the foundation of plmDCA and EVcouplings; full DCA needs Potts-model inference but APC-corrected MI runs in pure numpy and is informative on its own.

Skeleton (full implementation: `examples/mi_apc.py`):
```python
def mi_matrix_apc(alignment):
    # Compute pairwise MI across columns -> column_means -> APC = outer(means) / overall_mean
    # Return mi - apc
    ...
```

For production-grade contact prediction, switch to plmDCA (Ekeberg et al 2013 Phys Rev E) or EVcouplings (Hopf et al 2017 Nat Biotechnol). APC-corrected MI scales to a few hundred columns; deeper analyses need approximate likelihood methods.

**APC over-correction on short alignments:** APC subtracts each column-pair's product of column-average MIs. For alignments with <100 columns, the column-averages are noisy estimates dominated by their constituent column-pairs, and APC removes signal proportional to noise. Empirically on Pfam alignments <100 columns, APC-corrected MI underperforms raw MI for contact prediction (Cocco et al 2018 Rep Prog Phys review). Apply APC only when L > 100 and Neff/L > 1; below this, raw MI plus a phylogenetic-distance threshold is more reliable. plmDCA and EVcouplings auto-skip APC when input depth is insufficient.

## A2M / A3M Conventions

A2M (HMMER) and A3M (HHsuite, ColabFold) encode insert vs match columns via case (uppercase = match column residue, lowercase = insert). A3M does not pad inserts across sequences and must be reformatted to A2M before loading as a rectangular MSA. See `alignment/alignment-io` A2M / A3M Conventions section for the full character table, BioPython load pattern, and the `reformat.pl` reference-sequence pitfall.

## Streaming Large Alignments

For Pfam-scale streaming (multi-gigabyte Stockholm or A3M databases that exceed RAM), use `pyhmmer.easel.MSAFile` with `compute_weights(method='pb')` for in-flight Henikoff weighting. See `alignment/alignment-io` Streaming Large Stockholm Databases section for the full code pattern.

## Quick Reference: Common Operations

| Task | Code |
|------|------|
| Get column | `alignment[:, col_idx]` |
| Get sequence | `alignment[seq_idx]` |
| Column count | `alignment.get_alignment_length()` |
| Sequence count | `len(alignment)` |
| Find gaps | `str(record.seq).count('-')` |
| Consensus | Use custom `consensus_sequence()` function |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `IndexError` | Column index out of range | Check `get_alignment_length()` |
| Unequal sequence lengths | Invalid MSA | Ensure all sequences same length |
| Empty Counter | All gaps in column | Handle gap-only columns |

## Related Skills

- alignment/multiple-alignment - Run MSA tools (MAFFT, MUSCLE5, ClustalOmega) to generate alignments
- alignment/alignment-io - Read/write alignment files in various formats
- alignment/pairwise-alignment - Create pairwise alignments
- alignment/msa-statistics - Calculate conservation metrics
- alignment/alignment-trimming - ClipKIT, trimAl, BMGE, Divvier modes and decision trees
- alignment/structural-alignment - Twilight-zone alternative when sequence MSA is unreliable
- phylogenetics/modern-tree-inference - Build trees from processed alignments

## References

- Henikoff S, Henikoff JG. 1994. Position-based sequence weights. JMB 243:574-578.
- Dunn SD, Wahl LM, Gloor GB. 2008. Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction. Bioinf 24:333-340.
- Ekeberg M, Lovkvist C, Lan Y, Weigt M, Aurell E. 2013. Improved contact prediction in proteins: using pseudolikelihoods to infer Potts models. Phys Rev E 87:012707.
- Hopf TA, Ingraham JB, Poelwijk FJ, Scharfe CPI, Springer M, Sander C, Marks DS. 2017. Mutation effects predicted from sequence co-variation. Nat Biotechnol 35:128-135.
- Eddy SR. 2011. Accelerated profile HMM searches. PLOS CB 7:e1002195.
- Larralde M et al. 2023. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinf 39:btad214.
- Cocco S et al. 2018. Inverse statistical physics of protein sequences: a key issues review. Rep Prog Phys 81:032601.
- Marks DS et al. 2011. Protein 3D structure computed from evolutionary sequence variation. PLOS One 6:e28766.
- Jumper J et al. 2021. Highly accurate protein structure prediction with AlphaFold. Nature 596:583-589.
- Simmons MP, Ochoterena H. 2000. Gaps as characters in sequence-based phylogenetic analyses. Syst Biol 49:369-381.
- Mueller K. 2006. Incorporating information from length-mutational events into phylogenetic analysis. Cladistics 22:368-381.
- Dwivedi B, Gadagkar SR. 2009. Phylogenetic inference under varying proportions of indel-induced alignment gaps. BMC Evol Biol 9:211.
- Velankar S et al. 2013. SIFTS: structure integration with function, taxonomy and sequences resource. NAR 41:D483-D489.
