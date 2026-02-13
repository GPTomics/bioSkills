---
name: bio-rna-structure-secondary-structure-prediction
description: Predicts RNA secondary structures using minimum free energy folding and partition function analysis with ViennaRNA (RNAfold, RNAalifold, RNAcofold). Computes base-pair probabilities, centroid structures, and consensus structures from alignments. Use when predicting RNA folding, evaluating structural stability, or comparing structures across homologs.
tool_type: cli
primary_tool: ViennaRNA
---

## Version Compatibility

Reference examples tested with: Infernal 1.1+, matplotlib 3.8+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Secondary Structure Prediction

**"Predict the secondary structure of my RNA sequence"** → Compute minimum free energy (MFE) folding, base-pair probabilities via partition function, and consensus structures from alignments using thermodynamic models.
- CLI: `RNAfold` for single-sequence MFE/partition folding
- CLI: `RNAalifold` for consensus structure from alignment
- CLI: `RNAcofold` for RNA-RNA interaction structure

Predict RNA secondary structures using thermodynamic models. ViennaRNA provides MFE folding, partition function analysis, consensus structure prediction from alignments, and RNA-RNA interaction prediction.

## RNAfold: Single Sequence Folding

### MFE Structure

```bash
# Basic MFE folding (reads sequence from stdin or file)
echo "GGGAAACCC" | RNAfold

# With partition function (-p) and base-pair probabilities
echo "GGGAAACCC" | RNAfold -p

# Output PostScript dot plot and structure plot
echo ">myRNA" > input.fa
echo "GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA" >> input.fa
RNAfold -p --noPS < input.fa  # Suppress PostScript files
```

### Key RNAfold Options

| Option | Description |
|--------|-------------|
| `-p` | Compute partition function and base-pair probabilities |
| `--MEA` | Compute maximum expected accuracy structure |
| `-d2` | Dangling end energies on both sides of helices (default) |
| `-T 37` | Temperature in Celsius (default: 37) |
| `--noLP` | No lonely pairs (isolated base pairs) |
| `--noPS` | Suppress PostScript output files |
| `-C` | Read structure constraints from input |
| `--shape` | Incorporate SHAPE reactivity data |

### Constrained Folding

```bash
# Force specific positions paired/unpaired
# Constraint notation: '.' = unconstrained, 'x' = unpaired, '(' ')' = forced pair
echo -e ">constrained\nGGGCUAUUAGCUCAGUUGGUUAGAGCGCACC\n...xxxx.........................." | RNAfold -C
```

## RNAalifold: Consensus Structure from Alignment

Predicts a consensus structure from a multiple sequence alignment, combining thermodynamic stability with covariation evidence.

```bash
# Input: Stockholm or ClustalW alignment format
RNAalifold --aln alignment.sto

# With covariation weighting and partition function
RNAalifold --cfactor 0.6 --nfactor 0.5 -p alignment.sto

# RIBOSUM scoring for better covariation detection
RNAalifold --ribosum_scoring alignment.sto
```

| Option | Description |
|--------|-------------|
| `--cfactor` | Covariation weight (default: 1.0, lower = more thermodynamic) |
| `--nfactor` | Non-compatible penalty (default: 1.0) |
| `--ribosum_scoring` | Use RIBOSUM matrices for covariation |
| `-p` | Partition function for consensus |

## RNAcofold: RNA-RNA Interaction

Predicts the hybridization structure of two RNA molecules.

```bash
# Two sequences separated by '&'
echo "GCGCGC&GCGCGC" | RNAcofold

# With partition function
echo "GCGCGC&GCGCGC" | RNAcofold -p
```

## LinearFold: Fast Folding for Long Sequences

For sequences longer than ~5,000 nt, LinearFold provides O(n) time complexity instead of O(n^3).

```bash
# LinearFold (if installed separately)
echo "GGGAAACCC" | linearfold

# ViennaRNA also supports --maxBPspan for long sequences
RNAfold --maxBPspan 300 < long_sequence.fa
```

## ViennaRNA Python API

```python
import RNA

sequence = 'GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA'

# MFE folding
structure, mfe = RNA.fold(sequence)
print(f'Structure: {structure}')
print(f'MFE: {mfe:.2f} kcal/mol')

# Partition function and base-pair probabilities
fc = RNA.fold_compound(sequence)
structure_pf, pf_energy = fc.pf()
print(f'Ensemble energy: {pf_energy:.2f} kcal/mol')

# Base-pair probability matrix
bpp = fc.bpp()

# Centroid structure (most representative of the ensemble)
centroid, centroid_dist = fc.centroid()
print(f'Centroid: {centroid}')
print(f'Distance to ensemble: {centroid_dist:.2f}')

# MEA structure (maximum expected accuracy)
mea_struct, mea_val = fc.MEA()
print(f'MEA structure: {mea_struct}')
print(f'MEA value: {mea_val:.2f}')
```

### Folding with Constraints (Python)

```python
import RNA

sequence = 'GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA'

md = RNA.md()
md.uniq_ML = 1  # Unique multiloop decomposition

fc = RNA.fold_compound(sequence, md)

# Force position 10 unpaired (0-indexed)
fc.hc_add_up(10, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

# Force positions 1-3 paired with 70-72
fc.hc_add_bp(1, 72, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)

structure, mfe = fc.mfe()
print(f'Constrained: {structure} ({mfe:.2f} kcal/mol)')
```

### SHAPE-Constrained Folding (Python)

```python
import RNA

sequence = 'GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA'

fc = RNA.fold_compound(sequence)

# SHAPE reactivities: negative values = no data
# Deigan et al. (2009) parameters: m=1.8, b=-0.6 (default for SHAPE)
reactivities = [-999] + [0.1, 0.05, 0.8, 0.9, 0.2, 0.1, 0.3]  # 1-indexed, -999 = missing
fc.sc_add_SHAPE_deigan(reactivities, 1.8, -0.6)

structure, mfe = fc.mfe()
print(f'SHAPE-guided: {structure} ({mfe:.2f} kcal/mol)')
```

## Structure Comparison

```python
import RNA

struct1 = '(((....)))'
struct2 = '(((....).))'

# Base-pair distance
bp_dist = RNA.bp_distance(struct1, struct2)
print(f'Base-pair distance: {bp_dist}')

# Tree edit distance (more sophisticated)
tree1 = RNA.make_tree(RNA.expand_Full(struct1))
tree2 = RNA.make_tree(RNA.expand_Full(struct2))
tree_dist = RNA.tree_edit_distance(tree1, tree2)
print(f'Tree edit distance: {tree_dist}')
```

## Structure Formats

| Format | Description | Example |
|--------|-------------|---------|
| Dot-bracket | Parentheses for pairs, dots for unpaired | `(((...)))` |
| CT (connect) | Tab-delimited: index, base, prev, next, pair, index | Standard for mfold |
| BPSEQ | Three columns: position, nucleotide, pair partner (0=unpaired) | Used by comparative databases |
| WUSS | Extended dot-bracket with pseudoknot notation | `<<..AA..>>..aa` |

### Format Conversion

```python
import RNA

sequence = 'GGGAAACCC'
structure = '(((...)))'

# Dot-bracket to base-pair list
pt = RNA.ptable(structure)
pairs = [(i, pt[i]) for i in range(1, len(pt)) if pt[i] > i]
print(f'Base pairs: {pairs}')

# Dot-bracket to BPSEQ
for i in range(1, len(sequence) + 1):
    print(f'{i} {sequence[i-1]} {pt[i]}')
```

## Visualization

### Forna (web-based)

```python
# Generate JSON for forna viewer (http://rna.tbi.univie.ac.at/forna/)
import json

forna_data = {
    'sequence': 'GGGCUAUUAGCUCAGUUGGUUAGAGCGCACC',
    'structure': '((((....((((......))))....))))'
}
print(json.dumps(forna_data))
```

### R2DT (standardized 2D layouts)

```bash
# R2DT provides template-based 2D layouts for known RNA families
# Requires Docker
docker run -v $(pwd):/data rnacentral/r2dt draw /data/input.fa /data/output/
```

### Matplotlib Dot Plot

```python
import RNA
import matplotlib.pyplot as plt
import numpy as np

sequence = 'GGGCUAUUAGCUCAGUUGGUUAGAGCGCACCCCUGAUAAGGGUGAGGUCGCUGAUUCGAAUUCAGCAUAGCCCA'
fc = RNA.fold_compound(sequence)
fc.pf()
bpp = fc.bpp()

n = len(sequence)
matrix = np.zeros((n, n))
for i in range(1, n + 1):
    for j in range(i + 1, n + 1):
        matrix[i-1][j-1] = bpp[i][j]

fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(matrix, cmap='YlOrRd', origin='lower', vmin=0, vmax=1)
ax.set_xlabel('Position')
ax.set_ylabel('Position')
ax.set_title('Base-pair probability matrix')
plt.tight_layout()
plt.savefig('bpp_dotplot.png', dpi=150)
```

## Quality Thresholds

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| MFE z-score | < -2.0 | Sequence folds significantly better than shuffled controls |
| Ensemble diversity | < 5.0 | Low diversity indicates a well-defined structure |
| Base-pair probability | > 0.9 | High confidence for individual pairs |
| Covariation score | > 0.0 | Positive covariation supports predicted pair |

## Related Skills

- ncrna-search - Classify structured RNAs by family using Infernal/Rfam
- structure-probing - Use experimental SHAPE/DMS data to constrain predictions
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation
- sequence-manipulation/sequence-properties - Sequence composition analysis
