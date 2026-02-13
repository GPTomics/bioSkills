---
name: bio-rna-structure-structure-probing
description: Analyzes experimental RNA structure probing data from SHAPE-MaP and DMS-MaPseq experiments using ShapeMapper2. Converts mutation rates to per-nucleotide reactivity profiles that constrain structure prediction. Use when processing SHAPE-MaP or DMS-MaPseq sequencing data to obtain experimental RNA structure information.
tool_type: cli
primary_tool: ShapeMapper2
---

## Version Compatibility

Reference examples tested with: STAR 2.7.11+, eggNOG-mapper 2.1+, matplotlib 3.8+, numpy 1.26+, pandas 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structure Probing

**"Process my SHAPE-MaP experiment to get RNA reactivity profiles"** → Convert mutation rates from SHAPE-MaP or DMS-MaPseq sequencing data into per-nucleotide reactivity profiles, then use reactivities as constraints for thermodynamic structure prediction.
- CLI: `shapemapper` (ShapeMapper2) for end-to-end SHAPE-MaP processing
- CLI: `RNAfold --shape` (ViennaRNA) for SHAPE-constrained folding

Analyze experimental RNA structure probing data (SHAPE-MaP, DMS-MaPseq) to obtain per-nucleotide reactivity profiles. High reactivity indicates flexible/single-stranded nucleotides; low reactivity indicates base-paired/structured positions. Reactivities constrain thermodynamic folding for more accurate structure prediction.

## Platform Note

ShapeMapper2 is Linux-only. On macOS, use Docker or Singularity:

```bash
# Docker
docker pull shapemapper2/shapemapper2
docker run -v $(pwd):/data shapemapper2/shapemapper2 shapemapper \
    --target /data/target.fa --modified --R1 /data/mod_R1.fq.gz --R2 /data/mod_R2.fq.gz \
    --untreated --R1 /data/unmod_R1.fq.gz --R2 /data/unmod_R2.fq.gz \
    --out /data/results

# Singularity
singularity pull shapemapper2.sif docker://shapemapper2/shapemapper2
singularity exec -B $(pwd):/data shapemapper2.sif shapemapper \
    --target /data/target.fa --modified --R1 /data/mod_R1.fq.gz ...
```

## ShapeMapper2 Pipeline

### Basic SHAPE-MaP Analysis

```bash
shapemapper \
    --target target_rna.fa \
    --name my_rna \
    --modified --R1 modified_R1.fastq.gz --R2 modified_R2.fastq.gz \
    --untreated --R1 untreated_R1.fastq.gz --R2 untreated_R2.fastq.gz \
    --out results/ \
    --nproc 8 \
    --min-depth 5000
```

### With Denatured Control

A denatured control normalizes for sequence-dependent mutation biases.

```bash
shapemapper \
    --target target_rna.fa \
    --name my_rna \
    --modified --R1 modified_R1.fastq.gz --R2 modified_R2.fastq.gz \
    --untreated --R1 untreated_R1.fastq.gz --R2 untreated_R2.fastq.gz \
    --denatured --R1 denatured_R1.fastq.gz --R2 denatured_R2.fastq.gz \
    --out results/ \
    --nproc 8
```

### Amplicon Mode (Targeted)

```bash
# For PCR-amplified targets
shapemapper \
    --target target_rna.fa \
    --name my_rna \
    --amplicon \
    --modified --R1 modified_R1.fastq.gz --R2 modified_R2.fastq.gz \
    --untreated --R1 untreated_R1.fastq.gz --R2 untreated_R2.fastq.gz \
    --out results/ \
    --nproc 8
```

### Key ShapeMapper2 Options

| Option | Description |
|--------|-------------|
| `--target` | FASTA file with target RNA sequence(s) |
| `--name` | Name for output files |
| `--modified` | Modified sample FASTQ files |
| `--untreated` | Untreated control FASTQ files |
| `--denatured` | Denatured control FASTQ files |
| `--amplicon` | Amplicon sequencing mode |
| `--nproc` | Number of processors |
| `--min-depth` | Minimum read depth per nucleotide (default: 5000) |
| `--min-qual-to-count` | Minimum base quality (default: 20) |
| `--overwrite` | Overwrite existing output |
| `--star-aligner` | Use STAR instead of Bowtie2 for alignment |

### ShapeMapper2 Output Files

```
results/
├── my_rna_profile.txt          # Per-nucleotide reactivity profile
├── my_rna_shapemapper_log.txt  # Run log and QC metrics
├── my_rna_histograms/          # Mutation rate histograms
├── my_rna_profile.pdf          # Reactivity bar plot
└── my_rna_map.shape            # SHAPE reactivities for RNAfold
```

## Reactivity Interpretation

### SHAPE Reactivity Scale

| Value | Interpretation |
|-------|---------------|
| < 0.3 | Low reactivity, likely base-paired |
| 0.3 - 0.7 | Moderate, partially flexible or dynamic |
| > 0.7 | High reactivity, likely single-stranded |
| -999 | No data (insufficient read depth) |

### Quality Filters

| Metric | Threshold | Rationale |
|--------|-----------|-----------|
| Read depth | >= 5,000 | Minimum for reliable mutation rate estimation |
| Effective depth | >= 1,000 | After quality filtering |
| Modification rate (untreated) | < 0.5% | High background suggests RNA damage or mapping artifact |
| Modification rate (modified) | 1-10% | Too low = undermodified; too high = degraded |

## Structure-Guided Folding with SHAPE Data

### RNAfold with SHAPE Constraints

```bash
# Use ShapeMapper2 .shape output directly with RNAfold
RNAfold --shape=results/my_rna_map.shape < target_rna.fa

# Adjust SHAPE slope/intercept for different reagents
# Default (1NAI/NMIA): m=1.8, b=-0.6
RNAfold --shape=results/my_rna_map.shape --shapeMethod="Dm1.8b-0.6" < target_rna.fa

# Convert to hard constraints at extreme positions
# Reactivity > 0.7 forced unpaired, < 0.3 allowed to pair
RNAfold --shape=results/my_rna_map.shape --shapeConversion="S" < target_rna.fa
```

### Python: SHAPE-Constrained Folding

```python
import RNA
import pandas as pd

def load_shape_profile(profile_file):
    '''Load ShapeMapper2 reactivity profile.'''
    df = pd.read_csv(profile_file, sep='\t')
    reactivities = df['Reactivity_profile'].tolist()
    sequence = ''.join(df['Nucleotide'].tolist())
    return sequence, reactivities


def fold_with_shape(sequence, reactivities, m=1.8, b=-0.6):
    '''
    Fold RNA constrained by SHAPE reactivities.

    m, b: Deigan et al. (2009) parameters.
    m=1.8, b=-0.6: standard for SHAPE (1M7, NMIA, NAI).
    m=1.1, b=-0.3: suggested for DMS-MaPseq (A/C only).
    '''
    fc = RNA.fold_compound(sequence)

    # Prepend -999 for 0-indexed padding (ViennaRNA expects 1-indexed)
    shape_data = [-999.0] + [r if r != -999 else -999.0 for r in reactivities]
    fc.sc_add_SHAPE_deigan(shape_data, m, b)

    structure, mfe = fc.mfe()

    # Also compute partition function for confidence
    fc2 = RNA.fold_compound(sequence)
    fc2.sc_add_SHAPE_deigan(shape_data, m, b)
    _, pf_energy = fc2.pf()
    centroid, _ = fc2.centroid()

    return {
        'mfe_structure': structure,
        'mfe_energy': mfe,
        'centroid_structure': centroid,
        'ensemble_energy': pf_energy
    }


def shape_agreement(structure, reactivities, low_threshold=0.3, high_threshold=0.7):
    '''
    Assess agreement between SHAPE data and predicted structure.

    Returns fraction of nucleotides where reactivity agrees with structure
    (low reactivity = paired, high reactivity = unpaired).
    '''
    agree, total = 0, 0
    for i, (char, react) in enumerate(zip(structure, reactivities)):
        if react == -999:
            continue
        total += 1
        paired = char in '()'
        if paired and react < low_threshold:
            agree += 1
        elif not paired and react > high_threshold:
            agree += 1
    return agree / total if total > 0 else 0.0
```

## DMS-MaPseq Analysis

DMS methylates unpaired A and C residues. SEISMIC-RNA (successor to DREEM) or rf-count process DMS-MaPseq data.

### SEISMIC-RNA

```bash
# Install
pip install seismic-rna

# Align reads to target
seismic align target.fa reads_R1.fq.gz reads_R2.fq.gz --out-dir seismic_out

# Relate mutations to reference
seismic relate seismic_out/align --out-dir seismic_out

# Compute mutation rates per position
seismic mask seismic_out/relate --out-dir seismic_out

# Cluster structural conformations (if RNA adopts multiple states)
seismic cluster seismic_out/mask --out-dir seismic_out --max-clusters 3
```

### rf-count (RNAframework)

```bash
# Count mutations with rf-count
rf-count -t target.fa -r modified.bam -rc untreated.bam -o rf_results/

# Normalize reactivities with rf-norm
rf-norm -t target.fa -i rf_results/ -o rf_norm/ -sm 3 -nm 2
```

### DMS vs SHAPE Differences

| Feature | SHAPE (1M7/NAI) | DMS-MaPseq |
|---------|-----------------|------------|
| Reactive nucleotides | All four (A, C, G, U) | Primarily A, C |
| Resolution | All positions | ~50% of positions |
| In-cell probing | NAI-N3 for in vivo | Standard DMS works in vivo |
| Reagent cost | Moderate | Low |
| Analysis tools | ShapeMapper2 | SEISMIC-RNA, rf-count |

## Visualization

### Plot Reactivity Profile

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def plot_reactivity_profile(profile_file, output_file='reactivity_profile.png', title=None):
    '''Plot SHAPE/DMS reactivity bar chart with confidence coloring.'''
    df = pd.read_csv(profile_file, sep='\t')
    positions = df.index + 1
    reactivities = df['Reactivity_profile'].values

    valid = reactivities != -999
    colors = np.where(reactivities > 0.7, '#FF0000',   # high = red = unpaired
             np.where(reactivities > 0.3, '#FFA500',   # moderate = orange
                                          '#000000'))   # low = black = paired

    fig, ax = plt.subplots(figsize=(max(8, len(positions) * 0.08), 4))
    ax.bar(positions[valid], reactivities[valid], color=colors[valid], width=1.0, edgecolor='none')
    ax.axhline(y=0.3, color='gray', linestyle='--', linewidth=0.5)
    ax.axhline(y=0.7, color='gray', linestyle='--', linewidth=0.5)
    ax.set_xlabel('Nucleotide position')
    ax.set_ylabel('Reactivity')
    ax.set_title(title or 'SHAPE reactivity profile')
    ax.set_xlim(0, len(positions) + 1)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f'Saved reactivity profile to {output_file}')
```

### Arc Diagram with Reactivities

```python
import RNA
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


def plot_arc_diagram(sequence, structure, reactivities=None, output_file='arc_diagram.png'):
    '''Draw arc diagram colored by SHAPE reactivity.'''
    n = len(sequence)
    pt = RNA.ptable(structure)

    fig, ax = plt.subplots(figsize=(max(8, n * 0.1), 4))

    if reactivities is not None:
        valid_react = [r for r in reactivities if r != -999 and r >= 0]
        vmax = np.percentile(valid_react, 95) if valid_react else 1.0
        cmap = plt.cm.YlOrRd
        for i in range(n):
            r = reactivities[i]
            color = 'gray' if r == -999 else cmap(min(r / vmax, 1.0))
            ax.plot(i + 1, 0, 'o', color=color, markersize=3)
    else:
        for i in range(n):
            ax.plot(i + 1, 0, 'o', color='black', markersize=3)

    for i in range(1, n + 1):
        j = pt[i]
        if j > i:
            center = (i + j) / 2
            radius = (j - i) / 2
            arc = patches.Arc((center, 0), j - i, j - i, angle=0,
                              theta1=0, theta2=180, color='steelblue', linewidth=0.5)
            ax.add_patch(arc)

    ax.set_xlim(0, n + 1)
    ax.set_ylim(-1, n / 2 + 1)
    ax.set_xlabel('Position')
    ax.set_title('Arc diagram')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()
```

## Related Skills

- secondary-structure-prediction - Unconstrained and SHAPE-constrained RNA folding
- clip-seq/binding-site-annotation - Protein-RNA binding site annotation
- epitranscriptomics/m6a-peak-calling - RNA modification detection
