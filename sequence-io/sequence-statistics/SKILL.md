---
name: bio-sequence-statistics
description: Calculate assembly and sequence statistics (N50/L50, auN, NG50/NGA50, length distribution, GC content with ambiguity handling, summary reports) using Biopython. Use when analyzing sequence datasets, generating QC reports, or comparing genome assemblies.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Sequence Statistics

**"Calculate N50 and other assembly statistics"** -> Compute sequence count, length distribution, N50/L50, auN, GC content, and nucleotide composition for FASTA datasets.
- Python: `SeqIO.parse()`, `gc_fraction()` (BioPython)

Calculate comprehensive statistics for sequence datasets using Biopython.

## The Governing Principle

N50 measures CONTIGUITY, not correctness. A misassembled scaffold that wrongly joins distant regions can post a large N50 while being biologically wrong; N50 says nothing about base accuracy or join correctness. Treat contiguity metrics as one axis of assembly quality, alongside completeness (BUSCO) and correctness (read-backed validation, reference alignment).

Two further traps shape every reported number:
- N50 is a discontinuous, threshold-based statistic. Near the 50% crossing, contig lengths can differ by megabases, so a tiny change can jump N50 by megabases. Prefer auN (smooth, threshold-free) for robust comparison.
- N50 uses ASSEMBLY size as the denominator, so it is not comparable across assemblies of the same genome. Use NG50 (GENOME size denominator) to compare assemblies on a common baseline.

## Required Imports

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import statistics
```

## N50, L50, and Nx Statistics

**Goal:** Report the length at which half the assembled bases reside in equal-or-longer contigs (N50), and how many contigs that takes (L50).

**Approach:** N50 is the minimal length `x` such that contigs of length `>= x` together cover `>= 50%` of total assembly length. Sort lengths DESCENDING, take the cumulative sum, and return the length at which the cumulative sum first reaches or crosses 50%. L50 is the COUNT of contigs at that crossing. Three details silently break naive implementations: the sort must be descending, the crossing test must be `>=` (not strict `>`), and the denominator must be the assembly total (not a genome estimate).

**Reference (BioPython 1.83+):**
```python
def n50_l50(lengths):
    '''Return (N50 length, L50 count) for a list of contig lengths.'''
    sorted_lengths = sorted(lengths, reverse=True)
    half = sum(sorted_lengths) / 2
    cumsum = 0
    for count, length in enumerate(sorted_lengths, start=1):
        cumsum += length
        if cumsum >= half:
            return length, count
    return 0, 0

lengths = [len(r.seq) for r in SeqIO.parse('assembly.fasta', 'fasta')]
n50, l50 = n50_l50(lengths)
print(f'N50: {n50:,} bp  L50: {l50} contigs')
```

### Any Nx Value (N75, N90)

```python
def calculate_nx(lengths, x):
    '''Nx where x is a percentage (50 for N50, 90 for N90).'''
    sorted_lengths = sorted(lengths, reverse=True)
    threshold = sum(sorted_lengths) * (x / 100)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= threshold:
            return length
    return 0

lengths = [len(r.seq) for r in SeqIO.parse('assembly.fasta', 'fasta')]
print(f'N50: {calculate_nx(lengths, 50):,} bp')
print(f'N90: {calculate_nx(lengths, 90):,} bp')
```

## auN: Robust Contiguity (Preferred for Comparison)

**Goal:** Replace the discontinuous N50 with a smooth, threshold-free contiguity score that responds to every join.

**Approach:** auN (Heng Li, 2020) is the area under the Nx curve, equivalently a length-weighted average length: each contig contributes its own length weighted by the fraction of the assembly it represents. Connecting any two contigs always raises auN, even when N50 stays unchanged (joining two contigs both above, or both below, the N50 contig leaves N50 fixed). No single straddling contig arbitrarily sets the score.

Formula: auN = sum_i(L_i^2) / sum_j(L_j)

**Reference (BioPython 1.83+):**
```python
def calculate_aun(lengths):
    '''auN = sum(L_i^2) / sum(L_j); a length-weighted mean length.'''
    total = sum(lengths)
    return sum(length * length for length in lengths) / total if total else 0

lengths = [len(r.seq) for r in SeqIO.parse('assembly.fasta', 'fasta')]
print(f'auN: {calculate_aun(lengths):,.0f} bp')
```

## NG50, NGx, NA50, NGA50: Cross-Assembly and Misassembly-Aware

NG50/NGx use the GENOME size as the denominator instead of the assembly size, so two assemblies of the same genome share one baseline and become directly comparable. They require a known or estimated genome size (QUAST `--est-ref-size` or a reference). NA50/NGA50 are computed on alignment blocks broken at misassembly breakpoints; NGA50 markedly below NG50 signals misassemblies.

**Reference (BioPython 1.83+):**
```python
def calculate_ngx(lengths, genome_size, x=50):
    '''NGx uses genome_size (not assembly size) as the denominator.'''
    sorted_lengths = sorted(lengths, reverse=True)
    threshold = genome_size * (x / 100)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= threshold:
            return length
    return 0  # assembly never covers x% of the genome

lengths = [len(r.seq) for r in SeqIO.parse('assembly.fasta', 'fasta')]
print(f'NG50: {calculate_ngx(lengths, genome_size=3_100_000_000):,} bp')
```

## Length Distribution

Median contig length is near-useless for assemblies: it is dominated by the many tiny contigs and sits among fragments, ignoring where the sequence mass lives. N50 and auN are mass-weighted precisely to answer "in contigs of what size does the bulk of the genome reside?" Report median for read-length QC, not for assembly contiguity.

```python
lengths = [len(r.seq) for r in SeqIO.parse('sequences.fasta', 'fasta')]
print(f'Count: {len(lengths)}  Total: {sum(lengths):,} bp')
print(f'Min: {min(lengths):,}  Max: {max(lengths):,}  Mean: {statistics.mean(lengths):,.1f} bp')
```

### Length Histogram Data

```python
from collections import Counter

lengths = [len(r.seq) for r in SeqIO.parse('sequences.fasta', 'fasta')]
bin_size = 100  # 100-bp length bins
histogram = Counter((l // bin_size) * bin_size for l in lengths)
for length_bin in sorted(histogram):
    print(f'{length_bin}-{length_bin + bin_size}: {histogram[length_bin]}')
```

## GC Content: Choose the Ambiguity Mode Explicitly

`gc_fraction(seq, ambiguous=...)` returns a FRACTION 0-1 (the old `Bio.SeqUtils.GC()` returned a PERCENT 0-100 and was REMOVED in 1.82; swapping names without rescaling is a silent 100x error). The `ambiguous=` argument changes the answer, so set it on purpose:

| Mode | Numerator | Denominator | `GCGCNNNN` |
|------|-----------|-------------|------------|
| `remove` (default) | G, C, S | only unambiguous + S/W (N excluded) | 1.0 |
| `ignore` | G, C, S | full length (N dilutes GC) | 0.5 |
| `weighted` | G, C, S + each ambiguous code x its expected GC | full length | 0.75 |

A naive `(G + C) / len` silently equals the `ignore` mode, under-reporting GC whenever N is present. `remove` reports GC among called bases; `ignore` reports GC over the full sequence including gaps/Ns; `weighted` apportions each IUPAC code its expected GC.

```python
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

seq = Seq('GCGCNNNN')
gc_fraction(seq, ambiguous='remove')    # 1.0  - N dropped from both
gc_fraction(seq, ambiguous='ignore')    # 0.5  - N counted in denominator
gc_fraction(seq, ambiguous='weighted')  # 0.75 - N contributes 0.5 each
```

### Per-Sequence GC Distribution

```python
gc_values = [gc_fraction(r.seq, ambiguous='remove') for r in SeqIO.parse('sequences.fasta', 'fasta')]
print(f'Mean GC: {statistics.mean(gc_values):.1%}')
print(f'Median GC: {statistics.median(gc_values):.1%}')
print(f'Range: {min(gc_values):.1%} - {max(gc_values):.1%}')
```

## Comprehensive Summary Report

**Goal:** Generate a complete QC summary (counts, lengths, N50/L50, auN, GC) for any FASTA file in one pass.

**Approach:** Load all records once, compute length and GC arrays, derive N50/L50 from the cumulative sorted lengths and auN from the squared-length sum, and package into a dictionary.

**Reference (BioPython 1.83+):**
```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import statistics

def sequence_summary(fasta_file):
    records = list(SeqIO.parse(fasta_file, 'fasta'))
    lengths = [len(r.seq) for r in records]
    gc_values = [gc_fraction(r.seq, ambiguous='remove') for r in records]

    sorted_lengths = sorted(lengths, reverse=True)
    total_bp = sum(lengths)
    half = total_bp / 2

    cumsum, n50, l50 = 0, 0, 0
    for count, length in enumerate(sorted_lengths, start=1):
        cumsum += length
        if cumsum >= half:
            n50, l50 = length, count
            break

    aun = sum(length * length for length in lengths) / total_bp if total_bp else 0

    return {
        'file': fasta_file, 'sequences': len(records), 'total_bp': total_bp,
        'min_length': min(lengths), 'max_length': max(lengths),
        'mean_length': statistics.mean(lengths), 'median_length': statistics.median(lengths),
        'n50': n50, 'l50': l50, 'aun': aun,
        'gc_mean': statistics.mean(gc_values),
        'gc_std': statistics.stdev(gc_values) if len(gc_values) > 1 else 0,
    }

stats = sequence_summary('assembly.fasta')
print(f'Sequences: {stats["sequences"]:,}  Total: {stats["total_bp"]:,} bp')
print(f'N50: {stats["n50"]:,} bp (L50: {stats["l50"]})  auN: {stats["aun"]:,.0f} bp')
print(f'GC: {stats["gc_mean"]:.1%} (+/- {stats["gc_std"]:.1%})')
```

## Compare Multiple Assemblies

**Goal:** Build a side-by-side table of key metrics across assembly files.

**Approach:** Run `sequence_summary` on each file and format the results into an aligned table; auN is the most reliable single column for ranking contiguity.

**Reference (BioPython 1.83+):**
```python
from pathlib import Path

files = sorted(Path('assemblies/').glob('*.fasta'))
print(f'{"File":<30} {"Seqs":>8} {"Total bp":>15} {"N50":>12} {"auN":>12}')
print('-' * 80)
for fasta_file in files:
    s = sequence_summary(str(fasta_file))
    print(f'{fasta_file.name:<30} {s["sequences"]:>8,} {s["total_bp"]:>15,} {s["n50"]:>12,} {s["aun"]:>12,.0f}')
```

## Nucleotide Composition

```python
from collections import Counter

def nucleotide_composition(fasta_file):
    counts = Counter()
    for record in SeqIO.parse(fasta_file, 'fasta'):
        counts.update(str(record.seq).upper())
    total = sum(counts.values())
    return {base: count / total for base, count in counts.items()}

comp = nucleotide_composition('sequences.fasta')
for base in ['A', 'T', 'G', 'C', 'N']:
    if base in comp:
        print(f'{base}: {comp[base]:.2%}')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| N50 looks far too small | Sorted ascending instead of descending | Sort lengths with `reverse=True` before the cumulative sum |
| N50 off by one contig near the crossing | Strict `>` test misses the exact-50% case | Use `>= 50%` (minimal length that reaches/exceeds half) |
| N50 not comparable between assemblies | Used assembly size as denominator | Use NG50 with the genome size for cross-assembly comparison |
| GC values off by 100x | Treated `gc_fraction` (0-1) like old `GC()` (0-100) | Multiply by 100 only for display; never mix the two |
| GC silently low when Ns present | Default `remove` vs naive `(G+C)/len` (= `ignore`) | Pass `ambiguous=` explicitly to match intent |
| Huge N50 on a wrong assembly | N50 measures contiguity, not correctness | Pair with BUSCO completeness and read-backed/reference validation; prefer auN |

## References

- Li H (2020). "auN: a new metric to measure assembly contiguity." Technical blog post, https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity (auN = area under the Nx curve; smooth, threshold-free contiguity).

## Related Skills

- read-sequences - Parse sequences for statistics calculation
- batch-processing - Calculate stats across multiple files
- fastq-quality - Quality score statistics for FASTQ files
- sequence-manipulation/sequence-properties - Per-sequence GC content and properties
- alignment-files/bam-statistics - samtools stats/flagstat for alignment statistics
