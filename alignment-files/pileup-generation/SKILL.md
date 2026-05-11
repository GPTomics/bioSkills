---
name: bio-pileup-generation
description: Generate pileup data for variant calling using samtools mpileup and pysam. Use when preparing data for variant calling, analyzing per-position read data, or calculating allele frequencies.
tool_type: cli
primary_tool: samtools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Pileup Generation

Generate pileup data for variant calling and position-level analysis.

**"Generate pileup from BAM"** → Produce per-position read summaries showing depth, bases, and qualities.
- CLI: `samtools mpileup -f ref.fa input.bam`
- Python: `bam.pileup(chrom, start, end)` (pysam)

**"Count alleles at a position"** → Extract per-base read support at a specific genomic coordinate.
- Python: iterate `pileup_column.pileups` and count bases (pysam)

## What is Pileup?

Pileup shows all reads covering each position in the reference, used for:
- Variant calling (with bcftools)
- Coverage analysis
- Allele frequency calculation
- SNP/indel detection

## samtools mpileup vs bcftools mpileup (Deprecation)

`samtools mpileup -g/-u` (BCF output for variant calling) has been **deprecated since samtools 1.9** -- the genotype-likelihood code now lives in `bcftools mpileup`, which keeps mpileup logic versioned alongside `bcftools call` and avoids version-skew bugs.

| Use case | Recommended tool |
|----------|------------------|
| Quick allele counts at known sites | `samtools mpileup` or pysam pileup |
| Germline variant calling (small genomes, simple cohorts) | `bcftools mpileup` -> `bcftools call` |
| Germline WGS / WES production | DeepVariant or HaplotypeCaller (not mpileup) |
| Somatic SNV/indel | Mutect2 / VarDict / VarScan2 (direct from BAM) |
| Long-read small variants | clair3 / DeepVariant ONT (direct from BAM) |
| Long-read SV | Sniffles / cuteSV (direct from BAM) |
| Ultra-low-frequency (ctDNA / MRD) | fgbio consensus -> `bcftools call` or hot-spot Mutect2 |
| Per-position allele counts (custom) | pysam pileup |

`samtools mpileup` (without `-g`) is still the standard tool for human-readable per-position read summaries.

### Basic Pileup
```bash
samtools mpileup -f reference.fa input.bam > pileup.txt
```

### Pileup Specific Region
```bash
samtools mpileup -f reference.fa -r chr1:1000000-2000000 input.bam
```

### Regions from BED
```bash
samtools mpileup -f reference.fa -l targets.bed input.bam
```

### Multiple BAM Files
```bash
samtools mpileup -f reference.fa sample1.bam sample2.bam sample3.bam > pileup.txt
```

## Output Format

Text pileup format (6 columns per sample):
```
chr1    1000    A    15    ...............    FFFFFFFFFFF
chr1    1001    T    12    ............      FFFFFFFFFFFF
```

| Column | Description |
|--------|-------------|
| 1 | Chromosome |
| 2 | Position (1-based) |
| 3 | Reference base |
| 4 | Read depth |
| 5 | Read bases |
| 6 | Base qualities |

### Read Bases Encoding

| Symbol | Meaning |
|--------|---------|
| `.` | Match on forward strand |
| `,` | Match on reverse strand |
| `ACGT` | Mismatch (uppercase = forward) |
| `acgt` | Mismatch (lowercase = reverse) |
| `^Q` | Start of read (Q = MAPQ as ASCII) |
| `$` | End of read |
| `+NNN` | Insertion of N bases |
| `-NNN` | Deletion of N bases |
| `*` | Deleted base |
| `>` / `<` | Reference skip (intron) |

## Quality Filtering Options

### Minimum Mapping Quality
```bash
samtools mpileup -f reference.fa -q 20 input.bam
```

### Minimum Base Quality
```bash
samtools mpileup -f reference.fa -Q 20 input.bam
```

### Combined Quality Filters
```bash
samtools mpileup -f reference.fa -q 20 -Q 20 input.bam
```

### Maximum Depth (Critical Trap)
```bash
# samtools mpileup default -d 8000 silently truncates targeted / mt-DNA / amplicon / UMI-deduped data
# bcftools mpileup default -d 250 is far lower; both must be set explicitly when piping
samtools mpileup -f reference.fa -d 0 input.bam        # no cap
samtools mpileup -f reference.fa -d 1000000 input.bam  # explicit high cap

# WRONG -- samtools 8000 cap, then bcftools 250 cap re-applied
samtools mpileup -f ref.fa in.bam | bcftools call -mv

# RIGHT -- single tool, explicit -d
bcftools mpileup -d 1000000 -f ref.fa in.bam | bcftools call -mv
```

## BAQ: Base Alignment Quality (Critical Default)

When `-f ref.fa` is passed, BAQ is enabled by default. BAQ Phred-scales the probability that a base is misaligned (HMM realignment over a small window) and reduces base quality near indels. Tradeoffs: ~30% slower; suppresses FP SNVs near indels; hurts indel detection sensitivity.

| Flag | Behavior |
|------|----------|
| (default with `-f`) | BAQ on (computed from CIGAR if MD missing) |
| `-B` / `--no-BAQ` | Disable BAQ -- raw qualities |
| `-E` / `--redo-BAQ` | Force recompute (after BQSR; if MD stale) |

**BAQ ON for:** short-read germline SNV (BWA, Bowtie2, HISAT2), short-read somatic SNV.

**BAQ OFF (`-B`) for:** long-read variant calling (ONT, PacBio HiFi), SV calling, RNA-seq near splice junctions, viral / amplicon, ultra-deep ctDNA from consensus reads (consensus quality already inflated), aDNA (qualities pre-rescaled by mapDamage).

`-A` (count anomalous read pairs / orphans) is required for amplicon -- amplicon reads are by design not properly paired.

`-aa` (output all positions, including zero-coverage) is required for ARTIC SARS-CoV-2 consensus generation.

### Library-Typed Flags Cheat Sheet

| Library | Flags |
|---------|-------|
| Short-read germline WGS (BWA) | `-q 20 -Q 20 -d 0` (BAQ on default) |
| Short-read tumor WGS | `-q 1 -Q 13 -d 0 -B` (low MAPQ kept; BAQ off) |
| Amplicon viral (ARTIC) | `-aa -A -d 600000 -B -Q 20` |
| Capture / exome | `-q 20 -Q 20 -d 250` |
| Long-read ONT R10.4+ | `-q 30 -Q 0 -B -d 0 --max-BQ 50` |
| PacBio HiFi | `-q 20 -Q 0 -B -d 0` |
| RNA-seq variants | `-q 20 -Q 20 -B -d 0` |
| Forensic / aDNA | `-q 0 -Q 0 -A -d 0 -B` |

## Variant Calling Pipeline (Modern: bcftools mpileup)

**Goal:** Call variants from alignment data using the pileup-based approach.

**Approach:** Use `bcftools mpileup` (not `samtools mpileup -g`) so genotype-likelihood code is co-versioned with `bcftools call`. Apply quality and depth caps explicitly; annotate FORMAT fields needed for downstream filtering.

### Modern Germline Calling
```bash
bcftools mpileup -f reference.fa -d 1000000 -q 20 -Q 20 \
    --annotate FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
    input.bam | \
  bcftools call -mv -Oz -o variants.vcf.gz
bcftools index -t variants.vcf.gz
```

### Multi-Sample Joint Calling
```bash
bcftools mpileup -f reference.fa --threads 4 -d 250 -q 20 -Q 20 \
    -a FORMAT/AD,FORMAT/DP s1.bam s2.bam s3.bam | \
  bcftools call -mv --threads 4 -Oz -o joint.vcf.gz
```

For somatic / low-VAF, prefer Mutect2 / Strelka2 / DeepVariant -- materially better than mpileup-based callers.

### Overlap Detection Defaults

When fragment length < 2 * read_length, R1 and R2 overlap. `bcftools mpileup` detects overlap and counts once by default; `samtools mpileup` does not. Doubled support inflates somatic VAFs at sites covered by overlapping pairs (especially in cfDNA / FFPE). Prefer `bcftools mpileup` or use `samtools mpileup --ignore-overlaps`.

## pysam Python Alternative

### Basic Pileup
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000):
        print(f'{pileup_column.reference_name}:{pileup_column.pos} depth={pileup_column.n}')
```

### Access Reads at Position
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1000001, truncate=True):
        print(f'Position: {pileup_column.pos}')
        print(f'Depth: {pileup_column.n}')

        for pileup_read in pileup_column.pileups:
            if pileup_read.is_del:
                print('  Deletion')
            elif pileup_read.is_refskip:
                print('  Reference skip')
            else:
                qpos = pileup_read.query_position
                base = pileup_read.alignment.query_sequence[qpos]
                qual = pileup_read.alignment.query_qualities[qpos]
                print(f'  {base} (Q{qual})')
```

### Count Alleles at Position
```python
import pysam
from collections import Counter

def allele_counts(bam_path, chrom, pos):
    counts = Counter()

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, pos, pos + 1, truncate=True):
            if pileup_column.pos != pos:
                continue

            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del:
                    counts['DEL'] += 1
                elif pileup_read.is_refskip:
                    continue
                else:
                    qpos = pileup_read.query_position
                    base = pileup_read.alignment.query_sequence[qpos]
                    counts[base.upper()] += 1

    return dict(counts)

counts = allele_counts('input.bam', 'chr1', 1000000)
print(counts)  # {'A': 45, 'G': 5}
```

### Calculate Allele Frequency
```python
import pysam
from collections import Counter

def allele_frequency(bam_path, chrom, pos, min_qual=20):
    counts = Counter()

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, pos, pos + 1, truncate=True,
                                         min_base_quality=min_qual):
            if pileup_column.pos != pos:
                continue

            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                qpos = pileup_read.query_position
                base = pileup_read.alignment.query_sequence[qpos]
                counts[base.upper()] += 1

    total = sum(counts.values())
    if total == 0:
        return {}

    return {base: count / total for base, count in counts.items()}

freq = allele_frequency('input.bam', 'chr1', 1000000)
for base, f in sorted(freq.items(), key=lambda x: -x[1]):
    print(f'{base}: {f:.1%}')
```

### Pileup with Quality Filtering
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000,
                                     truncate=True,
                                     min_mapping_quality=20,
                                     min_base_quality=20):
        print(f'{pileup_column.pos}: {pileup_column.n}')
```

### Generate Pileup Text
```python
import pysam

def pileup_text(bam_path, ref_path, chrom, start, end):
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        with pysam.FastaFile(ref_path) as ref:
            for pileup_column in bam.pileup(chrom, start, end, truncate=True):
                pos = pileup_column.pos
                ref_base = ref.fetch(chrom, pos, pos + 1)
                depth = pileup_column.n

                bases = []
                for pileup_read in pileup_column.pileups:
                    if pileup_read.is_del:
                        bases.append('*')
                    elif pileup_read.is_refskip:
                        bases.append('>')
                    else:
                        qpos = pileup_read.query_position
                        base = pileup_read.alignment.query_sequence[qpos]
                        if base.upper() == ref_base.upper():
                            bases.append('.' if not pileup_read.alignment.is_reverse else ',')
                        else:
                            bases.append(base.upper() if not pileup_read.alignment.is_reverse else base.lower())

                print(f'{chrom}\t{pos+1}\t{ref_base}\t{depth}\t{"".join(bases)}')

pileup_text('input.bam', 'reference.fa', 'chr1', 1000000, 1000100)
```

## Pileup Options Summary

| Option | Description | Common pitfall |
|--------|-------------|----------------|
| `-f FILE` | Reference FASTA | Triggers BAQ ON by default |
| `-r REGION` | Restrict to region | |
| `-l FILE` | BED file of regions | |
| `-q INT` | Min mapping quality | Aligner-dependent semantics |
| `-Q INT` | Min base quality | `-Q 0` with default overlap detection has subtle behavior |
| `-d INT` | Max depth | **Default 8000 silently truncates**; bcftools mpileup default is 250 |
| `-B` | Disable BAQ | Often correct for long reads, SV, viral, consensus |
| `-A` | Count anomalous pairs | Required for amplicon |
| `-aa` | Output all positions | Required for consensus generation |
| `--ignore-overlaps` | Disable mate-overlap correction | Rarely correct |
| `--max-BQ INT` | Cap BQ (default 60) | Useful for ONT (Q values inflated) |
| `-g` (DEPRECATED) | Old BCF output | Use `bcftools mpileup` instead |

## Quick Reference

| Task | Command |
|------|---------|
| Basic pileup | `samtools mpileup -f ref.fa in.bam` |
| Quality filter | `samtools mpileup -f ref.fa -q 20 -Q 20 in.bam` |
| Region | `samtools mpileup -f ref.fa -r chr1:1-1000 in.bam` |
| To bcftools | `bcftools mpileup -f ref.fa -d 1000000 in.bam \| bcftools call -mv` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `No FASTA reference` | Missing -f option | Add `-f reference.fa` |
| `Reference mismatch` | Wrong reference | Use same reference as alignment |
| Out of memory | High coverage region | Use `-d` to cap depth |

## Related Skills

- alignment-filtering - Filter BAM before pileup
- reference-operations - Index reference for pileup; M5 cross-check
- bam-statistics - mosdepth, depth tool selection
- variant-calling/variant-calling - Full variant calling workflows
- variant-calling/vcf-basics - VCF/BCF I/O
- variant-calling/joint-calling - Multi-sample joint calling
