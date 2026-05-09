---
name: bio-alignment-filtering
description: Filter alignments by flags, mapping quality, and regions using samtools view and pysam. Use when extracting specific reads, removing low-quality alignments, or subsetting to target regions.
tool_type: cli
primary_tool: samtools
---

## Version Compatibility

Reference examples tested with: pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment Filtering

**"Filter my BAM file to keep only high-quality reads"** → Select reads by FLAG bits, mapping quality, and genomic regions using samtools view or pysam.
- CLI: `samtools view` with `-F`/`-f`/`-q`/`-L` flags (samtools)
- Python: `pysam.AlignmentFile` iteration with attribute filters (pysam)

Filter alignments by flags, quality, and regions using samtools and pysam.

## Filter Flags

| Option | Description |
|--------|-------------|
| `-f FLAG` | Include reads with ALL bits set |
| `-F FLAG` | Exclude reads with ANY bits set |
| `-G FLAG` | Exclude reads with ALL bits set |
| `-q MAPQ` | Minimum mapping quality |
| `-L BED` | Include reads overlapping regions |

## Common FLAG Values

| Flag | Hex | Meaning |
|------|-----|---------|
| 1 | 0x1 | Paired |
| 2 | 0x2 | Proper pair |
| 4 | 0x4 | Unmapped |
| 8 | 0x8 | Mate unmapped |
| 16 | 0x10 | Reverse strand |
| 32 | 0x20 | Mate reverse strand |
| 64 | 0x40 | First in pair (read1) |
| 128 | 0x80 | Second in pair (read2) |
| 256 | 0x100 | Secondary alignment |
| 512 | 0x200 | Failed QC |
| 1024 | 0x400 | Duplicate |
| 2048 | 0x800 | Supplementary |

## Filter by FLAG

### Keep Only Mapped Reads
```bash
samtools view -F 4 -o mapped.bam input.bam
```

### Keep Only Unmapped Reads
```bash
samtools view -f 4 -o unmapped.bam input.bam
```

### Keep Only Properly Paired
```bash
samtools view -f 2 -o proper.bam input.bam
```

### Remove Duplicates
```bash
samtools view -F 1024 -o nodup.bam input.bam
```

### Remove Secondary and Supplementary
```bash
samtools view -F 2304 -o primary.bam input.bam
```

### Keep Only Primary Alignments
```bash
samtools view -F 256 -F 2048 -o primary.bam input.bam
# Or combined: -F 2304
```

### Keep Read1 Only
```bash
samtools view -f 64 -o read1.bam input.bam
```

### Keep Read2 Only
```bash
samtools view -f 128 -o read2.bam input.bam
```

### Forward Strand Only
```bash
samtools view -F 16 -o forward.bam input.bam
```

### Reverse Strand Only
```bash
samtools view -f 16 -o reverse.bam input.bam
```

## Filter by Mapping Quality

### Minimum MAPQ
```bash
samtools view -q 30 -o highqual.bam input.bam
```

### MAPQ and Mapped
```bash
samtools view -F 4 -q 30 -o filtered.bam input.bam
```

### Aligner-Aware MAPQ Thresholds

MAPQ scales differ by aligner; the same `-q 30` filter does different things. See sam-bam-basics for the full MAPQ-by-aligner table. Filtering recommendations:

| Aligner | "Drop ambiguous" | "High confidence" |
|---------|------------------|-------------------|
| BWA-MEM / BWA-MEM2 | `-q 1` | `-q 30` (or `-q 60` for unique only) |
| Bowtie2 | `-q 1` | `-q 23` (~99% accurate per Heng Li) |
| **STAR** | `-q 255` | `-q 255` (255 is the unique-mapped sentinel; -q 60 drops everything) |
| HISAT2 | `-q 1` | `-q 60` |
| minimap2 (DNA, long-read) | `-q 1` | `-q 60` |
| pbmm2 (PacBio) | `-q 1` | `-q 60` |

For Phred-scaled aligners (BWA, minimap2), MAPQ Q maps to ~10^(-Q/10) probability of wrong mapping. For STAR, the values 0/1/2/3/255 are sentinels, not probabilities.

### Drop Ambiguous Across Aligners (Universal)
```bash
samtools view -q 1 in.bam   # exclude MAPQ=0; works for all aligners
```

## Filter by Region

### Single Region
```bash
samtools view -o region.bam input.bam chr1:1000000-2000000
```

### Multiple Regions
```bash
samtools view -o regions.bam input.bam chr1:1000-2000 chr2:3000-4000
```

### Regions from BED File
```bash
samtools view -L targets.bed -o targets.bam input.bam
```

### Combine Region and Quality
```bash
samtools view -q 30 -L targets.bed -o filtered.bam input.bam
```

## Combined Filters

### Standard Quality Filter

**Goal:** Produce a clean BAM containing only primary, mapped, non-duplicate reads with high mapping confidence.

**Approach:** Combine FLAG exclusion (-F for unmapped + secondary + duplicate + supplementary) with a MAPQ threshold.

**Reference (samtools 1.19+):**
```bash
samtools view -F 3332 -q 30 -o filtered.bam input.bam
# 3332 = 4 (unmapped) + 256 (secondary) + 1024 (duplicate) + 2048 (supplementary)
```

### Variant Calling Prep -- Assay-Aware

**Goal:** Choose a filter that matches what the downstream caller expects. Stripping supplementary alignments breaks SV callers; requiring proper-pair drops valid spliced RNA-seq reads.

| Assay / caller | Recommended filter | Why |
|----------------|-------------------|-----|
| Germline WGS short-variant (HaplotypeCaller, DeepVariant) | `-f 2 -F 3328 -q 20` | Primary, no dup, proper pair, MAPQ>=20 |
| Somatic short-variant (Mutect2, Strelka2) | `-F 3328 -q 1` | Drop only MAPQ=0; somatic callers handle low MAPQ; chimeric reads at SVs may carry real somatic SNVs |
| Long-read short-variant (clair3, DeepVariant ONT) | `-F 3328 -q 5` | Long-read MAPQ scale is lower |
| Long-read SV (Sniffles, cuteSV) | `-F 1024` only | **Keep supplementary** -- SA tag is the SV signal |
| Short-read SV (Manta, GRIDSS, Delly, SvABA) | `-F 1024` only | Same -- supplementary required |
| ChIP-seq peak calling | `-F 1804 -q 30` | Drop dup + secondary + supp + unmapped + mate-unmapped + QC-fail |
| ATAC-seq | `-F 1804 -q 30 -f 2` | Same plus proper pair |
| RNA-seq quantification (STAR) | `-q 255` | Unique only (STAR sentinel) |
| RNA-seq quantification (HISAT2) | `-F 256 -q 60` | Different aligner semantics |
| RNA-seq variant (after `SplitNCigarReads`) | `-F 3328 -q 20` | Standard germline after split-N-trim |
| Panel / amplicon | After `samtools ampliconclip`; `-F 1024 -q 20` | Primer overlap makes proper-pair unreliable |
| ctDNA / cfDNA (UMI) | After fgbio consensus; do not pre-filter raw | |

**Reference (samtools 1.19+):**
```bash
# Short-variant germline
samtools view -f 2 -F 3328 -q 20 -o clean.bam input.bam
# 3328 = 256 (secondary) + 1024 (duplicate) + 2048 (supplementary)

# SV calling: KEEP supplementary
samtools view -F 1024 -o sv_input.bam input.bam   # NOT -F 2304 or -F 3328

# ChIP-seq / ATAC-seq common filter
samtools view -F 1804 -q 30 -o filtered.bam input.bam
# 1804 = 4 + 8 + 256 + 512 + 1024 = unmapped + mate-unmapped + secondary + QC-fail + duplicate
```

**Cost of getting this wrong:** filtering `-F 2304` or `-F 3328` before SV calling produces zero SV calls -- a single-flag mistake that silently invalidates the analysis.

## Subsample Reads (Deterministic, Pair-Consistent)

`samtools view -s SEED.FRAC` -- integer is the hash seed; fractional is the keep fraction. The hash is on QNAME, so:
1. Mate consistency: read1 and read2 are kept or dropped together.
2. Reproducibility: same seed + same fraction returns the same reads.
3. **Sequential downsampling requires different seeds.** `-s 1.5` then `-s 1.25` keeps a nested 5/8 of the original (not 12.5%). Use different integer seeds for independent samples.

```bash
# 10% with seed 42 (always the same reads; pair-consistent)
samtools view -s 42.1 -b -o subset.bam input.bam

# Sequential cuts with INDEPENDENT seeds
samtools view -s 1.5 -b in.bam > half1.bam
samtools view -s 2.25 -b half1.bam > quarter.bam   # 12.5% of original

# Coverage-matching to a target read count
total=$(samtools view -c -F 2304 input.bam)
target=10000000
frac=$(awk -v t=$target -v n=$total 'BEGIN{printf "%.6f", t/n}')
samtools view -s "1.${frac#*.}" -b -o matched.bam input.bam

# Tumor-normal coverage matching (pull tumor down to normal)
normal_reads=$(samtools view -c -F 2308 normal.bam)
tumor_reads=$(samtools view -c -F 2308 tumor.bam)
if [ "$tumor_reads" -gt "$normal_reads" ]; then
    frac=$(awk -v n=$normal_reads -v t=$tumor_reads 'BEGIN{printf "%.6f", n/t}')
    samtools view -s "1.${frac#*.}" -b -o tumor_matched.bam tumor.bam
fi
```

A subsampled BAM without an integer seed (`-s 0.1`) is non-reproducible -- production pipelines should reject it.

## Expression Filtering

`samtools view -e EXPR` (or `--expr`, since samtools 1.16) supports arbitrary expression filtering on tags, FLAG, MAPQ, RNAME, CIGAR, etc. Powerful for filtering by `NM`, `AS`, `NH`, `cs`, etc. that the FLAG-based filters cannot reach:
```bash
# Reads with >=2 mismatches (NM tag)
samtools view -e '[NM] >= 2' in.bam

# Soft clip on the left, on chr1
samtools view -e 'cigar=~"^[0-9]+S" && rname=="chr1"' in.bam

# Combine with FLAG and MAPQ
samtools view -F 2308 -q 30 -e '[NM] <= 5 && [AS] >= 100' in.bam

# Drop reads with low mapped fraction (samtools-internal helpers)
samtools view -e 'sclen / qlen < 0.2' in.bam
```

Note: in samtools 1.16+, `![NM]` is true only if NM is missing (was buggy in earlier versions); NULL values from missing tags propagate through arithmetic.

## Filter by Read Group
```bash
samtools view -r library_A in.bam              # single read group
samtools view -R rg_list.txt in.bam            # multiple via file (one ID per line)
```

## pysam Python Alternative

### Basic Filtering
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('filtered.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if read.is_unmapped:
                continue
            if read.mapping_quality < 30:
                continue
            if read.is_duplicate:
                continue
            outfile.write(read)
```

### Filter with Function

**Goal:** Apply a multi-criteria quality filter to produce clean alignments for downstream analysis.

**Approach:** Define a predicate checking mapped status, primary alignment, duplicate flag, and MAPQ; stream reads through it.

**Reference (pysam 0.22+):**
```python
import pysam

def passes_filter(read):
    if read.is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if read.is_duplicate:
        return False
    if read.mapping_quality < 30:
        return False
    return True

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('filtered.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if passes_filter(read):
                outfile.write(read)
```

### Filter by Region
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('region.bam', 'wb', header=infile.header) as outfile:
        for read in infile.fetch('chr1', 1000000, 2000000):
            outfile.write(read)
```

### Filter from BED File

**Goal:** Extract only reads overlapping target regions defined in a BED file.

**Approach:** Parse BED into a list of (chrom, start, end) tuples, then fetch reads from each region and write to output.

**Reference (pysam 0.22+):**
```python
import pysam

def read_bed(bed_path):
    regions = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions

regions = read_bed('targets.bed')

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('targets.bam', 'wb', header=infile.header) as outfile:
        for chrom, start, end in regions:
            for read in infile.fetch(chrom, start, end):
                outfile.write(read)
```

### Subsample (Pair-Consistent)

Hash on QNAME so mates stay together (a fresh `random.random()` per read drops mates inconsistently and breaks paired-end tools):
```python
import pysam
import zlib

fraction = 0.1
seed = 42
threshold = int(0xffffffff * fraction)

def template_hash(qname, seed):
    return zlib.crc32(qname.encode()) ^ seed

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('subset.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if template_hash(read.query_name, seed) <= threshold:
                outfile.write(read)
```

## Quick Reference

| Task | samtools command |
|------|------------------|
| Mapped only | `view -F 4` |
| Unmapped only | `view -f 4` |
| Properly paired | `view -f 2` |
| Primary only | `view -F 2304` |
| No duplicates | `view -F 1024` |
| High MAPQ | `view -q 30` |
| Region | `view file.bam chr1:1-1000` |
| BED regions | `view -L file.bed` |
| Subsample 10% (reproducible) | `view -s 42.1` |
| Standard filter | `view -F 3332 -q 30` |

## Common Filter Combinations

| Purpose | Flags |
|---------|-------|
| Clean reads | `-F 3332 -q 30` (mapped, primary, no dups, high qual) |
| Variant calling | `-f 2 -F 3328 -q 20` (proper pair, primary, no dups) |
| Coverage analysis | `-F 1284 -q 1` (mapped, primary, no dups) |
| Count unique | `-F 2304` (primary only) |

Flag breakdowns:
- 2304 = 256 + 2048 (secondary + supplementary)
- 3328 = 256 + 1024 + 2048 (secondary + duplicate + supplementary)
- 3332 = 4 + 256 + 1024 + 2048 (unmapped + secondary + duplicate + supplementary)
- 1284 = 4 + 256 + 1024 (unmapped + secondary + duplicate)

## Related Skills

- sam-bam-basics - FLAG semantics, MAPQ-by-aligner, secondary vs supplementary
- alignment-sorting - Sort before/after filtering
- alignment-indexing - Required for region filtering
- alignment-amplicon-clipping - Primer clipping for amplicon panels
- duplicate-handling - Mark duplicates before filtering
- bam-statistics - Check filter effects
