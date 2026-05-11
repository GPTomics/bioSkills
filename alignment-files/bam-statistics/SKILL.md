---
name: bio-bam-statistics
description: Generate alignment statistics using samtools flagstat, stats, depth, coverage, and mosdepth. Use when assessing alignment quality, calculating coverage, or generating QC reports.
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

# BAM Statistics

**"Get alignment statistics and coverage from my BAM file"** → Generate read counts, mapping rates, per-chromosome statistics, depth profiles, and coverage summaries.
- CLI: `samtools flagstat`, `samtools stats`, `samtools depth`, `samtools coverage` (samtools)
- Python: `pysam.AlignmentFile` with `pileup()` and `get_index_statistics()` (pysam)

Generate alignment statistics using samtools and pysam.

## Quick Summary Commands

| Question | Best tool | Why |
|----------|-----------|-----|
| Quick read counts by FLAG category | `samtools flagstat` | Fast; counts secondary+supp in totals |
| Per-chromosome counts | `samtools idxstats` | Fast (needs index); **counts secondary+supp** |
| Insert size, MAPQ, error, GC | `samtools stats -r ref.fa` | Comprehensive; feeds MultiQC |
| Per-position depth (small region) | `samtools depth` or pysam pileup | Slow on full genome |
| Per-position depth (genome-wide) | **`mosdepth`** | 3-10x faster than `samtools depth` |
| Per-region coverage (BED) | `mosdepth --by regions.bed` | Production default |
| Coverage histogram / cumulative | `mosdepth -t 4 --no-per-base` | Single-pass histogram |
| Breadth at depth thresholds | `mosdepth --thresholds 1,10,30,100` | Standard exome QC |
| Targeted enrichment QC | `picard CollectHsMetrics` | PCT_OFF_BAIT, FOLD_80_BASE_PENALTY, AT/GC dropout |
| Cross-sample contamination | `verifybamid2`, `somalier` | FREEMIX < 0.01 expected |

### What Each Tool Counts (and Doesn't)

| Counting category | flagstat | stats | idxstats |
|-------------------|----------|-------|----------|
| Primary alignments | `in total` minus supp | `raw total sequences` | mapped column |
| Secondary | `secondary` line | filtered out | counted in mapped |
| Supplementary | `supplementary` line | filtered out | counted in mapped |
| Mapping rate denominator | total *including* supp | primary only | mapped+unmapped |

For long-read data where one read produces many supplementary alignments, the senior cross-check:
```
input_read_count = flagstat_total - secondary - supplementary
                 = stats_raw_total_sequences
```
Reports of "the file has 1.2M reads" where the input was actually 800k with 400k supplementary chimeric splits trace to flagstat misinterpretation.

## samtools flagstat

Fast summary of alignment flags.

```bash
samtools flagstat input.bam
```

Output:
```
10000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
50000 + 0 supplementary
0 + 0 duplicates
9800000 + 0 mapped (98.00% : N/A)
9950000 + 0 paired in sequencing
4975000 + 0 read1
4975000 + 0 read2
9700000 + 0 properly paired (97.49% : N/A)
9750000 + 0 with itself and mate mapped
100000 + 0 singletons (1.01% : N/A)
25000 + 0 with mate mapped to a different chr
10000 + 0 with mate mapped to a different chr (mapQ>=5)
```

### Multi-threaded
```bash
samtools flagstat -@ 4 input.bam
```

### Output to File
```bash
samtools flagstat input.bam > flagstat.txt
```

## samtools idxstats

Per-chromosome read counts (requires index).

```bash
samtools idxstats input.bam
```

Output format: `chrom length mapped unmapped`
```
chr1    248956422    5000000    1000
chr2    242193529    4800000    800
chrM    16569        50000      100
*       0            0          150000
```

### Parse idxstats
```bash
# Total mapped reads
samtools idxstats input.bam | awk '{sum += $3} END {print sum}'

# Mitochondrial percentage
samtools idxstats input.bam | awk '
    /^chrM/ {mt = $3}
    {total += $3}
    END {print mt/total*100 "% mitochondrial"}'
```

## samtools stats

Comprehensive statistics including insert size, base quality, and more.

```bash
samtools stats input.bam > stats.txt
```

### View Summary Numbers
```bash
samtools stats input.bam | grep "^SN"
```

Key summary fields:
- `raw total sequences` - Total reads
- `reads mapped` - Mapped reads
- `reads mapped and paired` - Properly paired
- `insert size average` - Mean insert size
- `insert size standard deviation` - Insert size spread
- `average length` - Mean read length
- `error rate` - Mismatch rate

### Generate Plots (with plot-bamstats)
```bash
samtools stats input.bam > stats.txt
plot-bamstats -p plots/ stats.txt
```

### Stats for Specific Region
```bash
samtools stats input.bam chr1:1000000-2000000 > region_stats.txt
```

## samtools depth

Per-position read depth.

### Basic Depth
```bash
samtools depth input.bam > depth.txt
```

Output: `chrom position depth`

### Depth at Specific Positions
```bash
samtools depth -r chr1:1000-2000 input.bam
```

### Include Zero-Depth Positions
```bash
samtools depth -a input.bam > depth_with_zeros.txt
```

### Maximum Depth Cap (Critical Trap)
```bash
# samtools <1.13: silently caps at 8000 -- pipelines pinned to old versions are silently wrong
# samtools 1.13+: no cap by default
# Always set explicitly to be portable:
samtools depth -d 0 input.bam
```

Pipelines that routinely break the 8000 cap: targeted oncology hotspots (5000-50000x), mitochondrial DNA (small genome, large read share), amplicon viral (ARTIC: 1000-100000x per amplicon), UMI-deduped capture (14000-17000x post-collapse), highly expressed transcripts (rRNA, mt-RNA).

### Overlapping Pair Correction
```bash
# When fragment length < 2 * read_length, R1 and R2 overlap.
# Default samtools depth double-counts overlap; -s deducts:
samtools depth -s input.bam
```
Without `-s`, doubled support inflates somatic VAFs at sites covered by overlapping pairs (especially in fragmented samples: FFPE, cfDNA). `mosdepth` does not double-count overlap; `bcftools mpileup` corrects by default; `samtools mpileup` does not.

### mosdepth (Modern Default)
```bash
mosdepth -t 4 sample input.bam                                                       # genome-wide per-base
mosdepth -t 4 --by exome.bed --thresholds 1,10,20,30,100 --no-per-base sample input.bam   # exome QC
mosdepth -t 4 --quantize 0:1:10:100: sample input.bam                                # CNV-style bands
mosdepth -t 4 -f ref.fa sample input.cram                                            # CRAM with reference
```
`mosdepth` filters dup/secondary/supp by default; configurable via `--flag`. Memory ~ 4 bytes x longest chrom (1 GB for human chr1, 12+ GB for axolotl). Does not honor base quality; use `samtools depth -q INT` if needed.

### Depth from BED Regions
```bash
samtools depth -b regions.bed input.bam
```

### Calculate Mean Depth
```bash
samtools depth input.bam | awk '{sum += $3; n++} END {print sum/n}'
```

## samtools coverage

Per-chromosome or per-region coverage statistics (faster than depth).

```bash
samtools coverage input.bam
```

Output columns:
- `#rname` - Reference name
- `startpos` - Start position
- `endpos` - End position
- `numreads` - Number of reads
- `covbases` - Bases with coverage
- `coverage` - Percentage of bases covered
- `meandepth` - Mean depth
- `meanbaseq` - Mean base quality
- `meanmapq` - Mean mapping quality

### Coverage for Specific Region
```bash
samtools coverage -r chr1:1000000-2000000 input.bam
```

### Coverage from BED
```bash
samtools coverage -b regions.bed input.bam
```

### Histogram Output
```bash
samtools coverage -m input.bam
```

## pysam Python Alternative

### Count Reads
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    total = mapped = paired = proper = 0
    for read in bam:
        total += 1
        if not read.is_unmapped:
            mapped += 1
        if read.is_paired:
            paired += 1
        if read.is_proper_pair:
            proper += 1

    print(f'Total: {total}')
    print(f'Mapped: {mapped} ({mapped/total*100:.1f}%)')
    print(f'Properly paired: {proper} ({proper/paired*100:.1f}%)')
```

### Per-Chromosome Counts
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for stat in bam.get_index_statistics():
        print(f'{stat.contig}: {stat.mapped} mapped, {stat.unmapped} unmapped')
```

### Calculate Depth at Position
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup in bam.pileup('chr1', 1000000, 1000001):
        print(f'Position {pileup.pos}: depth {pileup.n}')
```

### Mean Depth in Region
```python
import pysam

def mean_depth(bam_path, chrom, start, end):
    depths = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup in bam.pileup(chrom, start, end, truncate=True):
            depths.append(pileup.n)

    if depths:
        return sum(depths) / len(depths)
    return 0

depth = mean_depth('input.bam', 'chr1', 1000000, 2000000)
print(f'Mean depth: {depth:.1f}x')
```

### Coverage Statistics

**Goal:** Compute coverage breadth and depth for a genomic region from a BAM file.

**Approach:** Iterate pileup columns in the region, count covered positions and accumulate depth, then derive percentages and means.

**Reference (pysam 0.22+):**
```python
import pysam

def coverage_stats(bam_path, chrom, start, end):
    covered = 0
    total_depth = 0

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup in bam.pileup(chrom, start, end, truncate=True):
            covered += 1
            total_depth += pileup.n

    length = end - start
    pct_covered = covered / length * 100
    mean_depth = total_depth / length if length > 0 else 0

    return {
        'length': length,
        'covered_bases': covered,
        'pct_covered': pct_covered,
        'mean_depth': mean_depth
    }

stats = coverage_stats('input.bam', 'chr1', 1000000, 2000000)
print(f'Coverage: {stats["pct_covered"]:.1f}%')
print(f'Mean depth: {stats["mean_depth"]:.1f}x')
```

### Insert Size Distribution

**Goal:** Compute the insert size distribution to assess library preparation quality.

**Approach:** Iterate properly paired read1 records, accumulate template lengths into a Counter, then compute summary statistics.

**Reference (pysam 0.22+):**
```python
import pysam
from collections import Counter

insert_sizes = Counter()

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        if read.is_proper_pair and read.is_read1 and read.template_length > 0:
            insert_sizes[read.template_length] += 1

sizes = list(insert_sizes.keys())
mean_insert = sum(s * c for s, c in insert_sizes.items()) / sum(insert_sizes.values())
print(f'Mean insert size: {mean_insert:.0f}')
print(f'Min: {min(sizes)}, Max: {max(sizes)}')
```

## Quick Reference

| Task | Command |
|------|---------|
| Quick counts | `samtools flagstat input.bam` |
| Per-chrom counts | `samtools idxstats input.bam` |
| Full stats | `samtools stats input.bam` |
| Coverage summary | `samtools coverage input.bam` |
| Per-position depth | `samtools depth input.bam` |
| Mean depth | `samtools depth input.bam \| awk '{sum+=$3;n++}END{print sum/n}'` |

## QC Thresholds Are Assay-Specific

A single "mapping rate > 95%" rule rejects valid ATAC, ChIP, RNA-seq, metagenomics, and aDNA samples. The threshold question is "is this rate normal for this assay?" not "is this rate above 95%?"

| Metric | WGS PCR-free | WGS PCR | WES | Targeted panel | Deep panel (UMI) | RNA-seq | scRNA (10x) | ATAC | ChIP | Long-read | aDNA |
|--------|--------------|---------|-----|----------------|------------------|---------|-------------|------|------|-----------|------|
| Mapping rate | >99% | >98% | >95% | >95% | >95% | >90% | >70% | >50% | >60% | >95% | 1-50% |
| Duplicate rate | <5% | 5-15% | 20-50% | 20-50% | 50-90% pre-consensus | (skip) | (use UMI) | 10-30% | 5-30% | n/a | 20-60% |
| Proper pair rate | >95% | >95% | >85% | >80% | >80% | >70% | n/a | >50% | >70% | n/a | >60% |
| Mean MAPQ | bimodal at 0/60 | bimodal | bimodal | bimodal | bimodal | bimodal incl 255 (STAR) | 0/3/255 | 30-55 | 30-55 | 30-50 | 20-40 |
| Mt fraction | 0.1-2% | 0.1-2% | <1% | <0.1% | <0.1% | varies | varies | **<10% target** | <2% | n/a | varies |

Mean MAPQ is misleading; the distribution is bimodal (0 and aligner-max). The fraction at MAPQ >= 30 is more informative:
```bash
samtools view -c -F 2308 -q 30 in.bam   # primary, mapped, MAPQ>=30
samtools view -c -F 2308 in.bam          # primary, mapped (denominator)
# For STAR/STARsolo, use -q 255 instead of -q 30 (255 is the unique-mapping sentinel)
```

## What Flagstat Does Not Reveal

A 99% flagstat mapping rate does NOT mean the data is usable. Common false-positive scenarios:

1. **Adapter readthrough**: short fragments (insert < 2 * read_length) sequence into adapter; aligners soft-clip the adapter portion and flag the read as MAPPED. Detect:
   ```bash
   samtools stats input.bam | grep "bases soft-clipped"   # >5% suggests adapter contamination
   ```
2. **Off-target enrichment** (capture/WES): detect via `picard CollectHsMetrics` PCT_OFF_BAIT or PCT_SELECTED_BASES.
3. **Low-complexity pile-up**: telomere/centromere reads mass at MAPQ-0; counted as mapped but useless. Detect via MAPQ distribution.
4. **Cross-sample contamination**: detect via `verifybamid2` or `somalier` (FREEMIX > 1% degrades somatic calling; > 5% breaks germline calling).
5. **Wrong reference build**: a BAM aligned to GRCh37 viewed against GRCh38 looks fine to flagstat but produces nonsense pileups. Compare `@SQ M5:` from BAM header with `samtools dict ref.fa` -- see alignment-validation.

## Insert Size Caveats

`samtools stats` reports the IS section only for FR-oriented properly paired reads. So:
- Mate-pair libraries (RF orientation): IS section *empty* -- proper-pair flag not set for RF
- ATAC-seq: bimodal/multimodal expected (nucleosome ladder ~50/~180/~340 bp). Unimodal suggests poor transposition.
- RNA-seq: TLEN includes intron span -- mean meaningless
- Bisulfite (PBAT): orientation reversed; samtools may not flag proper pair

## Related Skills

- sam-bam-basics - View alignment files; aligner-aware MAPQ semantics
- alignment-indexing - idxstats requires index; secondary+supp counted
- alignment-validation - Insert size by library, contamination, sample-swap detection
- duplicate-handling - Library-aware duplicate rate expectations
- alignment-filtering - Filter before stats
- sequence-io/sequence-statistics - FASTA/FASTQ statistics
