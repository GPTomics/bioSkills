---
name: bio-alignment-sorting
description: Sort alignment files by coordinate or read name using samtools and pysam. Use when preparing BAM files for indexing, variant calling, or paired-end analysis.
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

# Alignment Sorting

Sort alignment files by coordinate or read name using samtools and pysam.

**"Sort a BAM file"** → Reorder reads by genomic coordinate (for indexing/variant calling) or by name (for paired-end processing).
- CLI: `samtools sort -o sorted.bam input.bam`
- Python: `pysam.sort('-o', 'sorted.bam', 'input.bam')`

## Sort Orders

| Order | Flag | Use Case |
|-------|------|----------|
| Coordinate | default | Indexing, visualization, variant calling |
| Name | `-n` | Paired-end processing, fixmate, markdup |
| Tag | `-t TAG` | Sort by specific tag value |

## samtools sort

### Sort by Coordinate (Default)
```bash
samtools sort -o sorted.bam input.bam
```

### Sort by Read Name
```bash
samtools sort -n -o namesorted.bam input.bam
```

### Multi-threaded Sorting
```bash
samtools sort -@ 8 -o sorted.bam input.bam
```

### Control Memory Usage
```bash
samtools sort -m 4G -@ 4 -o sorted.bam input.bam
```

### Set Temporary Directory
```bash
samtools sort -T /tmp/sort_tmp -o sorted.bam input.bam
```

### Specify Output Format
```bash
# Output as BAM (default)
samtools sort -O bam -o sorted.bam input.bam

# Output as CRAM
samtools sort -O cram --reference ref.fa -o sorted.cram input.bam
```

### Sort by Tag
```bash
# Sort by cell barcode (10x Genomics)
samtools sort -t CB -o sorted_by_barcode.bam input.bam
```

### Pipe from Aligner
```bash
bwa mem ref.fa reads.fq | samtools sort -o aligned.bam
```

## samtools collate vs sort -n

| Tool | Algorithm | Speed | Memory | Output guarantee |
|------|-----------|-------|--------|------------------|
| `sort -n` | Full lexicographic sort by QNAME | Slowest | Spills to `-T` | Strict total order by name |
| `collate` | Hash-bucket grouping | ~3-10x faster | Bounded | Mates adjacent; between-mate order undefined |

Use `collate` when extracting paired FASTQ, re-aligning, or streaming through markdup. Use `sort -n` only when a tool requires true lexicographic name order (e.g. RSEM, Salmon alignment-mode).

```bash
# Fast paired FASTQ extraction
samtools collate -O -u in.bam tmp_prefix | \
    samtools fastq -1 R1.fq.gz -2 R2.fq.gz -0 /dev/null -s /dev/null -n -

# Markdup pre-processing (collate beats sort -n here)
samtools collate -O -u in.bam tmp_prefix | \
    samtools fixmate -m -u - - | \
    samtools sort -u - | \
    samtools markdup - out.bam
```

### Sort Order Required by Downstream Tool

| Operation | Required sort |
|-----------|---------------|
| `samtools index` | coordinate (hard requirement) |
| `samtools fixmate -m` | name (or collate; needs mates adjacent) |
| `samtools markdup` | coordinate (after fixmate) |
| GATK MarkDuplicatesSpark | coordinate or queryname |
| `samtools mpileup` / `bcftools mpileup` | coordinate |
| GATK HaplotypeCaller, Mutect2 | coordinate |
| featureCounts / HTSeq | coordinate or name (`-p` for paired) |
| umi_tools dedup | coordinate (with index) |
| fgbio GroupReadsByUmi | queryname (hard requirement) |
| fgbio CallMolecularConsensusReads | TemplateCoordinate (`fgbio SortBam`) |
| Sniffles, cuteSV, Manta, Delly | coordinate (need SA tags) |
| Salmon alignment-mode | name |
| RSEM (with STAR `--quantMode TranscriptomeSAM`) | name (hard requirement) |

## Check Sort Order

### From Header
```bash
samtools view -H input.bam | grep "^@HD"
# SO:coordinate = coordinate sorted
# SO:queryname = name sorted
# SO:unsorted = not sorted
```

### Verify Sorted
```bash
# Check if coordinate sorted (returns 0 if sorted)
samtools view input.bam | awk '$4 < prev {exit 1} {prev=$4}'
```

## pysam Python Alternative

### Sort with pysam
```python
import pysam

pysam.sort('-o', 'sorted.bam', 'input.bam')
```

### Sort by Name
```python
pysam.sort('-n', '-o', 'namesorted.bam', 'input.bam')
```

### Sort with Options
```python
pysam.sort('-@', '4', '-m', '2G', '-o', 'sorted.bam', 'input.bam')
```

### Avoid In-Python Sorting

Do not load BAM records into a list and call `sorted()`. `pysam.sort()` calls samtools' external-merge sort which spills to disk; loading reads into memory blows up around ~30M reads (~10 GB human BAM). Always delegate to `pysam.sort()`:
```python
import pysam

pysam.sort('-@', '4', '-m', '2G', '-T', '/tmp/sortpfx',
           '-o', 'sorted.bam', 'input.bam')
```

### Check Sort Order in pysam
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    hd = bam.header.get('HD', {})
    sort_order = hd.get('SO', 'unknown')
    print(f'Sort order: {sort_order}')
```

### Stream Sort from Aligner
For streaming from aligners, use shell pipes (simpler and more reliable):
```python
import subprocess

subprocess.run(
    'bwa mem ref.fa reads.fq | samtools sort -o aligned.bam',
    shell=True, check=True
)
```

## samtools merge

Combine multiple BAM files into one. `samtools merge` does NOT validate sort-order consistency across inputs; mismatched inputs silently produce a malformed output.

### Verify Sort Order Consistency First
```bash
for f in *.bam; do samtools view -H "$f" | head -1; done | sort -u
# Should print exactly ONE line, e.g. "@HD VN:1.6 SO:coordinate"
```

### Safe Merge (dedup @PG/@CO and @RG)
```bash
# -c dedups @PG and @CO; -p dedups @RG ID collisions
samtools merge -c -p -@ 8 merged.bam sample1.bam sample2.bam sample3.bam
```

When merging BAMs from different lanes / machines / aligners, RG IDs may collide. `-c` and `-p` deduplicate header records, but RG IDs that genuinely refer to different lane-level read groups must be made unique upstream (`samtools addreplacerg`) before merge -- otherwise GATK BQSR (which keys models by RGID/PU) silently produces wrong recalibration.

### Merge with Threads / from File List
```bash
samtools merge -@ 4 merged.bam sample1.bam sample2.bam sample3.bam
samtools merge -b files.txt merged.bam   # one BAM path per line
```

### Force Overwrite
```bash
samtools merge -f merged.bam sample1.bam sample2.bam
```

### Merge Specific Region
```bash
samtools merge -R chr1:1000000-2000000 merged_region.bam sample1.bam sample2.bam
```

### pysam Merge
```python
import pysam

pysam.merge('-c', '-p', '-f', 'merged.bam', 'sample1.bam', 'sample2.bam', 'sample3.bam')
```

## Common Workflows

**Goal:** Combine sorting with other alignment processing steps into efficient pipelines.

**Approach:** Pipe aligner output directly into `samtools sort` to avoid writing unsorted intermediates, then index for downstream access.

### Align and Sort
```bash
bwa mem -t 8 ref.fa R1.fq R2.fq | samtools sort -@ 4 -o aligned.bam
samtools index aligned.bam
```

### Re-sort by Name for Duplicate Marking
```bash
# Full workflow: sort by name, fixmate, sort by coord, markdup
samtools sort -n -o namesorted.bam input.bam
samtools fixmate -m namesorted.bam fixmate.bam
samtools sort -o sorted.bam fixmate.bam
samtools markdup sorted.bam marked.bam
```

### Convert Name-sorted to Coordinate-sorted
```bash
samtools sort -o coord_sorted.bam name_sorted.bam
samtools index coord_sorted.bam
```

### Extract FASTQ from Sorted BAM
```bash
# Collate first to group pairs
samtools collate -u -O input.bam /tmp/collate | \
    samtools fastq -1 R1.fq -2 R2.fq -0 /dev/null -s /dev/null -
```

## Performance Tips

| Parameter | Effect |
|-----------|--------|
| `-@ N` | Use N additional threads |
| `-m SIZE` | Memory per thread (e.g., 4G) |
| `-T PREFIX` | Temp file location (use fast SSD scratch) |
| `-l LEVEL` | Compression level (1-9, default 6) |

### Compression Level Decision

| Level | Use | Wall-time vs default | Size vs default |
|-------|-----|----------------------|------------------|
| `-l 0` / `-u` | Pipe between samtools tools | 0% (skips BGZF) | +200-400% |
| `-l 1` | Final output if disk is cheap | ~+10% | ~+30% |
| `-l 6` | Default | baseline | baseline |
| `-l 9` | Archival, write-once | ~+50-100% | ~-2-5% |

```bash
# WRONG -- pipe re-compresses then decompresses every step
samtools fixmate -m in.bam - | samtools sort -o out.bam

# RIGHT -- uncompressed (-u) between piped samtools commands
samtools fixmate -m -u in.bam - | samtools sort -o out.bam
```

### Optimal Settings for Large Files
```bash
# 8 threads, 2GB per thread, low compression for output written to fast disk
samtools sort -@ 8 -m 2G -l 1 -T /scratch/sortpfx -o sorted.bam input.bam
```

## Quick Reference

| Task | Command |
|------|---------|
| Sort by coordinate | `samtools sort -o out.bam in.bam` |
| Sort by name | `samtools sort -n -o out.bam in.bam` |
| Sort with threads | `samtools sort -@ 8 -o out.bam in.bam` |
| Collate pairs | `samtools collate -o out.bam in.bam` |
| Merge BAMs | `samtools merge out.bam in1.bam in2.bam` |
| Check sort order | `samtools view -H in.bam \| grep "^@HD"` |
| Sort + index | `samtools sort -o out.bam in.bam && samtools index out.bam` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `out of memory` | Insufficient RAM | Use `-m` to limit per-thread memory |
| `disk full` | Temp files filling disk | Use `-T` to specify different location |
| `truncated file` | Interrupted sort | Re-run sort from original |

## Related Skills

- sam-bam-basics - View and convert alignment files
- alignment-indexing - Index after coordinate sorting
- duplicate-handling - Requires name-sorted input for fixmate
- alignment-filtering - Filter before or after sorting
