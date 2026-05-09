---
name: bio-alignment-indexing
description: Create and use BAI/CSI indices for BAM/CRAM files using samtools and pysam. Use when enabling random access to alignment files or fetching specific genomic regions.
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

# Alignment Indexing

Create indices for random access to alignment files using samtools and pysam.

**"Index a BAM file"** → Create a .bai/.csi index enabling random access to genomic regions.
- CLI: `samtools index file.bam`
- Python: `pysam.index('file.bam')`

## Index Types

| Index | Extension | Max contig | Bin shift | When required |
|-------|-----------|-----------|-----------|---------------|
| BAI | `.bai` / `.bam.bai` | 2^29-1 = ~536 Mbp | fixed (16 kb) | Default for human, mouse, fly, fish |
| CSI | `.csi` / `.bam.csi` | 2^(min_shift + depth*3) | configurable via `-m` | **Required** for any contig >536 Mbp |
| CRAI | `.crai` / `.cram.crai` | chunk-based | n/a | CRAM only |
| TBI | `.tbi` | 2^29-1 | fixed | tabix VCF/BED -- same limit as BAI |

### Which Index for Which Genome

| Genome | Largest contig | Index |
|--------|---------------|-------|
| GRCh38 / GRCh37 (human) | 248 Mbp | BAI |
| GRCm39 (mouse) | 195 Mbp | BAI |
| GRCz11 (zebrafish), TAIR10 (Arabidopsis) | 78 Mbp / 30 Mbp | BAI |
| Wheat IWGSC (Triticum aestivum) | ~830 Mbp avg | **CSI** |
| Pine, fir, axolotl, sugar pine | multi-Gbp | **CSI with larger `-m`** |
| Long-read assembly with very large contigs | varies | check `cut -f2 ref.fa.fai \| sort -nr \| head -1` |

For polyploid plants and salamander-scale genomes, increase the bin shift:
```bash
# Default CSI matches BAI bin layout: 2^(14 + 5*3) = 2^29 ≈ 512 Mbp per contig
samtools index -c file.bam

# Larger min_shift for contigs >512 Mbp (wheat, axolotl, sugar pine)
samtools index -c -m 18 file.bam   # 2^(18+15) = 2^33 = ~8.5 Gbp per contig
```

**Index file precedence trap:** when both `.bai` and `.csi` exist, samtools uses `.bai`. After re-indexing to CSI for a long contig, delete the old `.bai` or operations fail confusingly.

## samtools index

### Create BAI Index
```bash
samtools index input.bam
# Creates input.bam.bai
```

### Create CSI Index
```bash
samtools index -c input.bam
# Creates input.bam.csi
```

### Specify Output Name
```bash
samtools index input.bam output.bai
```

### Multi-threaded Indexing
```bash
samtools index -@ 4 input.bam
```

### Index CRAM
```bash
samtools index input.cram
# Creates input.cram.crai
```

## Index Requirements

Indexing requires coordinate-sorted files:
```bash
# Check sort order
samtools view -H input.bam | grep "^@HD"
# Should show SO:coordinate

# Sort if needed, then index
samtools sort -o sorted.bam input.bam
samtools index sorted.bam
```

## Using Indices for Region Access

**Goal:** Extract reads overlapping specific genomic coordinates from an indexed BAM.

**Approach:** With the index present, `samtools view` or `pysam.fetch()` can jump directly to the relevant file offset instead of scanning the entire file.

### samtools view with Region
```bash
# Requires index file present
samtools view input.bam chr1:1000000-2000000
```

### Multiple Regions
```bash
samtools view input.bam chr1:1000-2000 chr2:3000-4000
```

### Regions from BED File
```bash
samtools view -L regions.bed input.bam
```

## pysam Python Alternative

### Create Index
```python
import pysam

pysam.index('input.bam')
# Creates input.bam.bai
```

### Create CSI Index
```python
pysam.index('input.bam', 'input.bam.csi', csi=True)
```

### Fetch with Index
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    # fetch() requires index
    for read in bam.fetch('chr1', 1000000, 2000000):
        print(read.query_name)
```

### Check if Indexed
```python
import pysam
from pathlib import Path

def is_indexed(bam_path):
    bam_path = Path(bam_path)
    return (bam_path.with_suffix('.bam.bai').exists() or
            Path(str(bam_path) + '.bai').exists() or
            bam_path.with_suffix('.bam.csi').exists())

if not is_indexed('input.bam'):
    pysam.index('input.bam')
```

### Fetch Multiple Regions
```python
regions = [('chr1', 1000, 2000), ('chr1', 5000, 6000), ('chr2', 1000, 2000)]

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for chrom, start, end in regions:
        count = sum(1 for _ in bam.fetch(chrom, start, end))
        print(f'{chrom}:{start}-{end}: {count} reads')
```

### Count Reads in Region
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    count = bam.count('chr1', 1000000, 2000000)
    print(f'Reads in region: {count}')
```

### Get Reads Covering Position
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000000, 1000001):
        if read.reference_start <= 1000000 < read.reference_end:
            print(f'{read.query_name} covers position 1000000')
```

## Index File Locations

samtools looks for indices in two locations:
```
input.bam.bai   # Standard location
input.bai       # Alternative location
```

For CRAM:
```
input.cram.crai
```

## idxstats - Index Statistics

### Get Per-Chromosome Counts
```bash
samtools idxstats input.bam
```

Output format:
```
chr1    248956422    5000000    0
chr2    242193529    4500000    0
*       0            0          10000
```

Columns: reference name, length, mapped reads, unmapped reads.

### What "mapped" Actually Counts (Caveat)

The mapped column counts **every alignment record with that RNAME, including secondary AND supplementary**. For long-read minimap2 output, where a single read can produce many supplementary chimeric alignments, idxstats overcounts input reads -- typically 1.5-3x.

For unique read counts, use primary-only:
```bash
samtools view -c -F 2304 input.bam chr1   # primary only
```

Cross-check unmapped consistency (a senior sanity check):
```bash
samtools idxstats file.bam | awk '{sum+=$4} END {print sum}'   # idxstats unmapped (sum across all rows; PE orphans get a contig RNAME)
samtools view -c -f 4 -F 2304 file.bam                         # primary unmapped (should match)
```

### Sum Total Mapped Reads
```bash
samtools idxstats input.bam | awk '{sum += $3} END {print sum}'
```

### pysam idxstats
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for stat in bam.get_index_statistics():
        print(f'{stat.contig}: {stat.mapped} mapped, {stat.unmapped} unmapped')
```

## FASTA Index (faidx)

Related but different - index reference FASTA for random access:

```bash
samtools faidx reference.fa
# Creates reference.fa.fai

# Fetch region from indexed FASTA
samtools faidx reference.fa chr1:1000-2000
```

### pysam FastaFile
```python
with pysam.FastaFile('reference.fa') as ref:
    seq = ref.fetch('chr1', 1000, 2000)
    print(seq)
```

## Quick Reference

| Task | samtools | pysam |
|------|----------|-------|
| Create BAI | `samtools index file.bam` | `pysam.index('file.bam')` |
| Create CSI | `samtools index -c file.bam` | `pysam.index('file.bam', csi=True)` |
| Fetch region | `samtools view file.bam chr1:1-1000` | `bam.fetch('chr1', 0, 1000)` |
| Count in region | `samtools view -c file.bam chr1:1-1000` | `bam.count('chr1', 0, 1000)` |
| Index stats | `samtools idxstats file.bam` | `bam.get_index_statistics()` |
| Index FASTA | `samtools faidx ref.fa` | Automatic with FastaFile |

## Index Staleness

If the BAM was modified after indexing, the index points to wrong file offsets and region queries return wrong (or zero) reads. Quick check:
```bash
if [ input.bam -nt input.bam.bai ]; then
    echo "Index older than BAM; re-indexing"
    samtools index input.bam
fi
```

## Contig-Naming Sanity Check

A leading cause of "my variant calling produced empty VCFs" tickets: querying `chrM` against a BAM that uses `MT` (or `chr1` vs `1`). Always inspect contig conventions before region queries:
```bash
samtools view -H input.bam | grep '^@SQ' | head -3
# Compare with reference dict:
samtools dict ref.fa | head -3
```

UCSC convention uses `chr1`/`chrM`; Ensembl/NCBI uses `1`/`MT`. The two are not interchangeable; tools fail with "contig not found" or silently return zero reads.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `random alignment retrieval only works for indexed BAM` | Missing index | Run `samtools index file.bam` |
| `file is not sorted` | Unsorted BAM | Sort first with `samtools sort` |
| `chromosome not found` | Wrong chromosome name | Check names with `samtools view -H` |
| Region query returns zero reads on a known-populated locus | Stale BAI / `chr` vs no-`chr` mismatch | Re-index; verify naming convention |
| BAI silently truncates reads on contigs >536 Mbp | Plant / amphibian / amplified genome | Use CSI: `samtools index -c file.bam` |

## Related Skills

- sam-bam-basics - View and convert alignment files
- alignment-sorting - Sort BAM files (required before indexing)
- alignment-filtering - Filter by regions using index
- bam-statistics - Use idxstats for quick counts
- sequence-io/read-sequences - Index FASTA with SeqIO.index_db()
