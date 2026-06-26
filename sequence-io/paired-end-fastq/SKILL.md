---
name: bio-paired-end-fastq
description: Handle paired-end FASTQ files (R1/R2) using Biopython while keeping mates synchronized. Use when working with Illumina paired reads, synchronizing pairs, filtering both mates together with orphan routing, interleaving/deinterleaving, or matching mates by read name.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Paired-End FASTQ

**"Work with my paired-end FASTQ files"** -> Iterate R1/R2 pairs in sync, filter both mates together (routing orphans out), interleave/deinterleave files, and match mates by read name.
- Python: `SeqIO.parse()` with `zip()` iteration (BioPython)

## The Governing Principle: R1 and R2 Are Parallel Streams

Aligners (bwa mem, bowtie2, STAR) consume R1 and R2 as two parallel streams and pair the i-th record of each file: same order, same count. They assume the k-th read in R1 is the mate of the k-th read in R2.

This makes independent per-mate processing the #1 paired-end correctness trap. Filtering or trimming ONE mate without the other DESYNCS the files:
- Best case: the aligner detects a read-name mismatch and crashes loudly.
- Worst case: it silently pairs the wrong R1 with the wrong R2 -> mismapping, corrupt insert sizes, no error at all.

Governing rule: never filter, trim, sort, or subsample one mate independently. Process both mates as a unit. When a read fails but its mate passes, route the survivor to a separate singleton/orphan file rather than leaving a gap that desyncs the stream. Proper paired trimmers (Trimmomatic PE, fastp, cutadapt `-p`) do exactly this: synchronized paired output plus separate orphan files.

## Required Import

```python
from Bio import SeqIO
```

## Read-Name Conventions: How Mates Are Matched

Two distinct naming layers exist. File naming (which file is R1 vs R2) is separate from read naming (how a tool decides two records are mates).

### File naming patterns
- `sample_R1.fastq` / `sample_R2.fastq`
- `sample_1.fastq` / `sample_2.fastq`
- `sample_R1_001.fastq` / `sample_R2_001.fastq` (Illumina bcl2fastq)
- `sample.R1.fastq.gz` / `sample.R2.fastq.gz`

### Read naming: mate matched by shared ID up to the first whitespace

| Era | Mate marker | Example | How mates match |
|-----|-------------|---------|-----------------|
| Pre-CASAVA 1.8 | SUFFIX `/1`, `/2` on the read name | `@HWUSI-EAS100R:6:73:941:1973#0/1` | Strip the trailing `/1`/`/2`; the rest is the shared ID |
| CASAVA 1.8+ | SECOND field after a SPACE | `@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG` | The text before the space is IDENTICAL for both mates; the `1:`/`2:` lives only after the space |

In CASAVA 1.8+ the second field is `<read>:<is_filtered>:<control>:<index>` -> `read`=1 or 2 (mate number), `is_filtered`=Y (failed chastity) or N, `control`=0 normally, `index`=barcode.

Most tools (and Biopython) take the read ID as everything up to the first whitespace. Biopython puts that token in `record.id` and the full header line in `record.description`. So for 1.8+ data, `r1.id == r2.id` directly; the `1:`/`2:` distinction is only visible in `record.description`. For pre-1.8 data, strip the `/1`/`/2` suffix before comparing.

A normalizer that handles both eras:

```python
def mate_key(record):
    return record.id.rsplit('/', 1)[0]
```

`record.id` already excludes anything after the first space, so this single rsplit covers both the space-format (1.8+) and the slash-suffix (pre-1.8) conventions.

## Iterate Pairs Together

### Basic Paired Iteration
```python
r1_records = SeqIO.parse('reads_R1.fastq', 'fastq')
r2_records = SeqIO.parse('reads_R2.fastq', 'fastq')

for r1, r2 in zip(r1_records, r2_records):
    print(f'R1: {r1.id}, R2: {r2.id}')
    print(f'Lengths: {len(r1.seq)}, {len(r2.seq)}')
```

`zip` stops at the shorter iterator. If R1 has 1000 reads and R2 has 998, `zip` silently processes 998 and drops the tail with no warning. Verify counts match (see Paired Statistics) before trusting a `zip` loop on files of unknown provenance.

### Verify Pair Matching
```python
def iterate_pairs(r1_file, r2_file, format='fastq'):
    r1_records = SeqIO.parse(r1_file, format)
    r2_records = SeqIO.parse(r2_file, format)

    for r1, r2 in zip(r1_records, r2_records):
        if mate_key(r1) != mate_key(r2):
            raise ValueError(f'Pair mismatch: {r1.id} vs {r2.id}')
        yield r1, r2

for r1, r2 in iterate_pairs('reads_R1.fastq', 'reads_R2.fastq'):
    process_pair(r1, r2)
```

## Filter Pairs Together (Synchronized, With Orphan Routing)

**Goal:** Quality-filter paired reads so that R1 and R2 stay in lockstep, and reads whose mate was discarded are routed to orphan files instead of silently desyncing the stream.

**Approach:** Stream both files together with `zip`. Evaluate both mates. If both pass, write to the paired outputs. If exactly one passes, write the survivor to its orphan file. This mirrors the four-output behavior of Trimmomatic PE (paired R1, paired R2, orphan R1, orphan R2).

**Reference (BioPython 1.83+):**
```python
def filter_pairs_synced(r1_in, r2_in, r1_out, r2_out, r1_orphan, r2_orphan, min_qual=25):
    '''Keep a pair only if both mates pass; route lone survivors to orphan files.'''
    r1_records = SeqIO.parse(r1_in, 'fastq')
    r2_records = SeqIO.parse(r2_in, 'fastq')

    counts = {'paired': 0, 'r1_orphan': 0, 'r2_orphan': 0}
    with open(r1_out, 'w') as p1, open(r2_out, 'w') as p2, \
         open(r1_orphan, 'w') as o1, open(r2_orphan, 'w') as o2:
        for r1, r2 in zip(r1_records, r2_records):
            r1_ok = sum(r1.letter_annotations['phred_quality']) / len(r1.seq) >= min_qual
            r2_ok = sum(r2.letter_annotations['phred_quality']) / len(r2.seq) >= min_qual
            if r1_ok and r2_ok:
                SeqIO.write(r1, p1, 'fastq')
                SeqIO.write(r2, p2, 'fastq')
                counts['paired'] += 1
            elif r1_ok:
                SeqIO.write(r1, o1, 'fastq')
                counts['r1_orphan'] += 1
            elif r2_ok:
                SeqIO.write(r2, o2, 'fastq')
                counts['r2_orphan'] += 1
    return counts
```

The paired outputs stay synchronized because every pass-both write goes to BOTH files in the same iteration. Mean quality is one criterion; swap the test for a length threshold, adapter check, or any predicate, but always apply it to both mates and route orphans the same way. `min_qual=25` is a common Q-score cutoff (Phred 25 ~= 0.3% error); tune per experiment.

## Interleave Pairs

Interleaved FASTQ holds both mates in one file alternating R1, R2, R1, R2 (record 2k = forward, 2k+1 = reverse). `bwa mem -p` reads this format. The strict alternation IS the pairing, so a single mate-less read shifts every downstream record by one and desyncs everything. Only interleave files that are known to be synchronized, and pull orphans out first.

### Create Interleaved File

**Goal:** Merge synchronized R1/R2 files into one interleaved file.

**Approach:** Zip both iterators and yield alternating records through a generator so nothing is materialized in memory.

**Reference (BioPython 1.83+):**
```python
def interleave_pairs(r1_file, r2_file, output_file, format='fastq'):
    r1_records = SeqIO.parse(r1_file, format)
    r2_records = SeqIO.parse(r2_file, format)

    def interleaved():
        for r1, r2 in zip(r1_records, r2_records):
            yield r1
            yield r2

    count = SeqIO.write(interleaved(), output_file, format)
    return count // 2  # Number of pairs

pairs = interleave_pairs('reads_R1.fastq', 'reads_R2.fastq', 'reads_interleaved.fastq')
```

## Deinterleave

### Split Interleaved to Paired Files

**Goal:** Recover separate R1/R2 files from an interleaved file, streaming to avoid loading everything.

**Approach:** Parse once; route even-indexed records to R1, odd-indexed to R2.

**Reference (BioPython 1.83+):**
```python
def deinterleave_streaming(interleaved_file, r1_file, r2_file, format='fastq'):
    records = SeqIO.parse(interleaved_file, format)

    pairs = 0
    with open(r1_file, 'w') as r1_h, open(r2_file, 'w') as r2_h:
        for i, record in enumerate(records):
            if i % 2 == 0:
                SeqIO.write(record, r1_h, format)
            else:
                SeqIO.write(record, r2_h, format)
                pairs += 1
    return pairs
```

Even/odd splitting only stays correct if the interleaved file has perfect alternation. If an upstream per-read filter left an orphan in the file, every record after it lands in the wrong output. Guard by checking `mate_key` equality between each even/odd pair after splitting, or deinterleave with a name check.

## Paired Statistics

### Count and Verify Pairs
```python
def paired_stats(r1_file, r2_file):
    r1_count = sum(1 for _ in SeqIO.parse(r1_file, 'fastq'))
    r2_count = sum(1 for _ in SeqIO.parse(r2_file, 'fastq'))

    if r1_count != r2_count:
        print(f'WARNING: Unequal counts! R1={r1_count}, R2={r2_count} -> files are desynced')
    else:
        print(f'Pairs: {r1_count}, total reads: {r1_count * 2}')
    return r1_count, r2_count
```

### Paired Quality Summary
```python
def paired_quality_summary(r1_file, r2_file):
    r1_quals, r2_quals = [], []
    for r1, r2 in zip(SeqIO.parse(r1_file, 'fastq'), SeqIO.parse(r2_file, 'fastq')):
        r1_quals.append(sum(r1.letter_annotations['phred_quality']) / len(r1.seq))
        r2_quals.append(sum(r2.letter_annotations['phred_quality']) / len(r2.seq))
    print(f'R1 mean quality: {sum(r1_quals)/len(r1_quals):.1f}')
    print(f'R2 mean quality: {sum(r2_quals)/len(r2_quals):.1f}')
```

R2 commonly shows lower mean quality than R1 (the reverse read is sequenced later in the run); a modest R1/R2 gap is expected, not a defect.

## Find Paired Files

### Auto-Detect R2 from R1
```python
from pathlib import Path

def find_r2(r1_path):
    r1_path = Path(r1_path)
    name = r1_path.name
    patterns = [('_R1', '_R2'), ('_R1_', '_R2_'), ('.R1.', '.R2.'), ('_1', '_2')]

    for p1, p2 in patterns:
        if p1 in name:
            r2_path = r1_path.parent / name.replace(p1, p2, 1)
            if r2_path.exists():
                return r2_path
    return None
```

Order matters: test the specific `_R1` patterns before the bare `_1`, otherwise `sample_R1.fastq` would match `_1` and produce `sample_R2.fastq` only by luck. `replace(..., 1)` replaces the first occurrence only, so a sample named `sample_R1_lane_R1.fastq` swaps just the first token.

## Compressed Paired Files
```python
import gzip

def iterate_gzipped_pairs(r1_gz, r2_gz):
    with gzip.open(r1_gz, 'rt') as r1_h, gzip.open(r2_gz, 'rt') as r2_h:
        for r1, r2 in zip(SeqIO.parse(r1_h, 'fastq'), SeqIO.parse(r2_h, 'fastq')):
            yield r1, r2
```

Use text mode `'rt'`, not `'rb'`, when handing a gzip handle to `SeqIO.parse`; the parser expects decoded text.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Aligner reports "mismatched read names" or "unpaired reads" | R1 and R2 desynced by filtering/trimming one mate independently | Always filter both mates together; route orphans to separate files (see synchronized filter) |
| Silent mismapping, nonsensical insert sizes, no error | A per-mate operation dropped reads from one file -> i-th records no longer mates | Re-pair from source; never trust outputs from independent per-mate filtering |
| Mates never recognized as pairs | Mixing pre-1.8 `/1` `/2` data with 1.8+ space-format ids, or comparing full descriptions instead of the pre-space ID | Match on `mate_key` (ID up to first space, `/1`/`/2` stripped), not the whole header |
| Deinterleave produces shifted/wrong pairs | An orphan in the interleaved file broke the strict R1,R2 alternation | Remove orphans before interleaving; verify each even/odd pair with `mate_key` after splitting |
| `zip` loop processes fewer reads than expected | R1 and R2 have unequal counts; `zip` stops at the shorter and silently drops the tail | Run `paired_stats` first; counts must be equal |
| Memory error on large files | `list(SeqIO.parse(...))` materializes every record | Stream with generators; for random access use `SeqIO.index`/`index_db` |

## Related Skills

- read-sequences - Parse individual FASTQ files and choose parse vs index
- fastq-quality - Phred encoding and quality interpretation before paired filtering
- filter-sequences - Single-file filtering criteria (apply to both mates here)
- compressed-files - gzip vs BGZF handling for paired files
- read-qc/quality-reports - FastQC/MultiQC per-mate quality assessment
- alignment-files/sam-bam-basics - After filtering, align paired reads with bwa mem; proper pairs in BAM
