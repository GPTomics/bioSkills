---
name: bio-fastq-quality
description: Work with FASTQ quality scores using Biopython - access Phred scores, filter and trim by quality, compute per-position profiles, and convert between Sanger/Phred+33, Solexa, and Illumina/Phred+64 encodings. Use when analyzing read quality, filtering or trimming low-quality bases, generating quality reports, or deciding which FASTQ quality encoding a file uses before parsing.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# FASTQ Quality Scores

**"Filter my FASTQ reads by quality score"** -> Access, analyze, and filter per-base quality scores, trim low-quality bases, and generate per-position quality profiles.
- Python: `SeqIO.parse()` with `record.letter_annotations['phred_quality']` (BioPython)
- CLI alternative: `pysam.FastxFile` (`.get_quality_array()` returns offset-removed Phred ints, but ALWAYS subtracts 33 - it cannot read Phred+64/Solexa correctly, so use it only on confirmed Phred+33 data; for legacy encodings stay on `SeqIO` with the explicit variant string)

## The Governing Principle: Never Guess the Offset

A FASTQ file does not record which quality encoding it uses. The same ASCII byte means different Phred scores under different encodings, and the offsets (33 vs 64) differ by exactly 31. Choosing the wrong format string has two failure modes:

- LOUD (safe): a quality character lies outside the chosen parser's legal range -> `ValueError` naming the wrong QualityIO parser. The run stops.
- SILENT (dangerous): every character lies in the ASCII overlap region legal for both encodings -> no error, and every score is off by exactly 31. Reading Phred+33 data as `fastq-illumina` makes all scores 31 too LOW; reading Phred+64 data as `fastq-sanger` makes them 31 too HIGH. QC, filtering, and trimming silently operate on garbage scores.

Auto-detection is provably ambiguous: ASCII >= 64 is legal in every variant, so a high-quality Sanger file (all Q >= 31) and a low-quality Illumina-1.3 file can be byte-identical in their quality lines. There is no reliable way to detect the encoding from content alone (Biopython docs: this "cannot be detected reliably automatically"). Determine the encoding from the sequencing instrument and run metadata, not by guessing. Scanning for the minimum ASCII byte can only RULE OUT encodings (see "Ruling Out Encodings" below), never confirm one.

## The Four FASTQ Encodings

| Variant | Score type | Offset | Format string | ASCII chars | Q range |
|---------|-----------|--------|---------------|-------------|---------|
| Sanger / Phred+33 | Phred | 33 | `'fastq'` / `'fastq-sanger'` | `!`(33)..`~`(126) | 0..93 |
| Solexa / Illumina 1.0 (Solexa+64) | Solexa odds | 64 | `'fastq-solexa'` | `;`(59)..`~`(126) | -5..62 |
| Illumina 1.3+ (Phred+64) | Phred | 64 | `'fastq-illumina'` | `@`(64)..`~`(126) | 0..62 |
| Illumina 1.5-1.7 (Phred+64, B-tail) | Phred | 64 | `'fastq-illumina'` | `B`(66)..`~`(126) | 2..62 |
| Illumina 1.8+ (Phred+33) | Phred | 33 | `'fastq'` / `'fastq-sanger'` | `!`(33)..~`J`(74) | 0..~41 |

Almost all data produced since 2011 is Phred+33 (`'fastq'`). Phred+64 and Solexa appear only in legacy datasets, but the cost of misreading them is silent corruption, so the encoding must be confirmed before parsing legacy files.

`'fastq'` is an alias for `'fastq-sanger'`; both are Phred+33. The wrong choice produces a LOUD `ValueError` only when an out-of-range character appears, and SILENT 31-shifted scores otherwise.

## Phred vs Solexa: Two Different Score Definitions

The Solexa encoding is not just a different offset; it uses a different score formula, which is why it needs a separate parser.

- Phred: Q = -10 * log10(P_error). Always >= 0.
- Solexa: Q = -10 * log10(P/(1 - P)) - an ODDS score. It goes NEGATIVE when P > 0.5 (floor -5), which is why Solexa quality strings include ASCII 59-63.

The two scales are asymptotically equal at high quality (rounded scores above ~Q10-13 are interchangeable) but diverge for poor-quality bases. The round trip Phred -> Solexa -> Phred is LOSSY in that low-quality region: Cock et al. (2010) note that Solexa scores 9 and 10 both map to Phred 10. Do not convert legacy Solexa data to Phred and back if the low-Q values matter.

## Accessing Quality Scores

Quality scores live in `record.letter_annotations['phred_quality']` as a list of ints. The attribute is `letter_annotations` (NOT `per_letter_annotations`, which does not exist). Solexa data parsed with `'fastq-solexa'` stores `record.letter_annotations['solexa_quality']` instead, and those values can be negative.

```python
from Bio import SeqIO

for record in SeqIO.parse('reads.fastq', 'fastq'):
    quals = record.letter_annotations['phred_quality']
    print(record.id, quals[:10])
```

`letter_annotations` is length-locked to `len(record.seq)`: assigning a list of the wrong length raises. To edit sequence and quality together, slice the record (slicing keeps qualities in sync) or build a fresh record.

| Phred Score | Error Probability | Accuracy |
|-------------|-------------------|----------|
| 10 | 1 in 10 | 90% |
| 20 | 1 in 100 | 99% |
| 30 | 1 in 1000 | 99.9% |
| 40 | 1 in 10000 | 99.99% |

## Code Patterns

### Calculate Average Quality per Read
```python
for record in SeqIO.parse('reads.fastq', 'fastq'):
    quals = record.letter_annotations['phred_quality']
    print(f'{record.id}: {sum(quals) / len(quals):.1f}')
```

### Filter Reads by Mean Quality
```python
def high_quality_reads(records, min_avg_qual=20):
    for record in records:
        quals = record.letter_annotations['phred_quality']
        if sum(quals) / len(quals) >= min_avg_qual:
            yield record

records = SeqIO.parse('reads.fastq', 'fastq')
SeqIO.write(high_quality_reads(records, 25), 'filtered.fastq', 'fastq')
```

### Filter by Minimum Quality at Any Position
```python
def all_bases_above(records, min_qual=20):
    for record in records:
        if min(record.letter_annotations['phred_quality']) >= min_qual:
            yield record
```

### Trim Low-Quality 3' End

**Goal:** Drop trailing bases below a quality cutoff while keeping qualities aligned to the trimmed sequence.

**Approach:** Walk inward from the 3' end to the first base that meets the cutoff, then slice the record; slicing a SeqRecord trims `letter_annotations` in step with the sequence.

**Reference (BioPython 1.83+):**
```python
def trim_low_quality(record, min_qual=20):
    quals = record.letter_annotations['phred_quality']
    trim_pos = len(quals)
    for i in range(len(quals) - 1, -1, -1):
        if quals[i] >= min_qual:
            trim_pos = i + 1
            break
    return record[:trim_pos]

records = SeqIO.parse('reads.fastq', 'fastq')
SeqIO.write((trim_low_quality(r) for r in records), 'trimmed.fastq', 'fastq')
```

### Sliding Window Quality Trim

**Goal:** Truncate a read at the first position where average quality in a sliding window drops below a threshold (the Trimmomatic SLIDINGWINDOW model).

**Approach:** Slide a fixed-size window across the quality list; when the window mean falls below the cutoff, slice the record at that position.

**Reference (BioPython 1.83+):**
```python
def sliding_window_trim(record, window_size=5, min_avg_qual=20):
    quals = record.letter_annotations['phred_quality']
    for i in range(len(quals) - window_size + 1):
        if sum(quals[i:i + window_size]) / window_size < min_avg_qual:
            return record[:i] if i > 0 else None
    return record
```

### Per-Position Quality Profile

**Goal:** Compute mean quality at each read position to spot systematic drops (typically 3' degradation).

**Approach:** Accumulate scores by position across reads, then average each position. NovaSeq binning (see below) makes per-position values cluster at a few discrete levels - expected, not a defect.

**Reference (BioPython 1.83+):**
```python
from collections import defaultdict

position_quals = defaultdict(list)
for record in SeqIO.parse('reads.fastq', 'fastq'):
    for i, q in enumerate(record.letter_annotations['phred_quality']):
        position_quals[i].append(q)

for pos in sorted(position_quals)[:20]:
    quals = position_quals[pos]
    print(f'Position {pos}: mean={sum(quals) / len(quals):.1f}')
```

### Count Reads by Quality Threshold
```python
thresholds = [20, 25, 30, 35]
counts = {t: 0 for t in thresholds}
for record in SeqIO.parse('reads.fastq', 'fastq'):
    avg = sum(record.letter_annotations['phred_quality']) / len(record.seq)
    for t in thresholds:
        if avg >= t:
            counts[t] += 1
```

## The Illumina 1.5-1.7 B-Tail

In Illumina 1.5-1.7 (Phred+64) data, Q0 and Q1 are reserved, and ASCII `B` (Q2) at the 3' end is a Read Segment Quality Control Indicator, NOT a real Q2 measurement. A run of trailing `B`s marks a region the instrument deemed unreliable. A trimmer that treats `B` as literal Q2 keeps those junk bases instead of removing them. When trimming legacy Phred+64 data, drop trailing `B`/Q2 runs as flags rather than scores.

## NovaSeq / NextSeq Quality Binning

Modern Illumina instruments quantize quality on-instrument (RTA software, baked into the BCL), so the binned values arrive in the FASTQ - they are not introduced downstream. NovaSeq 6000 (RTA3) emits only four values: Q2, Q12, Q23, Q37. NovaSeq X / X Plus (RTA4, XLEAP-SBS chemistry) uses a different, software-version-dependent bin set whose high bins shifted to roughly Q9/Q24/Q40 (the exact ranges depend on the Control/RTA software version), so its spike values differ from the 6000 - confirm them against the run's instrument and software version rather than assuming the 6000 set. Consequences:

- Per-base quality histograms collapse to spikes at the bin values. This is expected; it is not a data problem.
- Mean quality stays meaningful (each bin approximates the mean of its input range).
- GATK BQSR interacts with binning: with only four input levels, recalibration tables are coarse and corrections are blunter than on unbinned data.

## Converting Between Encodings

`SeqIO.convert` (or parse + write) re-encodes legacy data to standard Phred+33. Specify the SOURCE encoding explicitly; an out-of-range character raises, but overlap-region characters convert silently with the wrong offset if the source is mislabeled.

```python
from Bio import SeqIO

SeqIO.convert('old_illumina.fastq', 'fastq-illumina', 'standard.fastq', 'fastq')
SeqIO.convert('solexa.fastq', 'fastq-solexa', 'standard.fastq', 'fastq')
```

Per-score conversion helpers return floats:

```python
from Bio.SeqIO.QualityIO import phred_quality_from_solexa, solexa_quality_from_phred

phred_quality_from_solexa(10)   # Solexa -> Phred (float)
solexa_quality_from_phred(30)   # Phred -> Solexa (float)
```

Writing `'fastq-solexa'` from a Phred-only record forces a lossy on-the-fly conversion and emits a `BiopythonWarning` when `max(qualities) >= 62.5`. There is no clean Phred-to-Solexa write path that avoids the lossy step, so keep modern data in Phred+33.

## Ruling Out Encodings (Heuristic Only)

The minimum ASCII byte present can EXCLUDE encodings but cannot confirm one: an ASCII >= 64 file is consistent with all four variants. Use this only to narrow candidates, then confirm against instrument metadata.

```python
def candidate_encodings(filepath, sample_size=1000):
    '''Narrow FASTQ encoding candidates from the minimum quality byte. Confirm with run metadata.'''
    min_byte = 126
    count = 0
    with open(filepath) as handle:
        for i, line in enumerate(handle):
            if i % 4 == 3:
                for char in line.strip():
                    min_byte = min(min_byte, ord(char))
                count += 1
                if count >= sample_size:
                    break
    if min_byte < 59:
        return ['fastq']                        # only Phred+33 reaches below ASCII 59
    if min_byte < 64:
        return ['fastq-solexa']                 # ASCII 59-63 unique to Solexa+64
    return ['fastq', 'fastq-solexa', 'fastq-illumina']  # ambiguous - metadata decides
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `ValueError: ... not in correct range (...right QualityIO parser?)` | Wrong format string; a char is out of the chosen parser's range | Use the encoding the instrument produced (`'fastq'`, `'fastq-illumina'`, or `'fastq-solexa'`) |
| Scores look uniformly ~31 too high or too low; QC silently off | Overlap-region 31-shift from a mislabeled offset | Confirm encoding from metadata; never guess. Phred+33 read as `fastq-illumina` is 31 low; Phred+64 read as `fastq-sanger` is 31 high |
| `AttributeError`/`KeyError` on `per_letter_annotations` | That attribute does not exist | Use `record.letter_annotations['phred_quality']` (or `['solexa_quality']` for Solexa) |
| `KeyError: 'phred_quality'` on Solexa data | Parsed with `'fastq-solexa'`, which stores `'solexa_quality'` | Read `['solexa_quality']`, or convert to Phred on write |
| Trailing `B`/Q2 bases survive trimming | Illumina 1.5-1.7 B-tail treated as real Q2 | Strip trailing `B` runs as QC flags, not scores |
| Quality histogram shows discrete spikes | NovaSeq 4-level binning (Q2/Q12/Q23/Q37) | Expected on binned instruments; not a data problem |

## References

Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM (2010). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research 38(6):1767-1771.

Ewing B, Green P (1998). Base-calling of automated sequencer traces using phred. II. Error probabilities. Genome Research 8(3):186-194.

Ewing B, Hillier L, Wendl MC, Green P (1998). Base-calling of automated sequencer traces using phred. I. Accuracy assessment. Genome Research 8(3):175-185.

## Related Skills

- read-sequences - Parse FASTQ records and choose parse vs index for large files
- filter-sequences - Filter reads by length and content alongside quality
- paired-end-fastq - Keep R1/R2 synchronized when filtering paired reads
- sequence-statistics - Summary statistics across read sets
- read-qc/quality-reports - FastQC-style aggregate quality reports
- alignment-files/sam-bam-basics - Align filtered reads; quality scores carry into BAM
