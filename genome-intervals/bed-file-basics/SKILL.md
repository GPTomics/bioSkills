---
name: bio-genome-intervals-bed-file-basics
description: Handles BED-format genomic intervals (BED3 through BED12, narrowPeak/broadPeak) and the coordinate-system substrate the whole interval category rests on, with bedtools (CLI) and pybedtools/pyranges/pandas (Python). Covers the 0-based half-open vs 1-based-closed convention boundary and the start-1/end-unchanged conversion, the silent failures (chrom-name mismatch, CRLF, lexicographic-vs-version sort under -sorted), genome/chrom.sizes generation, sorting contracts, BED12 block invariants, validation, makewindows, cross-assembly liftover (liftOver/CrossMap), and BED<->VCF/BAM/FASTA conversion. Use when reading, creating, validating, sorting, lifting between genome builds, or converting interval files, preparing inputs for bedtools/tabix/bigBed, or debugging an off-by-one or empty-overlap result.
tool_type: mixed
primary_tool: bedtools
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, pybedtools 0.10+, pyranges 0.x (the `pyranges1` rewrite ships as a separate package), samtools 1.19+, pandas 2.2+, UCSC liftOver / CrossMap 0.7+ (the bare `CrossMap` entry point replaced `CrossMap.py` at 0.7.0).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

pyranges has a major-version API split: pyranges 0.x and the `pyranges1` rewrite differ in method names and DataFrame access; both keep `Chromosome/Start/End` columns and 0-based half-open coordinates. Check `import pyranges; pyranges.__version__` before chaining methods. Operations that need chromosome lengths (`slop`, `complement`, `shuffle`, `makewindows -g`, `-sorted` ordering) require a genome/chrom.sizes file. If code throws an error, introspect the installed tool and adapt rather than retrying.

# BED File Basics

**"Work with this interval file without shifting everything by one base"** -> Establish the coordinate convention from the format, read/create/validate/sort the intervals, and convert across format boundaries with the correct base shift.
- CLI: `bedtools sort -i in.bed`, `bedtools getfasta`, `bedtools makewindows -g genome.txt -w 10000`, `sort -k1,1 -k2,2n`
- Python: `pybedtools.BedTool('in.bed')`, `pr.read_bed('in.bed')` (pyranges), `pd.read_csv(sep='\t', comment='#')`

## The Single Most Important Modern Insight -- A Coordinate Is a Bare Integer With No Self-Describing Convention

A `start` column is just an `int`. Nothing in the file says whether it is 0-based or 1-based, so the convention lives in the analyst's head, keyed off the **format**, not the data. BED is 0-based half-open `[start, end)`; GTF/GFF, SAM, VCF, and wiggle are 1-based fully-closed `[start, end]`. Three load-bearing consequences:

1. **The conversion is `start - 1, end unchanged` -- and it throws no error if wrong.** A botched convention shift still parses, still runs, and silently shifts every answer by one base: gene bodies 1 bp short, boundary SNPs flipping in/out, exact-edge intersections toggling. The end is numerically identical between BED-half-open and GFF-closed because GFF's last *included* base and BED's first *excluded* position are the same boundary. The symmetry instinct (subtract 1 from both) is the classic bug. The reflex: convert `start_bed = start_1based - 1` (end unchanged) and **test the round-trip on a 1 bp feature** -- BED `chr1 5 6` == GFF `chr1 6 6`, both one base. Length is `end - start` in BED (no `+1`); `end - start + 1` in GTF.

2. **The two truly silent file failures.** (a) **Chrom-name mismatch** (`chr1` vs `1`, `chrM` vs `MT`): intersecting a `chr`-prefixed file against a bare-numeral one yields a perfectly valid **empty** result -- "no overlap" looks like biology, not a bug. Confirm shared naming (`cut -f1 a.bed | sort -u`) before any cross-file op. (b) **CRLF line endings** from Excel/Windows glue `\r` onto the last field (`end` becomes `"100\r"`); the tell is "works for some tools, breaks for others." `cat -A` shows `^M$`; fix with `dos2unix`. Never open a BED in Excel -- it date-mangles `SEPT9`->`9-Sep` and float-truncates large coordinates.

3. **bedtools `-sorted` assumes both inputs share the SAME chromosome order.** With a lexicographic (`chr1, chr10, chr2`) vs version (`chr1, chr2, chr10`) sort mismatch, modern bedtools (>=~2.25) detects the inconsistency and **errors out** (exit 1, `chromomsome sort ordering ... is inconsistent`); older versions silently swept past and **dropped chr10-chr22**. Either sort both files with the identical command, or pass `-g genome.txt` (derived from the same reference FASTA) to pin the expected order. For one-off work, omit `-sorted` (the in-memory path tolerates any order). The mismatch that stays SILENT on every version is a chromosome-NAME difference (`chr1` vs `1`), which returns an empty result with no error.

## Tool Taxonomy

| Tool | Role | Mechanism | When |
|------|------|-----------|------|
| bedtools | CLI interval algebra (Quinlan 2010 *Bioinformatics* 26:841) | streaming, sorted-input reference implementation | shell pipelines, large files, reproducible one-liners |
| pybedtools | Python wrapper over bedtools (Dale 2011 *Bioinformatics* 27:3423) | shells out per op; BedTool objects + iterators | inside a Python analysis; chaining with pandas |
| pyranges | pure-Python interval engine (Stovner 2020 *Bioinformatics* 36:918) | vectorized PyRanges/pandas, no bedtools binary | large in-memory joins, dataframe-native workflows |
| pandas | flat tabular read | `read_csv(sep='\t')`; knows NO coordinate semantics | quick filter/inspect; the analyst enforces 0-based + sort manually |
| UCSC bedToBigBed / tabix | indexed/compressed BED for random access | requires `sort -k1,1 -k2,2n`, no track lines | browser tracks, region queries on huge files |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Quick create/sort/filter on the command line | bedtools + coreutils `sort -k1,1 -k2,2n` | no Python overhead; reproducible |
| Inside a pandas/Python pipeline | pybedtools or pyranges | stays in-process; pyranges if no bedtools binary |
| Convert VCF/GTF/SAM positions to BED | subtract 1 from start, end unchanged | the convention boundary; test on a 1 bp feature |
| Empty intersect / "no overlap found" | check chrom naming (`chr1` vs `1`) FIRST | the most common silent null result |
| Using `-sorted` for speed/RAM | sort both files identically, or pass `-g genome.txt` | modern bedtools errors on a lexicographic-vs-version mismatch; old versions dropped chroms silently |
| Need chromosome lengths (slop/complement/windows) | generate genome.txt from the SAME FASTA | a stale/generic chrom.sizes rots slop/complement |
| Set operations on these intervals | -> interval-arithmetic | this skill is the format/coordinate substrate |
| Parse a GTF/GFF gene model | -> gtf-gff-handling | 1-based, parent/child hierarchy, not a flat BED |
| Peaks not yet called | -> chip-seq/peak-calling or atac-seq/atac-peak-calling | this category operates on existing intervals |
| Convert between assemblies (hg19<->hg38) | liftOver/CrossMap, report unmapped | a different problem from convention shifts |

## BED Columns (BED3 -> BED12)

The first 3 fields are required; the rest are optional but **positional** (cannot supply field N without 1..N-1), and the field count must be identical on every line.

```
BED3   chrom  start  end
BED4   + name
BED5   + score (int 0-1000; '.' allowed)
BED6   + strand (+/-/.)            # the common stranded-interval form
BED12  + thickStart thickEnd itemRgb blockCount blockSizes blockStarts   # transcript/exon models
```

narrowPeak is **BED6+4** (`signalValue pValue qValue peak`); broadPeak is **BED6+3** (drops `peak`). The `peak` column is a **0-based offset from chromStart** (absolute summit = `chromStart + peak`), `-1` if none; `pValue`/`qValue` are `-log10` scaled with `-1` meaning "not assigned", NOT p=0.1.

## Create and Read BED Files

```python
import pybedtools
import pandas as pd

intervals = [('chr1', 100, 200, 'peak1', 100, '+'), ('chr1', 300, 400, 'peak2', 200, '-')]
bed = pybedtools.BedTool(intervals)                                       # from list of tuples
bed = pybedtools.BedTool.from_dataframe(pd.read_csv('peaks.tsv', sep='\t'))  # from a DataFrame
bed.saveas('peaks.bed')

for iv in pybedtools.BedTool('peaks.bed'):
    print(iv.chrom, iv.start, iv.end, len(iv))   # start/end are ints; len(iv) == end - start
df = pybedtools.BedTool('peaks.bed').to_dataframe(names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
```

pandas reads BED as a flat table but knows nothing about coordinates: `pd.read_csv('in.bed', sep='\t', header=None, comment='#')` -- the analyst enforces 0-based and sorts manually.

## Generate the Genome / chrom.sizes File

**Goal:** Produce the authoritative chromosome-length file that `slop`, `complement`, `shuffle`, `makewindows -g`, and `-sorted` ordering depend on.

**Approach:** Index the exact reference FASTA the rest of the pipeline used and take its first two columns -- never download a generic `hg38.chrom.sizes` and hope it matches the assembly the BAM was aligned to.

```bash
samtools faidx ref.fa
cut -f1,2 ref.fa.fai > genome.txt   # chrom<TAB>size; bedtools reads only the first 2 cols, accepts the .fai directly
```

## Convert Across Format Boundaries

**Goal:** Move positions between BED and VCF/SAM/GTF/BAM without an off-by-one.

**Approach:** Subtract 1 from the 1-based start (end unchanged) going TO BED; for BAM, let bedtools do the conversion; verify a known single-base landmark afterwards.

```bash
grep -v '^#' in.vcf | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' > variants.bed   # VCF POS (1-based) -> BED; start-1
bedtools bamtobed -i in.bam > alignments.bed                                       # bedtools handles the convention
bedtools bamtobed -i in.bam -split > spliced.bed                                   # split spliced reads into blocks
bedtools getfasta -fi ref.fa -bed in.bed -name -fo out.fa                          # extract sequence (uses the .fai)
```

## Sort, Validate, and Make Windows

```bash
bedtools sort -i in.bed > sorted.bed                 # lexicographic, == sort -k1,1 -k2,2n
bedtools sort -i in.bed -faidx names.txt > sorted.bed  # reorder to an arbitrary reference contig order
awk -F'\t' '{print NF}' in.bed | sort -u             # field count consistent? (one value expected)
awk -F'\t' '$2 < 0 || $2 >= $3' in.bed               # negative or inverted intervals (bedtools also errors on these)
cut -f1 in.bed | sort -u                             # chromosome names (compare against the partner file)
bedtools makewindows -g genome.txt -w 10000 -s 5000 -i winnum > windows.bed   # 10 kb sliding windows, step 5 kb, numbered
```

## Cross-Assembly Liftover (a different problem from convention conversion)

**Goal:** Move coordinates between genome builds (hg19<->hg38, mm10<->mm39) - which is remapping to a different reference, NOT the 0-based/1-based convention shift above.

**Approach:** Map intervals through a chain file with UCSC `liftOver` (BED) or CrossMap (BED/VCF/GFF/BAM/bigWig), and ALWAYS inspect the unmapped file - regions that fail to map (assembly gaps, rearrangements, split/merged contigs) are dropped, and silently ignoring them biases everything downstream.

```bash
liftOver in.hg19.bed hg19ToHg38.over.chain.gz out.hg38.bed unmapped.bed   # UCSC; chain from UCSC goldenPath
wc -l unmapped.bed                                                         # NEVER skip: dropped regions are not random
CrossMap bed GRCh37_to_GRCh38.chain.gz in.bed > out.bed                    # CrossMap also does vcf/gff/bam/bigwig
```

A coordinate is meaningless without its assembly just as it is meaningless without its convention - record the build (and the chain provenance) alongside the file. Liftover is many-to-one and one-to-none in places; never assume a 1:1 round-trip.

## BED12 Block Invariants

BED12 is a referentially-integral structure, not a flat table. `blockStarts` are **offsets from chromStart**, not absolute coordinates; the **first blockStart must be 0**; `blockStarts[last] + blockSizes[last]` must equal `chromEnd - chromStart`; blocks are ascending and non-overlapping. `thickStart/thickEnd` (the CDS) are absolute and independent of the blocks -- a non-coding feature sets both to `chromStart`. Validate by reconstructing absolute exon coordinates (`chromStart + blockStarts[i]`) and confirming they fall inside `[chromStart, chromEnd]`; `bedtools bed12tobed6` explodes a model into per-exon BED6 as a quick sanity reconstruction.

## Per-Method Failure Modes

### Off-by-one at a convention boundary
**Trigger:** treating a GTF/VCF 1-based `start` as a BED `chromStart`. **Mechanism:** the convention is not stored in the data. **Symptom:** every feature 1 bp short; boundary intersections flip; no error. **Fix:** `start - 1, end unchanged`; test the round-trip on a 1 bp feature.

### Chrom-name mismatch
**Trigger:** intersecting `chr`-prefixed against bare-numeral files. **Mechanism:** chrom strings never match. **Symptom:** valid empty/zero result presented as biology. **Fix:** harmonize naming across all files and genome.txt before any cross-file op.

### Lexicographic-vs-version sort under `-sorted`
**Trigger:** two inputs sorted in different chromosome orders. **Mechanism:** sweep-line assumes one shared order and concludes chroms passed each other. **Symptom:** modern bedtools errors out (`chromomsome sort ordering ... is inconsistent`, exit 1); pre-2.25 silently dropped chr10-chr22 and favored low chroms. **Fix:** sort both identically or pass `-g genome.txt`; omit `-sorted` for one-off work.

### CRLF line endings
**Trigger:** file authored/round-tripped through Windows/Excel. **Mechanism:** `\r` glues onto the last field. **Symptom:** "not an integer" or cryptic last-column errors; works in some tools. **Fix:** `dos2unix` or `sed 's/\r$//'`; never edit BED in Excel.

### Stale / generic genome file
**Trigger:** a downloaded chrom.sizes that does not match the aligned FASTA. **Mechanism:** missing contigs or wrong lengths. **Symptom:** `slop`/`complement` run off real ends or drop contigs; `-sorted` order breaks. **Fix:** regenerate from the exact reference FASTA.

### pyranges 0.x idioms on pyranges1
**Trigger:** running 0.x method/attribute names against the `pyranges1` rewrite. **Mechanism:** the rewrite renamed methods and changed DataFrame access. **Symptom:** AttributeError. **Fix:** check `pyranges.__version__` and use the matching API.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Required minimum 3 columns (BED3) | UCSC/hts-specs BEDv1 | chrom/start/end; field count identical on every line |
| BED `score` integer 0-1000 | UCSC spec | maps to browser gray shade; analysis tools tolerate out-of-range, browser clamps |
| narrowPeak/broadPeak `-1` in pValue/qValue/peak | ENCODE narrowPeak.as | means "not assigned", NOT a real value (e.g. p=0.1 or summit at offset 0) |
| Convention shift: start - 1, end unchanged | BED 0-based vs GFF/VCF 1-based | only the start moves; ends coincide at the shared boundary |
| makewindows window/step (e.g. -w 10000 -s 5000) | analysis choice | resolution vs file size; state the size used, do not default silently |
| tabix BED needs `-p bed` (or `-0`) | tabix default is 1-based | `-p bed` sets chrom/start/end cols and 0-based interpretation |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Empty intersect / "no overlap" | chrom naming mismatch (`chr1` vs `1`, `chrM` vs `MT`) | harmonize naming across all files + genome.txt |
| Every feature 1 bp short | converted 1-based->BED without the start-1 (or subtracted from both) | `start - 1, end unchanged`; verify on a 1 bp feature |
| Every coordinate shifted one column right | UCSC MySQL/`SELECT *` table dump prepends a `bin` column (hierarchical binning index, Kent 2002) | drop field 1 before feeding bedtools: `cut -f2-` |
| `-sorted` errors or (old bedtools) favors low chroms | lexicographic vs version sort mismatch | sort both identically or pass `-g genome.txt`; or drop `-sorted` |
| "not an integer" on the last column | CRLF line endings | `dos2unix` / `sed 's/\r$//'` |
| Negative start / past-chrom-end after slop | missing or stale genome file | regenerate genome.txt from the aligned FASTA |
| `bedToBigBed`/`tabix` errors | unsorted input or track lines present | `sort -k1,1 -k2,2n`; strip `track`/`browser`/`#` lines |
| pyranges AttributeError | 0.x vs pyranges1 API mismatch | check `pyranges.__version__`, use matching method names |

## References

- Quinlan AR, Hall IM. 2010. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26:841-842.
- Dale RK, Pedersen BS, Quinlan AR. 2011. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. *Bioinformatics* 27:3423-3424.
- Stovner EB, Sætrom P. 2020. PyRanges: efficient comparison of genomic intervals in Python. *Bioinformatics* 36:918-919.
- Kent WJ, Sugnet CW, Furey TS, et al. 2002. The Human Genome Browser at UCSC. *Genome Res* 12:996-1006.
- The Browser Extensible Data (BED) format, hts-specs BEDv1. samtools.github.io/hts-specs/BEDv1.pdf.
- UCSC Genome Browser FAQ: Data File Formats. genome.ucsc.edu/FAQ/FAQformat.html.

## Related Skills

- interval-arithmetic - Set operations (intersect/merge/subtract) on the intervals defined here
- gtf-gff-handling - 1-based annotation parsing and the parent/child gene-model hierarchy
- coverage-analysis - Per-base depth that becomes bedGraph intervals
- alignment-files/sam-bam-basics - BAM-to-BED conversion and the SAM-1-based vs BAM-0-based distinction
- variant-calling/vcf-basics - VCF POS (1-based) to BED conversion and indel left-anchoring
