---
name: bio-genome-intervals-interval-arithmetic
description: Performs set operations on genomic intervals - intersect (-wa/-wb/-wo/-wao/-loj/-c/-v/-u), subtract (-A), merge (-d, -c/-o), complement, cluster, multiinter, unionbedg, map, and groupby - with bedtools (CLI) and pybedtools/pyranges/bioframe (Python). Covers the sorted-input contract and the -sorted chromosome-order footgun, reciprocal/fractional overlap (-f/-F/-r/-e) and the A-vs-B asymmetry, -split for spliced/BED12/BAM features, and jaccard/fisher as mechanics only. Use when finding overlapping or unique regions between BED/peak/feature files, building consensus peaksets, removing blacklisted regions, transferring annotation values onto intervals, or computing interval-set similarity; route overlap-significance testing to overlap-significance.
tool_type: mixed
primary_tool: bedtools
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, pybedtools 0.10+, pyranges 0.1+ (or 1.0+ - see note), bioframe 0.7+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `bedtools --version` then `bedtools <subcommand> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

pyranges has a major-version API split: pyranges 0.x and the 1.0 rewrite (package `pyranges1`) differ in method names and return shapes. Verify with `import pyranges; pyranges.__version__` before pasting v0 idioms. If code throws an error, introspect the installed package and adapt rather than retrying.

# Interval Arithmetic

**"Which of my peaks overlap promoters, and how do I combine/subtract/annotate interval sets?"** -> Apply exact, deterministic set operations to sorted interval files, guarding the preconditions (prior `sort`, the `-sorted` chromosome-order contract, `-split`) that otherwise corrupt the answer.
- CLI: `bedtools intersect -a a.bed -b b.bed -u`, `bedtools merge`, `bedtools subtract`, `bedtools map -c 4 -o mean`
- Python: `a.intersect(b, u=True)`, `a.merge()` (pybedtools); `pr_a.overlap(pr_b)` (pyranges); `bf.overlap(df1, df2)` (bioframe)

## The Single Most Important Modern Insight -- The Arithmetic Is Exact; the Danger Is the Silent Preconditions

The set operations themselves are exact and deterministic - bedtools, pyranges, and bioframe compute identical geometry on the same 0-based half-open intervals. The bugs are never in the arithmetic; they hide in four preconditions that fail quietly, returning a plausible wrong answer with exit code 0:

1. **`merge`, `map`, `closest`, `groupby` require prior `sort`.** `merge` only collapses records that are adjacent *in file order* - on unsorted input, overlapping intervals survive un-merged and downstream counts are wrong, with no warning.
2. **`-sorted` requires sorted input in a shared chromosome order.** It swaps `intersect`'s in-memory interval tree for a low-memory chromosome sweep. Modern bedtools (>=~2.25) detects unsorted or differently-ordered `-sorted` input and errors out (exit 1: `... is not sorted` / `chromomsome sort ordering ... is inconsistent`); older versions silently swept past overlaps and under-reported. Pass `-g genome.txt` to pin the expected chromosome order (reproducible, and it catches the subtler missing-chromosome cases). The mismatch that stays SILENT on every version is a chromosome-NAME difference (`chr1` vs `1`), which returns an empty result with no error.
3. **`-split` changes whether the count is exons or the spanning envelope.** A BED12 record or spliced BAM read (CIGAR `N`) spans introns; without `-split` bedtools intersects the whole intron-spanning envelope, silently inflating RNA-seq overlaps. With `-split` it intersects only the blocks (exons).
4. **A raw overlap count is not association.** Long features, clustered features, and uneven coverage all inflate it; the number means nothing without a null. bedtools `fisher` is a weak analytic screen, not the answer - route rigorous significance to overlap-significance.

## Tool Taxonomy

| Tool | Role | Mechanism | When |
|------|------|-----------|------|
| bedtools | CLI interval algebra (reference implementation) | streaming sweep on sorted input; in-memory tree otherwise | shell pipelines, large files, reproducible one-liners |
| pybedtools | Python wrapper over bedtools | shells out to the bedtools binary; BedTool objects, flags as kwargs | inside a Python script; need exact bedtools parity; chaining with pandas |
| pyranges | pure-Python vectorized engine | native NumPy/pandas PyRanges; no bedtools dependency | large in-memory joins, no bedtools install, dataframe-native; mind the v0/v1 split |
| bioframe | functions on a plain pandas DataFrame | vectorized pandas merges; columns `chrom/start/end` | data already in pandas / the cooler-Hi-C ecosystem |

All three Python engines compute the same overlaps; the porting bugs are about default strand handling and return shape (pyranges `overlap` vs `join` vs `intersect`; bioframe `overlap` with `how=`), not geometry.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Quick overlap on the command line | `bedtools intersect -u` | no Python overhead; reproducible one-liner |
| Inside a pandas/Python pipeline | pybedtools or pyranges/bioframe | stays in-process; pyranges/bioframe need no bedtools binary |
| Whole-genome-scale intersect | `intersect -sorted -g genome.txt` | low-memory sweep; modern bedtools errors on a sort/order mismatch, `-g` pins the expected chromosome order |
| Spliced reads / BED12 vs exons | add `-split` | otherwise the intron-spanning envelope is intersected (RNA-seq inflation) |
| Are two SV/CNV calls the same event? | `-f 0.5 -r` (50% reciprocal) | one-sided fractions let a giant interval swallow a tiny one |
| Transfer/aggregate B values onto A | `bedtools map -c COL -o OP` | columnar alternative to `intersect -wo \| groupby` |
| Build consensus peakset from replicates | `cat \| sort \| merge -d N` | collapses replicate peaks within N bp |
| Multi-sample shared-region map | `multiinter` / `unionbedg` | presence/absence (intervals) or stacked signal matrix |
| Is the overlap more than chance? | -> overlap-significance | raw count is length/coverage-confounded; needs a permutation null |
| Peaks not yet called | -> chip-seq/peak-calling or atac-seq/atac-peak-calling | this category operates on existing intervals |

## Intersect - the Workhorse

The output-mode flags do not change *what overlaps*; they change *what gets printed* (the #1 source of confusion). Full flag semantics are in usage-guide.md.

```bash
bedtools intersect -a peaks.bed -b genes.bed -u            # whole A, once, if it overlaps >=1 B
bedtools intersect -a peaks.bed -b genes.bed -v            # A features with NO overlap (set difference)
bedtools intersect -a peaks.bed -b genes.bed -c            # per-A count of B hits (0 if none)
bedtools intersect -a peaks.bed -b genes.bed -wa -wb       # whole A + whole B, one line per pair ("join")
bedtools intersect -a peaks.bed -b genes.bed -loj          # left outer join: every A, NULL B if none
bedtools intersect -a peaks.bed -b genes.bed -wo           # A+B+bp-of-overlap, only A with overlap
bedtools intersect -a peaks.bed -b genes.bed -wao          # like -wo but A-with-no-overlap kept (B=., overlap=0)
```

```python
import pybedtools

a = pybedtools.BedTool('peaks.bed')
b = pybedtools.BedTool('genes.bed')
a.intersect(b, u=True)            # flags become kwargs
a.intersect(b, wa=True, wb=True)
a.intersect(b, c=True)
```

## Subtract, Merge, Complement, Cluster

```bash
bedtools subtract -a a.bed -b b.bed             # clip the overlapping portions out of A (A can fragment)
bedtools subtract -a a.bed -b b.bed -A          # drop the ENTIRE A feature if any part overlaps B
bedtools sort -i a.bed | bedtools merge -d 0    # collapse overlapping + book-ended; -d 0 is the default
bedtools sort -i a.bed | bedtools merge -c 4,5 -o distinct,sum   # summarize columns while merging
bedtools complement -i a.bed -g genome.txt      # the gaps: genome NOT covered by A (genome file required)
bedtools sort -i a.bed | bedtools cluster -d 0  # assign a cluster id to overlapping/adjacent features
```

Valid `-o` operations: `sum, min, max, absmin, absmax, mean, median, mode, antimode, stdev, sstdev, collapse, distinct, count, count_distinct, first, last`. `merge -d 0` merges overlapping and book-ended (touching) features but NOT a 1 bp gap; `-d 1` does.

## Map - Transfer Values, and Groupby - Aggregate

**Goal:** Summarize a column of overlapping B features onto each A interval (e.g. mean signal per gene).

**Approach:** For each sorted A interval, `map` collects overlapping B features and applies an aggregation `-o` to a B column `-c`; `groupby` is the single-file SQL-style aggregator after an `intersect -wo`.

```bash
bedtools map -a genes.bed -b scores.bedgraph -c 4 -o mean      # both inputs MUST be sorted
bedtools intersect -a genes.bed -b peaks.bed -wo \
  | bedtools groupby -g 1,2,3,4 -c 13 -o sum                   # group on A cols, sum the overlap-bp col
```

```python
import pybedtools

genes = pybedtools.BedTool('genes.bed').sort()
scores = pybedtools.BedTool('scores.bedgraph').sort()
genes.map(scores, c=4, o='mean')
```

## Multi-Sample: Multiinter and Unionbedg

```bash
bedtools multiinter -header -names s1 s2 s3 -i s1.bed s2.bed s3.bed   # which files cover each sub-interval
bedtools unionbedg -header -names s1 s2 s3 -i s1.bg s2.bg s3.bg       # stack bedGraph signal into a matrix
```

`multiinter` is the interval presence/absence map (build a consensus by filtering its `num`/`list` columns); `unionbedg` is its signal-track analog.

## Jaccard and Fisher - Mechanics Only

`jaccard` is a single similarity scalar `|A n B| / |A u B|` in [0,1], useful for all-vs-all dataset clustering - it is NOT a significance test (no p-value). `fisher` builds a 2x2 table and returns a Fisher p, but it estimates the in-neither cell from a mean-interval-size/genome-size heuristic, ignores genome structure, and is prone to inflation - treat it as a fast triage screen only.

```bash
bedtools jaccard -a a.bed -b b.bed -g genome.txt        # both sorted; reports jaccard + n_intersections
bedtools fisher  -a a.bed -b b.bed -g genome.txt         # weak analytic null; validate any low p by simulation
```

For a defensible enrichment p-value (size-preserving permutation in an accessible workspace, GAT/regioneR/LOLA/GREAT), route to overlap-significance.

## Per-Method Failure Modes

### Merge without sorting first
**Trigger:** `bedtools merge` (or `cluster`/`map`/`groupby`) on unsorted input. **Mechanism:** merge only collapses records adjacent in file order. **Symptom:** overlapping intervals survive un-merged; counts wrong, no error. **Fix:** `bedtools sort -i in.bed | bedtools merge`.

### `-sorted` on unsorted or differently-ordered input
**Trigger:** `intersect -sorted` on unsorted input or files in different chromosome orders. **Mechanism:** the sweep walks both files in lockstep assuming a shared order. **Symptom:** modern bedtools (>=~2.25) errors out (`... is not sorted` / `chromomsome sort ordering ... is inconsistent`, exit 1); pre-2.25 returned a silently smaller set. **Fix:** sort every input identically and pass `-g genome.txt` to pin the order; on an old bedtools, suspect this when a result is surprisingly small.

### Missing `-split` on spliced features
**Trigger:** intersecting BED12 / spliced BAM without `-split`. **Mechanism:** the intron-spanning envelope is treated as solid. **Symptom:** intronic positions "overlap" exons; RNA-seq overlap inflated/smeared. **Fix:** add `-split` whenever an operand is BED12 or a spliced alignment and exon-level truth is required.

### `-f` vs `-F` swapped, or default 1 bp overlap
**Trigger:** thresholding the wrong set, or no `-f` at all. **Mechanism:** `-f` is a fraction of A, `-F` a fraction of B (default `-f 1e-9` = any 1 bp); A and B play asymmetric roles. **Symptom:** a tiny peak "inside" a 2 Mb gene by one base; swapping `-a`/`-b` changes counts. **Fix:** threshold the *small* set; use `-r` for "same event" concordance.

### complement/shuffle without a genome file
**Trigger:** `complement` (or `closest`/`map` order assumptions) without `-g`. **Mechanism:** bedtools cannot know where chromosomes end. **Symptom:** error, or gaps/coordinates that run past chromosome ends. **Fix:** pass a correct `-g genome.txt` built from the same assembly.

### chrom-naming mismatch (`chr1` vs `1`)
**Trigger:** BED uses `chr1`, genome/other file uses `1`. **Mechanism:** chromosomes never match. **Symptom:** empty/zero output, no error. **Fix:** harmonize naming across all inputs and the genome file.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Overlap fraction `-f` (state explicitly) | analysis choice | default `-f 1e-9` (1 bp) is rarely the biological question; threshold the small set |
| 50% reciprocal overlap (`-f 0.5 -r`) | SV/CNV field convention | "are these the same event"; one-sided lets a big interval swallow a small one |
| Merge `-d` (e.g. 100 bp for replicate consensus) | replicate-merge convention | collapses near-coincident replicate peaks; tune per assay/resolution |
| `merge -d 0` (default) | bedtools default | merges overlapping + book-ended, NOT a 1 bp gap (use `-d 1` for that) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Empty intersect output | chrom naming mismatch (`chr1` vs `1`) | harmonize naming across files + genome.txt |
| `merge` left overlaps behind | input not sorted | `sort` before `merge`/`cluster`/`map`/`groupby` |
| `-sorted` errors or (old bedtools) returns too few | unsorted or mismatched chromosome order | sort all inputs identically; add `-g genome.txt` to pin order |
| RNA-seq overlap looks inflated | missing `-split` on BED12/spliced BAM | add `-split` |
| Negative start / past-chromosome-end | wrong/missing `-g genome.txt` | pass a correct chrom-sizes file |
| pyranges AttributeError | 0.x vs 1.0 API mismatch | check `pyranges.__version__`; use matching method names |

## References

- Quinlan AR, Hall IM. 2010. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26:841-842.
- Dale RK, Pedersen BS, Quinlan AR. 2011. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. *Bioinformatics* 27:3423-3424.
- Stovner EB, Sætrom P. 2020. PyRanges: efficient comparison of genomic intervals in Python. *Bioinformatics* 36:918-919.
- Open2C, Abdennur N, Fudenberg G, Flyamer IM, Galitsyna AA, Goloborodko A, Imakaev M, Venev SV. 2024. Bioframe: operations on genomic intervals in pandas dataframes. *Bioinformatics* 40:btae088.

## Related Skills

- bed-file-basics - BED format, coordinate systems, and the conversions this skill depends on
- overlap-significance - Whether an overlap count exceeds a matched null (permutation, GAT/regioneR/LOLA/GREAT)
- proximity-operations - closest, window, flank, slop for adjacency rather than membership
- coverage-analysis - per-base depth and bedGraph signal feeding map/unionbedg
- gtf-gff-handling - exon/feature models whose `-split` behavior this skill depends on
- chip-seq/peak-calling - source of the peak BED files these operations consume
- atac-seq/consensus-peakset - replicate merge via merge/multiinter
