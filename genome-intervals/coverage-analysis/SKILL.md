---
name: bio-genome-intervals-coverage-analysis
description: Computes and interprets sequencing read depth and coverage over a genome, windows, or target regions with mosdepth (windowed depth, cumulative distribution, --quantize callable BEDs), bedtools genomecov/coverage (bedGraph tracks, per-target stats), samtools depth/coverage (per-base depth, per-contig depth+breadth). Covers the breadth-vs-mean distinction, the cumulative-coverage curve, evenness (CV/Fano/fold-80/Gini), what each tool silently counts (duplicates, secondary/supplementary, MAPQ, read span vs fragment, mate-overlap), the samtools-depth 8000-cap version trap, and the bedtools coverage -a/-b orientation flip. Use when assessing sequencing adequacy, building coverage tracks, computing breadth at a depth threshold, defining callable regions, or QCing target-capture uniformity.
tool_type: mixed
primary_tool: bedtools
---

## Version Compatibility

Reference examples tested with: bedtools 2.31+, mosdepth 0.3+, samtools 1.19+, pybedtools 0.10+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

`samtools depth` behaviour changed across versions: pre-1.13 capped depth at 8000 and truncated silently (`-d 0` = unlimited); 1.13+ rewrote the subcommand with NO cap and `-d/-m` deprecated/ignored. Always check `samtools --version` before trusting a max-depth number. If code throws an error, introspect the installed tool and adapt rather than retrying.

# Coverage Analysis

**"Is my sequencing deep enough to answer the question?"** -> Measure depth as a distribution over positions, then report median, breadth at a depth threshold, and an evenness number -- never the mean alone.
- CLI: `mosdepth --by 500 prefix in.bam` (windowed depth + cumulative dist), `samtools coverage in.bam` (per-contig depth+breadth), `bedtools genomecov -ibam in.bam -bga` (bedGraph track)
- Python: `pybedtools.BedTool('in.bam').genome_coverage(bga=True)` (pybedtools); parse `prefix.mosdepth.global.dist.txt` for the breadth curve

## The Single Most Important Modern Insight -- Mean Depth Is a Budget, Not a Result; Report Breadth Off a Cumulative Curve

"30x WGS" describes what was paid for, not what was achieved. Coverage is a **distribution over positions**, and the mean is its worst summary: it is dragged **up** by a fat right tail (repeats, rDNA, mitochondria, PCR pileups, segmental dups) while staying **blind** to a hard left wall of zeros and near-zeros (GC-extreme exons, poorly-mappable regions, capture dropout). Two libraries with identical mean 30x can differ completely -- one even and callable everywhere, one spiky with 20% of the target uncallable. The mean hides both failures. Four load-bearing moves:

1. **Report MEDIAN, not mean.** The median is robust to the right tail. When mean/median exceeds ~1.1-1.2 the distribution is skewed and the mean is overstating typical depth -- that gap is a free evenness diagnostic.
2. **Report a BREADTH / cumulative-coverage curve.** "% of target >= 1x, >= 10x, >= 20x, >= 30x" is the honest summary, because adequacy is a breadth statement: a base that was not covered deeply enough is uncallable no matter how deep the rest of the genome is. mosdepth's `*.mosdepth.global.dist.txt` IS this curve. The killer question for any "mean = 30x" claim is "breadth at 20x?".
3. **Quantify EVENNESS** (CV, Fano factor, Picard fold-80, or Gini) -- an even 30x and a spiky 30x are different experiments, and a spiky library cannot be rescued by sequencing deeper (extra reads follow the same biased distribution; the holes stay holes). Fix the library (PCR-free, better capture, UMIs), not the lane count.
4. **Say WHAT WAS COUNTED.** A depth number is meaningless until the recipe is stated: duplicates dropped (only if MARKED first)? secondary/supplementary included? MAPQ filter? read span or fragment? mate-overlap corrected? per-base or per-region? The tools disagree on every one of these by default.

## Tool Taxonomy

| Tool | Counts what (defaults) | Per-base or region | When |
|------|------------------------|--------------------|------|
| mosdepth | corrects mate-overlap by default (off under `--fast-mode`/`-x`); `-Q` MAPQ filter; emits cumulative dist + summary | windowed (`--by`), per-region, or callable bins (`--quantize`) | the modern fast default for WGS/WES/targeted; gives the breadth curve directly |
| samtools coverage | per-reference summary (added 1.10); `coverage` column = breadth %, `meandepth` = depth | per-contig | quick "is this contig actually covered?" -- spots high-mean/low-breadth pileups |
| samtools depth | drops UNMAP/SECONDARY/QCFAIL/DUP by default; `-Q`/`-q` filters; `-s` de-double-counts overlap; CRAM needs `--reference` | per-base | exact per-base depth over small regions; watch the 8000-cap version trap |
| bedtools genomecov | counts READ coverage by default (double-counts mate overlap); `-pc` = fragment; `-split` for spliced | per-base / bedGraph / histogram | bedGraph tracks, genome-wide depth histogram |
| bedtools coverage | per-A-interval stats from B reads; `-a`/`-b` flipped at v2.24.0 | per-region (or `-d` per-base) | per-target counts/breadth/mean over a BED |
| Picard CollectHsMetrics | capture-kit QC | per-target panel | exome/panel uniformity: on-target %, fold-80, PCT_TARGET_BASES_20X |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| WGS / WES breadth + adequacy | `mosdepth --by` then parse `*.global.dist.txt` | emits the cumulative curve + median directly; fast |
| Quick per-contig depth & breadth glance | `samtools coverage` | one line/contig; `coverage` col = breadth, `meandepth` = depth |
| Exact per-base depth, small region | `samtools depth -a -r chr:from-to` | per-base; add `-s` for short-insert; check version for 8000 cap |
| bedGraph coverage TRACK for a browser | `bedtools genomecov -ibam -bga` (or `-bg`) | `-bga` marks zero-coverage gaps; convert to bigWig -> bigwig-tracks |
| Per-target counts/breadth/mean over a BED | `bedtools coverage -a targets.bed -b in.bam` | A = targets, B = reads (post-v2.24.0); `-mean` for mean depth |
| Callable-region BED (NO/LOW/CALLABLE/HIGH) | `mosdepth --quantize 0:1:4:150:` | lightweight CallableLoci replacement at scale |
| Target-capture uniformity QC | -> Picard CollectHsMetrics (fold-80, on-target %) | the capture-QC standard; off-target loss + bait unevenness |
| Spliced/RNA-seq depth | add `-split` (genomecov/coverage) | without it an intron (N CIGAR) is counted as covered |
| Short-insert VAF (amplicon/cfDNA) | correct mate-overlap: `samtools depth -s` / `genomecov -pc` / mosdepth default | naive per-base double-counts the overlap, corrupting VAFs |
| Normalized cross-sample track | -> chip-seq/chipseq-visualization (deepTools bamCoverage) | library-size correction (RPGC/CPM/BPM) for comparison |
| Pileup/variant evidence from BAM | -> alignment-files/pileup-generation | depth is upstream of per-call DP/AD |

## mosdepth -- The Modern Default

**Goal:** Get the median depth and the full breadth curve for a BAM in one fast pass.

**Approach:** Run mosdepth windowed (or whole-genome), then read the cumulative distribution file -- it already holds breadth at every depth threshold; no histogram integration needed.

```bash
mosdepth --by 500 -Q 20 sample in.bam     # --by 500 = 500 bp windows; -Q 20 = drop MAPQ<20 (repeat coverage collapses, intentionally)
# Outputs: sample.mosdepth.summary.txt (mean/min/max per chrom + total)
#          sample.mosdepth.global.dist.txt (cumulative: chrom, depth, proportion >= depth)
#          sample.regions.bed.gz (per-window mean depth)
```

The `*.global.dist.txt` rows are `chrom  depth  proportion_of_bases_at_least_this_depth` -- the breadth curve directly. Read median as the depth where proportion crosses 0.5. `--fast-mode`/`-x` is ~2x faster but SILENTLY disables mate-overlap correction -- fine for a rough WGS glance, wrong for VAF-sensitive short-insert data.

**Goal:** Emit a callable-region BED (NO_COVERAGE / LOW / CALLABLE / HIGH) without GATK3.

**Approach:** Use `--quantize` to bin depth and merge adjacent equal-bin runs into a compact BED.

```bash
mosdepth --quantize 0:1:4:150: callable in.bam   # bins: [0,1)=NO_COVERAGE, [1,4)=LOW, [4,150)=CALLABLE, [150,inf)=HIGH
# 4 = min callable depth (tune to caller); 150 = excessive-depth ceiling (flags rDNA/artifact pileups)
zcat callable.quantized.bed.gz | head
```

## bedtools genomecov -- Tracks and the Histogram Default

```bash
bedtools genomecov -ibam in.bam -bga > cov.bedGraph   # -bga = bedGraph INCLUDING zero-coverage runs; -bg omits zeros
bedtools genomecov -ibam in.bam -pc -bg > frag.bedGraph # -pc = FRAGMENT coverage (mate overlap counted once); default counts reads (double-counts overlap)
bedtools genomecov -ibam in.bam -split -bg > rna.bedGraph # -split = skip N-CIGAR gaps (introns); MANDATORY for spliced RNA-seq
bedtools genomecov -ibam in.bam > hist.txt            # NO output flag = a 5-col HISTOGRAM, not a track
```

The bare default is a **histogram**, not a bedGraph -- 5 columns: `chrom  depth  bases_at_that_depth  chrom_size  fraction_of_chrom`, with a final `genome` block for the whole genome. Breadth/mean must be integrated from it yourself (sum `fraction` over `depth >= threshold`) -- which is exactly why mosdepth's ready-made dist file is preferred.

## bedtools coverage -- Per-Target Stats (mind the orientation)

```bash
bedtools coverage -a targets.bed -b in.bam > per_target.bed   # stats reported FOR each A interval
bedtools coverage -a targets.bed -b in.bam -mean > mean.bed   # -mean = mean depth per A interval
```

As of bedtools **v2.24.0** coverage is computed for the **`-a`** file (it was `-b` before) -- A = the regions stats are wanted for (targets), B = the reads. The default appends 4 columns to each A interval: (1) count of B features overlapping, (2) bases in A covered >=1x, (3) length of A, (4) fraction of A covered (col2/col3 = per-interval breadth). `-d` = per-base depth within each interval; `-hist` = depth histogram per interval plus an `all` summary; `-counts` = just the overlap count (faster).

## samtools depth / coverage

```bash
samtools coverage in.bam                          # per-contig: rname..numreads covbases coverage(=breadth%) meandepth meanbaseq meanmapq
samtools depth -a -Q 20 -r chr1:1-100000 in.bam   # -a = report zero-depth positions; -Q = min MAPQ; -r = region
samtools depth -s in.bam                          # -s = count overlapping mate pair only ONCE (short-insert de-double-count)
samtools depth -a --reference ref.fa in.cram      # CRAM REQUIRES --reference
```

In `samtools coverage` the column literally named `coverage` is **breadth** (% bases >=1x), and `meandepth` is depth -- a contig with `coverage=9.7` and `meandepth=3.5` is 3.5x over only 9.7% of the contig (a localized pileup), NOT "9.7x coverage". `samtools depth` drops UNMAP/SECONDARY/QCFAIL/DUP by default (so duplicates are excluded -- but only if they were MARKED first). Without `-a`/`-aa`, zero-depth positions are omitted, so a naive `sum/lines` mean over-counts by dropping the zeros.

## Per-Method Failure Modes

### Reporting mean depth as the result
**Trigger:** citing "mean = 30x" as adequacy. **Mechanism:** mean is inflated by the repeat/rDNA tail and blind to GC/mappability holes. **Symptom:** a genome with large uncallable gaps looks fine. **Fix:** report median + breadth at the caller's threshold (mosdepth dist).

### samtools depth silent 8000 cap
**Trigger:** pre-1.13 samtools on high-depth loci (rDNA, mito, amplicon, ctDNA). **Mechanism:** old default `-d/-m` capped depth at 8000 and truncated with no warning. **Symptom:** depth plateaus near 8000. **Fix:** `samtools --version`; on old builds add `-d 0`; 1.13+ has no cap (flag ignored). Note `mpileup` has its own separate 8000 default.

### Mate-overlap double-counting
**Trigger:** naive per-base depth on short-insert libraries (amplicon, cfDNA, FFPE). **Mechanism:** the two mates of a short fragment both cover the overlap, counted twice but not independent. **Symptom:** locally doubled depth, corrupted/inflated VAFs. **Fix:** `samtools depth -s`, `genomecov -pc` (fragment), or mosdepth default -- and do NOT use mosdepth `--fast-mode`/`-x`, which turns the correction off.

### Coverage off an un-deduped BAM
**Trigger:** depth on a BAM whose duplicates were never marked. **Mechanism:** dedup-aware tools drop the DUP flag, but nothing was flagged. **Symptom:** inflated depth at amplified (often GC-extreme) loci, fatter right tail. **Fix:** Picard MarkDuplicates / `samtools markdup` FIRST, then measure.

### MAPQ filtering and repeat coverage
**Trigger:** choosing a MAPQ threshold without considering repeats. **Mechanism:** repeats give low MAPQ; `-Q 20+` makes repeat coverage vanish, MAPQ 0 lets multimappers pile up or smear. **Symptom:** repeats read as either empty or noisy -- no neutral choice. **Fix:** mask repeats (ENCODE blacklist / mappability) and report breadth over the MAPPABLE genome, not the whole genome.

### bedtools coverage -a/-b backwards
**Trigger:** pre-v2.24.0 muscle memory / old tutorials. **Mechanism:** semantics flipped to report stats for `-a` at v2.24.0. **Symptom:** well-formed output describing per-read instead of per-target stats. **Fix:** A = targets, B = reads; sanity-check the row count equals the target count.

### genomecov default misread / no -split on RNA-seq
**Trigger:** expecting a bedGraph from bare `genomecov`, or omitting `-split` on spliced reads. **Mechanism:** bare default is a histogram; without `-split` an N-CIGAR intron is counted as covered. **Symptom:** misparsed histogram, or every spliced gene appears fully covered across introns. **Fix:** add `-bg`/`-bga` for a track; always `-split` for spliced data.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| WGS germline ~30x mean -> ~95% of genome >= 20x | field convention (approx) | het-SNP sensitivity plateaus ~30x; frame as breadth, not mean |
| WES germline ~100x on-target -> ~90-95% target >= 10-20x | field convention (approx; ACMG-style, lab-dependent) | capture unevenness + off-target loss eat the raw mean |
| Somatic bulk tumor ~60-100x+ | field convention (approx) | low-VAF subclones need depth ~ 1/VAF; impure tumors need more |
| ctDNA/UMI panels 1000s-50000x raw | field convention (approx) | raw depth != usable depth after UMI collapse; report effective depth |
| Long-read WGS ~20-30x (HiFi ~30x, ONT SV ~20x+) | moving convention (approx) | flatter GC bias + better repeat mappability reach more genome per x |
| mean/median > ~1.1-1.2 = skewed | distribution diagnostic | the tail is inflating the mean; investigate dups/repeats/rDNA |
| Fano factor = 1 (Poisson ideal); real >> 1 | Lander & Waterman 1988 | overdispersion = evenness problem; deeper sequencing won't fill holes |
| Picard fold-80 ~1.3-2 good, >3 poor | practitioner heuristic (Picard defines only the metric) | fold extra sequencing to lift 80% of targets to the mean |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Depth plateaus at ~8000 | pre-1.13 samtools default cap | `samtools --version`; add `-d 0`; upgrade to 1.13+ |
| Inflated VAFs in amplicon/cfDNA | mate-overlap double-counting | `samtools depth -s` / `genomecov -pc` / mosdepth (not `--fast-mode`) |
| "30x" but variants missing in some genes | GC-shallow / uncallable holes hidden by mean | report breadth at threshold; mask blacklist; check fold-80 |
| `samtools coverage` "coverage" looks tiny | it is breadth %, not depth | read `meandepth` for depth; `coverage` = % bases >=1x |
| genomecov gives a histogram not a track | no `-bg`/`-bga` flag | add `-bga` (with zeros) or `-bg` |
| Every spliced gene fully covered | missing `-split` on RNA-seq | add `-split` to genomecov/coverage |
| bedtools coverage stats look per-read | `-a`/`-b` backwards (pre-2.24 habit) | A = targets, B = reads |
| CRAM depth errors / empty | missing reference | `samtools depth --reference ref.fa` |

## References

- Quinlan AR, Hall IM. 2010. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26:841-842.
- Pedersen BS, Quinlan AR. 2018. Mosdepth: quick coverage calculation for genomes and exomes. *Bioinformatics* 34:867-868.
- Danecek P, Bonfield JK, Liddle J, et al. 2021. Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
- Aird D, Ross MG, Chen WS, et al. 2011. Analyzing and minimizing PCR amplification bias in Illumina sequencing libraries. *Genome Biol* 12:R18.
- Benjamini Y, Speed TP. 2012. Summarizing and correcting the GC content bias in high-throughput sequencing. *Nucleic Acids Res* 40:e72.
- Amemiya HM, Kundaje A, Boyle AP. 2019. The ENCODE blacklist: identification of problematic regions of the genome. *Sci Rep* 9:9354.
- Lander ES, Waterman MS. 1988. Genomic mapping by fingerprinting random clones: a mathematical analysis. *Genomics* 2:231-239.

## Related Skills

- bedgraph-handling - bedGraph tracks this skill emits, and their normalization
- bigwig-tracks - Convert the coverage bedGraph to an indexed bigWig for browsers
- interval-arithmetic - Intersect coverage/callable BEDs with target regions
- alignment-files/pileup-generation - Per-base pileup upstream of depth and per-call DP
- alignment-files/bam-statistics - flagstat/idxstats and dup rate that explain coverage confounders
- chip-seq/chipseq-visualization - deepTools normalized coverage tracks for cross-sample comparison
- data-visualization/genome-tracks - Render the coverage tracks built here
