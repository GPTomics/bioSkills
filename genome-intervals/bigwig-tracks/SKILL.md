---
name: bio-genome-intervals-bigwig-tracks
description: Reads, queries, and writes bigWig indexed binary signal tracks (coverage, fold-change, conservation, methylation-rate) with pyBigWig (Python) and the UCSC Kent tools (bedGraphToBigWig, bigWigToBedGraph, bigWigInfo, bigWigSummary, bigWigAverageOverBed) and deepTools (multiBigwigSummary, computeMatrix, bigwigCompare). Covers the central trap that a wide query returns a precomputed zoom-level summary (by default the mean, which annihilates narrow peaks) not per-base data, when exact=True/values() is mandatory, the NaN-not-zero gap-handling fork, choosing mean vs max vs sum vs coverage by biological question, and the sorted-bedGraph plus chrom.sizes build requirement. Use when extracting signal at regions, computing mean signal per gene/peak, building a browser track from bedGraph, comparing tracks, or building TSS/gene-body metaprofiles.
tool_type: mixed
primary_tool: pyBigWig
---

## Version Compatibility

Reference examples tested with: pyBigWig 0.3.22+, numpy 1.26+, ucsc-bedgraphtobigwig/ucsc-tools 469+, deeptools 3.5+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` (or `bigWigInfo` with no args for usage) then `<tool> --help` to confirm flags
- Python: `pip show pyBigWig` then `help(pyBigWig.bigWigFile.stats)` to check signatures

Building any bigWig needs a chrom.sizes file (`name<TAB>length`) and a coordinate-sorted bedGraph; pyBigWig's numpy return path requires numpy present at compile time. If code throws an error, introspect the installed tool and adapt rather than retrying.

# BigWig Tracks

**"Get the signal from my bigWig over these regions / build a browser track."** -> Query an indexed binary signal track, choosing the summary statistic and exactness that match the biological question, or build one from a sorted bedGraph + chrom.sizes.
- CLI: `bigWigAverageOverBed in.bw regions.bed out.tab`, `bigWigSummary in.bw chr s e N -type=max`, `bedGraphToBigWig in.sorted.bedGraph chrom.sizes out.bw`, `bigWigInfo in.bw`
- Python: `bw=pyBigWig.open('x.bw')`, `bw.stats(chr,s,e,type='max',exact=True)`, `bw.values(chr,s,e,numpy=True)`, `bw.intervals(chr,s,e)` (pyBigWig)

## The Single Most Important Modern Insight -- A Wide Query Returns a Zoom-Level Summary, Not the Underlying Data

bigWig is fast (Kent 2010) because it stores, alongside base-resolution values, a ladder of precomputed zoom levels holding per-bin sum/sumSquared/min/max/nBasesCovered. A B+ tree resolves the chromosome, an R-tree (cirTree) finds the data blocks in O(log n), per-block zlib keeps it ~10x smaller than bedGraph, and **the zoom ladder answers a wide region in near-constant time by reading a precomputed summary instead of the base data.** That speed is bought with two stacked approximations, both ON by default, and a third trap at the moment the signal is reduced to a single number:

1. **WHICH statistic.** Over one wide bin `mean` (the default) dilutes a narrow tall feature toward background: a 200 bp ChIP summit of 500 in a 1 Mb sea of 1 averages to ~1.1 -- indistinguishable from background, while `type='max'` returns 500. **Same file, same coordinates, opposite conclusions, decided by the named statistic.** `mean` is faithful for broad features (domains, gene-body coverage) and a lie for narrow ones. `max`=peak height, `sum`=total amount (scales with width), `coverage`=fraction of bases with any data (ignores magnitude), `std`=variability.
2. **WHERE the number comes from.** `exact=False` (the pyBigWig default, and what `bigWigSummary` and every zoomed-out browser do) computes from the nearest zoom level, not base data. Fine for exploration and broad features; **`exact=True` (or `values()`) is mandatory whenever a number enters a result** -- a per-region average in a table, a threshold call, anything a reviewer recomputes. Plausible-but-zoom-approximated is the worst failure: it does not error, it rounds the biology.
3. **NaN is NOT zero.** Uncovered positions are no-data, surfaced as `NaN` in `values()` and as gaps between `intervals()` runs -- never 0. On a region 30% covered at signal 10: `np.mean` -> **NaN** (poisons), `np.nanmean` -> **10** (covered-only, = `bigWigAverageOverBed` `mean` column), gaps-as-zero (`np.nan_to_num().mean()`, deepTools `--missingDataAsZero`, the `mean0` column) -> **3.0**. A >3x swing in the headline number, and which is correct is **biological**: coverage/read-depth tracks -> gaps are zero (`mean0`); rate/ratio tracks (methylation %, log2FC, conservation) -> gaps are undefined (`mean`/`nanmean`).

Name the biological question first; the statistic, the `exact` flag, and the gap-handling then follow deterministically. Left on default, all three conspire to hand back a fast, confident, wrong answer.

## Tool Taxonomy

| Tool | Role | Mechanism | When |
|------|------|-----------|------|
| pyBigWig | Python read/write | C-extension over libBigWig; `stats`/`values`/`intervals`/`addEntries` | inside a Python pipeline; custom per-region extraction; writing a bigWig |
| bigWigAverageOverBed | mean signal per BED region | one row per feature; `name,size,covered,sum,mean0,mean` | the right tool for "average signal per gene/peak"; gives both mean0 and mean |
| bigWigSummary | region -> N equal bins | reads zoom levels (like `exact=False`); `-type=mean/min/max/std/coverage` | quick binned profile at the command line |
| bigWigInfo | header/stats sanity check | version, zoom-level count, basesCovered, min/max/mean without parsing data | first thing to run on an unfamiliar file |
| bedGraphToBigWig / wigToBigWig | build bigWig | needs sorted input + chrom.sizes | converting a coverage bedGraph/WIG to a track |
| bigWigToBedGraph / bigWigToWig | bigWig -> text | `-chrom/-start/-end` for a sub-region | exact arithmetic; inspecting values as text |
| multiBigwigSummary | score matrix across many bigWigs | mean per bin over genome `bins` or a `BED-file` | track correlation/PCA (-> `plotCorrelation`/`plotPCA`) |
| computeMatrix | signal across many regions | `reference-point` (TSS/peak center) or `scale-regions` (gene body) | metaprofiles/heatmaps (-> `plotHeatmap`/`plotProfile`) |
| bigwigCompare | combine two bigWigs bin-by-bin | `--operation log2/ratio/subtract/...` `--pseudocount` | a log2(IP/input) or (treat-control) track |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Mean signal per gene/peak (one number per BED row) | `bigWigAverageOverBed` | purpose-built; pick `mean` (covered-only) vs `mean0` (gaps as zero) deliberately |
| Peak height / "is there a binding event here?" | `stats(type='max')` or `bigWigSummary -type=max` | `mean` dilutes a narrow peak to background |
| Total signal over an exon/gene (an amount) | `stats(type='sum')` or the `sum` column | extensive quantity; do not use `mean` for a total |
| A number going into a table/threshold | `stats(..., exact=True)` or `values()` | the default `exact=False` reads zoom levels, not base data |
| Per-base values for plotting/analysis | `values(numpy=True)` | one number per base; `nan` for gaps -> `np.nanmean`, never `np.mean` |
| "Is this region even assayed/mappable?" | `stats(type='coverage')` | fraction with data; a different axis from magnitude |
| Compare many tracks (correlation/PCA) | `multiBigwigSummary` -> `plotCorrelation`/`plotPCA` | mean per bin; bin size matters |
| Metaprofile/heatmap over TSS or gene bodies | `computeMatrix reference-point`/`scale-regions` -> `plotHeatmap`/`plotProfile` | match anchored-point vs whole-body mode |
| A ratio/difference track | `bigwigCompare --operation log2 --pseudocount` | pseudocount only meaningful for log2/ratio |
| Build a normalized coverage track from a BAM | -> chip-seq/chipseq-visualization or atac-seq/footprinting (deepTools `bamCoverage`) | generation is upstream; library-size normalization lives there |
| Render the track in a browser figure | -> data-visualization/genome-tracks | pyGenomeTracks/IGV; zoom-out IS the summary trap made visual |
| Discrete features (peaks/genes), not signal | bigBed, not bigWig | continuous-vs-interval; one interval per base defeats the format |

## Inspect a File Before Trusting It

```bash
bigWigInfo coverage.bw                  # version, # zoom levels, basesCovered, min/max/mean/std
bigWigInfo -chroms coverage.bw          # chrom names + lengths (the file carries its own chrom list)
```

A zero or low zoom-level count means a zoomed-out browser will read base data slowly (or, with `maxZooms=0`, IGV breaks). `basesCovered` far below the genome size means most positions are no-data (NaN), which makes the `mean`-vs-`mean0` choice below load-bearing.

## Mean Signal per Region (the most common task)

**Goal:** Compute one signal number per gene/peak, choosing covered-only vs gaps-as-zero by the track's biology.

**Approach:** Use the purpose-built `bigWigAverageOverBed` (BED needs a unique name column) and read the right output column -- `mean` (covered bases only) for rate/ratio tracks, `mean0` (uncovered counted as zero) for coverage/depth tracks.

```bash
# BED4+ with a UNIQUE name in column 4; output columns: name size covered sum mean0 mean
bigWigAverageOverBed coverage.bw genes.bed signal_per_gene.tab
# -> read $6 (mean, covered-only) for methylation/log2FC; $5 (mean0, gaps=0) for read depth
```

The pyBigWig equivalent, when the extraction is inside a Python pipeline -- note `exact=True` because these numbers enter a result, and an explicit gap decision:

```python
import pyBigWig
import numpy as np

bw = pyBigWig.open('coverage.bw')
GAPS_ARE_ZERO = False   # True for read-depth/coverage tracks; False for rate/ratio (methylation, log2FC)

def region_signal(chrom, start, end):
    v = bw.values(chrom, start, end, numpy=True)                         # per-base, nan for gaps
    if GAPS_ARE_ZERO:
        return float(np.nan_to_num(v).mean())                           # mean0: gaps counted as 0 (= bigWigAverageOverBed mean0)
    return np.nanmean(v) if not np.all(np.isnan(v)) else float('nan')    # covered-only (= bigWigAverageOverBed mean); stats(type='mean') is also covered-only, NOT mean0
```

## Peak Height vs Total vs Coverage (statistic = question)

```python
import pyBigWig
bw = pyBigWig.open('chip.bw')
region = ('chr1', 1_000_000, 2_000_000)

peak  = bw.stats(*region, type='max', exact=True)[0]        # binding-event height; mean would dilute it
total = bw.stats(*region, type='sum', exact=True)[0]        # total signal (amount; scales with width)
assayed = bw.stats(*region, type='coverage', exact=True)[0] # FRACTION of bases with any data (0..1), ignores magnitude
profile = bw.stats(*region, type='max', nBins=1000)         # 1000-bin max profile; nBins keeps narrow features visible
```

`stats()` returns a list of length `nBins` (default 1). `type` is one of `mean`(default)/`max`/`min`/`coverage`/`std`/`sum`. Use `max` with `nBins>1` to see narrow features across a wide window; a single-bin `mean` over a megabase buries every peak.

## Per-Base Values (NaN is not zero)

```python
import pyBigWig
import numpy as np
bw = pyBigWig.open('coverage.bw')

v = bw.values('chr1', 1_000_000, 1_001_000, numpy=True)   # list by default; numpy=True -> ndarray, nan for gaps
covered_mean = np.nanmean(v)                               # ignores gaps (= bigWigAverageOverBed mean)
depth_mean = np.nan_to_num(v).mean()                       # gaps counted as zero (= mean0); only for depth tracks
raw = bw.intervals('chr1', 1_000_000, 1_001_000)          # [(start,end,value),...] the unresampled stored runs
bw.close()
```

## Build a Valid bigWig

**Goal:** Turn a coverage bedGraph into an indexed, browser-ready bigWig.

**Approach:** Coordinate-sort the bedGraph, supply a chrom.sizes whose names and lengths match the bedGraph exactly, and run `bedGraphToBigWig` (which builds the index + zoom levels).

```bash
sort -k1,1 -k2,2n coverage.bedGraph > coverage.sorted.bedGraph   # bedGraphToBigWig REQUIRES sorted, non-overlapping input
cut -f1,2 reference.fa.fai > chrom.sizes                          # or fetchChromSizes hg38 > chrom.sizes
bedGraphToBigWig coverage.sorted.bedGraph chrom.sizes coverage.bw
```

Writing directly with pyBigWig -- `addHeader` (ordered chrom list) MUST precede `addEntries`, and entries must be added in sorted (chrom, start) order matching the header:

```python
import pyBigWig
bw = pyBigWig.open('out.bw', 'w')
bw.addHeader([('chr1', 248956422)])   # ordered (name,length); maxZooms default 10; maxZooms=0 disables zoom and breaks IGV
bw.addEntries(['chr1'], [0], ends=[100], values=[1.5])           # mode (a) variable intervals
# mode (b) variableStep: bw.addEntries('chr1', [0,100], values=[1.5,2.3], span=20)
# mode (c) fixedStep:    bw.addEntries('chr1', 0, values=[1.5,2.3], span=20, step=30)
bw.close()                            # close() builds the R-tree index + zoom ladder
```

## Compare and Profile Tracks (deepTools)

```bash
bigwigCompare -b1 treat.bw -b2 control.bw -o log2ratio.bw --operation log2 --pseudocount 1   # NOT --ratio (older flag name)
multiBigwigSummary BED-file -b a.bw b.bw -o scores.npz --BED regions.bed                      # then plotCorrelation/plotPCA
computeMatrix reference-point -S signal.bw -R tss.bed -b 2000 -a 2000 -o matrix.gz            # anchored on TSS
plotHeatmap -m matrix.gz -o heatmap.png
```

Both `bigwigCompare` and `multiBigwigSummary` use mean-per-bin, so the zoom-level dilution caveat above applies; `computeMatrix --missingDataAsZero` is the same NaN-vs-zero fork inside deepTools.

## Per-Method Failure Modes

### Wide mean read as the peak
**Trigger:** `bw.stats(chrom, start, end)` (default `type='mean'`) over a region wide relative to the feature. **Mechanism:** the mean dilutes a narrow tall peak toward background. **Symptom:** "no signal here" that a browser zoom-in contradicts. **Fix:** use `type='max'` (or `nBins>1`, or `values()`) for narrow features.

### exact=False leaks zoom approximations into a result
**Trigger:** shipping default-`exact` `stats()` numbers into a table/threshold. **Mechanism:** the value is computed from the nearest zoom level, not base data, at up to 16x coarser granularity. **Symptom:** plausible numbers a reviewer cannot reproduce. **Fix:** pass `exact=True` (or use `values()`/`bigWigAverageOverBed`) whenever a number enters a result.

### Averaging NaN as zero (or poisoning to NaN)
**Trigger:** `np.mean(bw.values(...))` over a track with gaps, or `mean0` on a rate track. **Mechanism:** `np.mean` poisons to NaN; `nan_to_num`/`mean0`/`--missingDataAsZero` averages real gaps as zeros. **Symptom:** a >3x swing or a NaN where a number was expected. **Fix:** decide biologically -- coverage track -> `mean0`/zero; rate/ratio track -> `mean`/`np.nanmean`.

### addEntries before addHeader (or out of order)
**Trigger:** writing entries before the header, or in non-sorted order. **Mechanism:** the chrom list and offsets must exist and be ordered before data is appended. **Symptom:** runtime error or a corrupt file. **Fix:** `addHeader([(chrom,length),...])` first, add entries in (chrom, start) order matching the header; `close()` to finalize.

### chrom.sizes / naming mismatch on build
**Trigger:** `bedGraphToBigWig` with a chrom.sizes from a different assembly or naming (`chr1` vs `1`). **Mechanism:** the builder validates intervals against chrom lengths. **Symptom:** `chromosome not found`, or silently dropped/truncated intervals. **Fix:** derive chrom.sizes from the same reference (`cut -f1,2 ref.fa.fai`); harmonize naming; sort first.

### Forcing peaks into a bigWig (or dense signal into bigBed)
**Trigger:** storing called peaks as a bigWig. **Mechanism:** bigWig is continuous signal; discrete features with per-feature metadata belong in bigBed. **Symptom:** lost boundaries/names, or an enormous one-interval-per-base file. **Fix:** signal -> bigWig; intervals/features -> bigBed.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `exact=True` when a number enters a result | pyBigWig design | default `exact=False` reads zoom levels (up to ~16x coarser than data); fine for exploration only |
| Zoom ladder: smallest bin ~16x mean interval size, each level 4x the previous | Kent 2010 convention | the resolution at which a wide query is answered; `bigWigInfo -zooms` shows the actual levels |
| Index < ~1% of data; ~10x smaller than bedGraph | Kent 2010 | order-of-magnitude; exact ratio is data-dependent (sparse vs dense) |
| bedGraph must be sorted `-k1,1 -k2,2n`, non-overlapping | bedGraphToBigWig requirement | signal is a function (one value per base); unsorted/overlapping input errors out |
| computeMatrix flank `-b/-a` 2000-3000 bp at TSS | metaprofile convention | captures promoter-proximal signal; widen for distal features; state the value used |
| bin size (e.g. 10-50 bp) | resolution vs file size | finer bins preserve narrow features but enlarge the file; state the bin when reading values back |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Region reads flat but browser shows a peak | wide `mean` query diluted the peak | use `type='max'`, more bins, or zoom to feature resolution |
| Per-region numbers a reviewer cannot reproduce | `exact=False` zoom approximation | re-extract with `exact=True` / `bigWigAverageOverBed` |
| `np.mean` returns NaN | gaps in the track (NaN, not 0) | `np.nanmean`, or `np.nan_to_num` if gaps are biologically zero |
| `mean` and `mean0` differ a lot in `bigWigAverageOverBed` | track is sparsely covered | pick the column by biology (depth -> mean0; rate -> mean) |
| `bedGraphToBigWig` errors / drops intervals | unsorted input or chrom-name/length mismatch | `sort -k1,1 -k2,2n`; match chrom.sizes to the reference |
| IGV will not render zoom-out | bigWig built with `maxZooms=0` | rebuild with zoom levels (default 10) |
| `bigwigCompare` rejects `--ratio` | flag renamed | use `--operation log2` |

## References

- Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. 2010. BigWig and BigBed: enabling browsing of large distributed datasets. *Bioinformatics* 26:2204-2207.
- Ramirez F, Ryan DP, Gruning B, Bhardwaj V, Kilpert F, Richter AS, Heyne S, Dundar F, Manke T. 2016. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Res* 44:W160-W165.
- pyBigWig (Devon Ryan / deepTools project) - C-extension wrapping libBigWig; no journal paper, see https://github.com/deeptools/pyBigWig
- UCSC Kent utilities (bedGraphToBigWig, bigWigInfo, bigWigSummary, bigWigAverageOverBed) - https://github.com/ucscGenomeBrowser/kent; cite Kent 2010 for the format.

## Related Skills

- bedgraph-handling - The text bedGraph this skill converts to/from, and exact-arithmetic alternative
- coverage-analysis - Generates the per-base depth/bedGraph that becomes a bigWig
- bed-file-basics - The region BED files passed to bigWigAverageOverBed/computeMatrix
- chip-seq/chipseq-visualization - Generates normalized tracks (bamCoverage) and renders computeMatrix metaprofiles
- atac-seq/footprinting - Consumes bigWig signal over motif sites
- data-visualization/genome-tracks - Renders the bigWig in a browser figure (where zoom-out is the summary trap made visual)
