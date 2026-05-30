---
name: bio-genome-intervals-bedgraph-handling
description: Generates, normalizes, and converts bedGraph signal tracks (4-column chrom/start/end/value, 0-based half-open) with bedtools genomecov, deepTools bamCoverage/bamCompare/bigwigCompare, bedtools unionbedg, and UCSC bedGraphToBigWig. Covers why a raw coverage bedGraph is not comparable across samples until normalized, the CPM/RPKM/BPM/RPGC normalization menu and the conserved-total assumption that makes them wrong under a global perturbation, the strict sorted-non-overlapping-chrom.sizes bedGraphToBigWig contract that silently corrupts a bigWig, effective-genome-size selection, and bin-size aliasing. Use when building or normalizing a coverage/signal track from a BAM, comparing tracks across samples or conditions, converting bedGraph to a browser-ready bigWig, or diagnosing a track that looks plausible but reports wrong heights.
tool_type: mixed
primary_tool: deeptools
---

## Version Compatibility

Reference examples tested with: deeptools 3.5+, bedtools 2.31+, ucsc-bedgraphtobigwig 445+, pyBigWig 0.3.22+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- Python: `pip show <package>` then `help(module.function)` to check signatures

bedGraphToBigWig has a hard, under-advertised input contract: the bedGraph must be `LC_COLLATE=C`-sorted by chrom then start, contain non-overlapping intervals, and ship with a chrom.sizes derived from the exact assembly the reads were aligned to. deepTools effective-genome-size tables are occasionally updated between releases - re-check the installed version's table. If code throws an error, introspect the installed tool and adapt the example to match the actual API rather than retrying.

# bedGraph Handling

**"Make me a coverage/signal track I can compare across samples and load in a browser"** -> Generate a per-bin signal track, normalize it onto a common scale (or decide a spike-in is required), then convert the text bedGraph to an indexed bigWig under the strict sort/overlap/chrom.sizes contract.
- CLI: `bamCoverage -b s.bam -o s.bw --normalizeUsing RPGC --effectiveGenomeSize <N>`; `bedtools genomecov -ibam s.bam -bga`; `LC_COLLATE=C sort -k1,1 -k2,2n in.bdg | bedGraphToBigWig /dev/stdin chrom.sizes out.bw`
- Python: `pyBigWig.open('s.bw')` to read/extract; `bw.intervals(chrom, start, end)` returns the bedGraph rows

## The Single Most Important Modern Insight -- A Raw Coverage bedGraph Is a Library-Size Artifact, and the Wrong Normalization Is Worse Than None

Column 4 of a raw coverage bedGraph is not biology - it is sequencing depth. Two libraries of *identical* biology sequenced to different depths produce different heights, so any cross-sample statement ("more signal at this promoter in treatment") on un-normalized tracks is a category error. The modern path skips the text intermediate entirely: **deepTools bamCoverage takes BAM -> normalized bigWig in one step**, because bigWig is indexed, binary, random-access and bedGraph is flat text. Three load-bearing moves:

1. **Every library-size normalization (CPM/RPKM/BPM/RPGC=1x) assumes total signal is conserved across samples.** They all just rescale each library to a common total (per-million reads, or to 1x genome coverage). That model is correct when signal only *redistributes* locally - the usual case - and **actively wrong** when the perturbation changes *global* levels (histone-mark KD, BET-bromodomain inhibitor, global pol-II collapse). A genuine 3-fold global increase becomes, after CPM/RPGC, *no change* - the extra signal is spread thin and rescaled away. The model is unfalsifiable from the normalized data: forcing both libraries to the same total defines away any global difference. **Library-size normalization assumes the very thing under measurement does not happen.**
2. **There is no computational rescue for a global change after the fact.** The only fix is an external ruler decided AT THE BENCH - a spike-in of fixed foreign chromatin per cell (ChIP-Rx, Orlando 2014; defined reference epigenome, Bonhoure 2014) - scaled by the spike-in reads, not the sample reads. The wet-lab decision had to be made before sequencing; with no spike-in, the global scale is unrecoverable. The mechanics live in chip-seq/spike-in-normalization; the *decision* (could this perturbation change global levels?) belongs here, up front.
3. **bedGraph is scratch; bigWig is the artifact.** The text bedGraph is the last human-readable checkpoint - `awk '$4 > 1000'` to find blacklist pileups, confirm the sort/overlap invariants - before opaque binary. Inspect it, then ship bigWig. Never distribute a bedGraph as a final product: it is unindexed, so a browser reads the whole file to render any region.

## Normalization Taxonomy

| Method | What it assumes | When to use | When WRONG |
|--------|-----------------|-------------|------------|
| None | nothing (raw counts) | single-sample inspection only | any cross-sample comparison - depth confounds it |
| CPM | total mapped reads is the right denominator; total signal conserved | depth-only normalization; quick cross-sample on a common assay | a few high-coverage bins dominate (composition skew); global change |
| RPKM | as CPM plus bin length matters; total signal conserved | legacy default; depth + bin-length normalized | composition skew; global change; superseded by BPM for tracks |
| BPM (TPM-analog) | sum over all bins fixed at 1e6; total signal conserved | composition-aware cross-sample default; robust to a few dominant bins | global change (still a conserved-total rescale) |
| RPGC (1x) | mean genome-wide coverage = 1x; correct effective-genome-size; total signal conserved | field-standard ChIP/ATAC browser viewing; most interpretable height | wrong effective-genome-size (linear scaling error); global change |
| spike-in (external) | spike-in amount is constant per cell (a ruler that does not move) | global-level change plausible or under test | nothing computational - requires a bench step before sequencing |

All five library-size methods share one axiom: **total signal is conserved**. The decision is not which library-size method, it is whether library-size normalization is legitimate at all (see Decision Tree).

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| One BAM -> browser track, local redistribution | `bamCoverage --normalizeUsing RPGC --effectiveGenomeSize <N>` | one-step BAM->normalized bigWig; RPGC is the interpretable ChIP/ATAC standard |
| Cross-sample, composition skew likely | `bamCoverage --normalizeUsing BPM` | bins-per-million fixes the per-bin sum; robust to dominant bins |
| Global-level change plausible (KD/KO of a chromatin modifier, BET inhibitor) | spike-in -> chip-seq/spike-in-normalization | library-size normalization erases the global change by construction |
| RNA-seq coverage track | `bamCoverage --filterRNAstrand` or `genomecov -bga -split` | `-split`/strand handling so spliced reads do not paint introns |
| ChIP/ATAC track | add `--extendReads` (and `--centerReads` for footprints) | a read is a fragment END; raw read-end coverage is double-humped and wrong |
| Treatment vs input from raw BAMs | `bamCompare -b1 chip.bam -b2 input.bam --operation log2` | normalizes depth THEN does the arithmetic |
| Two already-normalized bigWigs | `bigwigCompare --operation log2` | arithmetic only - feeding un-normalized tracks manufactures a fake change |
| Stack N samples into a value matrix | `bedtools unionbedg -header -names ...` | union interval partition; feed the matrix to R/Python for testing |
| Sample-relatedness QC | `multiBigwigSummary bins` -> `plotCorrelation`/`plotPCA` | genome-wide value matrix for correlation/PCA |
| Need exact per-base arithmetic (not a browser) | keep bedGraph (`genomecov -bga`) | bedGraph is exact text; bigWig is binned/lossy |
| Convert finished bedGraph -> bigWig | `LC_COLLATE=C sort` then `bedGraphToBigWig` + matched chrom.sizes | the strict contract; inspect the text first |

## Generate a Normalized Track with bamCoverage (the modern default)

**Goal:** Turn one BAM into a normalized, browser-ready bigWig in a single command.

**Approach:** Let bamCoverage bin, normalize, and write bigWig directly; pick the normalization from the taxonomy, supply the effective-genome-size for RPGC, extend reads for ChIP/ATAC, and exclude chrX/chrM (and any spike-in contigs) from the scale-factor calculation.

```bash
BIN_SIZE=25                  # bp; smaller = finer + noisier + bigger. Match to feature width (sharp TF/ATAC 10-25; broad marks 50-200)
EFFGENOME=2913022398         # GRCh38 non-N length (faCount); use ONLY if multimappers were kept (see Effective Genome Size)

bamCoverage -b sample.bam -o sample.bw \
  --binSize $BIN_SIZE --normalizeUsing RPGC --effectiveGenomeSize $EFFGENOME \
  --extendReads --ignoreForNormalization chrX chrM -p 8
```

Defaults to verify: `--binSize 50`, `--normalizeUsing None`, `--scaleFactor 1.0`, `--extendReads` off, `--centerReads` off. For single-end ChIP supply the fragment length (`--extendReads 200`); paired-end infers it. `--scaleFactor` with `--scaleFactorsMethod None` is the hook for a bench-derived spike-in factor. `--outFileFormat bedgraph` writes the text form when the raw numbers are needed.

## Generate with bedtools genomecov (text, flexible, no normalization)

`-bg` collapses equal-coverage runs but **omits zero-coverage regions**; `-bga` additionally tiles zeros (use when downstream tools need explicit 0s). `-split` is mandatory for RNA-seq so spliced reads do not paint introns. `-scale 1000000/<nreads>` is a crude manual RPM; deepTools is preferred for real normalization.

```bash
bedtools genomecov -ibam sample.bam -bga -split > sample.bedgraph
```

## Convert bedGraph -> bigWig (the silent-corruption trap)

**Goal:** Produce a valid bigWig from a finished bedGraph without shipping a file that loads but lies.

**Approach:** C-locale-sort, guarantee non-overlapping intervals, derive chrom.sizes from the exact aligned-to FASTA, inspect the text, then convert.

```bash
samtools faidx ref.fa && cut -f1,2 ref.fa.fai > chrom.sizes   # chrom.sizes from the SAME FASTA the reads aligned to
LC_COLLATE=C sort -k1,1 -k2,2n sample.bedgraph > sample.sorted.bedgraph   # C locale: locale-aware sort triggers "is not case-sensitive sorted"
bedGraphToBigWig sample.sorted.bedgraph chrom.sizes sample.bw
```

If concatenation/merging introduced overlaps, collapse with an explicit aggregation BEFORE converting - and note `max` vs `mean` vs `sum` are different signals, there is no safe default:

```bash
bedtools merge -i sample.sorted.bedgraph -d 0 -c 4 -o max > sample.nonoverlap.bedgraph
```

bigWig round-trips losslessly: `bigWigToBedGraph sample.bw out.bedgraph` (optionally `-chrom=chr1 -start=1000 -end=2000`).

## Multi-Sample Arithmetic

**Goal:** Compare two tracks (treatment/input, two conditions) without letting a depth difference masquerade as biology.

**Approach:** From raw BAMs use bamCompare, which normalizes depth THEN applies the operation; only use bigwigCompare on bigWigs that are *already* on a common scale.

```bash
bamCompare -b1 chip.bam -b2 input.bam -o log2ratio.bw \
  --operation log2 --pseudocount 1 --binSize 25 --scaleFactorsMethod readCount
```

`--operation` (NOT `--ratio`) chooses log2/ratio/subtract/add/mean/reciprocal_ratio/first/second; default `log2`. `--scaleFactorsMethod readCount` (the default) scales by library size; `--scaleFactorsMethod SES` (signal-extraction scaling, Diaz 2012) instead estimates the factor from the shared background bins and is more robust than readCount for SHARP/punctate marks and TF ChIP where enrichment is a small genomic fraction; it DEGRADES for broad marks (H3K27me3/H3K9me3) where the diffuse enrichment cannot be cleanly separated from background, so use readCount (or spike-in) there. `--pseudocount` (default 1) prevents divide-by-zero in log2/ratio but pulls low-coverage bins toward 0 - a log2 track's apparent dynamic range is partly a pseudocount+bin-size artifact, do not read fold-changes off a browser track as measured. `bigwigCompare --skipZeroOverZero` drops bins that are 0 in both rather than flooding the output with `log2(1)=0`. Stack many samples and QC relatedness:

```bash
bedtools unionbedg -i s1.bdg s2.bdg s3.bdg -header -names s1 s2 s3 > matrix.txt   # inputs must be coordinate-sorted
multiBigwigSummary bins -b s1.bw s2.bw s3.bw -o scores.npz && plotCorrelation -in scores.npz --corMethod spearman --whatToPlot heatmap -o corr.png
```

## Read/Extract Signal with pyBigWig

```python
import pyBigWig

bw = pyBigWig.open('sample.bw')
mean_over_region = bw.stats('chr1', 1_000_000, 1_010_000, type='mean')[0]   # binned summary, not per-base
rows = bw.intervals('chr1', 1_000_000, 1_010_000)   # the underlying bedGraph rows: (start, end, value)
bw.close()
```

`bw.stats()`/`bw.values()` return what the bin resolution preserved, not a faithful per-base record - coarse bins silently change the values read back.

## Effective Genome Size (the two-table trap)

`--effectiveGenomeSize` feeds the RPGC scale factor and depends on the read-filtering regime. deepTools ships two tables that answer different questions:

| Build | Non-N length (faCount; multimappers KEPT) |
|-------|--------------------------------------------|
| GRCh38 | 2,913,022,398 |
| GRCh37 | 2,864,785,220 |
| GRCm38 (mm10) | 2,652,783,500 |
| dm6 | 142,573,017 |
| WBcel235 (C. elegans) | 100,286,401 |

When reads were instead **filtered to unique alignments / a MAPQ filter applied** (the common ChIP/ATAC case), use the read-length-dependent unique-k-mer value: GRCh38 is 2,701,495,711 (50 bp), 2,805,636,231 (100 bp), 2,862,010,428 (150 bp). The two GRCh38 numbers differ ~7% at short read length. RPGC scales linearly in this value, so the error cancels for *within-study* ratios but surfaces as a spurious constant fold-difference on cross-study integration (a public track, a collaborator's bigWig, a track made last year at a different read length). For non-model organisms there is no table - estimate it (`faCount` for non-N length, or unique-k-mers on the assembly).

## Per-Method Failure Modes

### Comparing un-normalized tracks across samples
**Trigger:** browser-comparing or quantifying raw coverage bedGraphs/bigWigs. **Mechanism:** column 4 scales with library size. **Symptom:** the deeper library looks like it has "more signal" everywhere. **Fix:** normalize during bamCoverage; never compare `--normalizeUsing None` tracks.

### Conserved-total assumption under a global change
**Trigger:** CPM/RPKM/BPM/RPGC on a perturbation that shifts global levels (chromatin-modifier KD/KO, BET inhibitor). **Mechanism:** every library-size method forces total signal to a constant. **Symptom:** a real global increase reads as no change; tracks look identical. **Fix:** spike-in decided at the bench -> chip-seq/spike-in-normalization. No computational rescue exists.

### Unsorted / overlapping input -> corrupt bigWig
**Trigger:** `bedGraphToBigWig` on non-C-sorted or overlapping input, or chrom.sizes from the wrong assembly. **Mechanism:** the contract is enforced inconsistently - some violations error, others build a bigWig that loads and shows wrong heights or silently drops chromosomes. **Symptom:** `is not case-sensitive sorted`, `overlapping regions`, `end coordinate bigger than`, or a silently wrong/incomplete track. **Fix:** `LC_COLLATE=C sort`; `bedtools merge -c 4 -o max/mean/sum`; chrom.sizes from the exact aligned-to FASTA; harmonize `chr1` vs `1`.

### Effective-genome-size drift
**Trigger:** grabbing the round 2.9e9 GRCh38 value regardless of multimapper filtering, or reusing a value across read lengths/assemblies. **Mechanism:** RPGC scales linearly in the value; the two tables differ ~7%. **Symptom:** invisible within a study; a spurious constant fold-difference on cross-study integration. **Fix:** match the value to the read length AND filtering regime; estimate it for non-model organisms.

### Bin-size aliasing
**Trigger:** a bin wider than ~half the feature, or comparing tracks built at different binSizes. **Mechanism:** binSize is a low-pass filter chosen once; a feature narrower than ~2 bins is averaged down or straddles a boundary (a phase artifact - replicates disagree by bin alignment). **Symptom:** sharp peaks shrink or split; bin-for-bin ratios meaningless at boundaries. **Fix:** match binSize to feature width (sharp TF/ATAC 10-25 bp, broad marks 50-200 bp); compared tracks MUST share binSize. `--smoothLength` is cosmetic, it cannot recover discarded information.

### ChIP/ATAC track without extendReads
**Trigger:** bamCoverage/genomecov on ChIP/ATAC without `--extendReads`. **Mechanism:** a read marks a fragment END, not the fragment. **Symptom:** double-humped peaks with a central dip; biased boundaries and quantification. **Fix:** `--extendReads` (paired-end infers; single-end supply the fragment length). RNA-seq mirror trap: without `-split` spliced reads paint introns.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| binSize default 50 bp; sharp TF/ATAC 10-25 bp, broad marks 50-200 bp | deepTools default + feature-width matching | binSize is a low-pass filter; finer is noisier/bigger, coarser aliases sharp features |
| GRCh38 effGenome 2,913,022,398 (multimappers kept) | deepTools faCount table | non-N genome length for the RPGC denominator |
| GRCh38 effGenome 2.70-2.86e9 by read length (unique alignments) | deepTools unique-k-mer table | ~7% below the non-N value; use when MAPQ/uniqueness-filtered |
| pseudocount default 1 (log2/ratio) | deepTools default | prevents divide-by-zero; biases low-coverage bins toward 0 |
| single-end fragment length ~200 bp (`--extendReads 200`) | typical sonicated ChIP fragment | wrong value distorts peak width; paired-end infers it |
| compared tracks must share binSize | signal-processing constraint | different grids make bin-for-bin ratios meaningless |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `is not case-sensitive sorted` | locale-aware sort | `LC_COLLATE=C sort -k1,1 -k2,2n` (works on a login node, fails in the scheduler when `$LC_*` differ) |
| `overlapping regions in bedGraph file` | concatenated/merged tracks | `bedtools merge -c 4 -o max/mean/sum` (choose the aggregation deliberately) |
| `end coordinate N bigger than ...` | chrom.sizes from a different assembly/patch | derive chrom.sizes from the exact aligned-to FASTA (`samtools faidx` + `cut -f1,2`) |
| Whole chromosomes missing from the bigWig, no error | `chr1` vs `1` / `MT` vs `chrM` naming mismatch | harmonize naming across bedGraph and chrom.sizes |
| Track line breaks sort/conversion | `track type=bedGraph ...` header row | remove the track line before sort/bedGraphToBigWig |
| RPGC normalization fails | `--effectiveGenomeSize` not supplied | pass the correct value for the build, read length, and filtering |
| Spliced reads paint introns | no `-split` (genomecov) / wrong RNA mode | `genomecov -bga -split` or bamCoverage `--filterRNAstrand` |

## References

- Ramírez F, Ryan DP, Grüning B, et al. 2016. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Res* 44:W160-W165.
- Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. 2010. BigWig and BigBed: enabling browsing of large distributed datasets. *Bioinformatics* 26:2204-2207.
- Quinlan AR, Hall IM. 2010. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics* 26:841-842.
- Orlando DA, Chen MW, Brown VE, et al. 2014. Quantitative ChIP-Seq normalization reveals global modulation of the epigenome. *Cell Reports* 9:1163-1170.
- Bonhoure N, Bounova G, Bernasconi D, et al. 2014. Quantifying ChIP-seq data: a spiking method providing an internal reference for sample-to-sample normalization. *Genome Res* 24:1157-1168.
- Diaz A, Park K, Lim DA, Song JS. 2012. Normalization, bias correction, and peak calling for ChIP-seq. *Stat Appl Genet Mol Biol* 11:Article 9.

## Related Skills

- coverage-analysis - Per-base depth generation and distribution-vs-mean diagnostics feeding bedGraph tracks
- bigwig-tracks - Reading, extracting, and writing the bigWig deliverable this skill produces
- chip-seq/spike-in-normalization - The bench-decided external-reference scaling when a global change makes library-size normalization wrong
- chip-seq/chipseq-visualization - Render the normalized signal tracks built here
- atac-seq/footprinting - Consumes high-resolution coverage/bigWig signal over motif sites
- data-visualization/genome-tracks - Render the bedGraph/bigWig tracks for figures
