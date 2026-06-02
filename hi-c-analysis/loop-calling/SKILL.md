---
name: bio-hi-c-analysis-loop-calling
description: Detects focal chromatin loops (point interactions / corner-dots) in balanced Hi-C and Micro-C contact maps and aggregates/validates a loop set. Covers de-novo calling with cooltools dots (HiCCUPS-style 4-background local enrichment with lambda-chunked FDR), chromosight (template-correlation), and Mustache (scale-space blob detection); aggregate peak analysis (APA) via cooltools pileup for confirmation; the depth/resolution prerequisite (de-novo needs ~5-10kb resolution = hundreds of millions to billions of valid pairs); consensus across callers and convergent-CTCF support as validation; and differential loops via union anchors plus chromosight quantify. Use when calling chromatin loops or dots from a cooler, deciding whether a map is deep enough to call de-novo vs running APA on known CTCF/cohesin anchors, building an aggregate peak pileup, comparing loops across conditions, or validating loop calls. For HiChIP/PLAC-seq/PCHi-C protein-anchored data use FitHiChIP/MAPS, not dots.
tool_type: mixed
primary_tool: cooltools
---

## Version Compatibility

Reference examples tested with: cooltools 0.7+, cooler 0.10+, bioframe 0.7+, chromosight 1.6+, mustache 1.3+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The `.cool` must be BALANCED before calling loops -- `dots`/`pileup` read the `weight` column and raw counts are unsupported. An `.mcool` is multi-resolution; pass a single-resolution URI (`file.mcool::/resolutions/10000`), not the bare `.mcool`. cooltools changed signatures around 0.5 -> 0.7 (view_df/expected_value_col conventions); verify `help(cooltools.dots)` for the installed version. The `view_df` passed to `expected_cis` MUST be the same one passed to `dots`/`pileup`.

# Chromatin Loop Calling

**"Where are the focal loops (CTCF/cohesin corner-dots, E-P contacts) in my Hi-C map?"** -> Test each off-diagonal pixel for focal enrichment against its local background (on a balanced, expected-normalized matrix), control FDR, then validate the set by aggregation and orthogonal support.
- Python: `cooltools.dots(clr, expected=cooltools.expected_cis(clr, view_df=arms), view_df=arms)`
- CLI: `chromosight detect --pattern loops --min-dist 20000 --max-dist 2000000 sample.cool::/resolutions/5000 out`

## The Single Most Important Modern Insight -- Loop Calling Is Depth-Limited, Not Algorithm-Limited; the First Question Is "How Deep Is the Map?"

The caller choice is second-order. The dominant variable in whether loops are found at all is sequencing depth / map resolution. Rao 2014 needed ~4.9 BILLION contacts in GM12878 to reach 1kb bins and call ~10,000 loops; robust de-novo calling realistically wants 5-10kb resolution, which is hundreds of millions to billions of valid cis pairs. Below that, every caller returns near-nothing or noise, and tuning the FDR will not rescue it. So the workflow forks on depth before any tool is chosen:

1. **Deep map (>=~500M-1B valid pairs, 5-10kb resolution):** de-novo calling is licensed. Run cooltools `dots` (or chromosight / Mustache), then validate (see below).
2. **Shallow map:** do NOT de-novo call. Run **APA / pileup on a KNOWN anchor set** -- loops imported from a deep reference map, or anchor pairs built from CTCF/cohesin ChIP-seq peaks. This is the single most important practical reframe in the skill: shallow data can still *confirm and quantify* a hypothesized loop set even when it cannot *discover* one.

Two corollaries that follow directly:

**De-novo calling DISCOVERS; APA CONFIRMS -- never conflate them.** APA aggregates many putative loops to surface mean signal no individual loop could pass FDR for. An enriched APA center pixel proves "this SET of pairs is enriched on average"; it does NOT prove any single pair is a loop and it cannot discover new loops. Presenting an APA pileup as evidence that "these loops exist" is the classic abuse. And the APA score is meaningless without a **corner control** -- center pixel divided by an off-diagonal corner block of the flank is the on-vs-off measurement; the bare center value alone says nothing.

**Loops form between convergent CTCF motifs -- biology AND a validation filter.** Loops preferentially link two CTCF motifs in CONVERGENT orientation (Rao 2014 observation; de Wit 2015, Sanborn 2015 extrusion mechanism; proven by CTCF-site inversion experiments that kill or reroute the loop). A called corner-dot whose two anchors carry convergent CTCF motifs is high-confidence; one with no CTCF/anchor support on a shallow map is likely a false positive. Not all loops are CTCF loops (E-P and polycomb loops exist), so convergent-CTCF is a strong positive filter, not a universal requirement.

## Loop-Caller Taxonomy

| Tool | Philosophy | Mechanism | When |
|------|-----------|-----------|------|
| cooltools `dots` | local enrichment (CPU HiCCUPS) | pixel must beat 4 local backgrounds (donut/horizontal/vertical/lower-left); Poisson p; lambda-binned BH-FDR | cooler/.mcool pipelines, the modern default; pure-CPU |
| Juicer HiCCUPS | local enrichment (GPU original) | same 4-kernel model on `.hic`; CUDA-bound | `.hic`/Juicer ecosystem with a GPU available |
| chromosight `detect` | template correlation | Pearson correlation of a loop/border/stripe kernel vs each window | want loops AND borders AND stripes from one engine; Micro-C-friendly |
| Mustache | scale-space blobs | Difference-of-Gaussians across scales; multi-scale catches loops of different sizes | mixed loop sizes, kb-resolution Micro-C, recovers more E-P/ChIA-PET loops |
| SIP | image processing | Gaussian blur + regional-max + watershed | `.hic` image-based alternative |
| cooltools `pileup` (APA) | CONFIRMATION, not discovery | aggregate snippets centered on an anchor set; measure center vs corner | validate/quantify a loop SET; works on shallow maps |

Forcato 2017 (*Nat Methods* 14:679) is the canonical finding that loop callers show LOW pairwise overlap and poor replicate reproducibility -- far worse than TAD callers. Practical consequence: a loop called by only one tool is suspect. Trust comes from consensus across >=2 callers plus orthogonal support (convergent CTCF, ChIA-PET/HiChIP), not from any single tool's list length.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Shallow map (tens of M pairs) | APA/pileup on KNOWN anchors (CTCF/cohesin ChIP or reference loops) -- STOP de-novo | callers return noise below ~5-10kb resolution |
| Deep cooler/.mcool, CPU only | `cooltools.dots` (default) | pure-CPU HiCCUPS reimplementation on balanced cooler |
| Deep `.hic` with a GPU | Juicer HiCCUPS | CUDA original built for billion-contact `.hic` scans |
| Mixed loop sizes / kb Micro-C | Mustache | scale-space natively spans loop sizes; sub-5kb-friendly |
| Want stripes/borders too | chromosight (swap `--pattern`) | same template engine; stripes are a SEPARATE class, not loops |
| Validate a call set | `cooltools.pileup` -> APA score vs corner control | aggregate enrichment + visual QC of the dot |
| Confirm anchors are real loops | convergent-CTCF check -> chip-seq/peak-annotation, atac-seq/footprinting | extrusion loops carry convergent CTCF motifs |
| Annotate loop anchors | -> chip-seq/peak-annotation, atac-seq/enhancer-gene-linking | E-P / TF context lives there |
| Anchor-overlap enrichment p-value | -> genome-intervals/overlap-significance | turn an anchor-overlap count into a permutation test |
| Two conditions, loop strength shift | union anchors -> chromosight `quantify` per condition -> test delta (or `diff_mustache`) | NO bin-level DESeq for loops; quantify a fixed coordinate set |
| HiChIP / PLAC-seq / PCHi-C | FitHiChIP / MAPS / HiC-DC+ -> chip-seq/peak-calling | protein-anchored, coverage-biased; HiCCUPS null is wrong |

## De-Novo Loop Calling with cooltools dots

**Goal:** Discover focal loops genome-wide on a deep, balanced map with honest FDR control.

**Approach:** Build chromosome-arm regions, compute the distance-decay expected on those arms, then run `dots` -- which convolves four local-background kernels and runs Benjamini-Hochberg FDR independently within geometrically-spaced lambda-bins of locally-adjusted expected. The arms `view_df` must be identical for `expected_cis` and `dots`.

```python
import cooler, cooltools, bioframe

clr = cooler.Cooler('matrix.mcool::/resolutions/10000')   # 10kb: a realistic de-novo floor; finer needs more depth
arms = bioframe.make_viewframe(clr.chromsizes)   # or cooltools.lib.read_viewframe_from_file('hg38_arms.bed', clr) for per-arm
expected = cooltools.expected_cis(clr, view_df=arms, nproc=4)   # distance-matched background; same view as dots
loops = cooltools.dots(
    clr, expected=expected, view_df=arms,
    max_loci_separation=10_000_000,   # ignore pixels farther than 10Mb from the diagonal
    n_lambda_bins=40, lambda_bin_fdr=0.1,   # FDR run independently per geometric lambda-bin (HiCCUPS default)
    clustering_radius=20_000,   # merge called pixels within 20kb into one loop
    nproc=4,
)
```

**The four backgrounds, and why lower-left is the clever one.** A pixel must beat ALL four local-background kernels, not one. **Donut** = is it a focal enrichment at all. **Horizontal** and **vertical** = is it actually a STRIPE pixel masquerading as a dot (these kernels exist to NOT call architectural stripes as loops). **Lower-left** = is it just a TAD/contact-domain CORNER -- a domain corner is enriched vs the donut but NOT vs its lower-left neighborhood, so requiring the pixel to also beat lower-left separates a genuine point loop from a generic domain corner. Skipping lower-left inflates calls with domain corners.

**Lambda-chunking is why HiCCUPS FDR is honest.** Contact counts span orders of magnitude with genomic distance, so a single genome-wide BH-FDR would be dominated by the high-count near-diagonal regime and over-call. `dots` bins pixels by their locally-adjusted expected into geometrically-spaced lambda-bins (`n_lambda_bins=40`) and runs BH-FDR independently within each (`lambda_bin_fdr=0.1`), so low-count and high-count regimes are each thresholded correctly.

## Template-Matching with chromosight

```bash
# detect: <contact_map> <prefix> are positional and come LAST
chromosight detect --pattern loops --threads 8 \
  --min-dist 20000 --max-dist 2000000 --pearson 0.4 \
  sample.cool::/resolutions/5000 sample_loops
# output: sample_loops.tsv -> chrom1,start1,end1,chrom2,start2,end2,bin1,bin2,score,pvalue,qvalue
```

The score is a Pearson correlation (-1..1) between a loop kernel and each windowed submatrix. The same engine finds borders and stripes by swapping `--pattern` (loops, loops_small, borders, hairpins, centromeres, stripes_left, stripes_right) -- but stripes are a separate feature class, NOT loops. `--pearson` is the correlation cutoff; raise it for fewer, higher-confidence calls.

## Scale-Space with Mustache

```bash
mustache -f sample.mcool -r 5000 -o loops.tsv -pt 0.1 -st 0.88 -norm weight -p 8
# output: BIN1_CHR BIN1_START BIN1_END BIN2_CHR BIN2_START BIN2_END FDR DETECTION_SCALE
```

`-pt` is the FDR/p-value threshold (default 0.1), `-st` the sparsity filter (default 0.88), `-norm weight` for a balanced `.cool` (`KR` for `.hic`). Mustache spans loop sizes natively via Difference-of-Gaussians across scales, which is why it adapts to kb-resolution Micro-C better than fixed-kernel HiCCUPS.

## Aggregate Peak Analysis (APA) -- Confirm, Don't Discover

**Goal:** Quantify whether a loop SET is enriched on average and visually QC the call set (a clean aggregate dot = mostly real; a smeared/absent center = contaminated).

**Approach:** Compute expected, pile up observed/expected snippets centered on each anchor pair, average across the stack, then report the APA score = center pixel divided by an off-diagonal corner-control block. Pass `expected_df` so snippets are O/E and comparable across genomic separations.

```python
import numpy as np
import cooltools

expected = cooltools.expected_cis(clr, view_df=arms, nproc=4)
stack = cooltools.pileup(clr, loops, view_df=arms, expected_df=expected, flank=100_000, nproc=4)   # bedpe two-anchor features
apa = np.nanmean(stack, axis=0)   # pileup returns (n_snippets, D, D); average over axis 0 -> 2D aggregate

center = apa.shape[0] // 2
corner = 3   # 3x3 corner-control block (Rao 2014 lower-left convention)
apa_score = apa[center, center] / np.nanmean(apa[-corner:, :corner])   # center vs lower-left corner; >1 = enriched
```

## Differential Loops -- Union Anchors, Not a Bin-Level Tool

**Goal:** Find loops whose strength changes between conditions.

**Approach:** There is NO DESeq-for-loops. Build a UNION anchor set across conditions, then score each loop's strength per condition at a FIXED coordinate set (`chromosight quantify`, which is purpose-built for this, or APA per condition), then test the strength delta.

```bash
# quantify scores a FIXED coordinate set; arg order: <bed2d> <contact_map> <prefix>
chromosight quantify --pattern loops union_anchors.bed2d condA.cool condA_q
chromosight quantify --pattern loops union_anchors.bed2d condB.cool condB_q
# compare the per-loop score columns; or diff_mustache.py -f1 A -f2 B -pt 0.05 -pt2 0.1 -r 5000 -o diff
```

`diffHic`, `multiHiCcompare`, and `dcHiC` operate on BINS or COMPARTMENTS, not focal loops -- do NOT use them as a loop-differential tool. Cross-reference hic-differential for the bin/compartment regime.

## Per-Method Failure Modes

### De-novo calling on a shallow map
**Trigger:** running `dots`/chromosight/Mustache on tens of millions of pairs or >=25kb bins. **Mechanism:** focal signal is below the noise floor without depth. **Symptom:** zero or a handful of scattered, irreproducible calls. **Fix:** STOP de-novo; run APA on a known anchor set (CTCF/cohesin ChIP or reference loops).

### APA reported without a corner control
**Trigger:** quoting the aggregate center-pixel value as the loop "strength." **Mechanism:** without an off-diagonal corner the number has no on-vs-off baseline. **Symptom:** a "high" APA that reflects distance-decay, not looping. **Fix:** APA score = center / corner-control block (Rao 2014 lower-left convention).

### APA presented as proof loops exist
**Trigger:** showing a pileup to claim "these N loops are real." **Mechanism:** APA surfaces mean enrichment across a SET; it cannot validate any single loop or discover new ones. **Symptom:** confident per-loop claims backed only by an aggregate. **Fix:** treat APA as set-level confirmation; for per-loop confidence use consensus + convergent-CTCF.

### Trusting a single caller's list
**Trigger:** reporting "Mustache found N loops" with no cross-check. **Mechanism:** callers have low pairwise overlap (Forcato 2017). **Symptom:** a list that barely overlaps a second tool or replicate. **Fix:** intersect >=2 callers and require convergent-CTCF / ChIA-PET / HiChIP support.

### Calling domain corners as loops
**Trigger:** a caller without a lower-left background (or a custom kernel set). **Mechanism:** a TAD corner beats the donut but is not a point loop. **Symptom:** "loops" sitting exactly at TAD corners with no anchor support. **Fix:** use `dots` (it beats all four backgrounds); cross-check anchors.

### Calling stripe pixels as dots
**Trigger:** detecting on a map with strong architectural stripes. **Mechanism:** a stripe pixel is enriched vs the donut but lies on a horizontal/vertical band. **Symptom:** "loops" smeared along a row/column. **Fix:** the horizontal/vertical kernels suppress these in `dots`; treat stripes as a separate class (chromosight stripes_*).

### 10kb-tuned kernels on 1kb Micro-C
**Trigger:** default HiCCUPS donut/peak-width on sub-5kb Micro-C. **Mechanism:** kernels are sized for 5-10kb Hi-C. **Symptom:** blurred or missed fine E-P loops. **Fix:** shrink the kernels for sub-5kb, or use Mustache/chromosight which adapt more gracefully.

### Raw (unbalanced) matrix into dots
**Trigger:** `dots` on a cooler with no `weight` column. **Mechanism:** dots requires balancing weights + expected. **Symptom:** error or meaningless output. **Fix:** `cooler balance` first; confirm `clr.matrix(balance=True)` is not all-NaN.

### HiCCUPS-style calling on HiChIP/PLAC-seq
**Trigger:** running `dots` on cohesin/H3K27ac HiChIP or PLAC-seq. **Mechanism:** protein-anchored data is coverage-biased; the Hi-C null is wrong. **Symptom:** distorted FDR, wrong loop counts. **Fix:** use FitHiChIP/MAPS/HiC-DC+ against a protein-anchored background.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| De-novo loop resolution 5-10kb | Rao 2014 depth scaling | ~4.9B contacts reached 1kb / ~10k loops; coarser bins blur anchors, shallow maps cannot resolve them |
| `max_loci_separation` 2-10Mb | loop size distribution | most loops are <2Mb; 10Mb is the cooltools default ceiling on diagonal distance |
| `n_lambda_bins=40`, `lambda_bin_fdr=0.1` | cooltools/HiCCUPS default | geometric lambda-binning + per-bin BH-FDR keeps FDR honest across the count dynamic range |
| `clustering_radius=20_000` | cooltools default | merges adjacent called pixels into one loop call |
| chromosight `--pearson` ~0.4 loops | chromosight default | template-correlation cutoff; raise for higher-confidence, fewer calls |
| Mustache `-pt 0.1`, `-st 0.88` | Mustache defaults | p/FDR threshold and sparsity filter |
| Consensus across >=2 callers | Forcato 2017 low overlap | single-caller lists are unreliable; require intersection or orthogonal support |
| APA score = center / corner block | Rao 2014 | the corner is the on-vs-off control; the bare center is uninterpretable |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `dots` returns nothing / scattered junk | map too shallow or resolution too coarse | check depth; below ~5-10kb resolution run APA on known anchors instead |
| `clr.matrix(balance=True)` all NaN | cooler not balanced | `cooler balance` / `cooler.balance_cooler` before calling loops |
| `expected`/`dots` shape or view error | different `view_df` for expected vs dots | reuse the same `view_df` (arms) for `expected_cis` and `dots` |
| Empty / wrong-resolution result on `.mcool` | bare `.mcool` passed | use `file.mcool::/resolutions/<bp>` URI |
| Empty result, no error | chrom naming mismatch (`chr1` vs `1`) across cooler/anchors/peaks | harmonize chromosome naming everywhere |
| `AttributeError` on a cooltools function | pre-0.7 vs 0.7+ signature change | `help(cooltools.dots)`; adapt to the installed signature |
| APA center looks high but loops are weak | no corner control / `expected_df` omitted | pass `expected_df` and divide center by a corner block |

## References

- Rao SSP, Huntley MH, Durand NC, et al. 2014. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. *Cell* 159(7):1665-1680.
- Open2C, Abdennur N, Abraham S, Fudenberg G, et al. 2024. Cooltools: enabling high-resolution Hi-C analysis in Python. *PLoS Comput Biol* 20(5):e1012067.
- Matthey-Doret C, Baudry L, Breuer A, et al. 2020. Computer vision for pattern detection in chromosome contact maps (chromosight). *Nat Commun* 11:5795.
- Roayaei Ardakany A, Gezer HT, Lonardi S, Ay F. 2020. Mustache: multi-scale detection of chromatin loops from Hi-C and Micro-C maps using scale-space representation. *Genome Biol* 21:256.
- Rowley MJ, Poulet A, Nichols MH, et al. 2020. Analysis of Hi-C data using SIP effectively identifies loops in organisms from C. elegans to mammals. *Genome Res* 30(3):447-458.
- Forcato M, Nicoletti C, Pal K, et al. 2017. Comparison of computational methods for Hi-C data analysis. *Nat Methods* 14:679-685.
- de Wit E, Vos ESM, Holwerda SJB, et al. 2015. CTCF binding polarity determines chromatin looping. *Mol Cell* 60(4):676-684.
- Sanborn AL, Rao SSP, Huang SC, et al. 2015. Chromatin extrusion explains key features of loop and domain formation. *PNAS* 112(47):E6456-E6465.
- Rao SSP, Huang SC, Glenn St Hilaire B, et al. 2017. Cohesin loss eliminates all loop domains. *Cell* 171(2):305-320.
- Haarhuis JHI, van der Weide RH, Blomen VA, et al. 2017. The cohesin release factor WAPL restricts chromatin loop extension. *Cell* 169(4):693-707.
- Schwarzer W, Abdennur N, Goloborodko A, et al. 2017. Two independent modes of chromatin organization revealed by cohesin removal. *Nature* 551:51-56.
- Krietenstein N, Abraham S, Venev SV, et al. 2020. Ultrastructural details of mammalian chromosome architecture (Micro-C). *Mol Cell* 78(3):554-565.
- Hsieh THS, Cattoglio C, Slobodyanyuk E, et al. 2020. Resolving the 3D landscape of transcription-linked mammalian chromatin folding (Micro-C). *Mol Cell* 78(3):539-553.
- Bhattacharyya S, Chandra V, Vijayanand P, Ay F. 2019. Identification of significant chromatin contacts from HiChIP data by FitHiChIP. *Nat Commun* 10:4221.

## Related Skills

- hic-data-io - Load and access the cooler files this skill calls loops on
- matrix-operations - Balancing and expected/O/E that dots and pileup depend on
- hic-visualization - Render called loops and APA pileups on the heatmap
- hic-differential - Bin/compartment-level differential (the regime loops are NOT in)
- tad-detection - TAD corners vs point loops; the lower-left background separates them
- chip-seq/peak-calling - CTCF/cohesin peaks to anchor and validate loops; HiChIP peak context
- chip-seq/peak-annotation - Annotate loop anchors with TF/CTCF peaks
- atac-seq/enhancer-gene-linking - E-P contacts complementing loop calls
- atac-seq/footprinting - TF footprints at loop anchors
- genome-intervals/overlap-significance - Permutation test for anchor/feature enrichment
