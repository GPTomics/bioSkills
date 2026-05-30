---
name: bio-copy-number-copy-ratio-segmentation
description: Normalize read-depth copy-ratio profiles and segment them into copy-number regions using circular binary segmentation (CBS, DNAcopy), hidden Markov models, HaarSeg, and fused-lasso methods. Covers GC-content, mappability, and replication-timing (wave-artifact) bias correction, panel-of-normals/PCA denoising, diploid-baseline centering, and algorithm selection by sequencing depth and event size. Use when choosing a segmentation algorithm, correcting depth bias, diagnosing oversegmentation or a mis-centered baseline, tuning CBS or HMM parameters, or understanding why a downstream CNV caller produced fragmented or shifted segments.
tool_type: mixed
primary_tool: DNAcopy
---

## Version Compatibility

Reference examples tested with: R 4.3+ with DNAcopy 1.76+, Python 3.10+ with numpy 1.26+, pandas 2.2+; QDNAseq 1.38+ (optional, GC/mappability normalization).

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('DNAcopy')` then `?segment` to confirm arguments
- Python: `pip show numpy pandas`

If code throws an error, introspect the installed package and adapt the example. CBS lives in Bioconductor `DNAcopy`; HMM segmentation is provided by caller-specific backends (CNVkit uses `pomegranate`; HaarSeg has its own R/Python packages).

# Copy-Ratio Segmentation

**"Turn noisy per-bin depth into clean copy-number segments"** -> Two stages, both error-prone. First, normalize the depth profile so the only remaining variation is copy number (not GC, mappability, or replication timing). Second, partition the normalized profile into segments of constant copy number. The segmentation algorithm choice has a *predictable* bias signature, and the diploid-baseline choice can invert every call.

- R: `DNAcopy::segment` (CBS, the reference implementation)
- Python: HMM via `pomegranate`; HaarSeg via `haarseg`
- The output feeds every CNV caller (cnvkit-analysis, gatk-cnv, allele-specific-copy-number)

## Stage 1: Why Depth Is Biased Before It Is Copy Number

Raw read depth confounds copy number with three systematic biases:

| Bias | Cause | Correction |
|------|-------|------------|
| GC content | PCR efficiency and probe hybridization vary with GC | Loess fit of depth vs GC (QDNAseq), or matched normal |
| Mappability | Multi-mapping reads under-counted in repetitive regions | Mappability track filter/weight; exclude low-mappability bins |
| Replication timing | Late-replicating DNA is under-represented — the "wave artifact" | Matched normal or PoN; GC correction alone does NOT remove it |
| Capture efficiency | Per-probe hybridization varies 10-100x (hybrid capture) | Panel of normals — the dominant bias for exomes/panels |

The key postdoc-level point: **GC correction alone is insufficient.** The wave artifact in cancer WGS is driven by replication timing, a biological signal GC normalization cannot flatten. Only a matched normal or a panel of normals removes it. This is why a CNVkit flat reference (GC-only) produces systematic false focal calls and why GATK tangent normalization exists.

## Stage 2: Segmentation Algorithm Taxonomy

| Algorithm | Model | Strength | Fails when |
|-----------|-------|----------|------------|
| CBS (circular binary segmentation) | Recursive t-statistic breakpoint test | High precision; excellent on small focal segments | Low depth (~3x): recall drops to ~42% (worse under over-dispersed counts); ~2 orders slower; fragments across assembly gaps |
| HMM | Hidden CN states, emission + transition | Depth-robust; high recall at low coverage | Less precise on small focal segments (~5 kb: ~76% precision vs CBS ~96%, Poisson model); EM finds only local optima |
| HaarSeg | Wavelet (Haar) multiscale edge detection | Very fast; good for shallow WGS | Less precise breakpoints than CBS; threshold-sensitive |
| Fused lasso (flasso) | L1-penalized piecewise-constant fit | Smooth; tunable sparsity | Penalty hard to set; can over-smooth focal events |
| ASPCF | Allele-specific piecewise-constant fit | Joint logR+BAF segmentation (ASCAT) | Needs BAF; see allele-specific-copy-number |

**Quantitative benchmark (Zhang et al 2024, Brief Bioinform):** the cited precision/recall numbers (CBS ~42% recall at 3x; HMM ~81% recall at 3x; CBS ~96% precision vs HMM ~76% on 5 kb focal segments under a Poisson model) summarise that paper's reported direction of the trade-off. Verify the exact figures against the published tables before quoting them in print; the qualitative trade-off (depth-vs-event-size, CBS-vs-HMM) is robust across recent benchmarks but the precise percentages depend on the simulation model (Poisson vs over-dispersed negative-binomial). There is no universally correct choice.

## Decision Tree

| Scenario | Algorithm | Rationale |
|----------|-----------|-----------|
| Panel / exome, adequate depth, focal events matter | CBS | Precise on small segments |
| Shallow WGS (< ~5x), broad events | HMM or HaarSeg | CBS recall degrades at low depth |
| Heterogeneous / impure tumor | HMM (e.g. CNVkit `hmm-tumor`) | Broader state transitions absorb noise |
| Germline, near-diploid | HMM with diploid-tight priors | Priors stabilize calls near CN=2 |
| Allele-specific (need BAF) | ASPCF / FACETS joint CBS | See allele-specific-copy-number |
| Very large WGS, speed-critical | HaarSeg | Near-linear; CBS is ~100x slower |

## Bias Correction — GC Loess Normalization

**Goal:** Remove GC-content bias from a per-bin depth profile.

**Approach:** Fit a loess curve of depth versus GC content, divide each bin by its fitted value, log2-transform. This corrects GC but not replication timing — use a normal for that.

```python
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

def gc_correct(bins):
    '''GC-correct a per-bin depth profile. bins: columns chrom, start, depth, gc.
    Returns log2 copy ratio relative to the GC-corrected genome median.'''
    df = bins[(bins['depth'] > 0) & bins['gc'].between(0.3, 0.7)].copy()
    fitted = lowess(df['depth'], df['gc'], frac=0.3, return_sorted=False)
    df['corrected'] = df['depth'] / fitted
    df['log2'] = np.log2(df['corrected'] / df['corrected'].median())
    return df
```

For exomes and panels, a panel of normals (per-bin median of normals, or PCA denoising) is preferred over GC-only correction because it also removes capture and replication-timing bias.

## Segmentation — CBS with DNAcopy

**Goal:** Segment a normalized log2 profile into copy-number regions.

**Approach:** Build a CNA object, smooth single-bin outliers, run CBS, then merge adjacent segments whose means differ by less than a noise-scaled threshold (`sdundo`).

```r
library(DNAcopy)

# bins: data frame with chrom, maploc (bin midpoint), log2
cna <- CNA(genomdat = bins$log2, chrom = bins$chrom, maploc = bins$maploc,
           data.type = 'logratio', sampleid = 'tumor')
cna <- smooth.CNA(cna)                          # damp single-bin outliers

# alpha = breakpoint significance; undo.splits='sdundo' merges segments whose means
# are within undo.SD noise standard deviations -- the main guard against oversegmentation.
seg <- segment(cna, alpha = 0.01, undo.splits = 'sdundo', undo.SD = 2,
               verbose = 1)
write.table(seg$output, 'tumor.segments.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)
```

## Failure Modes

### Oversegmentation / hyperfragmentation

**Trigger:** CBS `alpha` too liberal, `undo.SD` too small, or a noisy (high-MAD) profile; FACETS `cval` too low.

**Mechanism:** The breakpoint test fires on noise; the profile shatters into many tiny segments that do not correspond to real copy-number changes.

**Symptom:** Hundreds of short segments; segment count scales with noise, not biology; downstream integer CN incoherent with per-bin medians.

**Fix:** Raise `alpha` toward 0.01 or stricter, increase `undo.SD` (e.g. 2-3), or denoise the input first (better PoN, drop low-coverage bins). Three signatures (Steele 2022) had to be discarded as oversegmentation artifacts — fragmentation propagates into every downstream analysis, including copy-number signatures.

### The diploid-baseline centering trap

**Trigger:** Centering the log2 profile on its median or mode in a hyper-aneuploid or whole-genome-doubled genome.

**Mechanism:** Centering assumes the commonest log2 value is diploid. In a WGD genome the commonest state is tetraploid; centering on it shifts the whole profile so true diploid regions read as deletions and amplifications read as neutral.

**Symptom:** Genome-wide gain or loss inconsistent with biology; segmentation is fine but every call has the wrong sign.

**Fix:** Do not depth-center aneuploid genomes. Anchor the diploid baseline with BAF/SNV data via an allele-specific caller, which estimates absolute ploidy. GISTIC and most callers require a correctly centered seg file as input.

### CBS recall degrades at low depth

**Trigger:** CBS on shallow data (< ~5x WGS, or low-coverage bins).

**Mechanism:** The two-sample t-statistic loses power when per-bin variance swamps the mean difference; CBS misses real breakpoints (recall ~42% at 3x under a Poisson model, worse under over-dispersed counts), while the segments it does call stay fairly precise.

**Symptom:** Real events absent from the segmentation; recall poor on a genome with known CNVs.

**Fix:** Use HMM (depth-robust, ~81% recall at 3x) or HaarSeg for shallow data; or increase bin size to raise per-bin counts before segmenting.

### CBS fragments across assembly gaps

**Trigger:** CBS run over a profile with centromere/telomere gaps not handled as chromosome breaks.

**Mechanism:** CBS treats gapped data as independent subsets; spurious breakpoints appear at gap edges.

**Symptom:** Segment boundaries clustered at centromeres; tiny artifactual segments flanking gaps.

**Fix:** Segment per chromosome arm, or supply gap-aware chromosome coordinates so CBS does not bridge gaps.

### HMM EM converges to a local optimum

**Trigger:** HMM with poor initial parameters or too few iterations.

**Mechanism:** Baum-Welch EM is not globally optimal; emission/transition parameters can settle in a local optimum, mis-assigning states.

**Symptom:** Reruns give different state assignments; CN states inconsistent with the visible profile.

**Fix:** Use informative priors (diploid-centered for germline, broader for tumor), run multiple initializations, and sanity-check state means against the per-bin distribution.

## Reconciliation: When Segmentations Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| CBS shatters where HMM gives clean broad segments | Low depth — CBS over-fits noise | Trust HMM; CBS needs more depth |
| HMM misses a focal event CBS finds | HMM window resolution too coarse | Trust CBS for focal; HMM blurs small events |
| Both agree on arms, differ on focal boundaries | Different breakpoint resolution | Arm calls are robust; treat focal boundaries as approximate |
| Segmentation differs run-to-run | HMM local optima, or unfixed random seed | Fix seeds; use multiple HMM initializations |

**Operational rule:** Match the algorithm to depth and event size — CBS for adequate-depth focal work, HMM/HaarSeg for shallow or broad. Confirm the diploid baseline against an allele-specific ploidy estimate before any sign-dependent interpretation. Report arm-level segments with confidence; treat focal boundaries as algorithm-dependent.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| CBS `alpha` | 0.01 | DNAcopy default; breakpoint significance |
| CBS `undo.SD` | 2-3 | Merges segments within N noise SD; guards oversegmentation |
| Depth where CBS recall degrades | < ~5x WGS | Zhang 2024; CBS recall falls sharply (HMM is depth-robust) |
| GC range kept for loess | 0.3-0.7 | Extreme-GC bins are unreliable; standard restriction |
| Bin size, shallow WGS CNV | ~500 kb - 1 Mb | Larger bins raise per-bin counts for stable segmentation |
| Sample MAD usable | < 0.5 | Above this, segmentation chases noise regardless of algorithm |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Hundreds of tiny segments | Oversegmentation (liberal alpha / low cval / noisy input) | Tighten alpha, raise undo.SD, denoise input |
| Whole genome wrong-signed | Baseline centered on a non-diploid mode | Anchor ploidy with BAF; do not depth-center aneuploid genomes |
| Real events missed at low depth | CBS recall degrades | Use HMM/HaarSeg or larger bins |
| Breakpoints clustered at centromeres | CBS bridging assembly gaps | Segment per arm; supply gap-aware coordinates |
| Segmentation not reproducible | HMM local optima / unfixed seed | Fix seeds; multiple initializations |
| Wave artifact remains after GC correction | Replication-timing bias, not GC | Use a matched normal or PoN |

## References

- Olshen AB et al 2004. Circular binary segmentation for the analysis of array-based DNA copy number data. Biostatistics 5:557
- Venkatraman ES, Olshen AB 2007. A faster circular binary segmentation algorithm. Bioinformatics 23:657
- Zhang Y, Liu W, Duan J 2024. On the core segmentation algorithms of copy number variation detection tools. Brief Bioinform 25:bbae022
- Ben-Yaacov E, Eldar YC 2008. A fast and flexible method for the segmentation of aCGH data (HaarSeg). Bioinformatics 24:i139
- Scheinin I et al 2014. DNA copy number analysis of fresh and FFPE specimens by shallow WGS (QDNAseq). Genome Res 24:2022

## Related Skills

- copy-number/cnvkit-analysis - Read-depth caller exposing CBS/HMM/HaarSeg choices
- copy-number/gatk-cnv - Tangent normalization and ModelSegments segmentation
- copy-number/allele-specific-copy-number - ASPCF joint logR+BAF segmentation
- copy-number/recurrent-cnv - Copy-number signatures sensitive to segmentation quality
- copy-number/cnv-visualization - Visual diagnosis of oversegmentation and baseline shift
- genome-intervals/coverage-analysis - Per-bin depth computation upstream of segmentation
