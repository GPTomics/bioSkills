---
name: bio-copy-number-subclonal-copy-number
description: Resolve subclonal copy number, whole-genome doubling, and copy-number tumor evolution from bulk sequencing with Battenberg, TITAN, and MEDICC2. Covers clonal versus subclonal copy-number states, haplotype phasing for subclonal resolution, cancer cell fraction, whole-genome-doubling detection and timing relative to mutations, mirrored subclonal allelic imbalance, and copy-number phylogenies. Use when a tumor is heterogeneous and bulk data shows non-integer copy number, when calling subclonal CNAs, detecting or timing whole-genome doubling, reconstructing copy-number evolution, or deciding between Battenberg and TITAN.
tool_type: mixed
primary_tool: battenberg
---

## Version Compatibility

Reference examples tested with: R 4.3+ with Battenberg 2.2.10+ and TitanCNA 1.40+, MEDICC2 1.0+, Python 3.10+; impute2/Beagle phasing reference panels.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('Battenberg')` / `'TitanCNA')` then `?function`
- CLI: `medicc2 --help`
- Battenberg is GitHub-only (`Wedge-lab/battenberg`) and needs a 1000 Genomes impute/phasing reference and allele-counter; confirm reference data is installed

Battenberg and TITAN both consume allele-specific data (logR + BAF at heterozygous SNPs); they cannot run on relative copy ratio alone.

# Subclonal Copy Number and Tumor Evolution

**"This copy number is non-integer — is it noise, or are there subclones"** -> A tumor is a mixture of cell populations. When a copy-number change is present in only some cancer cells, bulk sequencing averages it into a *non-integer* state. A long non-integer segment is not noise — it is a subclonal copy-number alteration, and resolving it reveals the tumor's clonal architecture.

- R: `Battenberg` (phased clonal + subclonal CN), `TitanCNA` (HMM mixture of cell populations)
- CLI: `medicc2` (whole-genome-doubling-aware copy-number phylogenies)
- Input: allele-specific data — see allele-specific-copy-number for the clonal layer

## Clonal vs Subclonal — What the Tools Output

| Concept | Meaning |
|---------|---------|
| Clonal CNA | Present in all cancer cells; one copy-number state per segment |
| Subclonal CNA | Present in a fraction of cancer cells; the segment needs two states plus a fraction |
| Cancer cell fraction (CCF) | Fraction of cancer cells carrying the event |
| Mirrored subclonal allelic imbalance | Different subclones lose opposite haplotypes of the same region |

Battenberg fits a clonal allele-specific profile (ASCAT internally), then where a segment fits poorly as a single integer state, it models it as a mixture of two states with a subclonal fraction. TITAN uses an HMM whose states span multiple clonal clusters, jointly estimating per-cluster cellular prevalence. Both need haplotype phasing — subclonal allelic imbalance is only resolvable when SNPs are phased.

## Tool Selection

| Tool | Model | Best for | Fails when |
|------|-------|----------|------------|
| Battenberg | Phased clonal fit + per-segment subclonal mixture | WGS, subclonal CN to ~3% of cells, clonal-evolution studies | Low depth/purity; heavy compute; needs phasing reference |
| TITAN | HMM mixture across clonal clusters | WGS/WES, joint CN+LOH+subclonal prevalence, few clusters | Many subclones; cluster number must be chosen and swept |
| MEDICC2 | WGD-aware minimum-event copy-number phylogeny | Multi-sample / multi-region evolution | Single sample (no tree to build) |
| ASCAT/FACETS | Clonal allele-specific only | When subclonal resolution is not needed | Treats subclonal segments as noisy clonal — see allele-specific-copy-number |

## Whole-Genome Doubling — Detection and Timing

Whole-genome doubling (WGD) is a discrete, common (~30% of advanced cancers) evolutionary event, and it must be called explicitly because it changes how every copy number is read.

- **Detection:** A tumor has undergone WGD if more than ~50% of the autosomal genome has a major (more frequent) allele copy number >= 2. WGD tumors have median ploidy ~3.3 versus ~2.1 for non-WGD.
- **Relative vs absolute:** Depth gives *relative* copy number; WGD calling needs *absolute* allele-specific copy number (BAF anchors ploidy). A depth-only profile cannot distinguish a WGD genome from a non-WGD genome — this is the identifiability problem of allele-specific-copy-number in another guise.
- **Timing:** WGD is timeable relative to point mutations. Mutations that arose before WGD are carried at multiple copies (mutation copy number ~2); mutations after WGD sit at one copy. This dates WGD within the tumor's mutational history.

## Calling Subclonal CN with Battenberg

**Goal:** Fit clonal and subclonal allele-specific copy number genome-wide.

**Approach:** Generate phased allele counts against a 1000 Genomes reference, run the Battenberg pipeline; segments that fit poorly as one integer state are split into a two-state subclonal mixture with a cellular fraction.

```r
library(Battenberg)

# Battenberg orchestrates allele counting, phasing, ASCAT clonal fit, and the
# subclonal mixture step. Reference data (1000G impute panel) must be installed.
battenberg(
    samplename          = 'tumour_id',
    normalname          = 'normal_id',
    sample_data_file    = 'tumour.bam',
    normal_data_file    = 'normal.bam',
    ismale              = TRUE,
    imputeinfofile      = 'impute_info.txt',
    g1000prefix         = '1000G_loci/1000genomesloci2012_chr',     # SNP loci data
    g1000allelesprefix  = '1000G_alleles/1000genomesAlleles2012_chr', # SNP alleles (WGS)
    problemloci         = 'probloci.txt',
    gccorrectprefix     = 'GC_correction_hg38_chr',
    repliccorrectprefix = 'RT_correction_hg38_chr',
    genomebuild         = 'hg38',                                   # default is hg19
    nthreads            = 8)
# Output *_subclones.txt: per segment, nMaj1/nMin1 (state 1) + frac1, and nMaj2/nMin2 +
# frac2 when the segment is subclonal (two states).
```

## Calling Subclonal CN with TITAN

**Goal:** Jointly infer copy number, LOH, and the cellular prevalence of clonal clusters.

**Approach:** TITAN needs both allele counts (het SNPs) and corrected read depth. Load the allele counts; correct tumour/normal read depth for GC and mappability bias; overlay the resulting logR onto the het positions and log-transform; filter; then run the EM and sweep the cluster number — model selection picks the best.

```r
library(TitanCNA)

# Allele counts at het SNPs.
data <- loadAlleleCounts('tumour.allelicCounts.tsv', genomeStyle = 'UCSC')

# Read-depth correction is mandatory: correctReadDepth needs tumour + normal coverage
# WIGs and GC + mappability WIGs. genomeStyle MUST match loadAlleleCounts above
# (default 'NCBI' vs 'UCSC') or getPositionOverlap matches no chromosomes and logR is NA.
cnData <- correctReadDepth('tumour.wig', 'normal.wig', 'gc.wig', 'map.wig',
                           genomeStyle = 'UCSC')
data$logR <- log(2 ^ getPositionOverlap(data$chr, data$posn, cnData))
data <- filterData(data, 1:24, minDepth = 10, maxDepth = 200, map = NULL)

params <- loadDefaultParameters(copyNumber = 8, numberClonalClusters = 2,
                                symmetric = TRUE, data = data)
conv <- runEMclonalCN(data, params, maxiter = 20, txnExpLen = 1e15)
results <- viterbiClonalCN(data, conv)
# Sweep numberClonalClusters (1..5) and compare model fit; the S_Dbw validity index
# or the model log-likelihood selects the cluster number.
```

## Failure Modes

### Subclonal call from insufficient depth or purity

**Trigger:** Calling subclonal CN on shallow WGS or a low-purity tumor.

**Mechanism:** A subclonal segment's signal is the clonal deviation scaled by the subclone's cell fraction — already small, and below the noise floor at low depth/purity.

**Symptom:** Many "subclonal" segments with implausibly low fractions; calls not reproducible across reruns or regions.

**Fix:** Battenberg's ~3%-of-cells sensitivity assumes adequate WGS depth and purity. For low-depth or low-purity samples, treat only clonal CN as reliable and report subclonal calls as exploratory.

### Mirrored subclonal allelic imbalance misread

**Trigger:** A region where different subclones lost opposite haplotypes.

**Mechanism:** Bulk BAF averages the two opposite losses toward 0.5, so the region can look balanced (clonal, no LOH) when it is in fact subclonally rearranged on both haplotypes.

**Symptom:** A segment called clonal-balanced that conflicts with multi-region or single-cell data; BAF near 0.5 with an odd logR.

**Fix:** Phasing (Battenberg) is required to detect mirrored subclonal allelic imbalance. Multi-region or single-cell data resolves it definitively; a single bulk sample can miss it.

### WGD not called — every copy number off by a factor

**Trigger:** Interpreting copy number without first establishing WGD status.

**Mechanism:** The likelihood surface has near-equal modes at ploidy P and 2P; missing a WGD halves all copy numbers and mis-times every mutation.

**Symptom:** Copy numbers and mutation copy numbers inconsistent; "subclonal" gains that are actually clonal post-WGD states.

**Fix:** Call WGD explicitly (>50% of autosomes at major CN >= 2) from absolute allele-specific copy number. Cross-check ploidy against the odd/even CN fraction and clonal-SNV multiplicity before any subclonal interpretation.

### Over-interpreting one subclonal segment as a subclone

**Trigger:** Declaring a distinct tumor subclone from a single subclonal copy-number segment.

**Mechanism:** A single segment at an intermediate fraction can arise from segmentation error, a mis-fit clonal state, or genuine subclonality — one segment cannot distinguish these.

**Symptom:** A "subclone" supported by exactly one segment; clonal architecture claims that do not replicate.

**Fix:** Require multiple concordant subclonal segments at a consistent cell fraction, ideally corroborated by SNV-based subclonal reconstruction (cancer cell fraction clustering) and multi-region sampling.

### Single-region sampling misses spatial subclones

**Trigger:** Inferring clonal architecture from one biopsy of a spatially heterogeneous tumor.

**Mechanism:** A subclone confined to an unsampled region is invisible; a single region cannot capture branching evolution.

**Symptom:** Apparently simple clonal architecture contradicted by a second biopsy.

**Fix:** For evolution and architecture claims, use multi-region sampling and a phylogeny method (MEDICC2 for copy-number trees). Single-region subclonal calls describe that region only.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Battenberg subclonal vs TITAN clonal | Different mixture models; cluster number | Sweep TITAN clusters; compare cell fractions |
| WGD called by one tool, not another | Integer-multiple ploidy ambiguity | Check odd/even CN fraction and SNV multiplicity |
| Many low-fraction subclonal segments | Depth/purity too low | Trust only clonal CN; flag subclonal as exploratory |
| Subclonal CN vs SNV-based CCF disagree | CN and SNV subclones need not coincide | Integrate both; they answer different questions |

**Operational rule:** Report subclonal copy number as confident only when (1) depth and purity support it, (2) WGD status is established from absolute allele-specific CN, (3) multiple concordant segments support a subclone at a consistent fraction, and (4) for evolution claims, multi-region data and a copy-number phylogeny are used. A single subclonal segment is a hypothesis, not a subclone.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Battenberg subclonal sensitivity | ~3% of cells | Nik-Zainal 2012; requires adequate WGS depth/purity |
| WGD definition | > 50% of autosomes at major CN >= 2 | Bielski 2018; the operational WGD call |
| WGD median ploidy | ~3.3 (WGD) vs ~2.1 (non-WGD) | Bielski 2018 pan-cancer |
| Pre-WGD mutation copy number | >= ~1.75 | Pre-doubling mutations carried at multiple copies |
| TITAN clonal clusters | sweep 1-5, select by fit | Few clusters resolvable from one bulk sample |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Battenberg install/run fails | GitHub-only; missing 1000G reference | Install from GitHub; set up the impute reference |
| Non-integer segments treated as noise | Subclonal CNA not modeled | Use Battenberg/TITAN, not a clonal-only caller |
| All copy numbers half/double expected | WGD not called | Establish WGD from absolute CN; check SNV multiplicity |
| Subclones not reproducible | Low depth/purity, single segment | Require depth, concordant segments, multi-region |
| Balanced region conflicts with other data | Mirrored subclonal allelic imbalance | Use phased (Battenberg) or single-cell data |
| TITAN cluster number arbitrary | Cluster count not swept | Sweep 1-5; select by model fit |

## References

- Nik-Zainal S et al 2012. The life history of 21 breast cancers (Battenberg). Cell 149:994
- Ha G et al 2014. TITAN: inference of copy number architectures in clonal cell populations from tumor whole-genome sequence data. Genome Res 24:1881
- Bielski CM et al 2018. Genome doubling shapes the evolution and prognosis of advanced cancers. Nat Genet 50:1189
- Dewhurst SM et al 2014. Tolerance of whole-genome doubling propagates chromosomal instability. Cancer Discov 4:175
- Kaufmann TL et al 2022. MEDICC2: whole-genome doubling aware copy-number phylogenies for cancer evolution. Genome Biol 23:241

## Related Skills

- copy-number/allele-specific-copy-number - Clonal allele-specific CN, purity, ploidy
- copy-number/copy-ratio-segmentation - Segmentation feeding subclonal callers
- copy-number/hrd-scoring - Whole-genome-doubling correction for LST
- copy-number/recurrent-cnv - Copy-number signatures including WGD and chromothripsis
- copy-number/cnv-visualization - Visualizing subclonal segments and BAF
- variant-calling/vcf-basics - SNV calls for cancer cell fraction and WGD timing
