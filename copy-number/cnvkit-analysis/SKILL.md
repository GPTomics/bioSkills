---
name: bio-copy-number-cnvkit-analysis
description: Detect somatic and germline copy number variants from targeted, exome, and whole-genome sequencing with CNVkit, a read-depth caller that combines on-target and off-target (antitarget) coverage. Covers panel-of-normals construction, flat-reference tumor-only calling, hybrid/amplicon/WGS modes, CBS vs HMM segmentation selection, purity-aware integer calling, and reconciliation against GATK and allele-specific callers. Use when calling CNVs from hybrid-capture panels or exomes, deciding whether CNVkit (depth-only) is the right tool versus an allele-specific caller, building a panel of normals, diagnosing flat-reference false positives, or interpreting log2 ratios into copy-number states.
tool_type: cli
primary_tool: cnvkit
---

## Version Compatibility

Reference examples tested with: CNVkit 0.9.10+, samtools 1.19+, bedtools 2.31+, Python 3.10+, R 4.3+ with DNAcopy 1.76+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `cnvkit.py version` then `cnvkit.py batch --help` to confirm flags
- Python: `pip show cnvkit` then `python -c "import cnvlib; help(cnvlib.read)"`
- R: `packageVersion('DNAcopy')` (CBS backend)

If a command throws an unrecognized-argument or AttributeError, introspect the installed version and adapt the example rather than retrying. CNVkit segmentation methods (`hmm`, `hmm-tumor`, `hmm-germline`) depend on `pomegranate`; CBS depends on Bioconductor `DNAcopy`.

# CNVkit Copy Number Analysis

**"Detect copy number variants from my exome / panel data"** -> Run a read-depth pipeline: normalize on-target and off-target coverage against a reference, segment the log2-ratio profile, and call gains/losses. CNVkit is a *depth-only* caller — it estimates **relative** copy number and cannot, on its own, resolve tumor purity, ploidy, or allele-specific state. Choosing CNVkit is a decision that the experiment does not require allelic resolution.

- CLI: `cnvkit.py batch tumor.bam --normal normal.bam --targets panel.bed --fasta ref.fa`
- Python API: `cnvlib.read('sample.cnr')` for downstream filtering

## Where CNVkit Sits — Caller Taxonomy

| Caller | Signal used | Output | Purity/ploidy aware | Fails when |
|--------|-------------|--------|---------------------|------------|
| CNVkit | Depth (on + off-target) | Relative log2, threshold-called CN | No (manual `--purity`) | Sample purity < ~40%; hyper-aneuploid genome breaks median centering; balanced events invisible |
| GATK gCNV / somatic CNV | Depth (PCA/tangent denoised) | Copy-ratio segments, +/-/0 call | No (somatic); ploidy prior (germline) | Recurrent CNV in the PoN normalized away; ModelSegments gives no integer ASCN |
| ASCAT / Sequenza / FACETS | Depth + B-allele frequency | Integer allele-specific CN, purity, ploidy | Yes (jointly fit) | Near-diploid genome cannot anchor purity; low het-SNP density |
| ExomeDepth / GATK gCNV (cohort) | Depth across a cohort | Germline CN genotype | Germline ploidy only | < ~30-100 technically matched samples; common CNV |

CNVkit's niche: a fast, single-sample depth caller for hybrid-capture panels and exomes where antitarget reads recover genome-wide resolution. Its limit: it answers "is this region gained or lost relative to baseline" — not "how many absolute copies, on which haplotype, in what fraction of cells." For tumor integer CN, purity, LOH, or whole-genome doubling, escalate to allele-specific-copy-number.

## Decision Tree by Scenario

| Scenario | Recommended CNVkit configuration | Why |
|----------|----------------------------------|-----|
| Hybrid-capture panel or exome, tumor-normal | `batch` hybrid mode, matched normal as reference | Antitargets recover off-target resolution; matched normal cancels capture bias |
| Exome cohort, pooled normals available | Build pooled PoN reference, then `batch --reference` | Pooled reference averages out per-normal noise; 5-20+ normals |
| Amplicon / multiplex-PCR panel | `batch --method amplicon` | No usable off-target reads; antitarget bins are pure noise — must be dropped |
| Whole-genome sequencing | `batch --method wgs` | Genome-wide fixed bins; no target/antitarget split |
| Tumor-only, no normal of any kind | `batch` with flat reference (omit `--normal`) | Last resort; expect GC/capture-bias false positives — see failure mode below |
| FFPE / low-input / impure tumor | Add `--drop-low-coverage`; segment with `hmm-tumor` | FFPE dropout produces zero-coverage bins that CBS reads as deletions |
| Need absolute CN, LOH, purity | Do not use CNVkit alone | Escalate to allele-specific-copy-number (ASCAT/Sequenza/FACETS/PureCN) |

## Core Pipeline — Tumor-Normal Pair

The `batch` command wraps target/antitarget generation, coverage, reference building, fix, and segment:

```bash
cnvkit.py batch tumor.bam \
    --normal normal.bam \
    --targets panel.bed \
    --annotate refFlat.txt \
    --fasta reference.fa \
    --access access-excludes.bed \
    --output-reference reference.cnn \
    --output-dir results/ \
    --drop-low-coverage \
    --diagram --scatter
```

`--access` restricts antitarget bins to mappable, non-gap genome (generate once with `cnvkit.py access reference.fa -o access.bed`). `--drop-low-coverage` is effectively mandatory for tumor, FFPE, or any sample with coverage dropout.

## Panel of Normals — The Reference Determines Call Quality

A reference built from pooled normals is the single largest quality lever. Process is: build the reference from normals once, then run every tumor against it.

```bash
# Build pooled reference from process-matched normals (same capture kit, same lab)
cnvkit.py batch --normal normal*.bam \
    --targets panel.bed --annotate refFlat.txt --fasta reference.fa \
    --access access.bed \
    --output-reference pooled_reference.cnn

# Run each tumor against the pre-built reference
cnvkit.py batch tumor*.bam --reference pooled_reference.cnn \
    --output-dir results/ --drop-low-coverage --scatter --diagram
```

## Step-by-Step Pipeline (Fine-Grained Control)

```bash
cnvkit.py target panel.bed --annotate refFlat.txt --split -o targets.bed
cnvkit.py antitarget panel.bed --access access.bed -o antitargets.bed
cnvkit.py coverage tumor.bam targets.bed -o tumor.targetcoverage.cnn
cnvkit.py coverage tumor.bam antitargets.bed -o tumor.antitargetcoverage.cnn
cnvkit.py reference normal*.{target,antitarget}coverage.cnn --fasta reference.fa -o reference.cnn
cnvkit.py fix tumor.targetcoverage.cnn tumor.antitargetcoverage.cnn reference.cnn -o tumor.cnr
cnvkit.py segment tumor.cnr -o tumor.cns --drop-low-coverage
cnvkit.py call tumor.cns -o tumor.call.cns
```

## Segmentation Method Selection

CNVkit's `segment` step is where the bias-variance trade-off is set. The default CBS is not always correct — see copy-ratio-segmentation for the full algorithm comparison.

```bash
cnvkit.py segment tumor.cnr -m cbs -o tumor.cns          # default; precise on focal events
cnvkit.py segment tumor.cnr -m hmm-tumor -o tumor.cns    # heterogeneous tumor, broad states
cnvkit.py segment tumor.cnr -m hmm-germline -o tumor.cns # germline, priors near diploid
cnvkit.py segment tumor.cnr -m haar -o tumor.cns         # fast, low-depth WGS
```

Rule of thumb: CBS for panels/exomes with adequate depth (precise on small segments); `hmm-tumor` for impure or heterogeneous tumors where CBS over-fragments; `haar` for shallow WGS where CBS recall degrades.

## Purity-Aware Integer Calling

`call` converts segmented log2 ratios to copy-number states. The `clonal` method rescales by tumor purity before rounding to integers — without it, an impure tumor's true CN=4 amplification rounds to CN=3 or CN=2.

```bash
# Threshold method (default): fixed log2 cutpoints, no purity correction
cnvkit.py call tumor.cns -o tumor.call.cns

# Clonal method: rescale by purity, then round to integer CN
cnvkit.py call tumor.cns -m clonal --purity 0.65 --ploidy 2 -o tumor.call.cns

# Overlay B-allele frequency from a SNV VCF (for LOH visualization, NOT joint ASCN)
cnvkit.py call tumor.cns -m clonal --purity 0.65 --vcf tumor.vcf.gz -o tumor.call.cns
```

CNVkit can read BAF from a VCF and report a `baf` column, but it segments log2 and BAF *separately* and does not jointly fit purity from them. For a true joint allele-specific model (ASPCF, FACETS joint segmentation), use allele-specific-copy-number.

## Failure Modes

### Flat reference (tumor-only) — systematic false focal calls

**Trigger:** No `--normal` and no pooled PoN; CNVkit builds a flat reference (uniform log2 0) from the FASTA.

**Mechanism:** A flat reference corrects only GC and (optionally) RepeatMasker content via the FASTA. It cannot correct capture efficiency, which varies 10-100x across probes and is the dominant bias in hybrid capture. Per-probe capture bias is then misread as copy number.

**Symptom:** Recurrent "CNVs" at the same loci across unrelated tumor-only samples; spiky `.cnr` profiles; high MAD; calls concentrated at probe boundaries.

**Fix:** Never rely on a flat reference for clinical or focal calls. Build a pooled PoN from >= 5 process-matched normals. If truly no normal exists, treat tumor-only CNVkit output as hypothesis-generating only and escalate to PureCN (allele-specific-copy-number), which models a normal database explicitly.

### Antitarget bins on amplicon panels — pure noise

**Trigger:** Running default hybrid mode on an amplicon (multiplex-PCR) panel.

**Mechanism:** Amplicon panels produce essentially no off-target reads. Antitarget bins then contain a handful of stray reads, giving wildly variable log2 that the segmenter chases.

**Symptom:** Huge antitarget bin spread; nonsensical genome-wide segments between the targeted genes.

**Fix:** Use `--method amplicon`, which drops antitargets entirely and calls only from on-target bins. Accept that resolution is limited to the targeted genes.

### Low tumor purity — the death zone below ~40%

**Trigger:** Tumor cellularity below ~40% (common in breast, lung adenocarcinoma, melanoma, low-cellularity biopsies).

**Mechanism:** Each somatic CN change is diluted by 2-copy normal DNA. A true single-copy loss at 30% purity produces log2 ~ -0.23 — inside the diploid threshold band.

**Symptom:** Genome looks near-flat; few or no calls; known driver amplifications (e.g. ERBB2, MYC) missed.

**Fix:** CNVkit cannot rescue this. Confirm purity with an allele-specific caller (BAF gives an orthogonal purity estimate). Below ~20% purity, no depth-based caller is reliable — report as indeterminate.

### Hyper-aneuploid / whole-genome-doubled genome — baseline miscalled

**Trigger:** Tumor with >50% of the genome altered, or whole-genome doubling.

**Mechanism:** `call --center median` (or mode) assumes the commonest log2 state is diploid. In a WGD genome the commonest state is tetraploid; centering on it shifts the whole profile and inverts gain/loss calls.

**Symptom:** Genome-wide pattern of calls inconsistent with known biology; "deletions" everywhere or "gains" everywhere.

**Fix:** Do not trust depth-only centering on aneuploid tumors. Anchor the diploid baseline with BAF/SNV data via an allele-specific caller, which estimates absolute ploidy directly.

### FFPE / low-input dropout read as homozygous deletions

**Trigger:** Degraded FFPE DNA or low input; some bins have near-zero coverage.

**Mechanism:** Zero-coverage bins produce extreme negative log2; CBS joins them into spurious homozygous-deletion segments.

**Symptom:** Scattered tiny "CN=0" segments, often at hard-to-capture (high-GC) loci.

**Fix:** Always pass `--drop-low-coverage` to `batch` and `segment`. Inspect MAD; if MAD > 0.5, the sample is too noisy for confident focal calls.

## Reconciliation: When CNVkit Disagrees With Another Caller

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| CNVkit calls focal events GATK misses | GATK PoN/tangent normalization absorbed the event | Trust CNVkit if the event is rare; suspect GATK if PoN contained tumors |
| GATK/ASCAT call broad arm events CNVkit flattens | CNVkit centered on a non-diploid mode | Re-center against an allele-specific ploidy estimate |
| CNVkit and ASCAT disagree on integer CN | CNVkit purity guess wrong, or genome is WGD | Trust ASCAT/FACETS — joint BAF+depth fit resolves purity/ploidy |
| Tumor-only CNVkit calls absent in matched-normal rerun | Germline CNV or capture bias misread as somatic | Re-run with the matched normal; germline CNVs are not somatic events |

**Operational rule for high-confidence reporting:** Treat a CNVkit call as confident only when (1) the reference was a pooled PoN of process-matched normals, (2) sample MAD < 0.5, (3) the segment spans multiple bins with consistent weight, and (4) for any clinically actionable focal event, it is confirmed by an orthogonal caller or by allele-specific data. Depth-only calls on tumors are screening-grade, not definitive.

## Quality Control

```bash
cnvkit.py metrics results/*.cnr -s results/*.cns      # MAD, spread, bivar per sample
cnvkit.py sex results/*.cnr                           # detect sex / sample swaps
cnvkit.py segmetrics tumor.cnr -s tumor.cns --ci --pi --bootstrap 100 -o tumor.segmetrics.cns
cnvkit.py genemetrics tumor.cnr -s tumor.cns -t 0.2 --ci --bootstrap 100 -o tumor.genemetrics.tsv
```

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Sample MAD (noise) | < 0.5 acceptable; < 0.3 good | CNVkit docs; MAD is the median absolute deviation of bin log2 |
| Panel of normals size | >= 5; 10-20 preferred | Talevich 2016; pooling averages per-normal capture noise |
| Default call thresholds (log2) | -1.1, -0.25, 0.2, 0.7 | CNVkit `call -t` defaults: CN 0 / 1 / 2 / 3 / 4+ boundaries |
| Purity floor for depth calling | ~40% reliable; ~20% absolute floor | Below ~40% segmentation fails (sCNAphase, Gusnanto 2012) |
| `genemetrics -t` gain/loss | 0.2 (default) | 2^0.2 ~ 15% copy-ratio change; tune up for impure samples |
| Antitarget avg bin size | auto; ~target size x fold-enrichment | CNVkit docs; off-target bins should hold comparable read counts |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Spiky `.cnr`, recurrent calls across samples | Flat reference; capture bias misread | Build a pooled PoN; never use flat reference for focal calls |
| Antitarget spread huge, nonsense segments | Hybrid mode on an amplicon panel | Use `--method amplicon` |
| Scattered CN=0 micro-segments | FFPE zero-coverage bins | Add `--drop-low-coverage` to `batch` and `segment` |
| Integer CN systematically too low | `call` without `-m clonal --purity` | Supply purity, or use an allele-specific caller |
| Whole genome called gain or loss | Centering on a non-diploid mode (WGD) | Anchor ploidy with BAF; do not depth-center aneuploid tumors |
| `pomegranate` ImportError on HMM | HMM backend not installed | `pip install pomegranate`; or use `-m cbs` |

## References

- Talevich E et al 2016. CNVkit: genome-wide copy number detection from targeted DNA sequencing. PLoS Comput Biol 12:e1004873
- Olshen AB et al 2004. Circular binary segmentation for the analysis of array-based DNA copy number data. Biostatistics 5:557 (CBS)
- Gusnanto A et al 2012. Correcting for cancer genome size and tumour cell content in whole-genome copy number. Bioinformatics 28:40
- Benjamini Y, Speed TP 2012. Summarizing and correcting the GC content bias in high-throughput sequencing. Nucleic Acids Res 40:e72

## Related Skills

- copy-number/copy-ratio-segmentation - CBS vs HMM choice, depth normalization, bias correction
- copy-number/allele-specific-copy-number - ASCAT/Sequenza/FACETS/PureCN for purity, ploidy, integer ASCN
- copy-number/gatk-cnv - GATK depth-based alternative; tangent normalization
- copy-number/cnv-annotation - Gene and clinical annotation of CNV calls
- copy-number/cnv-visualization - Profile plots, segmentation views, cohort heatmaps
- copy-number/recurrent-cnv - GISTIC2 cohort-level recurrent and driver CNV
- alignment-files/bam-statistics - QC of input BAMs before calling
- long-read-sequencing/structural-variants - Complementary breakpoint-resolved SV calling
