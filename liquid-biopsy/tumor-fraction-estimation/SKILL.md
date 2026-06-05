---
name: bio-tumor-fraction-estimation
description: Estimates tumor fraction (the genome-wide proportion of cfDNA molecules that are tumor-derived, the cfDNA analogue of bulk-tumor purity) from shallow whole-genome sequencing with ichorCNA, an HMM over 1 Mb bins that jointly EM-estimates tumor fraction, ploidy, and subclonal prevalence over a normal/ploidy grid. Encodes the load-bearing reframes: tumor fraction is the quantity that travels across assays and is NOT mutation VAF (clonal-het VAF approximately TF/2), CNA-based estimation has a hard ~3 percent limit-of-detection floor, and near-diploid or copy-neutral-LOH genomes return a falsely low value. Selects the estimator by data type (sWGS to ichorCNA, deep panel to max-VAF, methylation to deconvolution, sub-3 percent to fragmentomics or methylation). Use when quantifying tumor burden from a liquid biopsy, picking a tumor-fraction estimator for a given assay, or reconciling a TF estimate against a panel VAF.
tool_type: r
primary_tool: ichorCNA
---

## Version Compatibility

Reference examples tested with: ichorCNA 0.6.0+ (GavinHaLab fork), HMMcopy 1.40+, R 4.2+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: ichorCNA is NOT an importable R function `runIchorCNA()` — it is a command-line script invoked as `Rscript scripts/runIchorCNA.R` with `optparse` flags, preceded by HMMcopy `readCounter` to build the WIG. Use the GavinHaLab fork (v0.6.0, 22 Nov 2024) for new work; the original broadinstitute/ichorCNA holds the wiki. The flag `--repTimeWig` does not exist — do not invent it.

# Tumor Fraction Estimation

**"Estimate tumor fraction from my cfDNA sample"** -> Estimate the genome-wide proportion of cfDNA molecules that are tumor-derived, mutation-agnostic, from copy-number signal.
- CLI: `readCounter` (HMMcopy) to bin the BAM, then `Rscript scripts/runIchorCNA.R` for the HMM
- R: parse `.params.txt` (tumor fraction = 1 - n); a Python subprocess wrapper is a thin alternative

## The Single Most Important Modern Insight -- tumor fraction is the quantity that travels across assays; it is NOT VAF and it is blind below ~3 percent

Tumor fraction (TF) is the fraction of cfDNA *molecules* that are tumor-derived — the cfDNA analogue of bulk-tumor purity. It is the burden metric that is comparable across assays and over time, which is exactly why it is the right unit to report. The two errors that dominate cfDNA work are unit confusion and floor confusion. Unit confusion: TF is NOT mutation VAF — a clonal heterozygous SNV in a diploid region sits at VAF approximately TF/2, so reporting a max-VAF as "tumor fraction" halves the true burden (and the factor changes entirely under LOH, amplification, or subclonality). Floor confusion: ichorCNA derives TF from copy-number deflection averaged over hundreds of 1 Mb bins, and that signal has a hard ~3 percent limit of detection. Below it the depth shift is smaller than per-bin sampling noise; a near-diploid or copy-neutral-LOH tumor returns a *falsely low* TF even at high true burden because it carries no depth signal. A low ichorCNA value is "low burden" only if the genome-wide plot is genuinely flat; otherwise it is uninformative, not negative.

## Estimator Landscape

| Estimator | Class | Input | Strength | Fails when |
|-----------|-------|-------|----------|------------|
| ichorCNA | CNA / depth (HMM) | sWGS 0.1-1x | Mutation-agnostic genome-wide burden; calibrated standard | TF < ~3%; near-diploid / copy-neutral-LOH genome |
| TitanCNA | CNA + allelic (B-allele) | deeper WGS with het-SNP depth | Resolves CNLOH via allelic imbalance | Needs informative het-SNP coverage (not 0.1x) |
| max-VAF / clonal-cluster MAF | Mutation / panel | deep targeted or WES | Sensitive to <0.1% VAF with UMI/duplex | Needs callable variants + CHIP filtering; CN-sensitive |
| Methylation deconvolution (CelFiE, CelFEER) | Methylation | WGBS/EM-seq or methyl panel | Dense per-molecule signal reaches below CNA floor | Needs a tumor-type methylation reference atlas |
| Fragmentomics (Griffin, DELFI) | Fragmentomic | sWGS | CN-independent corroboration at low TF | Quantifies "tumor signal," not a calibrated molecular fraction |

## Decision Tree by Data Type

| Data available | TF regime | Recommended | Why |
|---------------|-----------|-------------|-----|
| sWGS 0.1-1x, no known variants, aneuploid tumor | >= ~3% | ichorCNA | Mutation-agnostic genome-wide burden; the standard |
| sWGS, tumor type known to be near-diploid / quiet | any | mutation or methylation | ichorCNA underestimates with no depth signal |
| sWGS, TF suspected < 3% | < 3% | deep-panel max-VAF, methylation, or fragmentomics | Below the CNA floor (see fragment-analysis, methylation-based-detection) |
| Deep targeted / WES panel | down to <0.1% VAF | max-VAF excl. CHIP, or clonal-cluster MAF | Per-locus sensitivity; convert via TF approximately 2*VAF with CN care (see ctdna-mutation-detection) |
| Methylation (WGBS/EM-seq/panel) | very low | methylation deconvolution | Dense per-molecule signal; needs reference atlas |
| Targeted panel, want CN-based TF | >= few % | ichorCNA on off-target reads | Recovers genome-wide CN from off-target coverage |

Methodology evolves; verify current best practice against the live ichorCNA wiki and the relevant tool docs before committing to an estimator.

## ichorCNA Mechanics

ichorCNA is a hidden Markov model over copy-number states across 1 Mb bins. The emission per bin is the GC- and mappability-corrected log2 read-depth ratio (tumor vs a panel of normals). The HMM simultaneously segments the genome, calls large-scale CNAs (HOMD/DLOH/NEUT/GAIN/AMP/HLAMP up to `maxCN`), and by EM jointly estimates three global latent parameters: tumor fraction (via `n`), tumor ploidy (`phi`), and subclonal prevalence. The observed copy at a bin is a mixture: copy approximately 2*(1-TF) + TF*(tumor copy), and the sample ploidy identity is 2*(1-TF) + TF*tumor.ploidy. Because TF, ploidy, and per-bin tumor copy are all unknown, the same log-ratio can be explained by (low TF, large CN swing) or (high TF, small CN swing) — this ploidy/TF degeneracy is why ichorCNA fits over a grid of (`normal`, `ploidy`) start points and selects the maximum-likelihood solution.

### Bin the BAM and Run the HMM

**Goal:** Produce a calibrated tumor-fraction estimate plus genome-wide CN segments from a single sWGS BAM.

**Approach:** Bin coverage into 1 Mb WIG with HMMcopy `readCounter` (chromosome naming must match the BAM `@SQ` style), then run `runIchorCNA.R` with build-matched GC/map/centromere references and a protocol-matched panel of normals; read `.params.txt`.

```bash
# Step 1: 1 Mb bins. --chromosome style ('1' vs 'chr1') MUST match the BAM @SQ names.
readCounter --window 1000000 --quality 20 \
  --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
  tumor.bam > tumor.wig

# Step 2: the HMM. NOT an R function call -- it is a script with optparse flags.
Rscript scripts/runIchorCNA.R \
  --id tumor --WIG tumor.wig \
  --gcWig gc_hg38_1000kb.wig --mapWig map_hg38_1000kb.wig \
  --centromere GRCh38.centromere.txt \
  --normalPanel HD_ULP_PoN_1Mb_median.rds \
  --normal "c(0.5,0.6,0.7,0.8,0.9)" --ploidy "c(2,3)" --maxCN 7 \
  --estimateNormal TRUE --estimatePloidy TRUE --estimateScPrevalence TRUE \
  --scStates "c(1,3)" --txnE 0.9999999 --txnStrength 1e7 \
  --minMapScore 0.9 --genomeBuild hg38 --genomeStyle UCSC \
  --outDir ichor_out/
```

Key flags (verified defaults from `runIchorCNA.R`): `--maxCN 7` (lower to 3 for low-TF); `--normal "0.5"` and `--ploidy "2"` are grid *start* points, not fixed values (`--estimateNormal`/`--estimatePloidy` still estimate them; these are optparse `type=logical` flags so they need an explicit `TRUE`/`FALSE`, not a bare flag); `--txnE 0.9999999` and `--txnStrength 1e7` set the segment-length prior; `--minMapScore 0.9` drops low-mappability bins; `--gcWig`/`--mapWig`/`--centromere`/`--normalPanel` must all match the BAM's build and the 1 Mb bin size.

### Parse the Optimal Solution

**Goal:** Extract the calibrated tumor fraction, ploidy, and QC from ichorCNA output.

**Approach:** Read `.params.txt`; TF = 1 - n_est for the selected (max-loglik) solution; gate on the GC-Map MAD; inspect subclonal fractions and the genome-wide plot before trusting a borderline call.

```r
parse_ichor <- function(params_file) {
    p <- read.table(params_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    list(
        tumor_fraction = 1 - p$n_est[1],   # TF = 1 - normal fraction; selected solution is row 1
        ploidy = p$phi_est[1],
        loglik = p$loglik[1]
    )
}
```

The `.params.txt` also carries `Tumor Fraction` (= 1 - n), `Tumor Ploidy` (phi), `Fraction Genome Subclonal`, `Fraction CNA Subclonal`, and `GC-Map Correction MAD` (the data-noise QC). Companion outputs: `.cna.seg` (per-bin CN and log-ratio), `.seg` (IGV-compatible Viterbi segments), `.RData` (all grid solutions), and the genome-wide plot PDF — always inspect it for borderline calls because the ploidy/TF degeneracy can select a ploidy-3 alias of a ploidy-2 truth.

## Unit Confusion -- TF vs VAF vs ctDNA%

These three are routinely conflated; the relation is exact and copy-number-dependent. For a variant at local copy number `Cn` with mutant-copy multiplicity `m`:

```
VAF = (TF * m) / [ TF * Cn + 2 * (1 - TF) ]
```

For a clonal heterozygous SNV in a diploid region (`Cn`=2, `m`=1) this collapses to VAF approximately TF/2, equivalently **TF approximately 2*VAF**. The common errors:

- LOH variant (mutant on both copies, normal copy lost): `m`=`Cn`, so VAF -> TF, not TF/2 — treating it as TF/2 doubles the estimate.
- Amplified mutant allele inflates VAF above TF/2; a mutant on a deleted copy deflates it — so max-VAF over-estimates TF for amplified drivers and under-estimates for deleted ones.
- Subclonal variants carry an extra cancer-cell-fraction factor and understate TF.
- Germline heterozygous SNPs sit at VAF approximately 0.5 regardless of TF — never feed them into a TF-from-VAF calculation.
- CHIP (clonal hematopoiesis) variants are blood-derived, not tumor — exclude them from any max-VAF TF proxy.
- ctDNA% is loosely used for either TF or max-VAF; always pin down which a lab means.

Cross-check: for a clonal heterozygous driver in a diploid region, ichorCNA TF and 2*(panel VAF) should agree. TF >> 2*VAF implies a subclonal/deleted variant or a ploidy mis-call; TF << 2*VAF implies a near-diploid/CNLOH tumor or an amplified/LOH driver. Never average the two blindly (see ctdna-mutation-detection).

## Per-Method Failure Modes

### ~3 percent CNA floor
**Trigger:** TF below ~0.03 at 0.1x sWGS. **Mechanism:** the log2 deflection from a single-copy event is proportional to TF (~±0.02 at TF=0.03), smaller than per-bin sampling noise; only averaging over hundreds of bins recovers it. **Symptom:** TF collapses toward 0; replicate variability (MNSD) rises sharply. **Fix:** the floor scales with aneuploidy magnitude and coverage — it needs roughly one >100 Mb gain AND one >100 Mb loss; sequence deeper (>1-5x) or switch estimator class (fragment-analysis, methylation-based-detection).

### Near-diploid / copy-neutral-LOH
**Trigger:** quiet tumor type or CNLOH-rich genome. **Mechanism:** CNLOH has identical total coverage to diploid, indistinguishable on depth alone; ichorCNA is also tuned conservative and "may underestimate." **Symptom:** falsely low TF with a flat genome-wide plot. **Fix:** treat a flat low call as uninformative, not negative; escalate to a mutation/methylation assay; TitanCNA can use allelic imbalance if het-SNP depth exists.

### Mismatched panel of normals / references
**Trigger:** PoN, GC/map/centromere WIG, or build does not match the library prep, bin size, or genome build. **Mechanism:** the PoN models protocol-specific coverage bias; a mismatched PoN injects its own bias as spurious CN waviness. **Symptom:** wavy log-ratio, implausible TF. **Fix:** build/obtain a PoN from healthy-donor cfDNA on the exact protocol at the same bin size and build; keep hg19 vs hg38 and `1` vs `chr1` consistent end-to-end.

### Ploidy aliasing
**Trigger:** ploidy/TF degeneracy. **Mechanism:** the max-loglik solution is occasionally a ploidy-3 alias of a ploidy-2 truth. **Symptom:** doubled ploidy with halved TF. **Fix:** read all `.params.txt` solutions, inspect the plot; for low-TF samples force `--ploidy "c(2)"`.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Coverage 0.1-1x sWGS; 1 Mb bins | Adalsteinsson 2017; ichorCNA wiki | Finer bins add noise at 0.1x; ~0.1x is the calibrated ULP-WGS operating point |
| ~3% TF limit of detection at ~0.1x | Adalsteinsson 2017 (95% sens, 1125/1288 mixtures; 91% spec, 20/22 donors at 0.03 TF cutoff) | Below 0.03 the depth deflection falls under per-bin noise |
| 97.2-100% sensitivity to detect 3% TF (1x and 0.1x) | J Mol Diagn 2024 assay validation | Independent dilution/replicate validation; MNSD rises sharply below 3%, establishing 3% as the LOD |
| GC-Map Correction MAD < 0.15 good; > 0.3 distrust | ichorCNA FAQ | Residual post-correction noise; high MAD means the depth signal is unreliable |
| Manual-curation band 0.03-0.10 TF | ichorCNA wiki | Model can pick the wrong solution and tends to underestimate near the floor; inspect the plot |
| Low-TF recipe: `--normal "c(0.95,0.99,0.995,0.999)" --ploidy "c(2)" --maxCN 3 --estimateScPrevalence FALSE --scStates "c()"` | ichorCNA wiki | Seeds EM near TF 5/1/0.5/0.1%; ploidy and subclonality are unidentifiable when CN signal is weak |

## References

- Adalsteinsson VA, Ha G, Freeman SS, et al. 2017. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. *Nat Commun* 8(1):1324. — ichorCNA primary method; the ~3% LOD benchmark (95% sensitivity, 91% specificity at a 0.03 TF cutoff, ~0.1x).
- Assay Validation of Cell-Free DNA Shallow Whole-Genome Sequencing to Determine Tumor Fraction in Advanced Cancers. 2024. *J Mol Diagn* 26(5):413-422 (PMC11090203). — Independent validation: 97.2-100% sensitivity at 3% TF (1x and 0.1x); MNSD rising below 3% establishes 3% as the LOD.
- broadinstitute/ichorCNA and GavinHaLab/ichorCNA GitHub repositories and wiki (Usage, Output, Parameter-tuning, Create-Panel-of-Normals, FAQ). — `readCounter` command, `runIchorCNA.R` flag defaults, `.params.txt` fields, MAD QC thresholds, CNLOH/near-diploid underestimation, PoN construction.

## Related Skills

- cfdna-preprocessing - sWGS BAM input and minimal-processing path
- fragment-analysis - the estimator to use below the ~3% CNA floor
- ctdna-mutation-detection - max-VAF cross-check and the TF-vs-VAF reconciliation
- analytical-validation - the ~3% floor framed as a limit of detection
- copy-number/cnvkit-analysis - copy-number calling concepts
- copy-number/copy-ratio-segmentation - segmentation concepts
