---
name: bio-copy-number-allele-specific-copy-number
description: Infer integer allele-specific copy number, tumor purity, and ploidy from tumor sequencing by jointly modeling read depth (logR) and B-allele frequency (BAF) with ASCAT, Sequenza, FACETS, PURPLE, and PureCN (tumor-only). Covers the purity-ploidy identifiability problem, the diploid-baseline (dipLogR) anchor, major/minor copy number, loss of heterozygosity, sunrise/contour fit diagnostics, and reconciliation of conflicting fits. Use when tumor analysis needs absolute copy number rather than relative log2, when estimating purity and ploidy, calling LOH or copy-neutral LOH, resolving whole-genome doubling, running tumor-only allele-specific calling, or choosing among ASCAT, Sequenza, FACETS, and PureCN.
tool_type: mixed
primary_tool: ascat
---

## Version Compatibility

Reference examples tested with: ASCAT 3.1+, Sequenza 3.0+ (sequenza-utils 3.0+), FACETS 0.6+ (snp-pileup), PureCN 2.6+, R 4.3+, Python 3.10+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('ASCAT')` / `'sequenza'` / `'facets'` / `'PureCN'`, then `?function`
- CLI: `sequenza-utils --version`, `snp-pileup --help`

Sequenza 3.0 depends on the `copynumber` Bioconductor package, REMOVED from Bioconductor 3.18+ (2023). Install a maintained fork (`ShixiangWang/copynumber` or `igordot/copynumber`) before Sequenza will load. ASCAT's GC-correction function was renamed across 2.x->3.x (`ascat.GCcorrect` -> `ascat.correctLogR`) — verify against the installed version.

# Allele-Specific Copy Number

**"How many copies of each allele, in what fraction of cells, at what tumor purity"** -> Jointly model read depth and B-allele frequency to fit tumor purity, ploidy, and integer major/minor copy number per segment. Depth alone gives only *relative* copy ratio; depth + BAF gives *absolute* allele-specific copy number. This skill is required whenever the question involves LOH, absolute copy number, purity, ploidy, or whole-genome doubling — CNVkit and GATK somatic CNV cannot answer those.

- R: `ASCAT` (WGS, SNP array), `sequenza` (WES/WGS), `facets` (panel/WES/WGS), `PureCN` (tumor-only panel/WES)
- CLI: `purple` (Hartwig WGS pipeline, with AMBER + COBALT)

## The Identifiability Problem — Why This Is Hard

Purity and ploidy are **not identifiable from depth alone**. The same log-ratio profile is explained equally well by many (purity, ploidy) pairs: a homozygous deletion at 30% purity looks identical to a heterozygous deletion at 60% purity; an entire profile can be reinterpreted at 2x ploidy with halved purity. Every allele-specific caller breaks this degeneracy by adding BAF — allelic imbalance constrains which solution is real. The consequence: the likelihood surface is **multimodal**, the fit can lock onto an integer-multiple of the true ploidy, and a single point estimate must never be trusted without inspecting the fit diagnostic (ASCAT sunrise plot, Sequenza cellularity/ploidy contour, FACETS dipLogR). A real published example: the same tumor scored ploidy 4.27 by FACETS (WGS) and 2.42 by ASCAT (SNP array).

## Caller Taxonomy

| Tool | Input | Segmentation | Best for | Fails when |
|------|-------|--------------|----------|------------|
| ASCAT | SNP array or WGS logR+BAF | ASPCF (allele-specific PCF) | WGS, SNP6, large cohorts | Near-diploid genome with few aberrations cannot anchor purity -> defaults toward purity ~100% |
| Sequenza | Tumor-normal WES/WGS (seqz) | `copynumber` PCF | Exome, accessible install | Picks a near-diploid local optimum; needs manual review of alternative solutions |
| FACETS | Tumor-normal, snp-pileup | Joint logR+BAF CBS | Targeted panels, WES, clinical NGS | `cval` too low -> hyperfragmentation; EM locks onto an integer-multiple ploidy |
| PURPLE | WGS, AMBER+COBALT+SVs | Integrates SV breakpoints | WGS with matched SV calls | Targeted/WES (designed for WGS); needs the Hartwig tool stack |
| PureCN | Tumor-only WES/panel + PoN | Coverage + VCF, normal DB | No matched normal | Sparse hets; small panels; needs a well-built normal database |
| Battenberg | WGS logR+BAF, phased | ASCAT-based clonal + subclonal | Subclonal CN, clonal evolution | Heavy; needs phasing reference — see subclonal-copy-number |

## Decision Tree by Data Type

| Scenario | Recommended caller | Rationale |
|----------|--------------------|-----------|
| Tumor-normal WGS | ASCAT or PURPLE | PURPLE if SV calls available (resolves breakpoints); ASCAT otherwise |
| Tumor-normal WES | Sequenza or FACETS | Both joint logR+BAF; FACETS faster, Sequenza reports alternative solutions |
| Targeted panel (tumor-normal) | FACETS | Designed for panel het density; clinical-NGS standard |
| Tumor-only panel / WES | PureCN | Models a normal database; the standard tumor-only solution |
| SNP array (legacy) | ASCAT | ASCAT was built for SNP arrays |
| Subclonal CN / clonal evolution | Battenberg / TITAN | See subclonal-copy-number |
| Only relative gain/loss needed | CNVkit / GATK | Allele-specific machinery is unnecessary |

## FACETS — Tumor-Normal Allele-Specific CN

**Goal:** Fit purity, ploidy, and integer allele-specific CN for a panel or WES pair.

**Approach:** Pile up read counts at common SNPs with `snp-pileup`, then run the two-pass FACETS workflow — a high-`cval` purity run whose `dipLogR` seeds a low-`cval` sensitivity run for focal events.

```bash
# Step 1: pileup at dbSNP common sites (normal first, then tumor)
snp-pileup -g -q15 -Q20 -P100 -r25,0 dbsnp_common.vcf.gz \
    sample.snp_pileup.csv.gz normal.bam tumor.bam
```

```r
library(facets)
set.seed(1234)                                  # FACETS uses random initialization
rcmat <- readSnpMatrix('sample.snp_pileup.csv.gz')
xx <- preProcSample(rcmat)                      # gbuild default 'hg19'; pass gbuild='hg38' for GRCh38

# Pass 1: purity/ploidy at a coarse cval (panels ~150-300; WGS ~25-100)
oo1 <- procSample(xx, cval = 300)
fit1 <- emcncf(oo1)

# Pass 2: focal sensitivity, seeded by the diploid baseline from pass 1
oo2 <- procSample(xx, cval = 150, dipLogR = oo1$dipLogR)
fit2 <- emcncf(oo2)

cat('purity', fit2$purity, 'ploidy', fit2$ploidy, 'dipLogR', oo2$dipLogR, '\n')
# fit2$cncf has per-segment tcn.em (total CN) and lcn.em (minor CN); lcn.em == 0 -> LOH
plotSample(x = oo2, emfit = fit2)               # ALWAYS inspect this diagnostic plot
```

## Sequenza — Exome Tumor-Normal

**Goal:** Estimate cellularity/ploidy and allele-specific CN from a WES pair, with explicit alternative solutions.

**Approach:** Build a `seqz` file from the BAMs, bin it, then run the extract/fit/results chain; inspect the cellularity/ploidy contour and the reported alternative solutions.

```bash
sequenza-utils bam2seqz -n normal.bam -t tumor.bam --fasta ref.fa \
    -gc hg38.gc50.wig.gz -o sample.seqz.gz
sequenza-utils seqz_binning --seqz sample.seqz.gz -w 50 -o sample.bin.seqz.gz
```

```r
library(sequenza)
seqz <- sequenza.extract('sample.bin.seqz.gz')
CP <- sequenza.fit(seqz)                        # grid search over cellularity x ploidy
sequenza.results(seqz, CP, 'sampleID', out.dir = 'sequenza_out')
# Inspect *_CP_contours.pdf and *_alternative_solutions.txt before accepting the fit.
```

## ASCAT — WGS / SNP Array

**Goal:** Fit purity (rho), ploidy (psi), and allele-specific CN genome-wide.

**Approach:** Load logR/BAF, correct for GC (and optionally replication timing), segment with ASPCF, run the ASCAT fit, and read the sunrise plot.

```r
library(ASCAT)
ascat.bc <- ascat.loadData('Tumor_LogR.txt', 'Tumor_BAF.txt',
                           'Germline_LogR.txt', 'Germline_BAF.txt')
ascat.bc <- ascat.correctLogR(ascat.bc, GCcontentfile = 'GC_G1000.txt',
                              replictimingfile = 'RT_G1000.txt')   # RT optional
ascat.bc <- ascat.aspcf(ascat.bc)
ascat.output <- ascat.runAscat(ascat.bc, gamma = 1)   # gamma=1 for NGS; ~0.55 for arrays
# ascat.output$purity, $ploidy, $goodnessOfFit; $nA / $nB are major/minor CN per segment
# Inspect the sunrise plot: banding at multiples of ploidy signals an ambiguous fit.
```

## PureCN — Tumor-Only

**Goal:** Recover purity, ploidy, allele-specific CN, and LOH without a matched normal.

**Approach:** Build a normal database (PoN) once, then run `runAbsoluteCN` with the tumor coverage and a VCF; PureCN uses the normal DB and a mapping-bias model in place of a matched normal.

```r
library(PureCN)
ret <- runAbsoluteCN(
    tumor.coverage.file = 'tumor_coverage.txt.gz',
    vcf.file = 'tumor.vcf.gz',
    normalDB = readRDS('normalDB.rds'),       # built once from >= ~20 process-matched normals
    genome = 'hg38', sampleid = 'tumor',
    interval.file = 'baits_intervals.txt')
# ret$results[[1]]$purity / $ploidy; createCurationFile() flags fits needing manual review
```

## Failure Modes

### ASCAT defaults to ~100% purity on a near-diploid genome

**Trigger:** A tumor with very few copy-number aberrations and overall ploidy near 2.

**Mechanism:** ASCAT infers purity from the depth/BAF deviation of aberrant segments. With almost no aberrant segments there is nothing to anchor purity against, so the grid search drifts to the boundary.

**Symptom:** Reported purity ~1.0 (or implausibly high) with an almost flat profile; the sunrise plot is nearly featureless.

**Fix:** Treat purity as indeterminate, not 100%. Cross-check with an orthogonal estimate (SNV VAF mode for clonal mutations, pathology estimate). A genuinely quiet genome simply does not support a confident purity call.

### FACETS hyperfragmentation from too-low cval

**Trigger:** `cval` set too low for the data (e.g. panel data run at WGS-scale cval).

**Mechanism:** `cval` is the segmentation critical value; low values let the segmenter split on noise, shattering the profile into spurious micro-segments.

**Symptom:** Hundreds of tiny segments; `tcn.em`/`lcn.em` incoherent with `cnlr.median`; jagged `plotSample` output.

**Fix:** Use cval ~150-300 for panels/WES, ~25-100 for WGS. Run the two-pass workflow (coarse purity run -> dipLogR-seeded sensitivity run). If naive `tcn` and EM `tcn.em` disagree wildly, the fit is bad — re-tune cval.

### Integer-multiple ploidy flip

**Trigger:** Any allele-specific caller on a genome where the diploid baseline is ambiguous (few hets, low purity, or genuine WGD).

**Mechanism:** The likelihood surface has near-equal modes at ploidy P and 2P; the optimizer can select the wrong one, halving or doubling all copy numbers.

**Symptom:** Two callers disagree by a factor of ~2 in ploidy; "balanced" CN states that should be odd come out even (or vice versa); SNV multiplicities inconsistent with the called CN.

**Fix:** Inspect the fit diagnostic (sunrise/contour). Cross-check ploidy against the fraction of the genome at odd vs even CN and against clonal-SNV VAF. Prefer the solution consistent with known biology; if truly ambiguous, report both.

### Sequenza fails to load / picks a near-diploid optimum

**Trigger:** Fresh Sequenza install on Bioconductor 3.18+; or accepting `sequenza.fit`'s point estimate without review.

**Mechanism:** Sequenza depends on `copynumber`, removed from Bioconductor 3.18+. Separately, the LPP grid search can settle on a near-diploid local optimum when a higher-ploidy solution fits comparably.

**Symptom:** `copynumber` not available at load; or a ploidy ~2 call that conflicts with visible large-scale imbalance.

**Fix:** Install a maintained `copynumber` fork. Always inspect `*_CP_contours.pdf` and `*_alternative_solutions.txt`; if a non-diploid alternative fits nearly as well and matches the BAF pattern, prefer it.

### Low-purity death zone

**Trigger:** Tumor purity below ~40% (common in breast, lung adenocarcinoma, melanoma).

**Mechanism:** Allelic imbalance and depth deviation both shrink with purity; below ~40% the signal approaches the noise floor and segmentation fails.

**Symptom:** No confident fit; purity estimate unstable across reruns; flat BAF.

**Fix:** Below ~40% purity, allele-specific calling is unreliable; below ~20% it is not possible with bulk sequencing. Report indeterminate; consider deeper sequencing or microdissection.

## Reconciliation: When Callers Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Caller A ploidy ~= 2x caller B | Integer-multiple ploidy flip | Check odd/even CN fraction and SNV multiplicity; pick the biology-consistent fit |
| Purity differs widely, ploidy agrees | One caller hit a boundary on a quiet genome | Trust the caller whose diagnostic plot shows real structure |
| FACETS vs ASCAT integer CN differ | Different segmentation (CBS vs ASPCF) at boundaries | Compare segment edges; arm-level calls usually agree, focal may not |
| Tumor-only (PureCN) vs tumor-normal differ | Tumor-only has weaker purity constraint | Prefer the matched-normal fit when available |

**Operational rule:** Report an allele-specific fit as confident only when (1) the fit diagnostic (sunrise/contour/dipLogR) shows clear, non-degenerate structure, (2) purity is above ~40%, (3) ploidy is consistent with the odd/even CN fraction and with clonal-SNV multiplicity, and (4) for ambiguous cases, the alternative solutions have been reviewed. A bare purity/ploidy number with no diagnostic inspection is not a result.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Purity floor | ~40% reliable; ~20% absolute floor | Below ~40% segmentation fails (Gusnanto 2012; sCNAphase) |
| FACETS cval (panel/WES) | 150-300 | FACETS docs; lower -> hyperfragmentation |
| FACETS cval (WGS) | 25-100 | FACETS docs; scales with marker density |
| ASCAT gamma | 1.0 (NGS); ~0.55 (SNP array) | ASCAT docs; platform-specific logR shrinkage |
| LOH definition | minor CN (lcn) = 0 | Minor allele lost; total CN may still be >= 2 (CN-neutral LOH) |
| PureCN normal DB size | >= ~20 process-matched normals | PureCN docs; mapping-bias and coverage model |
| Het SNP density for stable BAF | thousands genome-wide / hundreds per arm | Sparse hets give noisy allele-fraction segmentation |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Sequenza: `copynumber` not found | Removed from Bioconductor 3.18+ | Install a maintained `copynumber` fork |
| `ascat.GCcorrect` not found | Renamed in ASCAT 3.x | Use `ascat.correctLogR` |
| FACETS profile shattered | cval too low | Raise cval; two-pass workflow |
| Purity reported ~1.0, flat genome | Near-diploid, unanchored | Report indeterminate; cross-check with SNV VAF |
| All CN halved or doubled vs expectation | Integer-multiple ploidy flip | Inspect fit diagnostic; check SNV multiplicity |
| PureCN unstable tumor-only fit | Sparse hets / weak normal DB | Larger normal DB; deeper sequencing; flag for curation |

## References

- Van Loo P et al 2010. Allele-specific copy number analysis of tumors. PNAS 107:16910 (ASCAT)
- Ross EM et al 2021. Allele-specific multi-sample copy number segmentation in ASCAT. Bioinformatics 37:1909
- Favero F et al 2015. Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data. Ann Oncol 26:64
- Shen R, Seshan VE 2016. FACETS: allele-specific copy number and clonal heterogeneity analysis. Nucleic Acids Res 44:e131
- Riester M et al 2016. PureCN: copy number calling and SNV classification using targeted short read sequencing. Source Code Biol Med 11:13
- Priestley P et al 2019. Pan-cancer whole-genome analyses of metastatic solid tumours (PURPLE). Nature 575:210

## Related Skills

- copy-number/copy-ratio-segmentation - logR normalization and segmentation feeding these callers
- copy-number/subclonal-copy-number - Battenberg/TITAN subclonal CN, whole-genome doubling
- copy-number/hrd-scoring - LOH/LST/TAI scars computed from allele-specific output
- copy-number/cnvkit-analysis - Relative depth-only calling (when allelic resolution is not needed)
- copy-number/gatk-cnv - GATK somatic CNV (relative; no purity/ploidy)
- variant-calling/vcf-basics - SNV VCFs supplying BAF and clonal-mutation cross-checks
