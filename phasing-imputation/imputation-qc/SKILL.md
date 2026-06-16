---
name: bio-phasing-imputation-imputation-qc
description: Assesses and filters phasing/imputation output - the quality metrics (Beagle DR2, Minimac R2 and EmpRsq, IMPUTE/GLIMPSE INFO), MAF-stratified filtering, true accuracy by masking, the differential-imputation confound, dosage-based downstream usage, and phasing switch-error QC. Covers why every routine quality score is an ESTIMATE of r2 from the posterior spread (not validation against truth), why it is confounded with MAF so a flat INFO>=0.3 cutoff is a hidden rare-variant filter, why concordance lies for rare variants while masked dosage-r2 by MAF is the gold standard, why separate case/control imputation manufactures false GWAS hits, and that the field name tells the tool (DR2=Beagle, R2=Minimac, INFO=GLIMPSE/IMPUTE). Use when filtering imputed variants before GWAS, validating accuracy, benchmarking phasing against trios, or diagnosing inflated association. Imputation is genotype-imputation; phasing is haplotype-phasing; panel ancestry is reference-panels; the test is population-genetics/association-testing.
tool_type: mixed
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, cyvcf2 0.31+, pandas 2.2+, numpy 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show cyvcf2 pandas numpy` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The imputation quality field is named differently by each engine: `DR2` (Beagle 5.x; `AR2` is legacy 4.x), `R2` (Minimac4; `ER2`/`EmpRsq` for typed sites), `INFO` (IMPUTE5 and GLIMPSE). `bcftools +fill-tags` computes AF/MAF/HWE but CANNOT produce an imputation Rsq - the quality number comes from the imputer only. Record the engine, panel, and build, because a quality number is only comparable within the same engine and panel.

# Imputation QC -- Filtering on a Self-Estimated Quality, Validating With Masked Truth

**"Filter my imputed genotypes by quality and check the accuracy"** -> Filter on the imputer's per-variant quality field with a MAF floor, and validate true accuracy by masking - because the routine quality score is the model grading its own posterior, it is confounded with allele frequency, and a flat cutoff silently deletes the rare variants of greatest interest.
- CLI: `bcftools view -e 'INFO/DR2<0.3 || INFO/AF<0.01 || INFO/AF>0.99' imputed.vcf.gz` (DR2 for Beagle; R2 for Minimac; INFO for GLIMPSE/IMPUTE; Minimac also emits a MAF tag, Beagle only AF)

Scope: QC of already-produced phasing/imputation output - what the numbers mean, which to trust, how to filter, and how poor QC propagates into false GWAS hits. Running imputation -> genotype-imputation. Phasing -> haplotype-phasing. Panel ancestry (which the metric cannot detect) -> reference-panels. The GWAS test on the filtered dosages -> population-genetics/association-testing. VCF field-parsing mechanics -> variant-calling/vcf-statistics. Read-backed phasing switch QC -> long-read-sequencing/haplotype-phasing.

## The Single Most Important Modern Insight -- The Quality Score Is the Model Grading Its Own Posterior, Confounded With MAF, So a Flat Cutoff Is a Hidden Rare-Variant Filter

Every routine imputation quality score (IMPUTE INFO, Minimac R2/Rsq, Beagle DR2) estimates the same quantity - the squared correlation r2 between the imputed dosage and the unobserved true genotype - from the posterior spread, without ever seeing truth (Marchini & Howie 2010 *Nat Rev Genet* 11:499). It is a self-report of confidence, not a measured accuracy. Three facts organize all of QC:

1. **The metric is structurally confounded with MAF.** The estimator is the fraction of the HWE-expected dosage variance 2p(1-p) the imputed dosages recover; the denominator collapses toward zero as the allele gets rarer, so the estimate is noisier and systematically lower for rare variants. A single flat INFO/R2 >= 0.3 cutoff therefore deletes rare variants far more aggressively than common ones - a hidden MAF filter applied without anyone deciding to. Filter MAF-stratified, or report accuracy per MAF bin.
2. **The only true accuracy comes from masking.** Hide genotypes that are actually known, re-impute them, and compute the real squared correlation against the held-out truth, binned by MAF (the dosage-r2-by-MAF curve). That is the gold standard; the per-variant R2 is a convenient proxy never validated against truth for that variant. Minimac's EmpRsq (masked, measured) vs Rsq (self-estimated) is the only place the two meet, and a gap between them is a QC alarm (panel/strand/ancestry mismatch). Concordance is disqualified for rare variants: a do-nothing imputer that always calls the major homozygote scores ~98% concordance at MAF 1% (0.99^2) while carrying zero information.
3. **The differential-imputation confound manufactures false hits, and every per-group QC metric passes.** Because imputation quality is a function of the data and panel, imputing cases and controls (or batches, or ancestries) separately makes the imputation error differ between groups, and a case-control test cannot distinguish that artifactual allele-frequency difference from real association. The fix is structural (impute together / harmonize), not a filter; the absence of a QC red flag never clears it.

## The Metrics, Precisely

All three estimate r2 (imputed dosage vs true genotype) from the posterior, never from observed truth. The variance-ratio intuition: `estimated r2 = Var(imputed dosage) / [2p(1-p)]`. Under perfect information the dosages equal the true 0/1/2 genotypes and recover the full HWE variance (ratio 1); under no information every dosage collapses to the mean 2p and the variance goes to 0 (ratio 0). The fraction of HWE variance recovered IS the estimate, which is also why it is noisier and lower at low MAF. That `Var(dosage)/2p(1-p)` form is specifically the Minimac/Beagle estimator (the spread of the point dosages across samples); IMPUTE INFO targets the same r2 but computes it differently, averaging each sample's within-individual posterior variance - the mechanistic reason the three numbers are not interchangeable across engines.

| Field | Engine | Meaning |
|-------|--------|---------|
| DR2 | Beagle 5.x | estimated squared correlation between estimated and true allele dose (AR2 is legacy 4.x) |
| R2 (Rsq) | Minimac4 | estimated r2 for all sites, from the dosage variance ratio |
| EmpRsq / ER2 (EmpR) | Minimac4 | empirical r2 from leave-one-out at TYPED sites (masked, measured); a negative EmpR flags a strand/allele flip |
| INFO | IMPUTE5, GLIMPSE | the IMPUTE information measure (posterior-variance ratio), same spirit |

The provenance matters: the field name tells the tool, the numbers are NOT comparable across engines (a Beagle DR2 of 0.8 is not a Minimac R2 of 0.8), and `bcftools +fill-tags` cannot produce any of them.

## The INFO x MAF Interaction and Thresholds

At low MAF there is little dosage variance to predict and few panel copies of the rare haplotype to anchor the estimate, so R2 is intrinsically noisier and downward-biased there - a property of the construction, not a bug. The standard cutoffs and their status:

| Threshold | Source / status | Rationale |
|-----------|-----------------|-----------|
| INFO/R2/DR2 >= 0.3 | the common GWAS cutoff; convention, NOT a theorem | WTCCC/early-IMPUTE-era practice, acknowledged as somewhat arbitrary; acts as a hidden MAF filter |
| INFO/R2 >= 0.8 | stricter, high-confidence analyses | the 0.3 and 0.8 pair are the two established values |
| MAF-stratified INFO | recommended remedy | a flat cutoff differentially deletes rare variants; filter or report per MAF bin |
| GP > 0.9 (or 0.8) hardcall threshold | when forced to hardcall | below this set missing; hardcalling at all loses power vs dosage regression |
| Meta-analysis N_eff = N * INFO | METAL/GWAMA/FinnGen convention | down-weight poorly-imputed variants; not a named theorem |
| EmpRsq vs Rsq large gap | QC alarm | self-estimate not matching measured accuracy = panel/strand/ancestry mismatch |

## True Accuracy by Masking

The accepted accuracy curve is aggregate r2 (squared Pearson correlation between imputed dosage and the masked-then-revealed true genotype) binned by MAF; it decreases monotonically as MAF falls. Workflow: mask typed genotypes (array sites or a held-out set), re-impute from the panel, compute squared correlation vs the withheld truth, bin by MAF. Leave-one-out (one variant at a time) is the per-variant version and is what produces Minimac's EmpRsq for typed sites. Never use concordance as the rare-variant accuracy metric (Ramnarine 2015 *PLoS One* 10:e0137601); it is dominated by the major-homozygote class and inflates rare-variant accuracy.

### MAF-stratified quality summary

**Goal:** Reveal the hidden-MAF-filter effect by reporting imputation quality per MAF bin instead of as one global mean, so the rare-variant tail a flat cutoff would delete is visible.

**Approach:** Parse the imputed VCF for the engine's quality field and allele frequency with cyvcf2, derive MAF, bin it, and report mean quality and the fraction passing a candidate cutoff per bin.

```python
import numpy as np
import pandas as pd
from cyvcf2 import VCF

def quality_by_maf(vcf_path, qual_key='DR2', cutoff=0.3):
    rows = []
    for v in VCF(vcf_path):
        q = v.INFO.get(qual_key)
        af = v.INFO.get('AF')
        if q is None or af is None:
            continue
        af = af[0] if isinstance(af, tuple) else af
        rows.append((min(af, 1 - af), float(q)))
    df = pd.DataFrame(rows, columns=['maf', 'qual'])
    bins = [0, 0.001, 0.01, 0.05, 0.5]   # rare-to-common; the rare bins are where a flat cutoff bites
    df['maf_bin'] = pd.cut(df['maf'], bins=bins)
    summary = df.groupby('maf_bin', observed=True).agg(n=('qual', 'size'), mean_qual=('qual', 'mean'), frac_pass=('qual', lambda q: (q >= cutoff).mean()))
    return summary

quality_by_maf('imputed.vcf.gz', qual_key='DR2', cutoff=0.3)
```

## Phasing QC -- Switch Error Rate

Switch error rate (SER) = switch errors / opportunities, an opportunity being each consecutive heterozygous-site pair, scored against trio/duo/benchmark truth. A flip error is two switches one site apart (an isolated mis-assignment, not a long-range switch); Hamming distance counts overall haplotype differences and inflates on block swaps. Trio-based SER is biased upward by genotype error, which is why SHAPEIT5's `switch` reports SER and a genotyping-error rate jointly. Magnitudes are always MAC- and N-stratified (sub-0.5% common-variant SER for modern tools; single-digit-percent and rising as MAC approaches 1). Tools: `whatshap compare`, `vcftools --diff-switch-error`, SHAPEIT5 `switch`.

## Dosage-Based Downstream Usage

Hardcall thresholding discards imputation uncertainty and sets low-confidence calls missing, losing power versus regression on the expected dosage, especially at low MAF (Huang 2014 *PLoS One* 9:e110679). Carry dosages: PLINK2 (`--vcf file dosage=DS`, `dosage=HDS` for Minimac4 phased, `.pgen`), SNPTEST (`-method expected`/`score`/`em`), REGENIE (BGEN v1.2/PGEN), and BOLT-LMM (BGEN v1.2) all accept dosages. In meta-analysis, INFO enters as a per-study filter and as an effective-N weight.

## Per-Method Failure Modes

### Flat INFO cutoff as a silent rare-variant filter
**Trigger:** applying one INFO/R2 >= 0.3 across all frequencies. **Mechanism:** the metric is confounded with MAF, so the cutoff removes a far higher fraction of rare than common variants. **Symptom:** the rare-variant tail vanishes with no record that a frequency filter was applied. **Fix:** filter MAF-stratified or report accuracy per MAF bin; pair any INFO cutoff with an explicit MAF floor stated in the methods.

### Concordance reported as rare-variant accuracy
**Trigger:** quoting genotype concordance for rare variants. **Mechanism:** concordance is dominated by the major-homozygote class; a do-nothing imputer scores ~98% at MAF 1%. **Symptom:** uniformly high concordance hiding catastrophic rare-variant failure. **Fix:** use masked dosage-r2 binned by MAF; reserve concordance for sanity checks on common variants.

### Cases and controls imputed separately
**Trigger:** imputing batches/groups independently. **Mechanism:** batch-differential imputation quality creates an artifactual allele-frequency difference indistinguishable from association. **Symptom:** genome-wide-significant hits that fail to replicate; every per-group QC passes. **Fix:** impute all samples together, or harmonize panel/version and verify quality does not differ by batch -> reference-panels.

### Filtering on the wrong quality field
**Trigger:** `bcftools view -i 'INFO/R2>0.3'` on a Beagle VCF (which has DR2, not R2). **Mechanism:** the field name is engine-specific. **Symptom:** the filter silently passes everything or errors on a missing tag. **Fix:** use DR2 for Beagle, R2 for Minimac, INFO for GLIMPSE/IMPUTE; confirm with `bcftools view -h`.

### Trusting Rsq where EmpRsq diverges
**Trigger:** reporting a high Rsq while the masked EmpRsq is much lower. **Mechanism:** the self-estimate assumes the model (panel, strand, ancestry) is right; a gap means it is not. **Symptom:** confident Rsq on systematically wrong imputation; a negative EmpR is an outright strand flip. **Fix:** treat the Rsq-EmpRsq gap as an alarm; check strand/build/ancestry -> reference-panels.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| INFO/R2/DR2 >= 0.3 (common), >= 0.8 (strict) | convention (WTCCC/early-IMPUTE era) | the common GWAS cutoffs; not derived, and a hidden MAF filter |
| Always pair the quality cutoff with a MAF floor | Magi 2012 *Genet Epidemiol* 36:785 | rare + low-R2 is the classic false-positive generator |
| Accuracy = masked dosage-r2 binned by MAF | Ramnarine 2015 *PLoS One* 10:e0137601 | concordance inflates rare-variant accuracy; r2 by MAF bin is the gold standard |
| Hardcall GP > 0.9 only when forced | convention | hardcalling loses power vs dosage regression |
| Impute cases and controls together | best-practice consensus | separate imputation manufactures batch-driven false positives |
| Switch-error magnitudes are MAC/N-stratified | Hofmeister 2023 *Nat Genet* 55:1243 | no single universal SER threshold; qualify by MAC bin and validation type |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Filter on INFO/R2 passes all Beagle variants | wrong field name (Beagle uses DR2) | use DR2; confirm with `bcftools view -h` |
| Rare-variant signal disappears after QC | flat INFO cutoff as a hidden MAF filter | filter MAF-stratified; state the MAF floor |
| Uniformly high "accuracy" for rare variants | concordance metric | use masked dosage-r2 by MAF bin |
| Genome-wide-significant hits do not replicate | cases/controls imputed separately, or hardcalled | impute together; regress on dosages |
| `bcftools +fill-tags` did not add an Rsq | fill-tags computes AF/MAF/HWE, not imputation quality | the quality field comes from the imputer |
| Negative Minimac EmpR at a site | strand/allele flip | re-align strand to the panel -> reference-panels |

## References

- Marchini J, Howie B. 2010. Genotype imputation for genome-wide association studies. *Nat Rev Genet* 11:499-511.
- Howie BN, Donnelly P, Marchini J. 2009. A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. *PLoS Genet* 5:e1000529.
- Das S, Forer L, Schonherr S, et al. 2016. Next-generation genotype imputation service and methods. *Nat Genet* 48:1284-1287.
- Browning BL, Zhou Y, Browning SR. 2018. A one-penny imputed genome from next-generation reference panels. *Am J Hum Genet* 103:338-348.
- Ramnarine S, Zhang J, Chen LS, et al. 2015. When does choice of accuracy measure alter imputation accuracy assessments? *PLoS One* 10:e0137601.
- Magi R, Asimit JL, Day-Williams AG, Zeggini E, Morris AP. 2012. Genome-wide association analysis of imputed rare variants: application to seven common complex diseases. *Genet Epidemiol* 36:785-796.
- Hofmeister RJ, Ribeiro DM, Rubinacci S, Delaneau O. 2023. Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank. *Nat Genet* 55:1243-1249.
- Huang KC, Sun W, Wu Y, et al. 2014. Association studies with imputed variants using expectation-maximization likelihood-ratio tests. *PLoS One* 9:e110679.
- Willer CJ, Li Y, Abecasis GR. 2010. METAL: fast and efficient meta-analysis of genomewide association scans. *Bioinformatics* 26:2190-2191.

## Related Skills

- genotype-imputation - Produces the dosages and the quality field this skill filters
- haplotype-phasing - Switch-error benchmarking of the phasing that precedes imputation
- reference-panels - Panel ancestry mismatch, which the quality metric cannot detect
- variant-calling/vcf-statistics - Generic VCF INFO/FORMAT field parsing
- population-genetics/association-testing - Consumes the filtered dosages
- long-read-sequencing/haplotype-phasing - Read-backed phasing switch QC
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate
