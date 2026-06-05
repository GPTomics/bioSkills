---
name: bio-analytical-validation
description: Treats a ctDNA assay as a molecule-counting experiment at the Poisson edge and builds its analytical-validation case the measurement-science way. Covers the genome-equivalent currency (~330 haploid copies/ng), the lambda = input_GE x VAF sampling ceiling (lambda>=3 for ~95% detection), the error-suppression ladder (raw NGS ~1e-3 -> single-strand UMI ~1e-4/1e-5 -> duplex <1e-7), the CLSI EP17 LoB/LoD/LoD95/LoQ framework, the per-locus-vs-panel-integrated LoD distinction that lets bespoke MRD reach ppm, contrived/SEQC2 reference standards, and honest LoD reporting conditioned on input mass + consensus depth + replicate detection rate. Use when stating or trusting a sensitivity claim, designing a dilution-series validation, deciding how many genome equivalents are needed at a target VAF, choosing a single-locus vs panel-integrated LoD, or auditing a "detects 0.1% VAF" claim.
tool_type: python
primary_tool: scipy
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Analytical Validation and Detection Limits

**"What is the real limit of detection of my ctDNA assay, and can I trust the number I am about to report?"** -> Quantify the Poisson sampling ceiling, the error-suppression floor, and the LoB/LoD/LoQ that together define a defensible sensitivity claim.
- Python: `scipy.stats.poisson` for detection-probability math, `scipy.stats.norm` for CLSI LoB/LoD, `statsmodels` Probit/Logit for a dilution-series LoD95 fit.

## The Single Most Important Modern Insight -- LoD Is Set by Genome Equivalents Sampled and Error Suppression, NOT by Sequencing Depth or the Caller

A ctDNA assay is a molecule-counting experiment at the Poisson edge. The mutant signal is a fixed, tiny number of physical template molecules in the tube, and the limit of detection is governed by two ceilings: how many genome equivalents were sampled (Poisson), and how low the background error floor was driven (error suppression). 1 ng of human DNA is ~330 haploid genome equivalents; the expected mutant-molecule count is lambda = input_GE x VAF. A 0.1% variant on 1,000 GE (~3.0 ng) has lambda = 1, so e^-1 ~= 37% of the time the mutant template was never in the tube and a perfect sequencer detects nothing. Past the point where every input molecule has been read once (sampling saturation, visible as a deduplication plateau in UMI families), additional read depth re-sequences the same physical molecules and adds zero information. Reporting an LoD as a bare VAF -- with no input mass, no unique-molecule (consensus) depth, no replicate detection rate -- is reporting an undefined quantity.

The second ceiling is the per-base background error rate, which sets the VAF floor independently: a 0.1% variant cannot be distinguished from noise if the assay manufactures that base at 0.1%. Error suppression is a ladder (raw NGS ~1e-3 -> single-strand UMI consensus ~1e-4/1e-5 -> duplex <1e-7), and single-strand consensus does NOT remove template-resident damage (C->T deamination, G->T 8-oxoG) because every PCR copy of that strand inherits the lesion -- only duplex strand-concordance catches it. The achieved LoD is the *worse* of the two ceilings: error dominates above ~0.1% VAF for tumor-naive single-locus calling, sampling dominates below it. The escape hatch is integration -- a bespoke panel summing mutant molecules across 16-50 loci against summed background reaches single-ppm even though each locus alone is ~1e-3 to 1e-4 (per-locus vs panel-integrated LoD).

## Methods Landscape

| Concept | Definition | Source |
|---------|------------|--------|
| LoB (Limit of Blank) | Highest signal expected from an analyte-free blank (95th pct): LoB = mean_blank + 1.645*SD_blank; the false-positive anchor on true negatives | CLSI EP17-A2 |
| LoD (Limit of Detection) | Lowest level reliably distinguishable from LoB: LoD = LoB + 1.645*SD_low; a sample at LoD is detected ~95% of the time | CLSI EP17-A2 |
| LoD95 | The concentration/VAF where detection probability = 95%; a point on a probit/logistic detection curve, not a separate definition | CLSI EP17-A2; Newman 2016 |
| LoQ (Limit of Quantitation) | Lowest level measurable with stated precision (e.g. CV<=20%); LoQ >= LoD always, so a "VAF" near the floor is detectable but not trustworthy | CLSI EP17-A2 |
| Per-locus LoD | Single-variant LoD; sampling- and error-limited (~0.05-0.1% VAF typical) | Newman 2014/2016 |
| Panel-integrated LoD | Evidence summed across N tracked variants via a >=k-of-N positivity rule (binomial over per-locus Poisson detection), reaching single-ppm at 16-50 loci; ~sqrt(N) variance-averaging is only a loose lower bound | Reinert 2019 |
| Reference standards | Contrived defined-VAF cell-line admixtures fragmented to ~160 bp into normal cfDNA; SEQC2 Sample A / HCC1395 truth sets | Fang 2021 (SEQC2) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| "How many GE for 95% detection at VAF X?" | Solve lambda = input_GE x VAF >= 3, so input_GE >= 3/VAF | 1 - e^-3 = 0.95; ~30,000 GE (~91 ng at 330 GE/ng) for a single 1e-4 variant -- often more than one tube provides |
| "Why is more depth not helping?" | Report unique (consensus) molecular coverage, not raw depth; check the dedup plateau | Past sampling saturation, depth re-reads the same molecules; the ceiling is GE in the tube |
| Single hotspot vs bespoke panel for low VAF | Single locus: error/sampling-limited ~0.1%; need ppm -> integrate across 16-50 clonal loci | Per-locus Poisson/error floor is escaped only by summing independent detections (panel-integrated LoD) |
| Reporting an LoD | Condition on input mass (GE) + consensus depth + replicate detection rate (e.g. "LoD95 0.1% VAF at 30 ng / 2x duplex / 95% of 20 replicates") | A bare VAF omits the input mass, the unique depth, and per-locus vs integrated -- it is undefined |
| Estimating LoD95 from a dilution series | Probit (or logistic) regression of detection (0/1) on VAF; read off the 95% point with a CI | CLSI EP17-A2 detection-curve method; binary detection is a clean GLM target |
| Distinguishing detection from quantitation | Set LoD for yes/no calls; set LoQ (CV<=20%) separately for any reported VAF/TF | MRD calls are binary and can sit far below LoQ; a near-floor VAF number is not quantitative |
| Validating against truth | Contrived SEQC2 Sample A / HCC1395 admixtures, fragmented to cfDNA-like ~160 bp | Real low-VAF patient material is scarce/unverifiable; commutability with plasma is the caveat |

## Genome-Equivalent and Poisson Detection Calculator

**Goal:** Convert an input mass and target VAF into an expected mutant-molecule count and a detection probability, so a sensitivity claim is anchored to molecules rather than to a VAF alone.

**Approach:** Convert ng to haploid genome equivalents (~330/ng), set lambda = input_GE x VAF, and read the detection probability as a Poisson tail P(X >= k) = 1 - cdf(k-1, lambda); invert for the minimum GE that puts lambda at the >=3 sampling-detection threshold.

```python
import numpy as np
from scipy.stats import poisson

GE_PER_NG = 330  # haploid ~3.3 pg -> strict 1 ng / 3.3 pg = 303; 330 is the common diploid-6.6 pg/rounding convention

def genome_equivalents(input_ng):
    return input_ng * GE_PER_NG

def detection_probability(input_ng, vaf, min_mutant_molecules=1):
    '''P(at least min_mutant_molecules present) under Poisson(lambda = GE * VAF).'''
    lam = genome_equivalents(input_ng) * vaf
    return float(poisson.sf(min_mutant_molecules - 1, lam))

def ge_for_sampling_detection(vaf, target_lambda=3.0):
    '''GE needed so lambda >= 3 -> ~95% chance the mutant molecule is present at all.'''
    return target_lambda / vaf

# A 0.1% variant on 3.0 ng (~990 GE) has lambda ~= 1 -> ~63% detected, ~37% missed by sampling alone.
detection_probability(3.0, 0.001)          # ~0.63
ge_for_sampling_detection(1e-4)            # 30000 GE (~91 ng) for a single 0.01% variant
```

## LoB and LoD95 from a Dilution Series (CLSI EP17 style)

**Goal:** Estimate the VAF at which the assay detects 95% of the time, from a contrived dilution series, and anchor it to the blank-derived false-positive floor.

**Approach:** Compute LoB from blank replicates (mean + 1.645*SD, one-sided 95th pct), then fit a probit GLM of binary detection on log10(VAF) (CLSI EP17 fits on log concentration) across the dilution series and invert it for the 95% detection point. The series must bracket the 0.95 crossing — all-detected upper levels cause near-complete separation and an unstable slope.

```python
import numpy as np
import statsmodels.api as sm
from scipy.stats import norm

def limit_of_blank(blank_signals):
    '''LoB = mean + 1.645*SD; one-sided 95th percentile of analyte-free blanks.'''
    blank_signals = np.asarray(blank_signals, dtype=float)
    return blank_signals.mean() + 1.645 * blank_signals.std(ddof=1)

def lod95_probit(vaf_levels, detected):
    '''Probit fit of detection (0/1) on log10(VAF); returns the VAF where P(detect) = 0.95.'''
    log_vaf = np.log10(np.asarray(vaf_levels, dtype=float))
    y = np.asarray(detected, dtype=float)
    X = sm.add_constant(log_vaf)
    fit = sm.GLM(y, X, family=sm.families.Binomial(link=sm.families.links.Probit())).fit()
    intercept, slope = fit.params
    return 10 ** ((norm.ppf(0.95) - intercept) / slope)
```

## Per-Locus to Panel-Integrated LoD

**Goal:** Combine independent per-locus detection probabilities into the panel-level detection probability that a bespoke MRD assay actually achieves, and find the integrated LoD.

**Approach:** Treat each tracked locus as an independent Poisson sampler at the same tumor VAF; a panel positive call requires at least k loci detected, so the panel detection probability is the binomial-tail over the per-locus probabilities -- this is why summing 16-50 loci reaches ppm.

```python
import numpy as np
from scipy.stats import poisson, binom

def panel_detection_probability(input_ng, vaf, n_loci, min_loci_positive=2):
    '''P(>= min_loci_positive of n_loci detected); >=2-of-N is the Signatera-style positivity rule.'''
    per_locus = float(poisson.sf(0, input_ng * 330 * vaf))
    return float(binom.sf(min_loci_positive - 1, n_loci, per_locus))

def panel_integrated_lod95(input_ng, n_loci, min_loci_positive=2, grid=None):
    '''Lowest VAF on a log grid where the >=k-of-N panel call hits 95%.'''
    grid = np.logspace(-6, -2, 400) if grid is None else np.asarray(grid)
    probs = [panel_detection_probability(input_ng, v, n_loci, min_loci_positive) for v in grid]
    hits = grid[np.asarray(probs) >= 0.95]
    return float(hits.min()) if hits.size else float('nan')
```

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| ~330 haploid genome equivalents per ng cfDNA | Standard (haploid ~3.3 pg) | Converts input mass to the molecule count that actually sets sensitivity; strict 1 ng / 3.3 pg = 303, with 330 the common diploid-6.6 pg/rounding convention |
| lambda = input_GE x VAF; lambda >= 3 for ~95% sampling-detection | Poisson, 1 - e^-3 = 0.95 | Below lambda~3 the mutant template is often simply absent from the tube regardless of sequencing |
| Raw NGS error floor ~1e-3 | Schmitt 2012 context; field consensus | Sets the per-base VAF floor before any consensus; a global VAF cutoff above this is noise-limited |
| Single-strand UMI consensus ~1e-4 to 1e-5 | Newman 2014/2016 (CAPP-Seq/iDES) | Majority-vote within a UMI family erases PCR/sequencing error not shared across the family |
| Duplex sequencing <1e-7 (theory <1/1e9 nt) | Schmitt 2012 *PNAS* 109:14508 | Requires the variant on BOTH original strands; independent strand errors cannot agree |
| iDES adds ~3-15x over baseline; ctDNA to ~4e-5 | Newman 2016 *Nat Biotechnol* 34:547 | Position/trinucleotide background model subtracts stereotyped artifacts per locus |
| ichorCNA tumor-fraction floor ~3% | Adalsteinsson 2017 *Nat Commun* 8:1324 | Copy-number-based TF estimation; sWGS/ULP-WGS cannot resolve TF below ~3% -- an LoD, not a VAF |
| Bespoke panel reaches single-ppm by integrating 16-50 loci | Reinert 2019 *JAMA Oncol* 5:1124 | Per-locus ~1e-4 floor escaped by summing independent detections; >=2-of-N positivity rule |
| LoQ >= LoD (e.g. CV<=20% for quantitation) | CLSI EP17-A2 | Detection (binary) is easier than quantitation (continuous); near-floor VAFs are not trustworthy numbers |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| "Assay detects 0.1% VAF" with no input mass | VAF reported as a standalone sensitivity spec | Condition the LoD on input GE + consensus depth + replicate detection rate; 0.1% on 100 GE is noise |
| Buying more sequencing depth to improve sensitivity | Conflating read depth with molecule count | Past the dedup plateau the assay is sampling-saturated; add plasma volume / conversion efficiency, not depth |
| Per-locus LoD quoted as the panel LoD (or vice versa) | Ignoring integration across tracked loci | State which is reported; a 50-variant panel's integrated LoD is orders of magnitude below any single locus |
| VAF used as the sensitivity unit | Omitting the molecule count behind the fraction | Pair every VAF with input GE; lambda = GE x VAF is the quantity that determines detection |
| Single-strand UMI assumed to remove damage artifacts | Template-resident C->T/G->T inherited by every copy | Use duplex strand-concordance for sub-1e-5 claims; single-strand votes unanimously for the lesion |
| Reporting a near-floor VAF as a measured value | Confusing LoD (detect) with LoQ (quantify) | Quantitative VAF/TF only at/above LoQ (CV<=20%); below it report detected/not-detected |
| Global VAF cutoff across all loci | Background error is position/context-dependent | Use a per-locus background model (iDES-style); a flat threshold loses sensitivity and specificity |

## References

- Diehl F, Schmidt K, Choti MA, et al. 2008. Circulating mutant DNA to assess tumor dynamics. *Nat Med* 14:985-990. -- ctDNA half-life ~114 min; molecule-counting framing of tumor dynamics.
- Schmitt MW, Kennedy SR, Salk JJ, et al. 2012. Detection of ultra-rare mutations by next-generation sequencing. *PNAS* 109:14508-14513. -- Duplex sequencing; theoretical error floor <1 per 1e9 nt.
- Newman AM, Bratman SV, To J, et al. 2014. An ultrasensitive method for quantitating circulating tumor DNA with broad patient coverage. *Nat Med* 20:548-554. -- CAPP-Seq; UMI-consensus error suppression.
- Newman AM, Lovejoy AF, Klass DM, et al. 2016. Integrated digital error suppression for improved detection of circulating tumor DNA. *Nat Biotechnol* 34:547-555. -- iDES; ~3-15x gain; ctDNA to ~4e-5.
- Razavi P, Li BT, Brown DN, et al. 2019. High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants. *Nat Med* 25:1928-1937. -- CHIP as the dominant non-tumor signal in the LoB blank.
- Adalsteinsson VA, Ha G, Freeman SS, et al. 2017. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. *Nat Commun* 8:1324. -- ichorCNA; copy-number tumor-fraction floor ~3%.
- Reinert T, Henriksen TV, Christensen E, et al. 2019. Analysis of plasma cell-free DNA by ultradeep sequencing in patients with stages I to III colorectal cancer. *JAMA Oncol* 5:1124-1131. -- Signatera; 16-variant integration; >=2-of-N positivity.
- Fang LT, Zhu B, Zhao Y, et al.; SEQC2 Consortium. 2021. Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing. *Nat Biotechnol* 39:1151-1160. -- SEQC2 Sample A / HCC1395 contrived reference standards.
- CLSI EP17-A2. 2012. Evaluation of Detection Capability for Clinical Laboratory Measurement Procedures; Approved Guideline -- Second Edition. Clinical and Laboratory Standards Institute. -- Governing LoB/LoD/LoQ definitions.

## Related Skills

- ctdna-mutation-detection - applies these limits to low-VAF somatic calls
- longitudinal-monitoring - per-timepoint LoD and left-censoring of undetectable samples
- tumor-fraction-estimation - the ~3% CNA-based detection floor as an LoD
- experimental-design/multiple-testing - repeated-surveillance specificity and FDR
- clinical-biostatistics/power-and-sample-size - validation-study design
