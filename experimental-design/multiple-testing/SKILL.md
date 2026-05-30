---
name: bio-experimental-design-multiple-testing
description: Controls error rates across thousands of simultaneous tests in genomics discovery using false-discovery-rate methods (Benjamini-Hochberg 1995; Benjamini-Yekutieli 2001 for arbitrary dependence; Storey q-value with pi0 estimation; local FDR; independent filtering Bourgon 2010; covariate-weighted FDR via IHW Ignatiadis 2016), plus family-wise error control (Bonferroni, Holm) and the GWAS genome-wide threshold. Covers the FDR-versus-FWER choice as the discovery-versus-confirmatory distinction, the dependence assumptions behind BH (PRDS) versus BY, pi0 estimation, the independent-filtering and false-coverage-rate traps, and reproducibility ranking via IDR (Li 2011). Use when correcting p-values from genome-wide tests, choosing between BH/BY/q-value/Bonferroni, setting an FDR threshold, applying IHW or independent filtering, or interpreting q-values. For confirmatory trials with few pre-specified endpoints (closed testing, graphical/gatekeeping), see clinical-biostatistics/multiplicity-graphical.
tool_type: mixed
primary_tool: qvalue
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: qvalue 2.34+, IHW 1.30+, R stats (base) p.adjust, statsmodels 0.14+, scipy 1.12+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws an error, introspect the installed package and adapt to the actual API. Note: `statsmodels.stats.multitest.multipletests` defaults to `method='hs'` (Holm-Sidak, an FWER method), NOT Benjamini-Hochberg — always pass `method='fdr_bh'`/`'fdr_by'`/`'bonferroni'`/`'holm'` explicitly.

# Multiple Testing Correction

**"Correct p-values for testing thousands of features"** -> Choose an error rate appropriate to the regime (FDR for discovery, FWER for confirmatory), apply a procedure whose dependence assumptions match the data, and report the adjusted quantity with its interpretation.
- R: `p.adjust(p, method = 'BH')`, `qvalue::qvalue()`, `IHW::ihw()`
- Python: `statsmodels.stats.multitest.multipletests(p, method='fdr_bh')`

## The Single Most Important Modern Insight -- FDR vs FWER Is a Choice About Which Error Matters

The choice between false-discovery-rate and family-wise-error control is not a technicality; it is a statement about which kind of mistake is costly. In **discovery** (20,000 genes, thousands of peaks), tolerating a small, controlled fraction of false positives among the rejections buys enormous power — FDR is the right currency, and Bonferroni would discard nearly every true effect. In **confirmatory** work (a handful of pre-specified endpoints), a single false positive is unacceptable and FWER/closed testing is the standard (that regime lives in clinical-biostatistics/multiplicity-graphical). Two further levers buy back power that plain BH leaves on the table: estimating **pi0** (the proportion of true nulls) turns BH into the more powerful **q-value** (Storey 2002 *J R Stat Soc B* 64:479; Storey & Tibshirani 2003 *PNAS* 100:9440), and weighting hypotheses by an **independent informative covariate** recovers power via **IHW** (Ignatiadis 2016 *Nat Methods* 13:577). The dependence structure matters: BH controls FDR under independence or positive regression dependence (PRDS); under arbitrary or negative dependence use **BY** (Benjamini & Yekutieli 2001 *Ann Stat* 29:1165).

## Algorithmic Taxonomy

| Method | Controls | Dependence assumption | When to use | Tool |
|--------|----------|------------------------|-------------|------|
| Bonferroni | FWER | any | tiny families; confirmatory | `p.adjust(method='bonferroni')` |
| Holm | FWER | any | uniformly beats Bonferroni | `p.adjust(method='holm')` |
| Hochberg / Hommel | FWER | positive dependence | step-up FWER, more power | `p.adjust(method='hochberg'/'hommel')` |
| Benjamini-Hochberg | FDR | independence / PRDS | genome-wide discovery default | `p.adjust(method='BH')` |
| Benjamini-Yekutieli | FDR | arbitrary (incl. negative) | unknown/negative dependence | `p.adjust(method='BY')` |
| Storey q-value | pFDR | independence / weak dependence | many true positives (pi0 << 1) | `qvalue::qvalue` |
| Local FDR | posterior null prob | two-groups model | per-feature null probability | `qvalue` ($lfdr); `locfdr` |
| IHW | FDR | covariate independent of null p | informative covariate available | `IHW::ihw` |
| IDR | reproducibility | replicate ranks | thresholding by replicate consistency | `idr` (ENCODE) |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Genome-wide DE / peaks, discovery | BH or q-value at FDR 0.05 | controlled false-positive fraction; high power |
| Many true positives expected | q-value (estimates pi0) | more powerful than BH when pi0 << 1 |
| Strong/unknown/negative dependence | BY | BH guarantee needs PRDS |
| Informative covariate (mean expr, peak width) | IHW | data-driven weights recover power |
| Per-feature "is this one real?" | local FDR | posterior null probability, not tail average |
| Reporting CIs only on significant hits | FCR-adjusted intervals | naive selected CIs under-cover |
| Small confirmatory gene panel | Bonferroni/Holm | FWER appropriate; power loss acceptable |
| GWAS | genome-wide threshold ~5e-8 | ~1M effective independent tests |
| Confirmatory trial, few endpoints | -> clinical-biostatistics/multiplicity-graphical | closed testing / gatekeeping |
| Applying padj to a finished DE table | -> differential-expression/de-results | method choice here; application there |

## FDR -- Benjamini-Hochberg and the q-value

```r
# Benjamini-Hochberg adjusted p-values (the genome-wide default)
padj <- p.adjust(pvalues, method = 'BH')
sum(padj < 0.05)                                  # discoveries at FDR 5%

# Storey q-value: estimates pi0 (fraction of true nulls) for more power when pi0 << 1
library(qvalue)
qobj <- qvalue(pvalues)
qobj$pi0                                           # estimated proportion of true nulls
q   <- qobj$qvalues                                # min FDR at which each feature is called
lfdr <- qobj$lfdr                                  # local FDR: posterior P(null | statistic)
```

## Dependence -- When BH Is Not Enough (BY)

```r
# BH controls FDR under independence or positive regression dependence (PRDS).
# Under arbitrary or negative dependence, use Benjamini-Yekutieli (more conservative).
padj_by <- p.adjust(pvalues, method = 'BY')        # valid under any dependence structure
```

## Covariate-Weighted FDR -- IHW

```r
# Weight hypotheses by an INDEPENDENT informative covariate (e.g. mean expression),
# which must be independent of the p-value under the null. Recovers power vs plain BH.
library(IHW)
res <- ihw(pvalue ~ mean_expression, data = de_table, alpha = 0.05)
de_table$padj_ihw <- adj_pvalues(res)
rejections(res)
```

## Independent Filtering -- Power for Free, If the Filter Is Independent

Filtering out features before testing increases power **only if** the filter statistic is independent of the test statistic under the null (Bourgon, Gentleman & Huber 2010 *PNAS* 107:9546). Overall mean count is independent and is why DESeq2 filters low-count genes automatically; a pre-test on variance or a preliminary t-test is **not** independent and biases the FDR. The DE filtering itself is executed in differential-expression; this skill governs whether a proposed filter is legitimate.

## Python Equivalent (mind the default)

```python
from statsmodels.stats.multitest import multipletests
# DEFAULT method is 'hs' (Holm-Sidak, FWER) -- ALWAYS pass method explicitly.
rej, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')   # Benjamini-Hochberg
rej_by, padj_by, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_by')  # BY
```

## GWAS and the Family-Definition Problem

The genome-wide significance threshold of ~5e-8 is a Bonferroni-style bound for roughly one million effectively independent common-variant tests; Dudbridge & Gusnanto 2008 (*Genet Epidemiol* 32:227) derived ~7.2e-8 for European-ancestry data, near the standard 5e-8. The GWAS test machinery lives in population-genetics/association-testing. More broadly, **what counts as "the family"** of tests is an analyst decision and part of the garden of forking paths: correcting within one contrast, across all contrasts, or across a whole paper are different alpha budgets. Pre-specify the family before seeing results.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| q-value finds many more hits than BH | pi0 << 1 (many true positives) | q-value legitimately more powerful; report pi0 |
| BY far more conservative than BH | strong/negative dependence penalty | if dependence is positive, BH is justified; state the assumption |
| IHW and BH differ substantially | informative, null-independent covariate | IHW gain is real if independence holds; verify the covariate |
| Filtering changed the hit count | filter not independent of the test statistic | use a null-independent filter (mean count), not variance/preliminary test |
| Per-feature local FDR high but BH q low | tail-average vs per-feature interpretation | report both; local FDR answers "is THIS one real?" |

## Per-Method Failure Modes

### Bonferroni on a transcriptome
- **Trigger:** Bonferroni across 20,000 genes in a discovery study.
- **Mechanism:** FWER control is far too strict for discovery.
- **Symptom:** almost nothing significant; true effects discarded.
- **Fix:** BH or q-value at a target FDR.

### BH under arbitrary/negative dependence
- **Trigger:** BH on strongly/negatively correlated statistics.
- **Mechanism:** BH guarantee requires independence or PRDS (Benjamini-Yekutieli 2001).
- **Symptom:** realized FDR exceeds nominal.
- **Fix:** BY when dependence is unknown or negative.

### statsmodels default is not BH
- **Trigger:** `multipletests(p)` expecting Benjamini-Hochberg.
- **Mechanism:** default `method='hs'` (Holm-Sidak, FWER).
- **Symptom:** far fewer significant calls than expected.
- **Fix:** pass `method='fdr_bh'` explicitly.

### Non-independent filtering
- **Trigger:** filter on variance or a preliminary test before the main test.
- **Mechanism:** filter statistic correlated with the test statistic under the null (Bourgon 2010).
- **Symptom:** anti-conservative FDR.
- **Fix:** filter only on a null-independent statistic (overall mean count).

### Selected CIs without FCR adjustment
- **Trigger:** reporting unadjusted CIs only for significant features.
- **Mechanism:** selection induces under-coverage (false coverage rate).
- **Symptom:** intervals too narrow; replication misses.
- **Fix:** FCR-adjusted intervals for the selected set.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| FDR < 0.05 discovery default | Benjamini-Hochberg 1995 *JRSS-B* 57:289 | 5% of calls expected false |
| FDR < 0.10 exploratory | common practice | more leads at higher false fraction |
| q-value uses estimated pi0 | Storey 2002 *JRSS-B* 64:479 | power gain when pi0 << 1 |
| BH valid under independence/PRDS; else BY | Benjamini-Yekutieli 2001 *Ann Stat* 29:1165 | dependence governs validity |
| GWAS ~5e-8 (7.2e-8 derived) | Dudbridge-Gusnanto 2008 *Genet Epidemiol* 32:227 | ~1M effective tests |
| Filter must be null-independent | Bourgon 2010 *PNAS* 107:9546 | otherwise FDR is biased |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Almost nothing significant genome-wide | Bonferroni in a discovery study | BH or q-value |
| Realized FDR exceeds nominal | BH under negative dependence | BY |
| Far fewer hits than expected in Python | statsmodels default 'hs' | `method='fdr_bh'` |
| FDR biased after pre-filtering | non-independent filter | filter on mean count only |
| Replication misses "significant" effects | unadjusted selected CIs | FCR-adjusted intervals |

## References

- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B* 57:289-300.
- Benjamini Y, Yekutieli D. 2001. The control of the false discovery rate in multiple testing under dependency. *Ann Stat* 29:1165-1188.
- Storey JD. 2002. A direct approach to false discovery rates. *J R Stat Soc B* 64:479-498.
- Storey JD, Tibshirani R. 2003. Statistical significance for genomewide studies. *PNAS* 100:9440-9445.
- Efron B. 2008. Microarrays, empirical Bayes and the two-groups model. *Stat Sci* 23:1-22.
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107:9546-9551.
- Ignatiadis N, Klaus B, Zaugg JB, Huber W. 2016. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods* 13:577-580.
- Li Q, Brown JB, Huang H, Bickel PJ. 2011. Measuring reproducibility of high-throughput experiments. *Ann Appl Stat* 5:1752-1779.
- Dudbridge F, Gusnanto A. 2008. Estimation of significance thresholds for genomewide association scans. *Genet Epidemiol* 32:227-234.

## Related Skills

- power-analysis - The FDR target feeds the power/EDR calculation
- sample-size - Replicate number depends on the FDR threshold chosen here
- batch-design - Surrogate variables change the effective number of tests
- differential-expression/de-results - Where the padj column is applied to a DE table
- population-genetics/association-testing - GWAS genome-wide significance machinery
- pathway-analysis/go-enrichment - Correcting enrichment p-values
- clinical-biostatistics/multiplicity-graphical - Confirmatory FWER / closed testing for trials with few endpoints
