# drugZ Chemogenomic - Usage Guide

## Overview

Decision-grade analysis of CRISPR drug-modifier screens with drugZ (Li & Hart 2019 *Genome Medicine*). Uses bidirectional Z-scores comparing drug vs vehicle (NOT Day-0) to identify sensitizing (synthetic-lethal) genes and resistance-conferring (suppressor) genes. 2-3x more sensitive than MAGeCK / STARS / edgeR / RIGER on drug screens. Covers vehicle-vs-Day-0 reference choice, sumZ / normZ math, per-direction FDR, and reconciliation with MAGeCK MLE.

## Prerequisites

```bash
# drugZ (Python)
git clone https://github.com/hart-lab/drugz
cd drugz && pip install -e .
# OR (if available)
pip install drugz

# Helpers
pip install pandas numpy scipy statsmodels matplotlib
```

Required inputs:
- Count matrix with `sgRNA`, `GENE`, then sample columns (tab-separated)
- Vehicle (DMSO or carrier) sample columns
- Drug-treated sample columns
- Same time point for both arms (e.g., Day 14 vehicle vs Day 14 drug)

## Quick Start

Tell the AI agent what to do:
- "Run drugZ on PARPi screen: vehicle vs olaparib at 14 days; identify sensitizers (BRCA1/2) and suppressors (PARP1 paradox)"
- "Compare drugZ vs MAGeCK MLE on the same chemogenomic screen; expect drugZ to call 2-3x more hits"
- "Drug screen at 3 doses: run drugZ per dose; require consistent direction across doses for high-confidence hits"
- "Diagnose drugZ calling many essentials as sensitizers -- recommend `--remove-genes-file` with CEGv2"

## Example Prompts

### Standard Drug Screen

> "Run drugZ on counts.txt with vehicle samples Veh_r1, Veh_r2, Veh_r3 (passed via -c flag) and drug samples Drug_r1, Drug_r2, Drug_r3 (via -x flag). Pseudocount 5 (-p 5). Output ranked by fdr_synth (sensitizers) and fdr_supp (suppressors)."

> "Run drugZ on a PARPi (olaparib) chemogenomic screen. Expect sensitizers in DDR pathway (BRCA1, BRCA2, RAD51, FANCD2). Identify novel resistance genes (suppressors)."

### Reference Choice

> "Diagnose why my drugZ output has no sensitizers. Confirm: am I comparing drug vs vehicle, or drug vs Day 0? If Day-0, re-run with vehicle samples."

### Dose Response

> "For my 3-dose PARPi screen, run drugZ at low, mid, and high dose vs vehicle. Identify dose-consistent sensitizers (same gene appearing at all 3 doses with same direction)."

### Removing Essentials from Null

> "Run drugZ with `-r CEGv2.txt` containing CEGv2 essentials. Compare to default output; the difference is the artifact from essentiality inflating the null."

### Method Comparison

> "Compare drugZ vs MAGeCK MLE on the same drug screen. Compute Jaccard similarity at FDR 0.05. For drugZ-only hits, investigate per-sgRNA evidence."

> "MAGeCK and drugZ disagree on PARPi resistance hits. drugZ shows 12 hits at fdr_supp <0.05; MAGeCK shows 4 at pos|fdr <0.05. Reconcile."

### Diagnostics

> "drugZ hits are dominated by RPS / RPL / EIF essentials. Diagnose: am I running with Day-0 baseline by mistake, or do I need to exclude essentials from the null?"

> "Replicate-to-replicate drugZ runs give different top hits. Diagnose: insufficient sgRNAs per gene, or biological variation between replicates?"

## What the Agent Will Do

1. Verify experimental design: vehicle and drug arms at matched timepoint
2. Confirm replicate count: 3+ each arm
3. Run drugZ via CLI with vehicle as control, drug as treatment
4. Apply pseudocount (default 5; increase for low-count screens)
5. Optionally exclude essentials via `--remove-genes-file`
6. Output ranked hits: sensitizers (fdr_synth) and suppressors (fdr_supp)
7. For multi-dose: per-dose drugZ then consistency check
8. Cross-check against MAGeCK MLE on same data (consensus tier)
9. Annotate drug target as expected suppressor; report novel resistance genes
10. Report: sensitizer list, suppressor list, dose-consistency, comparison with MAGeCK

## Tips

- The single most common silent failure: comparing drug to Day-0 instead of drug to vehicle. Always use vehicle as control. Day-0 comparison conflates drug effect with normal-culture proliferation.
- drugZ is bidirectional. The sensitizer (synth) and suppressor (supp) lists are independent. Drug target genes often appear as suppressors (loss of drug target = resistance); this is correct biology.
- drugZ's 2-3x sensitivity advantage over MAGeCK comes from parametric Z-scoring of small effects; this is the right tool for low-effect chemogenomic screens. For high-effect screens (essentiality), MAGeCK and BAGEL2 work equally well.
- For multi-dose, multi-condition, or multi-cell-line screens, drugZ does not natively handle the additional dimensions. Run drugZ per-condition and require consistency.
- The `-r` option excludes specified genes from the null distribution. Use it to remove CEGv2 essentials when they dominate the sensitizer list -- this is a common pitfall in dropout-heavy drug screens.
- For consensus hit calling, run BOTH drugZ and MAGeCK MLE on the same data. Hits in both are high confidence; drugZ-only at low LFC need orthogonal validation.
- Drug-target genes often appear as suppressors -- e.g., PARP1 in PARPi screen, MEK in MEKi screen. This is expected biology, not a bug.
- For low-quality screens (replicate Pearson <0.85), drugZ's small-effect sensitivity becomes a liability -- false positives appear. Improve replicate concordance first.

## Decision Cheat Sheet

| Screen design | Method |
|---------------|--------|
| Vehicle vs drug, single dose, 3+ replicates | drugZ |
| Multi-dose response | drugZ per dose + consistency |
| Time course + drug | MAGeCK MLE with time + dose covariates |
| Drug × cell-line panel | drugZ per line + meta-analysis |
| Synergy / antagonism | MAGeCK MLE with interaction term |
| Small-effect screen | drugZ (preferred for chemogenomic) |
| Essentiality + drug | drugZ for drug, BAGEL2 for essentiality |

## Thresholds

| Threshold | Value | Rationale |
|-----------|-------|-----------|
| Sensitizer FDR | `fdr_synth < 0.05` | Li & Hart 2019 |
| Suppressor FDR | `fdr_supp < 0.05` | Same |
| High-confidence sensitizer | `fdr_synth < 0.01 AND normZ < -3` | Conservative |
| Replicate Pearson | >0.85 | Pre-drugZ QC |
| Min sgRNAs/gene | 4-6 (Brunello-style) | Stable Z |
| Pseudocount | 5 (default); higher for low-count | Li & Hart 2019 |

## Validation Checklist

- [ ] Vehicle (NOT Day-0) as control
- [ ] 3+ replicates each arm
- [ ] Library coverage validated upstream
- [ ] Replicate Pearson >0.85 within arms
- [ ] Pre-drugZ MAGeCK count run for QC
- [ ] CEGv2 essentials excluded if dominating sensitizer list
- [ ] Hits cross-validated with MAGeCK MLE
- [ ] Drug-target gene annotated as expected suppressor (if in list)
- [ ] Dose-consistency verified (if multi-dose)

## Related Skills

- crispr-screens/mageck-analysis - MAGeCK MLE alternative
- crispr-screens/bagel-essentiality - Tumor-suppressor-sensitive alternative
- crispr-screens/hit-calling - Cross-method decision tree
- crispr-screens/screen-qc - Pre-drugZ replicate concordance
- crispr-screens/library-design - 6+ sgRNAs/gene library
- crispr-screens/copy-number-correction - Pre-correction for cancer-line drug screens
- crispr-screens/base-editing-analysis - Variant-function drug screens
- pathway-analysis/go-enrichment - Functional analysis of drug-modifier hits
