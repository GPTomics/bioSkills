---
name: bio-immunoinformatics-immunogenicity-scoring
description: Rank and prioritize neoantigen/epitope candidates by likely T-cell response using NeoFox feature annotation, PRIME2.0, BigMHC-IM, the Łuksza/Balachandran fitness model (agretopicity + foreignness), and pVACtools tiering. Encodes the field's hard truths that immunogenicity is the least-solved layer (dedicated scores ~AUROC 0.6-0.7, modest PPV), that scores are valid only for RANKING within one patient (never absolute go/no-go or cross-patient), that DAI has anchor-inflation and WT-denominator traps, and that stacking weak correlated scores into one number is a red flag. Use when ordering a candidate list for a vaccine. Binding lives in mhc-binding-prediction; calling in neoantigen-prediction.
tool_type: python
primary_tool: NeoFox
---

## Version Compatibility

Reference examples tested with: NeoFox 1.0+, pVACtools 4.1+, pandas 2.2+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: NeoFox annotates ~16 published features (it does not rank candidates automatically); PRIME 2.x requires MixMHCpred v3.0+ on PATH; BigMHC has separate `-m el` and `-m im` heads. ImmunoBERT is a PRESENTATION model, not an immunogenicity predictor — do not use it here. PRIME2.0 is the *Cell Systems* 2023 paper (the *Cell Reports Medicine* 2021 paper is PRIME v1). Re-verify tool versions and the supported-allele lists before scoring.

# Immunogenicity Scoring

**"Rank my neoantigen candidates by how likely a T cell responds"** -> Annotate presentation + recognition features and order candidates within a patient; never assign an absolute immunogenicity verdict.
- Python: `NeoFox` to compute the published feature panel; `PRIME` / `BigMHC -m im` for recognition scores
- CLI: pVACtools aggregate-report tiering as the auditable, rule-based default ranking

## The Single Most Important Modern Insight -- this is the least-solved layer; rank within a patient, never threshold

Binding/presentation is genuinely good (AUROC high-0.9s); immunogenicity is not close. Predicting whether a displayed peptide provokes a T-cell response requires knowing whether a cognate TCR exists in this patient's repertoire, whether that clone survived thymic negative selection (escaped tolerance), and whether it activates in a suppressive tumor microenvironment — none observable from sequence. Dedicated immunogenicity tools land around AUROC 0.6-0.7 on their own test sets and worse on independent data; in TESLA the dedicated in-silico immunogenicity scores correlated poorly with validated immunogenicity, while presentation strength, binding stability, abundance/expression, agretopicity, and foreignness carried the signal. Two operational rules follow. First, immunogenicity scores are calibrated within a context (a tool, an allele, often a patient's HLA), so they are legitimate for ordering one patient's candidate list and illegitimate for absolute go/no-go or cross-patient/cross-allele comparison. Second, a confident single composite number is a red flag: stacking weak, correlated, IEDB-bias-trained scores into one value launders the bias at higher apparent precision. The honest deliverable is an ordered, feature-annotated shortlist with its uncertainty stated out loud.

## Why "Best Binder" Lost to "Best Quality"

The best-binder heuristic fails on a tolerance argument: a peptide that binds MHC superbly but closely resembles a self-peptide the thymus presented has had its cognate T cells deleted, so display does not help. A moderate binder that looks strikingly un-self may have a full, un-tolerized repertoire. The modern requirement is conjunctive — a useful neoantigen must be both PRESENTED (binding) AND FOREIGN enough (different from self) to have escaped tolerance. The Łuksza/Balachandran fitness model formalizes this: quality = amplitude (how much better the mutant is presented than its WT, a DAI-like term) x recognition potential R (resemblance to known immunogenic foreign epitopes). This is why agretopicity and foreignness, not raw affinity, recur in every validated analysis.

## Tool Taxonomy

| Tool | Citation | What it scores | Note |
|------|----------|----------------|------|
| NeoFox | Lang 2021 | ~16 features at once (DAI, foreignness, dissimilarity, PRIME, PHBR, ...) | Annotates, does NOT rank — the right division of labor |
| pVACtools tiering | Hundal 2020 | Rule-based tiers + within-tier sort | Auditable default; quarantines anchor/subclonal traps |
| PRIME2.0 | Gfeller 2023 | Class I immunogenicity (presentation x TCR-recognition) | Strong; needs MixMHCpred v3.0+ |
| BigMHC-IM | Albert 2023 | Class I immunogenicity (transfer-learned) | High precision; pan-allelic |
| IEDB immunogenicity | Calis 2013 | Class I (AA + position) | Weak, allele-pooled, no self-comparison; one feature only |
| DeepImmuno | Li 2021 | Class I CNN | 9/10mer only; limited alleles |
| fitness model (foreignness) | Łuksza 2017; Balachandran 2017 | Quality = amplitude x recognition | The conceptual backbone |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Default: rank a patient's candidates | NeoFox features -> pVACtools tiering -> human curation | Transparent features + auditable tiers, not a black-box score |
| Need a single recognition score | PRIME2.0 or BigMHC-IM | Best-validated class I; report alongside features, not alone |
| "Is this one immunogenic, yes/no?" | Reframe to ranking | No honest tool gives an absolute verdict |
| CD4 / class II immunogenicity | Flag as a frontier (TLimmuno2 etc.) | Class II immunogenicity is even less solved |
| Final shortlist for synthesis | Feature-annotated table + expression/clonality filters | Presentation + abundance carry most real signal (TESLA) |

## Annotate Features, Then Rank Within Patient

**Goal:** Order one patient's candidates without collapsing fragile features into a single over-trusted number.

**Approach:** Compute the feature panel (NeoFox), apply the non-negotiable expression/clonality filters first, then sort by presentation + abundance + quality features, keeping the features visible side by side for human curation. Cross-patient comparison is invalid.

```python
import pandas as pd

def rank_within_patient(df, expr_col='gene_expression', vaf_col='rna_vaf'):
    '''Filter (not score) on expression/clonality first, then order by presentation,
    abundance, and quality. Returns a feature-annotated table for human curation, not
    a verdict. Scores are within-patient only - never compare across patients/alleles.'''
    keep = df[(df[expr_col] >= 1.0) & (df[vaf_col] >= 0.25)].copy()
    sort_cols = ['presentation_rank', 'gene_expression', 'agretopicity', 'foreignness']
    ascending = [True, False, False, False]
    cols = [c for c in sort_cols if c in keep.columns]
    asc = [a for c, a in zip(sort_cols, ascending) if c in keep.columns]
    return keep.sort_values(cols, ascending=asc)
```

## Compute Agretopicity (DAI) Defensively

**Goal:** Use the mutant-vs-WT binding gain without falling into its two traps.

**Approach:** Agretopicity (ratio, IC50_WT / IC50_MT; the DAI family — Duan 2014 uses the difference form) rewards a mutant that binds while WT does not. Trap 1: an anchor-position mutation inflates it without changing the TCR-facing surface (quarantine via the Anchor tier). Trap 2: when WT binds very poorly, the denominator explodes and the ratio is dominated by prediction noise — a value of 200 on a barely-estimable WT is not 100x more meaningful than a value of 2.

```python
def defensive_dai(df, wt='wt_ic50', mt='mt_ic50', anchor='mutation_at_anchor', wt_cap=5000):
    '''Flag anchor-inflated and denominator-unstable DAI rather than trusting the number.'''
    out = df.copy()
    out['dai'] = out[wt] / out[mt]
    out['dai_anchor_artifact'] = out[anchor]                 # surface unchanged -> DAI is artifact
    out['dai_unstable'] = out[wt] > wt_cap                   # WT barely presented -> ratio is noise
    out['dai_trustworthy'] = ~out['dai_anchor_artifact'] & ~out['dai_unstable']
    return out
```

## Per-Method Failure Modes

### Treating a score as a verdict
**Trigger:** "score > X means immunogenic" or comparing scores across patients. **Mechanism:** scores are calibrated within tool/allele/patient. **Symptom:** false confidence; cross-patient mis-ranking. **Fix:** rank within a patient; state uncertainty; never threshold absolutely.

### The composite-score illusion
**Trigger:** summing/modeling DAI + foreignness + dissimilarity + hydrophobicity + PRIME into one number. **Mechanism:** components are weak, correlated (several measure "un-selfness"), and trained on ill-defined negatives. **Symptom:** an authoritative-looking 3-decimal number hiding fragile assumptions. **Fix:** keep features side by side; use auditable tiers; let a human weigh axes.

### DAI anchor inflation / denominator instability
**Trigger:** trusting a high DAI. **Mechanism:** anchor mutation changes binding not TCR surface; tiny WT binding blows up the ratio. **Symptom:** top-ranked candidates that are anchor artifacts or noise. **Fix:** inspect mutation position and actual WT binding; quarantine via Anchor tier.

### Negative-set blindness
**Trigger:** trusting a new tool's headline AUROC. **Mechanism:** IEDB "negatives" conflate proven-non-immunogenic with untested; redrawing realistic negatives collapses performance. **Symptom:** great benchmark, poor real-world PPV. **Fix:** ask how negatives were defined before reading the number.

### CD4/class II blind spot
**Trigger:** optimizing a vaccine purely on class I immunogenicity. **Mechanism:** CD4 help drives durable efficacy but class II immunogenicity is a frontier. **Symptom:** optimizing the better-measured half of a two-armed problem. **Fix:** flag class II as unproven; include CD4 epitopes via mhc-class-ii-prediction.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Dedicated immunogenicity AUROC ~0.6-0.7 | TESLA; tool benchmarks | The honest performance ceiling; weak prior, not verdict |
| Gene TPM >= 1, RNA VAF >= 0.25 | pVACtools defaults | Unexpressed/low-VAF peptides are not displayed (filter first) |
| Subclonal at DNA VAF <= purity/4 | pVACtools | Clonal targets beat subclonal (McGranahan 2016) |
| Presentation + abundance carry the signal | Wells 2020 (TESLA) | Most predictive power is upstream of recognition scores |
| Rank within patient only | Score calibration | Cross-patient/allele comparison is invalid |
| Agretopicity ratio (amplitude); DAI difference | Łuksza 2017; Duan 2014 | Inspect position + WT binding; anchor inflation and denominator instability |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Absolute "immunogenic: yes/no" claim | Thresholded a within-context score | Reframe as within-patient ranking |
| Over-trusted single composite | Stacked weak correlated scores | Keep features visible; audit with tiers |
| High-DAI artifacts at top | Anchor mutation / unstable WT denominator | Defensive DAI; Anchor tier |
| Used ImmunoBERT as immunogenicity | It is a presentation model | Use PRIME/BigMHC-IM/Calis for recognition |
| Great AUROC, poor validation | Ill-defined negative set | Interrogate negatives; demand functional validation |
| Class II candidates over-trusted | CD4 immunogenicity is a frontier | Flag uncertainty; treat as unproven |

## References

- Wells DK, van Buuren MM, Dang KK, et al. 2020. Key parameters of tumor epitope immunogenicity revealed through a consortium approach improve neoantigen prediction (TESLA). *Cell* 183(3):818-834.
- Łuksza M, Riaz N, Makarov V, et al. 2017. A neoantigen fitness model predicts tumour response to checkpoint blockade immunotherapy. *Nature* 551:517-520.
- Balachandran VP, Łuksza M, Zhao JN, et al. 2017. Identification of unique neoantigen qualities in long-term survivors of pancreatic cancer. *Nature* 551:512-516.
- Calis JJA, Maybeno M, Greenbaum JA, et al. 2013. Properties of MHC class I presented peptides that enhance immunogenicity. *PLoS Computational Biology* 9(10):e1003266.
- Schmidt J, Smith AR, Magnin M, et al. 2021. Prediction of neo-epitope immunogenicity reveals TCR recognition determinants (PRIME). *Cell Reports Medicine* 2(2):100194.
- Gfeller D, Schmidt J, Croce G, et al. 2023. Improved predictions of antigen presentation and TCR recognition with MixMHCpred2.2 and PRIME2.0. *Cell Systems* 14(1):72-83.
- Albert BA, Yang Y, Shao XM, et al. 2023. Deep neural networks predict class I MHC epitope presentation and transfer learn neoepitope immunogenicity (BigMHC). *Nature Machine Intelligence* 5(8):861-872.
- Duan F, Duitama J, Al Seesi S, et al. 2014. Genomic and bioinformatic profiling of mutational neoepitopes reveals new rules to predict anticancer immunogenicity (DAI). *Journal of Experimental Medicine* 211(11):2231-2248.
- Richman LP, Vonderheide RH, Rech AJ. 2019. Neoantigen dissimilarity to the self-proteome predicts immunogenicity and response to immune checkpoint blockade. *Cell Systems* 9(4):375-382.
- Lang F, Riesgo-Ferreiro P, Löwer M, Sahin U, Schrörs B. 2021. NeoFox: annotating neoantigen candidates with neoantigen features. *Bioinformatics* 37(22):4246-4247.
- Hundal J, Kiwala S, McMichael J, et al. 2020. pVACtools: a computational toolkit to identify and visualize cancer neoantigens. *Cancer Immunology Research* 8(3):409-420.

## Related Skills

- immunoinformatics/neoantigen-prediction - produces the candidate list this skill ranks
- immunoinformatics/mhc-binding-prediction - the presentation features that carry most of the signal
- immunoinformatics/mhc-class-ii-prediction - CD4 immunogenicity, the under-served frontier
- immunoinformatics/epitope-prediction - epitope candidates feeding the ranking
- clinical-databases/somatic-signatures - clonal neoantigen burden as an ICI-response correlate
