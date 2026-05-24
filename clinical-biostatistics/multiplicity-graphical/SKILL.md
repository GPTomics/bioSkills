---
name: bio-clinical-biostatistics-multiplicity-graphical
description: Implements multiplicity control for confirmatory clinical trials using graphical procedures (Bretz-Maurer-Hommel), gatekeeping (parallel, serial, mixed), Hochberg/Hommel/Holm with PRDS, and the closed-testing principle (Marcus-Peritz-Gabriel; Goeman 2021 admissibility). Covers FDA Multiple Endpoints Final Guidance (October 2022), graphical procedures via R gMCP, primary + key-secondary + subgroup hierarchies, and FWER vs FDR distinction. Use when designing the multiplicity strategy for confirmatory trials with multiple primary or key secondary endpoints.
tool_type: r
primary_tool: gMCP
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: R `gMCP` 0.8.16+, `graphicalMCP` 0.2+, `gatekeeping`, `multcomp`, `multxpert`; Python `statsmodels` 0.14+ for basic FDR/FWER methods.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Multiplicity Control for Confirmatory Trials

**"Design the multiplicity strategy for my trial"** -> Specify a closed-testing procedure (graphical, gatekeeping, hierarchical, or step-down Bonferroni-Holm) that controls family-wise error rate at the trial-wide level across primary endpoints, key secondary endpoints, and subgroup analyses, with provable strong FWER control.

## The Foundational Theorem -- Closed Testing Is Necessary

**Marcus, Peritz & Gabriel 1976 *Biometrika* 63:655:** a hypothesis H_I (I ⊆ {1,...,m}) is rejected iff every intersection hypothesis ∩_{J⊇I} H_J is rejected by a valid α-level local test. Strong FWER control holds for ANY choice of local tests.

**Goeman, Hemerik & Solari 2021 *Ann Stat* 49:1218** tightens this: closed testing is not merely sufficient — it is *necessary* for admissibility under FDP/FWER/k-FWER. **Every admissible multiplicity procedure is equivalent to some closed test.** Graphical procedures, gatekeepers, Hommel, fixed-sequence, fallback — all are closed tests in disguise.

**FWER vs FDR philosophical divide:**

- **FWER:** P(any false positive among m tests) — regulatory standard for confirmatory inference (agency wants to bound per-trial false-positive rate)
- **FDR:** Expected proportion of false discoveries among rejections — exploratory standard (genomics, fMRI, biomarker screens) where many true positives expected

Confirmatory clinical trials use FWER essentially universally.

## Algorithmic Taxonomy

| Procedure | Type | FWER control | Power profile | Use case |
|-----------|------|--------------|---------------|----------|
| Bonferroni | Single-step | Yes, any dependence | Conservative; loses 30-50% power vs Hommel under positive dependence | Very small m; worst-case dependence |
| Holm 1979 | Step-down | Yes, any dependence | Better than Bonferroni; uniformly dominates | Default for any dependence pattern |
| Hochberg 1988 | Step-up | Yes under PRDS (Sarkar 1998) | Better than Holm under PRDS | Positive correlation; verify PRDS |
| Hommel 1988 | Step-up via closed tests | Yes under PRDS | Uniformly dominates Hochberg by 1-3% | Whenever Hochberg is valid |
| Fixed-sequence (hierarchical) | Sequential | Yes, any dependence | Full alpha for first; subsequent zero if any fail | When clear priority ordering; "key secondary" labelling |
| Parallel gatekeeping (Dmitrienko 2003) | Multi-family | Yes | Family-by-family; secondary tested if any primary rejects | Primary family + secondary family |
| Serial gatekeeping | Sequential families | Yes | Strict: family k tested only if ALL of family k-1 reject | Co-primary + secondary tiers |
| Mixed gatekeeping (Dmitrienko-Tamhane 2008) | Combination | Yes | Combines closed-testing local procedures across families | Complex hierarchies |
| Graphical procedures (Bretz-Maurer 2009) | Closed-test as directed graph | Yes by construction | Flexible; allocate alpha to hypotheses via graph weights | Modern standard for confirmatory SAPs |
| Graphical + Simes/parametric (Bretz et al 2011) | Closed-test with non-Bonferroni local tests | Yes when Simes valid | Gains power under correlation | Complex co-primary + key secondary + subgroup hierarchies |
| Maurer-Bretz 2013 entangled graphs | Memory-augmented graphs | Yes by construction | Alpha propagation depends on origin | Parent-descendant constraints |
| Benjamini-Hochberg 1995 | FDR | FDR controlled at level q | Higher power than FWER | Exploratory only; NOT for confirmatory regulatory |

**Postdoc reading list:**

- Marcus R, Peritz E, Gabriel KR 1976 *Biometrika* 63:655 (closed testing — the foundation)
- Goeman JJ, Hemerik J, Solari A 2021 *Ann Stat* 49:1218 (closed testing necessary for admissibility)
- Holm S 1979 *Scand J Stat* 6:65 (step-down Bonferroni)
- Hochberg Y 1988 *Biometrika* 75:800 (step-up Simes)
- Hommel G 1988 *Biometrika* 75:383 (closed Simes; dominates Hochberg)
- Sarkar SK 1998/2008 *Ann Stat* (PRDS for Hochberg validity)
- Bretz F, Maurer W, Brannath W, Posch M 2009 *Stat Med* 28:586 (graphical procedures — foundational paper)
- Bretz F, Posch K, Glimm E, Klinglmueller F, Maurer W, Rohmeyer K 2011 *Biom J* 53:894 (Simes/parametric extensions)
- Maurer W, Bretz F 2013 *Stat Med* 32:1739 (entangled graphs / memory)
- Dmitrienko A, Offen WW, Westfall PH 2003 *Stat Med* 22:2387 (parallel gatekeeping)
- Dmitrienko A, Tamhane AC, Wiens B 2008 *Biom J* (mixed/multistage gatekeeping)
- FDA 2022 *Multiple Endpoints in Clinical Trials* Final Guidance (October 2022)
- Pocock SJ, Ariti CA, Collier TJ, Wang D 2012 *Eur Heart J* (win-ratio)

## Decision Tree by Scenario

| Scenario | Recommended procedure | Why |
|----------|----------------------|-----|
| 2 co-primary endpoints (both must succeed) | No alpha split needed; per-endpoint alpha-level test; cite FDA 2022 | Co-primary doesn't split alpha; inflates n via joint power |
| 2 multiple primary endpoints (any-wins) | Graphical procedure or Holm with weights | Alpha must be allocated; graphical is flexible |
| 1 primary + 2 key secondary endpoints | Hierarchical (serial gatekeeping) OR graphical with alpha propagation | Modern SAPs favour graphical |
| 1 primary + 3 secondary + 4 subgroup analyses | Graphical procedure via gMCP with pre-specified weights | Complex hierarchies benefit from graph visualisation |
| Primary endpoint + tipping-point sensitivity | No multiplicity adjustment needed for sensitivity | Sensitivity is "what if" not "another claim" |
| Many exploratory biomarker subgroups | Benjamini-Hochberg FDR | Exploratory; not for label claims |
| Win-ratio composite (cardiology) | Single test; no multiplicity | Composite captures multiple events in single hierarchy |
| Subgroup analysis (pre-specified) | Graphical alpha allocation; small budget (≤20%) per Dane 2019 | Confirmatory subgroup discovery requires explicit allocation |
| Adaptive trial with treatment arm dropping | Combination tests (Bauer-Köhne 1994) + closed testing | See clinical-biostatistics/adaptive-designs |
| Group-sequential with multiple endpoints | gsDesign or rpact with multivariate alpha spending | Hierarchical alpha across both time and endpoints |

## Bretz-Maurer Graphical Procedures -- The Modern Standard

**The Bretz-Maurer-Brannath-Posch 2009 *Stat Med* 28:586 framework recast weighted Bonferroni-Holm closed tests as directed weighted graphs:**

- Vertices = elementary null hypotheses with local weights summing to 1
- Directed edges = alpha-propagation rule (when a hypothesis is rejected, its weight redistributes to descendants per edge weights)
- The graph IS the procedure: a single visual fully specifies a closed-test procedure across primary, key secondary, and subgroup hierarchies

### gMCP R package

```r
library(gMCP)

# Construct a graph for primary + 2 key secondary endpoints
# Primary endpoint at full alpha; if rejected, alpha propagates equally to secondaries
hypotheses <- c('Primary', 'Sec1', 'Sec2')
weights <- c(1, 0, 0)  # initial alpha all on primary
# Transition matrix: rows = source, columns = target
# When Primary rejects, weight 0.5 goes to each secondary; when Sec1/Sec2 rejects, alpha returns
transitions <- matrix(c(
    0,    0.5,  0.5,
    0,    0,    1,
    0,    1,    0
), nrow = 3, byrow = TRUE, dimnames = list(hypotheses, hypotheses))

graph <- graphMCP(m = transitions, weights = weights, hnames = hypotheses)
# Note: in current gMCP, the graph constructor is `graphMCP(m=, weights=, hnames=)`;
# `matrix2graph()` appeared in older tutorials and is not the canonical exported API
# -- verify with `?graphMCP` / `?gMCP` in the installed gMCP release before scripting.
# Set p-values from the trial
p_vals <- c(Primary = 0.018, Sec1 = 0.042, Sec2 = 0.038)

# Run the graphical procedure at alpha = 0.025
result <- gMCP(graph, pvalues = p_vals, alpha = 0.025)
print(result)
# Hierarchical rejection: Primary rejects -> alpha propagates to secondaries -> ...
```

### Standard SAP graph patterns

| Pattern | Graph topology | Use |
|---------|----------------|-----|
| Pure hierarchical (fixed sequence) | H1 → H2 → H3 with weight 1 on each transition | Strict ordering |
| Holm graph (equal weights) | Each Hi → Hj with weight 1/(m-1) | No priority ordering |
| Primary + secondaries | Primary → Sec1 (0.5), Sec2 (0.5); Sec1 ↔ Sec2 (1) | Pivotal labeling claims |
| Co-primary chain | H1 → H2 with full weight if BOTH H1a, H1b reject | Co-primary + secondary |
| Subgroup branch | Primary → Subgroup_OS (0.2), Sec1 (0.4), Sec2 (0.4) | Discovery subgroup with budget |

### Bretz et al 2011 -- Simes and parametric extensions

When endpoints are positively correlated, replace the Bonferroni-based intersection test with Simes (for positive dependence) or parametric (using known correlation):

```r
library(gMCP)
# Use Simes-based local tests at each intersection
result_simes <- gMCP(graph, pvalues = p_vals, alpha = 0.025, test = 'Simes')
# Or parametric with estimated correlation matrix
result_param <- gMCP(graph, pvalues = p_vals, alpha = 0.025, corr = correlation_matrix)
```

### Maurer-Bretz 2013 entangled graphs

**Entangled graphs** add memory: the alpha propagation can depend on the *origin* of the alpha. This allows parent-descendant constraints that a single non-entangled graph cannot express. Example: secondary endpoint Sec1 receives alpha only from Primary, never from Sec2.

**Postdoc argument:** purists argue memory makes the procedure non-coherent in Gabriel's sense; Glimm/Maurer/Bretz argue it matches real-world inferential intent.

**Gabriel coherence in plain terms:** a coherent procedure rejects a hypothesis H consistently regardless of which superset of H is being tested. Non-entangled graphs are coherent: if H1 is rejected via path A, it would also be rejected via path B. **Entangled (memory-bearing) graphs sacrifice coherence:** the same H may be rejected when alpha arrives from one parent but not from another, because the propagation history changes the available alpha. The trade-off is operational power -- entangled graphs can encode "secondary X is meaningful only if primary Y rejects, not if primary Z rejects" inferential intent that flat coherent procedures cannot express. Choose based on whether the SAP needs path-dependent priority.

## Gatekeeping Procedures

### Serial gatekeeping (hierarchical)

Test H1 at full alpha; only if it rejects, test H2 at full alpha; etc. **Maximises power for H1** but H_k becomes inferentially worthless once any H_j (j<k) fails.

```r
# Hierarchical / serial: just a chain graph in gMCP
hyp <- c('H1', 'H2', 'H3', 'H4')
weights <- c(1, 0, 0, 0)
trans <- matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0),
                 nrow=4, dimnames=list(hyp, hyp))
graph <- graphMCP(m = trans, weights = weights, hnames = hyp)
```

**Pre-specification of order is critical** — based on clinical importance, NOT expected effect size. Ordering by expected effect is data-driven and inflates Type-I.

### Parallel gatekeeping (Dmitrienko 2003)

Secondary family is tested only if *at least one* primary rejects. Bonferroni-based parallel gatekeeper has stepwise representation (Guilbaud 2007 *Biom J* 49:917).

### Mixed / multistage (Dmitrienko-Tamhane 2008)

Permits using any closed-testing local procedure (e.g., Holm in family 1, Hommel in family 2) and combining via closure principle. R `gMCP::generalMixGatekeeping` or `Mediana`/`MultXpert` packages.

**Postdoc tradeoff:** parallel gatekeeping power loss vs collapsing endpoints into a composite (which avoids multiplicity but dilutes effect if components move in opposite directions); whether tree gatekeeping (Dmitrienko 2008 *Stat Med* 27:2114) over-engineers vs equivalent graphical procedure.

## Hochberg vs Hommel vs Holm

**Holm 1979** — step-down rejective Bonferroni; FWER controlled under any joint dependence. **Conservative but robust.**

**Hochberg 1988** — step-up using ordered Simes critical values; needs Simes inequality which requires PRDS (Sarkar 1998, 2008). Under PRDS, Hochberg uniformly dominates Holm.

**Hommel 1988** — also Simes-based but uses closed-testing tableau directly (not step-up shortcut). Uniformly more powerful than Hochberg (typically 1-3% gain).

### When Hochberg fails (anti-conservative)

**Hochberg becomes Type-I-inflated under negative dependence** — relevant when comparing endpoints mathematically constrained to move in opposite directions (LDL-C and HDL-C; complementary efficacy and safety endpoints).

**Sarkar critique:** when PRDS cannot be proven, fall back to Holm. The lost power is the price of robustness.

```python
# Python: statsmodels supports Holm, Hochberg, Hommel
from statsmodels.stats.multitest import multipletests

p_vals = [0.018, 0.042, 0.038, 0.015]
for method in ['holm', 'hochberg', 'hommel', 'bonferroni']:
    reject, adj_p, _, _ = multipletests(p_vals, alpha=0.05, method=method)
    print(f'{method}: reject={reject}, adjusted={adj_p}')
```

## FDA Multiple Endpoints Final Guidance (October 2022)

**Federal Register 2022-22882** finalises 2017 draft. Key changes vs draft:

- Explicit recognition of newer methods including win-ratio (Pocock 2012 *Eur Heart J*) and weighted composites
- Clearer language that "key secondary" endpoints are those for which sponsor wishes to make label claims and which must be in a Type-I-error-controlled hierarchy
- Appendix with worked graphical-procedure examples

### Categories

| Category | Approach | Note |
|----------|----------|------|
| Composite | Single test; no multiplicity | Win-ratio, DOOR/RADAR, time-to-first-event |
| Co-primary (all-win) | Each at full alpha; n inflated for joint power | Power = product of marginals |
| Multiple primary (any-wins) | Alpha must be split (Bonferroni or graphical) | More n required than co-primary if effects similar |
| Primary + key secondary | Hierarchical or graphical | Modern preference: graphical for flexibility |

**Winner's bias warning:** when post-hoc-selected endpoints are emphasised, bias-corrected effect estimates are recommended (same selection-bias issue as adaptive design).

## The "Almighty Primary Endpoint" Critique

**Dmitrienko-D'Agostino 2017 *Stat Med* 36:4341 editorial** + 2024 *Pharm Stat* discussion: the dogma of a single primary endpoint causes systematic Type-II error inflation in trials with broad multi-domain benefit (heart failure drugs with effects on mortality, hospitalisation, symptoms, biomarkers).

**Win-ratio (Pocock-Ariti-Collier-Wang 2012)** and **hierarchical composite (DOOR/RADAR, Evans 2015)** are responses — they preserve a single inferential test while letting multiple endpoints contribute.

**FDA counter-position** (Hung, O'Neill, Wang): without a designated primary, sponsors and regulators negotiate over secondary endpoints post hoc, destroying inferential meaning. Hence the FDA 2022 guidance reaffirms key-secondary hierarchies.

## Per-Method Failure Modes

### Hochberg under negative dependence

- **Trigger:** Endpoints constrained to move in opposite directions (LDL vs HDL; efficacy vs harm).
- **Mechanism:** Hochberg's PRDS assumption fails; Simes inequality doesn't hold; Type-I inflated.
- **Symptom:** Replication with Holm finds non-significant where Hochberg rejected.
- **Fix:** Switch to Holm (no PRDS assumption); cite Sarkar 1998.

### Fixed-sequence with wrong ordering

- **Trigger:** Ordering by expected effect size rather than clinical priority.
- **Mechanism:** Data-driven ordering inflates Type-I.
- **Symptom:** Reviewer asks for pre-specified ordering rationale.
- **Fix:** Pre-specify order by clinical priority in SAP; document rationale.

### Graphical procedure without pre-specified weights

- **Trigger:** Weights chosen at analysis time to favour observed results.
- **Mechanism:** Equivalent to post-hoc multiplicity tuning; inflates Type-I.
- **Symptom:** Multiple "what if" graph variants in CSR.
- **Fix:** Pre-specify graph and weights in SAP; document at protocol design.

### Bonferroni when graph would gain power

- **Trigger:** Default conservative choice when no thought put into structure.
- **Mechanism:** Loses 30-50% power vs Hommel/graphical when m ~ 10 correlated tests.
- **Symptom:** Underpowered trial reaches non-significance where graphical procedure would.
- **Fix:** Design a proper graph in `gMCP`; cite Bretz-Maurer 2009.

### FDR used for confirmatory primary

- **Trigger:** SAP specifies BH-FDR for primary multiplicity.
- **Mechanism:** FDR controls expected proportion of false discoveries, not P(any false positive).
- **Symptom:** Regulatory reviewer rejects as non-confirmatory.
- **Fix:** FWER (graphical, Holm, Hochberg, Hommel) for confirmatory; FDR for exploratory only.

### Subgroup analyses claimed without multiplicity

- **Trigger:** Trial reports 10 subgroups with one significant at α=0.05.
- **Mechanism:** ~40% probability of at least one false positive under global null.
- **Symptom:** Cherry-picked subgroup claim in submission.
- **Fix:** Pre-specified graphical alpha allocation OR explicit hypothesis-generating label; cite EMA 2019 subgroup guideline.

### Win-ratio reported without hierarchical priority

- **Trigger:** Win-ratio composite with unspecified component priority.
- **Mechanism:** Component prioritisation drives the result; arbitrary choice = data-dependent answer.
- **Symptom:** Two analysts get different results from same data depending on hierarchy.
- **Fix:** Pre-specify hierarchy in SAP; sensitivity over alternative hierarchies.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| FWER for confirmatory; FDR for exploratory | ICH E9; FDA 2022 Multiple Endpoints | Regulatory standard universally |
| Bonferroni: ~10 tests -> 30-50% power loss | Sarkar 1998 PRDS | Conservative under positive dependence |
| PRDS required for Hochberg validity | Sarkar 2008 *Ann Stat* | Otherwise Type-I inflated; fall back to Holm |
| Subgroup α budget <=20% of total | Dane 2019 EFSPI white paper | Discipline against subgroup fishing |
| Key secondary requires hierarchy in SAP | FDA 2022 Final | Labeling claims need Type-I-controlled test |
| Composite avoids multiplicity but dilutes effect | Pocock 2012 *Eur Heart J* | Win-ratio captures heterogeneity in single test |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Hochberg applied to negatively-dependent endpoints | PRDS not checked | Switch to Holm (cite Sarkar 1998) |
| Fixed-sequence ordering data-driven | Post-hoc selection | Pre-specify clinical priority in SAP |
| Bonferroni at 10 correlated endpoints | Default conservatism | Graphical procedure (gMCP); 30-50% power gain |
| FDR for confirmatory primary | Misunderstanding error rates | FWER mandatory for confirmatory; FDR exploratory only |
| Graphical procedure run with multiple weight schemes | Post-hoc graph tuning | Pre-specify single graph in SAP |
| Subgroups significant without multiplicity | Cherry-picking | Pre-specified allocation OR explicit hypothesis-generating label |
| Co-primary treated as multiple primary | Confused alpha allocation | Co-primary: no alpha split; inflate n. Multiple primary: split alpha |
| Win-ratio component priority unspecified | Data-driven choice | Pre-specify hierarchy with rationale; sensitivity over alternatives |
| `multipletests` default `method='hs'` (Holm-Sidak) | Common Python mistake | Always specify `method='holm'`, `'hommel'`, etc., explicitly |
| Sensitivity analysis listed as a "key secondary" requiring alpha | Confusion about role | Sensitivity is "what if" not "another claim"; no alpha needed |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why this multiplicity procedure?" | Closed testing per Marcus-Peritz-Gabriel; specific implementation is graphical (Bretz-Maurer 2009) with pre-specified weights in SAP |
| "Why Hommel not Holm?" | PRDS holds (positive correlation among endpoints); Hommel dominates Holm by 1-3% with no Type-I cost |
| "Why graph weights X, Y, Z?" | Clinical priority: primary > key secondary > exploratory; weights reflect labelling claim hierarchy |
| "Are these endpoints positively correlated?" | Sensitivity analyses provided: Bonferroni, Holm, Hochberg, Hommel results all in CSR appendix; concordant |
| "Where is alpha for the subgroup analysis?" | Pre-specified 20% of primary alpha allocated; cite Dane 2019 |
| "Why not just composite endpoint?" | Composite would dilute differential effect on mortality vs hospitalisation; key-secondary hierarchy preserves component-level claims |
| "PRDS check for Hochberg?" | Endpoints positively correlated via simulation under null; PRDS holds; Hochberg/Hommel valid |
| "Sensitivity in the hierarchy?" | No — sensitivity is "what if" and does not require alpha. Listed as supportive not key secondary. |

## References

- Bretz F, Maurer W, Brannath W, Posch M. 2009. A graphical approach to sequentially rejective multiple test procedures. *Stat Med* 28:586-604.
- Bretz F, Posch K, Glimm E, Klinglmueller F, Maurer W, Rohmeyer K. 2011. Graphical approaches for multiple comparison procedures using weighted Bonferroni, Simes, or parametric tests. *Biom J* 53:894-913.
- Burman CF, Sonesson C, Guilbaud O. 2009. A recycling framework for the construction of Bonferroni-based multiple tests. *Stat Med* 28:739-761.
- Dmitrienko A, Offen WW, Westfall PH. 2003. Gatekeeping strategies for clinical trials that do not require all primary effects to be significant. *Stat Med* 22:2387-2400.
- Dmitrienko A, Tamhane AC, Wiens BL. 2008. General multistage gatekeeping procedures. *Biom J* 50:667-677.
- FDA. 2022. Multiple Endpoints in Clinical Trials. Final Guidance.
- Goeman JJ, Hemerik J, Solari A. 2021. Only closed testing procedures are admissible for controlling false discovery proportions. *Ann Stat* 49:1218-1238.
- Guilbaud O. 2007. Bonferroni parallel gatekeeping -- transparent generalizations, adjusted p-values, and short proofs. *Biom J* 49:917-927.
- Hochberg Y. 1988. A sharper Bonferroni procedure for multiple tests of significance. *Biometrika* 75:800-802.
- Holm S. 1979. A simple sequentially rejective multiple test procedure. *Scand J Stat* 6:65-70.
- Hommel G. 1988. A stagewise rejective multiple test procedure based on a modified Bonferroni test. *Biometrika* 75:383-386.
- Marcus R, Peritz E, Gabriel KR. 1976. On closed testing procedures with special reference to ordered analysis of variance. *Biometrika* 63:655-660.
- Maurer W, Bretz F. 2013. A graphical approach to multiple testing in clinical trials with simultaneous testing of multiple hypotheses. *Stat Med* 32:1739-1753.
- Pocock SJ, Ariti CA, Collier TJ, Wang D. 2012. The win ratio: a new approach to the analysis of composite endpoints in clinical trials. *Eur Heart J* 33:176-182.
- Sarkar SK. 2008. Generalizing Simes' test and Hochberg's stepup procedure. *Ann Stat* 36:337-363.

## Related Skills

- clinical-biostatistics/trial-reporting - Multiplicity strategy reporting per CONSORT 2025
- clinical-biostatistics/subgroup-analysis - Subgroup multiplicity allocation
- clinical-biostatistics/power-and-sample-size - Power adjustment for co-primary endpoints
- clinical-biostatistics/adaptive-designs - Combination tests for adaptive multiplicity
- clinical-biostatistics/effect-measures - Reporting multiple effect measures post-multiplicity adjustment
- experimental-design/multiple-testing - General methods (FDR, FWER, q-values)
