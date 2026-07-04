---
name: bio-temporal-genomics-differential-rhythmicity
description: Compares how a rhythm CHANGES between conditions, genotypes, treatments, tissues, or ages (differential rhythmicity), classifying each feature as gain-of-rhythm, loss-of-rhythm, phase change, amplitude change, unchanged-rhythmic, or arrhythmic-in-both, and distinguishing differential EXPRESSION (condition main effect) from differential RHYTHMICITY (condition x time interaction). Uses model-based approaches that borrow strength across conditions - LimoRhyde (sin/cos interaction terms in a limma/edgeR/DESeq2 design), dryR (BIC model selection across >=2 conditions), compareRhythms (direct gain/loss/change/same classification), DODR, CircaCompare - instead of the detect-then-Venn anti-pattern that overestimates reprogramming. Use when testing whether rhythms differ between conditions/genotypes/tissues/ages, classifying gain/loss/phase/amplitude change, or separating differential expression from differential rhythmicity. Not for detecting rhythms in one condition (see temporal-genomics/circadian-rhythms).
tool_type: r
primary_tool: limorhyde
---

## Version Compatibility

Reference examples tested with: limorhyde 1.0+, limma 3.50+, compareRhythms 1.0+, dryR (GitHub naef-lab), DODR 0.99+, CircaCompare 0.2+, R 4.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: differential rhythmicity estimates an INTERACTION (condition x time), so it needs MORE power than single-condition detection - the same matched, evenly-sampled, replicated grid must exist in BOTH conditions, and unmatched timepoints across conditions break the interaction model.

# Differential Rhythmicity (comparing rhythms between conditions)

**"How does this gene's rhythm change in the knockout / high-fat diet / aged tissue?"** -> classify each feature's rhythm CHANGE relative to a reference condition into gain, loss, phase, or amplitude change, using a model that borrows strength across conditions.
- R: `limorhyde()` (sin/cos basis for a limma/edgeR/DESeq2 interaction model); `compareRhythms()` (direct gain/loss/change/same classification); `dryseq()` (dryR, BIC model selection); `DODR`, `circacompare()`, `diffCircadian` (targeted comparators)

This is the sibling of temporal-genomics/circadian-rhythms, which asks the single-condition question "is this feature rhythmic at 24h?" (cosinor/JTK/RAIN detection). Do not re-run detection here; this skill answers the categorically different question of how a rhythm DIFFERS between groups.

## The governing principle: differential rhythmicity is a distinct question, and detect-then-Venn overestimates it

The #1 error is the detect-then-Venn anti-pattern: defining "genes that lost rhythm in the KO" as "rhythmic in WT MINUS rhythmic in KO" from two independently thresholded rhythmicity lists. This systematically OVERESTIMATES reprogramming, because near p~0.05 the two lists differ mostly from threshold noise, not biology - a gene at q=0.04 in WT and q=0.06 in KO is called "lost" when nothing changed (Pelikan 2022 *FEBS J* 289:6605). The fix is to never intersect two lists: fit ONE model spanning both conditions and test the condition x time INTERACTION directly, so a single calibrated test asks "did the rhythm change?" with strength borrowed across conditions.

The second load-bearing distinction is differential EXPRESSION vs differential RHYTHMICITY. A gene whose mean level shifts between conditions but whose oscillation is unchanged is differentially EXPRESSED, not differentially rhythmic. In the sin/cos framing these are orthogonal: the condition MAIN effect = differential expression (mean shift, adjusting for time), the condition:time INTERACTION = differential rhythmicity (amplitude/phase change). Reporting a main-effect hit as "rhythm reprogramming" is a common and wrong conflation.

## The four canonical outcome classes (relative to a reference group)

| Class | What changed | Interaction signature |
|-------|--------------|-----------------------|
| Gain of rhythm | arrhythmic in reference, rhythmic in test | interaction significant; reference amplitude ~0 |
| Loss of rhythm | rhythmic in reference, arrhythmic in test | interaction significant; test amplitude ~0 |
| Phase change | rhythmic in both, peak time shifted | interaction significant; amplitudes similar, acrophases differ |
| Amplitude change | rhythmic in both, oscillation damped/amplified | interaction significant; same phase, amplitudes differ |
| Unchanged-rhythmic (same) | rhythmic in both, same amplitude+phase | interaction NOT significant; both rhythmic |
| Arrhythmic-in-both | flat in both | neither main effect nor interaction; excluded by amplitude filter |

Interpret phase and amplitude CHANGE only for features confidently rhythmic in AT LEAST ONE condition; a "phase shift" between two genes that are arrhythmic in both is noise. Amplitude-filter first (peak-to-trough or relative amplitude), then classify.

## Method selection (which to pick and why)

| Method | Pick when | Mechanism | Fails / caveat |
|--------|-----------|-----------|----------------|
| LimoRhyde + limma/edgeR/DESeq2 | 2 conditions; want to fold DR into a standard DE pipeline with covariates/batch; count or microarray data | `limorhyde()` adds sin/cos time columns; condition:time interaction = DR, condition main effect = DE | Two-condition framing; the interaction test flags THAT a rhythm changed, not which class - read per-condition amplitude/phase to classify |
| dryR (`dryseq`) | >=2 conditions; want a parsimonious per-gene MODEL assignment (shared vs independent rhythm parameters) | BIC model selection over a family of shared/independent-parameter models, tailored to RNA-seq noise | Model-selection categories depend on the model family and BIC penalty; needs enough timepoints for BIC to discriminate models |
| compareRhythms | want gain/loss/change/same DIRECTLY, built to replace the Venn approach; microarray or RNA-seq | wraps model-selection (`mod_sel`) or hypothesis tests (`dodr`/`limma`/`voom`/`deseq2`/`edger`/`cosinor`); classifies vs reference | Two groups only; `mod_sel` needs no DE package but `deseq2`/`edger`/`voom` do; a feature must clear `amp_cutoff` in >=1 group to be reported |
| DODR | direct two-condition differential-rhythmicity test on already-detected rhythmic features | robust/rank comparison of rhythm shape (amplitude, phase, signal-to-noise) between conditions | Tests differential rhythmicity given rhythmicity; pre-filter to features rhythmic in >=1 group first |
| CircaCompare | a FEW targeted genes; want explicit estimates + p-values for the mesor/amplitude/phase DIFFERENCE | non-linear regression fitting both curves jointly with difference parameters | Compares two groups only if BOTH are rhythmic (amplitude non-zero); not a genome-scale screen |
| diffCircadian | a few genes; want likelihood-ratio tests separating differential amplitude vs phase vs basal vs fit | likelihood-based tests (`LR_diff`) on two conditions | Two conditions; targeted rather than transcriptome-wide throughput |

Methodology here is evolving (LimoRhyde2 reframes around effect-size/posterior shrinkage; benchmarks disagree on the best classifier) - verify current best practice against the latest tool docs before committing to one method, and prefer a screen (LimoRhyde/dryR/compareRhythms) followed by targeted confirmation (CircaCompare/diffCircadian) on hits.

## LimoRhyde + limma interaction test

**Goal:** Rank features by differential rhythmicity between two conditions while separately quantifying differential expression, in a single linear model.

**Approach:** Decompose measured time into a sin/cos basis with `limorhyde()`, fit `condition*(time_cos+time_sin)`, moderated-F-test the two interaction coefficients for DR and the condition main effect for DE.

```r
library(limorhyde); library(limma)

# limorhyde() decomposes measured time into a cosinor basis; prefix 'time_' names them time_cos, time_sin.
# period = 24h circadian; time is MEASURED (ZT/CT), not inferred pseudotime.
meta <- cbind(meta, limorhyde(meta$time, 'time_', period = 24))

# Differential RHYTHMICITY = condition:time interaction; differential EXPRESSION = condition main effect.
design <- model.matrix(~ condition * (time_cos + time_sin), data = meta)
fit <- eBayes(lmFit(expr, design))   # expr: features x samples, columns aligned to meta rows

dr_cols <- grep('conditionKO:time_', colnames(design), value = TRUE)  # the two interaction coefficients
dr <- topTable(fit, coef = dr_cols, number = Inf, sort.by = 'F')      # BH-adjusted adj.P.Val ranks DR
de <- topTable(fit, coef = 'conditionKO', number = Inf, sort.by = 'p')  # condition main effect = DE
```

The interaction F-test says a rhythm changed; it does not name the class. To assign gain/loss/phase/amplitude, fit a per-condition cosinor (or read `limorhyde2` posterior estimates) and compare amplitudes and acrophases between conditions for the significant features. For count data, run the same design through `voom`+limma, edgeR, or DESeq2 with the sin/cos and interaction columns.

## compareRhythms direct classification

**Goal:** Assign every feature directly to gain / loss / change / same relative to a reference condition, without intersecting two detection lists.

**Approach:** Pass a features x samples matrix plus an `exp_design` data.frame (numeric `time`, 2-level factor `group`), choose a `method`, and read the returned category per feature.

```r
library(compareRhythms)

# exp_design: one row per sample; numeric 'time', factor 'group' with EXACTLY 2 levels (reference first).
# data: numeric matrix, rows = features (rownames = ids), columns = samples matching exp_design rows.
# method='mod_sel' = BIC model selection (no DE package); 'deseq2'/'edger'/'voom' for RNA-seq counts;
# 'limma' for log-microarray; 'dodr'/'cosinor' also available. amp_cutoff = peak-to-trough floor (>=1 group).
res <- compareRhythms(data, exp_design = exp_design, period = 24,
                      method = 'mod_sel', amp_cutoff = 0.5, criterion = 'bic')
# res: data.frame with id + category (gain / loss / change / same, relative to the reference group).
```

For >2 conditions, dryR does BIC model selection across all conditions at once: `dryseq(counts, group, time)` assigns each gene a rhythm-parameter-sharing model and returns per-condition amplitude/phase/mean. Prefer it over pairwise interaction tests when the design has three or more groups and a parsimonious classification is wanted.

## The confounds a differential-rhythmicity claim must address

**Reduced bulk amplitude can be loss of SYNCHRONY, not loss of per-cell rhythm.** A bulk/tissue readout is the sum over many single-cell oscillators; if cells DESYNCHRONIZE (dephase) between conditions, the ensemble amplitude damps toward zero even though every cell still oscillates. Bulk "amplitude change" or "loss of rhythm" therefore has three indistinguishable causes: true loss of cell-autonomous rhythmicity, loss of inter-cell synchrony, or reduced single-cell amplitude. Report it as "reduced ensemble amplitude" and use single-cell or live-imaging assays to separate the causes; see single-cell/preprocessing for cell-level analysis.

**A rhythm change under light-dark may be a driven change, not a clock change.** Under an entraining LD cycle, an apparent rhythm can be masked (driven directly by light/feeding/temperature). A between-condition difference (e.g. a feeding-time or lighting manipulation) can shift the DRIVEN component without touching the endogenous clock. Only free-running (DD/constant) conditions license "the clock rewired"; under LD, a differential-rhythmicity hit may reflect a change in the environmental drive. Declare the light regime and use ZT (entrained) vs CT (free-running) accordingly.

## Design constraints (shared with detection, stricter here)

- Matched sampling GRID across conditions: the SAME timepoints, evenly spaced, in every condition. Unmatched timepoints break the interaction model (the condition:time terms become non-estimable or confounded) - this is the constraint most often violated when two datasets are compared post hoc.
- >=2 full cycles and >=6 (ideally 8-12) samples/cycle in EACH condition (the genome-scale rhythm-analysis design guidelines of Hughes et al. 2017 apply per condition). One cycle cannot separate a rhythm from a trend, so it certainly cannot compare rhythms.
- Replicates per timepoint in each condition. DR estimates an interaction (a difference of differences), so it is hungrier for power than detection; a design just adequate to detect a rhythm is usually too thin to confidently call a rhythm CHANGE.
- The amplitude filter is not optional. Classify phase/amplitude change only for features confidently rhythmic in at least one condition; significance without an effect-size floor over-detects (Laloum 2020 *PLoS Comput Biol* 16:e1007666), and the effect is worse for a difference test.

## Common Errors (trap -> fix)

| Trap | Fix |
|------|-----|
| Defining "lost rhythm" as (rhythmic in WT) minus (rhythmic in KO) via two thresholded lists | Fit one model across both conditions and test the condition:time INTERACTION (LimoRhyde/dryR/compareRhythms); the Venn approach overestimates reprogramming (Pelikan 2022) |
| Calling a condition MAIN-effect (mean-shift) hit "differential rhythmicity" | Main effect = differential EXPRESSION; differential RHYTHMICITY is the condition:time INTERACTION. Report them separately |
| Unmatched or unevenly-spaced timepoints across conditions | Use a matched, evenly-sampled grid in every condition; unmatched times make the interaction terms non-estimable or confound them with condition |
| Interpreting a "phase shift" for features arrhythmic in both conditions | Amplitude-filter FIRST; classify phase/amplitude change only for features confidently rhythmic in >=1 condition |
| Calling reduced BULK amplitude "loss of rhythm" | Ensemble amplitude damps from cell DESYNCHRONY too; report "reduced ensemble amplitude" and separate with single-cell or imaging assays |
| Claiming the clock "rewired" from LD (entrained) data | An LD rhythm change can be a driven/masking change (light/feeding); endogenous rewiring needs free-running (DD) conditions |
| Running DR on a design too sparse to even DETECT a rhythm | DR needs MORE power than detection (it estimates an interaction); ensure >=2 cycles, >=6/cycle, replicates in EACH condition before comparing |
| Reading the interaction F-test as the CLASS | The interaction says a rhythm changed, not which class; fit per-condition amplitude/phase (or LimoRhyde2 posteriors) to assign gain/loss/phase/amplitude |
| Using compareRhythms `deseq2`/`edger`/`voom` on already-normalized log data | Those methods expect RAW counts; use `mod_sel`/`limma`/`cosinor` for normalized or microarray data |

## Related Skills

temporal-genomics/circadian-rhythms - Single-condition rhythm DETECTION and parameter estimation (cosinor/JTK/RAIN); run it first, this skill compares its results across conditions
differential-expression/timeseries-de - Temporal differential expression (a monotone trend or between-timepoint change), which is differential EXPRESSION over time, not differential rhythmicity
temporal-genomics/temporal-clustering - Group differentially-rhythmic genes by the shape of their change
single-cell/preprocessing - Entry to cell-level analysis, to separate reduced ensemble amplitude (desynchrony) from true loss of per-cell rhythm

## References

- Singer JM, Hughey JJ. 2019. LimoRhyde: a flexible approach for differential analysis of rhythmic transcriptome data. J Biol Rhythms 34(1):5-18. doi:10.1177/0748730418813785
- Weger BD, Gobet C, David FPA, et al. 2021. Systematic analysis of differential rhythmic liver gene expression mediated by the circadian clock and feeding rhythms (dryR). PNAS 118(3):e2015803118. doi:10.1073/pnas.2015803118
- Pelikan A, Herzel H, Kramer A, Ananthasubramaniam B. 2022. Venn diagram analysis overestimates the extent of circadian rhythm reprogramming (compareRhythms). FEBS J 289(21):6605-6621. doi:10.1111/febs.16095
- Thaben PF, Westermark PO. 2016. Differential rhythmicity: detecting altered rhythmicity in biological data (DODR). Bioinformatics 32(18):2800-2808. doi:10.1093/bioinformatics/btw309
- Parsons R, Parsons R, Garner N, Oster H, Rawashdeh O. 2020. CircaCompare: a method to estimate and statistically support differences in mesor, amplitude and phase, between circadian rhythms. Bioinformatics 36(4):1208-1212. doi:10.1093/bioinformatics/btz730
- Ding H, Meng L, Liu AC, et al. 2021. Likelihood-based tests for detecting circadian rhythmicity and differential circadian patterns in transcriptomic applications (diffCircadian). Brief Bioinform 22(6):bbab224. doi:10.1093/bib/bbab224
- Hughes ME, Abruzzi KC, Allada R, et al. 2017. Guidelines for genome-scale analysis of biological rhythms. J Biol Rhythms 32(5):380-393. doi:10.1177/0748730417728663
- Laloum D, Robinson-Rechavi M. 2020. Methods detecting rhythmic gene expression are biologically relevant only for strong signal. PLoS Comput Biol 16(3):e1007666. doi:10.1371/journal.pcbi.1007666
