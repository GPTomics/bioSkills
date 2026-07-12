---
name: bio-differential-splicing
description: Detects differential alternative splicing between conditions using rMATS-turbo (binomial LRT on junction counts), leafcutter (Dirichlet-multinomial GLM on intron clusters), MAJIQ V3 deltapsi/HET (Bayesian posterior on LSVs), SUPPA2 (empirical-null on TPM-derived PSI), or Shiba (junction-imbalance-corrected, 2025 SOTA at low coverage). Reports FDR-corrected significance and delta PSI effect sizes. Tools differ in statistical model, annotation dependence, calibration regime, and replicate-count requirements. Use when comparing splicing patterns between treatment groups, tissues, or disease states.
tool_type: mixed
primary_tool: rMATS-turbo
---

## Version Compatibility

Reference examples tested with: rMATS-turbo 4.3+, SUPPA2 2.4+, leafcutter 0.2.9+, MAJIQ 3.0+, Shiba 0.5+, STAR 2.7.11+, regtools 1.0+, pandas 2.2+, R 4.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Splicing

Detect splicing changes between conditions. Tool choice is a decision about **statistical model**, **annotation dependence**, and **calibration regime** under the specific experimental design — not a preference. Wrong tool for the design produces uncalibrated FDR or systematic effect-size bias.

## Statistical Model Taxonomy

| Tool | Model | Test statistic | Min reps per group | Calibration regime | Fails when |
|------|-------|-----------------|---------------------|---------------------|------------|
| rMATS-turbo | Binomial counts with hierarchical PSI variance | LRT on \|ΔPSI\| > `cutoff` (default 0.0001) | n>=3 | Well-calibrated at n>=3 with adequate junction reads | Junction read imbalance; very low coverage; uncorrected for confounders |
| leafcutter | Dirichlet-multinomial GLM at cluster level | LRT on group factor | n>=2 (n>=3 preferred) | Strong at n>=3; novel-junction-friendly | Undersampled clusters (DM dispersion unstable); cluster topology arbitrariness |
| MAJIQ deltapsi | Beta-binomial bootstrap -> posterior over PSI per LSV | P(\|ΔPSI\| > T) threshold (T=0.2) | n>=3 | Replicate-structured n=3 vs n=3 | Cohorts where between-sample variability dominates between-group |
| MAJIQ HET | Same model, heterogeneity-aware | Per-LSV permutation-based test | n>=10 | n>=10 vs n>=10 cohort designs | Tightly-controlled small replicate experiments |
| SUPPA2 (empirical) | Empirical null from between-replicate ΔPSI | ECDF on \|ΔPSI\| conditioned on TPM | n>=4 | n>=4 vs n>=4 with paired-end deep sequencing | n<=3 vs n<=3 (sparse null collapses) |
| SUPPA2 (classical) | Wilcoxon rank-sum on PSI distributions | Wilcoxon p-value | n>=2 | Small samples; non-parametric backup | Cassette events with tight PSI distributions |
| Shiba (2025) | Beta-binomial with explicit junction-imbalance correction | LRT | n>=2 | n=2-3 vs n=2-3 | Established benchmarks limited (new tool) |
| LeafcutterMD | Dirichlet-multinomial outlier mode | Per-sample p-value | n=1 vs cohort >=20 | Single-patient vs cohort | Too few controls (<20) |
| FRASER 2.0 | Beta-binomial autoencoder on Intron Jaccard Index | Per-sample p-value with delta cutoff | n=1 vs cohort >=20 | n>=20 control cohort, single-patient query | See `outlier-splicing-detection` for this regime |

The first decision is which **regime** the design falls into: between-group with replicates, heterogeneous cohort, or single-sample-vs-cohort. Within each regime, tool choice is much smaller (1-2 options).

Comprehensive 2023-2026 benchmarks: Olofsson 2023 *Biochem Biophys Res Commun*; Tran 2025 *WIREs RNA*; Kubota 2025 *NAR*. Methodology evolves — verify benchmarks and tool docs before reporting. Default 2026 recommendation: run **two complementary tools** (rMATS + leafcutter) and require concordance for high-confidence calls.

## Decision Tree by Experimental Design

| Scenario | Recommended tool | Why | Threshold |
|----------|------------------|-----|-----------|
| Standard n=3 vs n=3, GENCODE-annotated | rMATS-turbo + leafcutter (concordance) | Two algorithmic families; concordant hits = high-confidence | FDR<0.05, \|ΔPSI\|>0.10 |
| n=2 vs n=2 small pilot | Shiba | Junction-imbalance correction matters most at low coverage | FDR<0.10, \|ΔPSI\|>0.10 |
| n=10+ vs n=10+ heterogeneous (clinical, GTEx-style) | MAJIQ V3 HET | HET designed for between-sample heterogeneity | P(\|ΔPSI\|>0.2)>0.95 |
| Single rare-disease patient vs panel of n>=20 | FRASER 2.0 (see outlier-splicing-detection) | Outlier detection statistical model is fundamentally different | padj<0.05, \|delta-jaccard\|>=0.1 |
| Time-course / multi-condition design | Custom DEXSeq or limma on PSI matrix | rMATS/leafcutter primarily 2-group | FDR<0.05 on time:group interaction |
| Paired tumor-normal | rMATS with `--paired-stats` | Paired test reduces inter-patient variance | FDR<0.05, paired \|ΔPSI\|>0.10 |
| Cancer with spliceosomal mutation (SF3B1, U2AF1) | leafcutter or MAJIQ denovo | Cryptic events not in annotation | FDR<0.05; check 3'ss shifts in IGV |
| TDP-43 loss / ALS post-mortem | leafcutter denovo | Cryptic exons not in annotation | FDR<0.05; expect UNC13A, STMN2 |
| Non-model organism without GENCODE-grade annotation | leafcutter | Annotation-free | FDR<0.05, \|ΔPSI\|>0.10 |
| Long-read available | rMATS-long, FLAIR diffSplice | See long-read-splicing | Tool-specific |

## rMATS-turbo Differential Analysis

**Goal:** Detect statistically significant differential splicing between two groups from BAMs.

**Approach:** Run rMATS-turbo without `--statoff`, then filter by FDR + ΔPSI + per-replicate coverage.

```bash
rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --libType fr-firststrand \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp \
    --novelSS \
    --cstat 0.05
```

`--cstat 0.05` tests `|ΔPSI| > 0.05`; raise to 0.10 for stricter discovery. `--novelSS` enables novel-junction discovery (recommended with STAR 2-pass). For paired designs, add `--paired-stats`.

```python
import pandas as pd
import numpy as np

se = pd.read_csv('rmats_output/SE.MATS.JC.txt', sep='\t')

def min_per_rep(s):
    return s.str.split(',').apply(lambda x: min(int(v) for v in x))

se['min_inc'] = min_per_rep(se['IJC_SAMPLE_1']).combine(min_per_rep(se['IJC_SAMPLE_2']), min)
se['min_skip'] = min_per_rep(se['SJC_SAMPLE_1']).combine(min_per_rep(se['SJC_SAMPLE_2']), min)

significant = se[
    (se['FDR'] < 0.05) &
    (se['IncLevelDifference'].abs() > 0.10) &
    ((se['min_inc'] + se['min_skip']) >= 10)
].copy()

significant['score'] = -np.log10(significant['FDR']) * significant['IncLevelDifference'].abs()
top = significant.nlargest(50, 'score')
```

## leafcutter Differential Intron Usage

**Goal:** Detect differential intron-cluster usage annotation-free, capturing novel junctions and complex multi-junction events.

**Approach:** Extract junctions with regtools, cluster introns by shared splice sites, run cluster-level Dirichlet-multinomial test.

```bash
for bam in *.bam; do
    regtools junctions extract -a 8 -m 50 -s XS "$bam" -o "${bam%.bam}.junc"
done
ls *.junc > juncfiles.txt

python leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -o leafcutter \
    -m 50 \
    -l 500000
```

```r
library(leafcutter)

groups <- data.frame(
    sample = c('s1', 's2', 's3', 's4', 's5', 's6'),
    group = c('control', 'control', 'control', 'treatment', 'treatment', 'treatment')
)
write.table(groups, 'groups.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

system('leafcutter_ds.R --num_threads 4 --exon_file gencode_exons.txt.gz \
    leafcutter_perind_numers.counts.gz groups.txt -o ds_results')

cluster_sig <- read.table('ds_results_cluster_significance.txt', header = TRUE, sep = '\t')
intron_effects <- read.table('ds_results_effect_sizes.txt', header = TRUE, sep = '\t')

sig_clusters <- subset(cluster_sig, p.adjust < 0.05)
```

**LeafCutter2** (Buen Abad Najar 2025 *bioRxiv*) extends leafcutter with NMD-aware classification of unproductive splicing — useful when AS-NMD coupling is the question.

## MAJIQ V3 Differential Analysis

**Goal:** Detect differential LSVs with full posterior distributions over ΔPSI; ideal for complex multi-junction events and heterogeneous cohorts.

**Approach:** Build splice graph -> compute coverage per group -> run deltapsi (replicate-structured) or heterogen (cohort-style).

```bash
majiq build annotation.gff3 -c settings.ini -j 8 -o build_output

majiq deltapsi \
    -grp1 build_output/ctrl1.majiq build_output/ctrl2.majiq build_output/ctrl3.majiq \
    -grp2 build_output/trt1.majiq build_output/trt2.majiq build_output/trt3.majiq \
    -n control treatment \
    -o deltapsi_output \
    --minreads 10 --minpos 3 \
    -j 8

majiq heterogen \
    -grp1 build_output/het_ctrl{1..20}.majiq \
    -grp2 build_output/het_trt{1..20}.majiq \
    -n control treatment \
    -o heterogen_output \
    -j 8

voila view -p 5000 -j 8 build_output/splicegraph.zarr deltapsi_output/control_treatment.deltapsi.voila -o voila_html
```

MAJIQ V3 (Aicher, Slaff, Jewell, Barash *bioRxiv* 2024; public release 2025) uses Zarr storage (`splicegraph.zarr`); V2's SQLite splicegraph is deprecated. MAJIQ reports posterior probability `P(|ΔPSI| > 0.2)`; thresholds are interpreted differently from FDR. Use HET for n>=10 vs n>=10 cohort designs (clinical, GTEx-style); deltapsi for tightly controlled n=3 vs n=3.

## SUPPA2 Differential Analysis

**Goal:** Quick differential splicing from existing transcript quantifications, useful as a sanity check or pilot.

**Approach:** Generate per-condition PSI files from Salmon TPM, then run `diffSplice` with empirical or classical p-values.

```bash
suppa.py generateEvents -i annotation.gtf -o events -f ioe -e SE SS MX RI

for ev in SE A5 A3 MX RI; do
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e ctrl_tpm.tsv -o ctrl_${ev}
    suppa.py psiPerEvent -i events_${ev}_strict.ioe -e trt_tpm.tsv -o trt_${ev}

    suppa.py diffSplice \
        -m empirical \
        -gc \
        -i events_${ev}_strict.ioe \
        -p ctrl_${ev}.psi trt_${ev}.psi \
        -e ctrl_tpm.tsv trt_tpm.tsv \
        -o diff_${ev}
done
```

For n<=3 designs, switch `-m classical` (Wilcoxon). Empirical null requires sufficient between-replicate observations to construct.

## Shiba for Low-Coverage / Few-Replicate Designs

**Goal:** Detect differential splicing with explicit junction-imbalance correction — addresses a known false-positive source for rMATS-style methods.

**Approach:** Shiba is a Snakemake-based pipeline configured via YAML. Install via bioconda, write a config file describing groups + BAMs, then run with snakemake.

```bash
conda install -c bioconda shiba

# Edit config.yaml with reference GTF, BAM groups, output dir, thresholds
# Then run the Snakemake workflow:
snakemake -s snakeshiba.smk \
    --configfile config.yaml \
    --cores 8 \
    --use-singularity \
    --singularity-args "--bind $HOME:$HOME"
```

Shiba (Kubota 2025 *NAR*) reportedly outperforms rMATS at n=2 vs n=2 by correcting differential mappability between inclusion and skipping junctions; community calibration still emerging. See https://sika-zheng-lab.github.io/Shiba/ for the full config.yaml schema.

## Per-Tool Failure Modes

### rMATS: Confounder-Blind LRT

**Trigger:** Sequencing batch, RIN, library prep date, or sex correlates with the comparison of interest.

**Mechanism:** rMATS' default LRT does not natively accept covariates the way DESeq2 does; the `--paired-stats` flag handles paired designs but not arbitrary covariates.

**Symptom:** Many "significant" hits driven by batch rather than biological condition; PCA on PSI matrix shows samples clustering by batch rather than group.

**Fix:** Either (a) include batch as a stratification (run rMATS within each batch), (b) regress PSI matrix against batch in R, then test residuals, or (c) switch to leafcutter which accepts a `confounders` argument in `differential_splicing`.

### leafcutter: Cluster Mis-Topology

**Trigger:** A cluster spans a complex topology (cassette + alternative donor in same cluster).

**Mechanism:** Cluster-level p-value reports "something in this cluster differs" but doesn't indicate which intron drove the change; downstream analysis needs per-intron effect sizes.

**Symptom:** Significant cluster, multiple introns with different ΔPSI directions, ambiguous biological interpretation.

**Fix:** Inspect cluster in leafviz; report per-intron effect sizes from `ds_results_effect_sizes.txt`; map cluster topology to canonical SE/A5SS/A3SS via flanking exon coordinates manually.

### MAJIQ HET: Power vs Type-1 Tradeoff

**Trigger:** HET module on n=5-10 cohorts (between regimes).

**Mechanism:** HET assumes between-sample variability dominates; for moderate-replicate designs (n=5-10), HET is conservative and deltapsi is more powerful.

**Symptom:** HET reports few hits in n=5-10 designs; deltapsi on the same data reports many.

**Fix:** Use deltapsi for n=3-5; reserve HET for n>=10 with explicit cohort heterogeneity.

### SUPPA2 Empirical: Sparse Null at Low Replicate

**Trigger:** n<=3 vs n<=3 with `--method empirical`.

**Mechanism:** Empirical null is ECDF of |ΔPSI| from between-replicate comparisons within each group, binned by transcript expression. Few replicates -> few null observations -> wide confidence on null distribution.

**Symptom:** Inflated FDR (15-30%); "significant" hits don't replicate or validate.

**Fix:** Use `-m classical` (Wilcoxon) for n<=3 vs n<=3; or switch tool entirely (leafcutter, Shiba).

## Reconciliation: When Tools Disagree

The two most common short-read tools answer slightly different questions: rMATS classifies on annotated event templates; leafcutter classifies on observed cluster usage. Disagreement is informative.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| rMATS sig, leafcutter not sig | rMATS junction imbalance OR rMATS event hits annotation that leafcutter clustered differently | Inspect locus in IGV; check Shiba on the same locus |
| leafcutter sig, rMATS not sig | Novel junction not in rMATS annotation; rMATS `--novelSS` may have missed it | Verify `--novelSS` was on; rerun if not |
| Both sig, opposite ΔPSI direction | Event class mismatch (rMATS calls SE positive, leafcutter sees A5SS shift in same cluster) | Manually map cluster topology to event class |
| Both sig, same direction | High-confidence call | Report; cross-validate with sashimi-plot |
| All tools null but biology suggests change | Underpowered design or wrong regime | Increase replicates; check whether outlier-splicing-detection regime applies |

**Operational rule:** for high-confidence reporting, require concordant detection in two tools from different algorithmic families (event-based + cluster-based, or LSV + isoform-based). Document both calls and any explainable disagreements.

## rMATS Output Columns Reference

| Column | Meaning |
|--------|---------|
| IJC_SAMPLE_1 / SJC_SAMPLE_1 | Comma-delimited inclusion / skipping junction counts per replicate, group 1 |
| IJC_SAMPLE_2 / SJC_SAMPLE_2 | Same for group 2 |
| IncFormLen / SkipFormLen | Effective lengths normalizing PSI for differential mapping opportunity |
| upstreamES/EE, downstreamES/EE | Flanking exon coordinates (genomic order; strand-agnostic in column meaning) |
| exonStart_0base / exonEnd | Cassette exon coordinates (0-based half-open) |
| PValue | LRT p-value of \|ΔPSI\| > cutoff |
| FDR | BH-adjusted PValue within event class |
| IncLevel1, IncLevel2 | Comma-delimited per-replicate PSI values |
| IncLevelDifference | mean(IncLevel1) - mean(IncLevel2); sign matches --b1 - --b2 order |

## Replicate Count and Power

| Design | Recommended tools | Expected power for ΔPSI=0.2 |
|--------|-------------------|------------------------------|
| n=2 vs n=2 | leafcutter or Shiba; **avoid SUPPA2** | Marginal; many real effects missed |
| n=3 vs n=3 | rMATS-turbo + leafcutter | Adequate at moderate coverage; standard |
| n=5 vs n=5 | rMATS or leafcutter, MAJIQ deltapsi | Good; recommended for publication |
| n=10+ vs n=10+ heterogeneous | MAJIQ-HET | Designed for this scale |
| Single patient vs n=20+ controls | leafcutterMD or FRASER2 | Outlier regime; see outlier-splicing-detection |

For an effect-size of |ΔPSI|=0.10 (typical biological signal), power generally requires n>=4 and >=20 junction reads per replicate. Below this, expect to miss most real changes.

## Significance and Effect-Size Thresholds

| Stringency | \|ΔPSI\| | FDR | Use case |
|------------|----------|-----|----------|
| Lenient | > 0.05 | < 0.10 | Discovery, exploratory, hypothesis generation |
| Standard | > 0.10 | < 0.05 | Publication; default reporting threshold |
| Stringent | > 0.20 | < 0.01 | Validation cohort, follow-up targets |

For MAJIQ: posterior probability `P(|ΔPSI| > 0.2) >= 0.95` is roughly equivalent to standard stringency. Always document tool, threshold, and rationale.

**Biologically meaningful ΔPSI varies by context:**
- A poison exon shift of |ΔPSI|=0.10 can halve functional protein (huge biology, modest number).
- A stoichiometric isoform shift of |ΔPSI|=0.10 may be physiologically silent.
- Therapeutic ASO target: SMA nusinersen aims for ΔPSI~+0.30 in SMN2 exon 7.

## Confounder Handling

**rMATS** does not natively accept arbitrary covariates. Workarounds:
1. **Stratification**: run rMATS within each batch separately and meta-analyze.
2. **PSI residuals (logit-transformed)**: PSI is bounded [0,1]; raw linear regression near the boundaries is biased. Logit-transform first, regress on confounders, then test residuals.
3. **Switch to leafcutter** (R function accepts `confounders` matrix; CLI accepts confounders as additional columns in the groups file).

```python
import numpy as np
import statsmodels.formula.api as smf

# logit-transform PSI before residualization (PSI is bounded [0,1])
eps = 1e-3
psi['logit_psi'] = np.log((psi['psi'].clip(eps, 1 - eps)) / (1 - psi['psi'].clip(eps, 1 - eps)))
psi['psi_resid'] = smf.ols('logit_psi ~ batch + RIN', data=psi).fit().resid
# then test psi_resid by group via Wilcoxon
```

**leafcutter** accepts confounders two ways:
- **R function**: `differential_splicing(counts, x, confounders=numeric_matrix)` accepts a numeric covariate matrix
- **CLI script**: `leafcutter_ds.R` reads confounders from **additional columns in the groups file** (3rd, 4th, ... columns), NOT from a `--confounders` flag

**MAJIQ** does not accept arbitrary confounders; use stratification or switch tool.

**Always check confounding before reporting:** PCA on PSI matrix; if PC1 separates by batch rather than group, the comparison is confounded.

## Multi-Group / Multi-Factor Designs

| Design | Approach |
|--------|----------|
| 3 groups (e.g. drug A, drug B, control) | Pairwise rMATS or leafcutter; OR limma/DESeq2 on logit-PSI matrix |
| Time-course (e.g. 0h, 6h, 24h) | DEXSeq on event counts with time as factor; or limma::lmFit on PSI matrix |
| 2x2 factorial (genotype × treatment) | DEXSeq with interaction term; rMATS pairwise on interaction subsets |
| Continuous covariate (dose, age) | limma::lmFit on logit-PSI ~ covariate |

For complex designs, custom regression on the PSI matrix is more flexible than rMATS/leafcutter pairwise.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `rMATS: numpy.AxisError` | rMATS version mismatch with numpy >=2.0 | Pin numpy<2.0 or update rMATS-turbo to >=4.3 |
| `leafcutter: zero variance in cluster` | Cluster has all-zero counts in a group | Pre-filter with `--min_samples_per_intron 5 --min_samples_per_group 3` |
| `MAJIQ: out of memory` | Default settings on >50-sample cohort | Use `--mem-profile` flag; chunk samples; consider HET for large cohorts |
| `SUPPA2: no events with sufficient coverage` | Salmon/kallisto TPM filter too strict upstream | Lower upstream TPM threshold; verify event annotations |
| `voila: missing splicegraph.zarr` (V3) or `splicegraph.sql` (V2; deprecated) | Forgot to keep build output directory | Re-run `majiq build`; output must persist for VOILA |
| `regtools: too many open files` | Many BAMs in one batch | `ulimit -n 4096` or batch in groups |

## Result Prioritization

**Goal:** Rank events by combined statistical and biological significance for follow-up.

**Approach:** Composite score combining FDR and effect size, then enrich for biology (RBP binding, NMD sensitivity, conservation, disease relevance).

```python
import pandas as pd
import numpy as np

sig['score'] = -np.log10(sig['FDR']) * sig['IncLevelDifference'].abs()
sig['exon_length'] = sig['exonEnd'] - sig['exonStart_0base']
sig['nmd_likely'] = (sig['exon_length'] % 3 != 0)
top_events = sig.nlargest(50, 'score')
```

Cross-reference top hits with:
- **eCLIP/ENCODE RBP target databases** (POSTAR3, oRNAment, RBP2GO) -> candidate trans-regulators
- **Disease-specific signatures**: SF3B1 cryptic 3'ss for MDS/CLL/UM; TDP-43 cryptic exons (UNC13A, STMN2) for ALS/FTD
- **Conservation**: VastDB cross-species PSI for evolutionary support
- **Splice-site predictions**: SpliceAI scores for the involved sites (see splice-variant-prediction)

## Common Pitfalls

- **Junction read imbalance** (cassette exon flanks have unequal mapping opportunity) inflates rMATS false positives; Shiba explicitly corrects this.
- **Comparing tool outputs naively** — MAJIQ posteriors and rMATS FDR are different scales; use threshold equivalents (P>0.95 ~ FDR<0.05 in many regimes) but confirm with simulation when reporting.
- **Forgetting NMD direction** — increased PSI of a poison exon decreases protein. Always check whether the alternative form is PTC-introducing using ORF-aware annotation.
- **Cryptic splicing in TDP-43 loss / SF3B1-mutant samples** — annotation-bound tools (rMATS, SUPPA2) miss these; need leafcutter or MAJIQ with denovo mode.
- **Forgetting strand** — wrong `--libType` halves usable junctions. Confirm with RSeQC `infer_experiment.py`.
- **Reporting one tool's call as ground truth** — discordance between rMATS and leafcutter is informative, not a problem to hide.
- **Skipping confounder check** — always run PCA on PSI matrix before final reporting.
- **Using empirical SUPPA2 at n<=3** — calibration collapses; use classical mode or different tool.

## Related Skills

- splicing-quantification - PSI estimation per event; foundational
- splicing-qc - Run BEFORE differential to verify library, depth, strandedness; avoid downstream surprises
- isoform-switching - DTU framework with NMD/ORF/domain consequences; complementary to event-level
- sashimi-plots - Visualize differential events for QC and reporting
- outlier-splicing-detection - Single-sample-vs-cohort regime (FRASER2/DROP); use when not 2-group
- splice-variant-prediction - SpliceAI / Pangolin for variant-driven mechanistic explanation of differential events
- long-read-splicing - Differential analysis from full-length isoforms; use when short-read insufficient
- read-alignment/star-alignment - STAR 2-pass cohort-style required upstream

## References

- Shen et al 2014 *PNAS* - rMATS original
- Wang et al 2024 *Nat Protoc* - rMATS-turbo
- Li et al 2018 *Nat Genet* - leafcutter (Dirichlet-multinomial GLM)
- Buen Abad Najar et al 2025 *bioRxiv* - LeafCutter2 (NMD-aware unproductive splicing)
- Vaquero-Garcia et al 2016 *eLife* - MAJIQ LSV framework
- Vaquero-Garcia et al 2023 *Nat Commun* - MAJIQ-HET heterogeneity module
- Aicher, Slaff, Jewell, Barash 2024 *bioRxiv* - MAJIQ V3
- Trincado et al 2018 *Genome Biol* - SUPPA2
- Kubota et al 2025 *NAR* - Shiba (junction-imbalance correction)
- Olofsson et al 2023 *Biochem Biophys Res Commun* 653:31-37 - benchmark across tools
- Tran et al 2025 *WIREs RNA* - methodology review
- Brown et al 2022 *Nature* - UNC13A cryptic exon (TDP-43 / ALS)
- Klim et al 2019 *Nat Neurosci* - STMN2 cryptic splicing (ALS)
- Darman et al 2015 *Cell Rep* - SF3B1 cryptic 3'ss
