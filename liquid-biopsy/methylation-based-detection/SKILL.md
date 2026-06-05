---
name: bio-methylation-based-detection
description: Detects cancer and infers tissue-of-origin from cfDNA methylation by choosing conversion chemistry (bisulfite vs EM-seq vs TAPS vs cfMeDIP), calling read-level methylation haplotypes rather than averaged beta values, and deconvolving a hematopoietic-dominated cfDNA mixture against a methylation atlas via NNLS/quadratic programming. Encodes the GRAIL/CCGA thesis that thousands of tissue-specific markers make methylation outperform sparse mutations for multi-cancer early detection (MCED) and localization, and that single concordantly-methylated fragments give ppm-level sensitivity. Uses MethylDackel for extraction (mbias-then-extract), MEDIPS/QSEA for enrichment data, scipy.optimize.nnls for deconvolution. Use when building an MCED or methylation-MRD assay, picking a conversion chemistry for low-input plasma, or deconvolving tissue-of-origin from cfDNA.
tool_type: mixed
primary_tool: MethylDackel
---

## Version Compatibility

Reference examples tested with: MethylDackel 0.6+, Bismark 0.24+, numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: MethylDackel `extract` bedGraph column order is fixed (chrom / start / end / methylation-% rounded to integer / count-methylated / count-unmethylated); always run `MethylDackel mbias` first and feed its suggested `--OT/--OB` trimming into `extract`. cfMeDIP data are coverage, not conversion — do not feed them into per-CpG bisulfite pipelines.

# Methylation-Based Detection

**"Detect cancer and find where it came from using cfDNA methylation"** -> Call per-CpG (and read-level) methylation from converted plasma DNA, then deconvolve tissue-of-origin against a reference atlas.
- CLI: `MethylDackel extract` for per-CpG methylation from bisulfite/EM-seq BAMs
- CLI: `MethylDackel mbias` to choose strand-specific trimming before extraction
- R: MEDIPS / QSEA for cfMeDIP enrichment (coverage, not conversion)
- Python: `scipy.optimize.nnls` for atlas-based tissue deconvolution

## The Single Most Important Modern Insight -- methylation beats mutations, and read-level haplotypes beat averaged beta

Methylation is the right altitude for multi-cancer early detection (MCED) and tissue-of-origin (TOO) in a single assay because the genome carries thousands of stable, cell-type-specific differentially methylated regions, whereas somatic mutations are sparse, recurrent only at a few driver loci, and carry no tissue label. On the same cfDNA inside CCGA, the methylation assay outperformed WGS-SNV/CNV approaches, which is why GRAIL down-selected to a targeted methylation panel (Liu 2020). One panel answers both "is there cancer?" and "where is it?" — mutations answer neither well.

The sensitivity engine is read-level, not site-level. Averaging beta across reads at a CpG discards phasing. A single tumor-derived fragment that is concordantly methylated across the k CpGs of a methylation haplotype block (Guo 2017) has a background probability of roughly p^k of arising from the hematopoietic ocean; for a block of 5-8 CpGs that is small enough that ONE such fragment is strong evidence, independent of tumor fraction. Per-site beta dilutes that signal into sampling noise and clonal-hematopoiesis variance and has essentially no power at parts-per-million tumor fraction. The correct primitive for detection is molecule counting over haplotype blocks, not site averaging.

## Conversion Chemistry Tradeoffs

The conversion step is chosen on the worst possible substrate — already-fragmented, low-input plasma DNA — so destructiveness is load-bearing, not a footnote.

| Method | Destructiveness | Min input | Base resolution | 5mC readout | Key bias / caveat |
|--------|-----------------|-----------|-----------------|-------------|-------------------|
| Bisulfite (WGBS/targeted) | Severe — depurinates/fragments, >90% loss possible | High (degradation eats low input) | Yes | C->T after conversion; remaining C = methylated | Complexity collapse on cfDNA; incomplete conversion -> false methylation; GC/coverage bias |
| EM-seq (Vaisvila 2021) | Much gentler — enzymatic, no chemical fragmentation | Picograms demonstrated | Yes | Same C->T readout as bisulfite, milder | Reads 5mC+5hmC together unless separated; APOBEC over/under-deamination edge cases |
| TAPS (Liu 2019) | Non-destructive — mild | Low / cfDNA-friendly | Yes | Direct: 5mC/5hmC -> T, unmethylated C untouched | Only a few % of Cs convert -> preserves complexity, lower seq cost; needs TET + pyridine borane |
| cfMeDIP-seq (Shen 2018) | No conversion (antibody enrichment) | Very low (>=5-10 ng) | No — region/enrichment-level only | Antibody pulls down methylated fragments | CpG-density bias; no single-CpG quantitation; needs MEDIPS/QSEA density modeling |

Bisulfite is gold standard for cell-line gDNA, not for low-input fragmented plasma; EM-seq and TAPS exist precisely to recover ctDNA molecules bisulfite destroys. Neither bisulfite nor EM-seq separates 5mC from 5hmC without added oxBS/TAB steps; TAPS variants (TAPSbeta, CAPS) can split the marks. cfMeDIP coverage is enrichment, not quantitation — density bias must be modeled before any absolute-methylation claim.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| MCED + tissue-of-origin in one assay | Targeted methylation panel + atlas deconvolution | Thousands of tissue markers carry both cancer and organ signal (Liu 2020; Loyfer 2023) |
| Ultra-low-input plasma, genome-wide | cfMeDIP-seq | Antibody enrichment works at ng-to-low input where conversion destroys the library (Shen 2018) |
| Need base resolution at low input | EM-seq or TAPS, not bisulfite | Gentle/non-destructive conversion preserves complexity bisulfite collapses (Vaisvila 2021; Liu 2019) |
| MRD / ppm-level detection | Read-level haplotype counting over pre-defined blocks | One concordant fragment is decisive; averaged beta has no power at low tumor fraction (Guo 2017) |
| Absolute methylation level from enrichment data | QSEA (Bayesian density + CNV + TMM) | Converts cfMeDIP coverage to BS-comparable values; MEDIPS gives differential coverage only (Lienhard 2017) |

Methodology evolves; verify current atlas versions and panel-marker coverage against live tool docs before committing — atlas markers are platform-specific and do not transfer across assays.

## Extract Per-CpG Methylation with MethylDackel

**Goal:** Produce per-CpG methylation calls from a bisulfite/EM-seq cfDNA BAM, with end-repair artifacts trimmed.

**Approach:** Run `mbias` first to read the suggested strand-specific trimming, then `extract` with `--mergeContext` and that trimming so each CpG is one row; parse the fixed 6-column bedGraph.

```bash
# Step 1: choose trimming. mbias prints a suggestion like --OT 2,0,0,98 and writes M-bias SVGs.
MethylDackel mbias ref.fa sample.bam sample_mbias

# Step 2: extract per-CpG (one row per CpG) with the suggested trimming.
MethylDackel extract ref.fa sample.bam -o sample --mergeContext --minDepth 1 --OT 2,0,0,98
# Output sample_CpG.bedGraph columns (fixed order):
#   chrom  start  end  methylation%(integer, rounded)  count_methylated  count_unmethylated
```

Add `--CHG --CHH` only to audit non-CpG methylation (a conversion-failure check); CpG is the default context. Low per-CpG depth gates can erase cfDNA signal — prefer region/molecule aggregation over a high `--minDepth`.

## Tissue-of-Origin Deconvolution Against an Atlas

**Goal:** Attribute cfDNA to its cell types of origin and surface a solid-tissue coefficient elevated above the hematopoietic baseline.

**Approach:** Model the observed methylation vector m ~ A*w with atlas A (rows = markers, cols = cell types), solve for non-negative mixing fractions w with NNLS over atlas-covered markers, then renormalize so the fractions sum to one.

```python
from scipy.optimize import nnls

def deconvolve_tissue(sample_beta, atlas):
    'sample_beta: Series indexed by marker; atlas: DataFrame markers x cell_types.'
    markers = sample_beta.index.intersection(atlas.index)
    w, _ = nnls(atlas.loc[markers].values, sample_beta.loc[markers].values)
    w = w / w.sum()
    return dict(zip(atlas.columns, w))
```

The simplex constraint (w >= 0, sum w = 1) is mandatory — unconstrained regression gives nonsense fractions (Moss 2018). Use Loyfer 2023's fragment-level WGBS atlas (39 cell types from 205 healthy samples) where the assay covers its markers; a generic atlas does not transfer, because WGBS-fragment markers differ from 450K/EPIC probes and from capture-panel coverage.

## Region-Level DMR Discovery

**Goal:** Define a discriminating panel by finding regions (not single CpGs) that separate cancer from normal cfDNA.

**Approach:** Aggregate per-CpG beta into pre-defined regions/blocks, test cancer vs normal per region, and control FDR with Benjamini-Hochberg specified explicitly — naming method='fdr_bh' rather than relying on the statsmodels default ('hs', Holm-Sidak).

```python
from scipy import stats
from statsmodels.stats.multitest import multipletests

def region_dmrs(cancer, normal, region_col='region'):
    'cancer/normal: long DataFrames with [region_col, beta]; one row per sample-region.'
    out = []
    for region, c in cancer.groupby(region_col)['beta']:
        n = normal.loc[normal[region_col] == region, 'beta'].dropna()
        c = c.dropna()
        if len(c) < 3 or len(n) < 3:
            continue
        _, p = stats.mannwhitneyu(c, n, alternative='two-sided')
        out.append((region, c.mean() - n.mean(), p))
    import pandas as pd
    res = pd.DataFrame(out, columns=['region', 'delta_beta', 'pvalue'])
    res['fdr'] = multipletests(res['pvalue'], method='fdr_bh')[1]
    return res.sort_values('fdr')
```

## cfMeDIP Enrichment Analysis

**Goal:** Get density-corrected differential methylation from antibody-enrichment coverage rather than conversion data.

**Approach:** Use MEDIPS (Lienhard 2014) for CpG-density-corrected differential coverage, or QSEA (Lienhard 2017) when absolute, BS-comparable methylation levels are needed — QSEA adds a Bayesian CpG-density model, CNV correction, and TMM effective-library-size normalization. Both are R/Bioconductor; do not pass cfMeDIP coverage through MethylDackel.

## Per-Method Failure Modes

### Averaged beta wastes the read-level signal
**Trigger:** Reporting region beta means / per-CpG DMR t-tests for an MCED or MRD assay. **Mechanism:** Averaging across reads discards fragment-level concordance, the exact signal that lets one tumor fragment be called. **Symptom:** No power at low tumor fraction despite deep coverage. **Fix:** Count concordantly-methylated molecules over haplotype blocks; reserve beta for discovery and QC.

### Per-CpG t-tests + naive BH are the wrong altitude
**Trigger:** Genome-wide per-CpG Welch t-tests with plain Benjamini-Hochberg. **Mechanism:** ~28M correlated CpGs violate BH independence (anticonservative) and per-site estimates are coverage-starved. **Symptom:** Inflated "DMR" lists that do not replicate. **Fix:** Region/block methods (dmrseq/metilene/methylKit windows) with permutation or correlation-aware FDR.

### Bisulfite degradation lowers complexity
**Trigger:** WGBS on low-input plasma. **Mechanism:** Chemical depurination/fragmentation destroys input, collapsing unique molecules. **Symptom:** Low library complexity, duplicate-heavy, lost ctDNA molecules. **Fix:** EM-seq or TAPS; track conversion completeness via CHH methylation.

### WBC background swamps the tumor coefficient
**Trigger:** Deconvolving with a mis-specified or unmatched hematopoietic reference. **Mechanism:** >90% of cfDNA is leukocyte/megakaryocyte-derived; reference error leaks variance into the small tumor term. **Symptom:** False or unstable TOO; a methylation analog of CHIP (clonal hematopoiesis/age/inflammation shifts the WBC methylome). **Fix:** Age/condition-matched background, fine-grained atlas, treat tumor as a small residual.

### cfMeDIP density bias / no base resolution
**Trigger:** Reading cfMeDIP coverage as methylation level, or running it through a per-CpG pipeline. **Mechanism:** Antibody enriches CpG-dense regions; there is no single-CpG quantitation. **Symptom:** Apparent hypermethylation tracking CpG density, not biology. **Fix:** Model density with MEDIPS coupling factor or QSEA's Bayesian model.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Specificity 99.3%, sensitivity 54.9%, TOO 93% (among detected) | Liu 2020 *Ann Oncol* 31(6):745 | CCGA2 targeted-methylation operating point; screening fixes high specificity and accepts modest sensitivity |
| Specificity 99.5%, sensitivity 51.5%; stage I 16.8% -> IV 90.1%; TOO 88.7% | Klein 2021 *Ann Oncol* 32(9):1167 | CCGA3 clinical validation; sensitivity is stage- and tumor-type-dominated, never quote it as uniform |
| cfMeDIP input >= 5-10 ng | Shen 2018 *Nature* 563:579 | Enrichment works below conversion-assay input floors but still needs ng-scale material |
| Deconvolution: w >= 0 and sum w = 1 | Moss 2018 *Nat Commun* 9:5068; Loyfer 2023 *Nature* 613:355 | Simplex constraint mandatory; unconstrained regression yields nonsense fractions |
| Avoid high per-CpG `--minDepth` for cfDNA | MethylDackel docs; community | Per-CpG depth gates erase low-coverage cfDNA signal; aggregate over regions/molecules instead |

## References

- Liu MC, et al. 2020. Sensitive and specific multi-cancer detection and localization using methylation signatures in cell-free DNA. *Annals of Oncology* 31(6):745-759. — CCGA2 targeted-methylation MCED; spec 99.3%, sens 54.9%, TOO 93%.
- Klein EA, et al. 2021. Clinical validation of a targeted methylation-based multi-cancer early detection test using an independent validation set. *Annals of Oncology* 32(9):1167-1177. — CCGA3 validation; stage-dependent sensitivity.
- Shen SY, et al. 2018. Sensitive tumour detection and classification using plasma cell-free DNA methylomes. *Nature* 563(7732):579-583. — cfMeDIP-seq.
- Moss J, et al. 2018. Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease. *Nature Communications* 9:5068. — NNLS atlas deconvolution.
- Loyfer N, et al. 2023. A DNA methylation atlas of normal human cell types. *Nature* 613(7943):355-364. — 39 cell types, 205 samples; fragment-level WGBS atlas.
- Liu Y, et al. 2019. Bisulfite-free direct detection of 5-methylcytosine and 5-hydroxymethylcytosine at base resolution (TAPS). *Nature Biotechnology* 37(4):424-429.
- Guo S, et al. 2017. Identification of methylation haplotype blocks aids in deconvolution of heterogeneous tissue samples and tumor tissue-of-origin mapping from plasma DNA. *Nature Genetics* 49(4):635-642. — methylation haplotype blocks.
- Vaisvila R, et al. 2021. Enzymatic methyl sequencing detects DNA methylation at single-base resolution from picograms of DNA (EM-seq). *Genome Research* 31(7):1280-1289. — page span verified from DOI 10.1101/gr.266551.120; confirm against the published PDF if citing the exact pages.
- Lienhard M, et al. 2014. MEDIPS: genome-wide differential coverage analysis of sequencing data derived from DNA enrichment experiments. *Bioinformatics* 30(2):284-286.
- Lienhard M, et al. 2017. QSEA — modelling of genome-wide DNA methylation from sequencing enrichment experiments. *Nucleic Acids Research* 45(6):e44.

## Related Skills

- cfdna-preprocessing - conversion chemistry and library choices upstream
- fragment-analysis - orthogonal genome-wide cfDNA signal
- analytical-validation - read-level detection framed as a limit-of-detection problem
- methylation-analysis/bismark-alignment - bisulfite read alignment
- methylation-analysis/dmr-detection - region-level differential methylation statistics
