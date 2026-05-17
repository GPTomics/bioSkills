---
name: bio-atac-seq-enhancer-gene-linking
description: Predict enhancer-gene regulatory connections from ATAC-seq using ABC, ENCODE-rE2G, HiChIP, or Cicero. Use when linking distal enhancers to target genes, choosing between contact-aware (ABC, ENCODE-rE2G), accessibility-only (Cicero), and orthogonal (HiChIP H3K27ac, EpiMap) approaches, validating predictions against CRISPRi-FlowFISH gold-standard, or building cell-type-specific regulatory maps for fine-mapping or therapeutic target discovery.
tool_type: mixed
primary_tool: ABC-Enhancer-Gene-Prediction
---

## Version Compatibility

Reference examples tested with: ABC-Enhancer-Gene-Prediction 0.2.2+ (Engreitz lab), ENCODE-rE2G v1.0+ (2024 release), Cicero 1.20+, GenomicInteractions 1.36+, FitHiChIP 9.1+, hicpro 3.1+, FAN-C 0.9+, MACS3 3.0+, samtools 1.19+, bedtools 2.31+.

Verify before use:
- CLI: `<tool> --version` then `<tool> --help` to confirm flags
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws unexpected errors, introspect the installed package and adapt rather than retrying.

# Enhancer-Gene Linking

**"Which gene does this distal accessible region regulate?"** -> Predict the enhancer's target gene using a model that combines accessibility activity, 3D contact frequency, and (optionally) sequence-based chromatin predictions. Output is a per-(enhancer, gene) score that can be thresholded for high-confidence calls.

- CLI: ABC pipeline (`run.neighborhoods.py`, `predict.py` from Engreitz lab)
- CLI: ENCODE-rE2G (Snakemake-based; ENCODE 4 standard 2024)
- R: Cicero (ATAC-only; covered in atac-seq/co-accessibility)
- CLI: FitHiChIP / hichipper for HiChIP H3K27ac loops
- Database: EpiMap (Boix 2021), GeneHancer, FANTOM5 (pre-computed reference)

ABC and ENCODE-rE2G are the canonical predictors when Hi-C/Micro-C data is available. Cicero is the ATAC-only fallback. CRISPRi-FlowFISH (Fulco 2019) is the gold-standard experimental validation.

## Algorithmic Taxonomy

| Method | Inputs | Mathematics | Strength | Fails when |
|--------|--------|-------------|----------|------------|
| ABC (Fulco 2019, Nasser 2021) | ATAC + H3K27ac + Hi-C/Micro-C | ABC = (Activity_E x Contact_E,G) / sum_e(Activity_e x Contact_e,G); threshold typically >= 0.02 | Mechanistically grounded; published gold-standard for human cell lines | Requires matched Hi-C / Micro-C; cell-type-specific; default contact uses average across 5 ENCODE cell types if Hi-C not available |
| ENCODE-rE2G (ENCODE 2024) | ATAC + H3K27ac + (Hi-C optional) | Logistic regression trained on CRISPRi-FlowFISH ground truth; uses ABC features + sequence features + distance | ENCODE 4 standard; pre-trained models for many cell types | Pre-trained models only available for ENCODE cell types; retraining requires CRISPRi data |
| Cicero (Pliner 2018) | scATAC peak-cell matrix | Graphical lasso on metacell co-accessibility | ATAC-only; works without Hi-C | Less concordant with Hi-C than ABC; cis-distance-limited; alpha-sensitive |
| HiChIP H3K27ac + FitHiChIP | H3K27ac HiChIP | Statistically significant loops at FDR < 0.05 | Direct experimental loop measurement; cell-type-specific; orthogonal to ATAC | Requires HiChIP wet-lab; only captures loops within HiChIP resolution (~10 kb) |
| Hi-C + HiCCUPS | Bulk Hi-C | Fold-enrichment loop calling | Most-validated 3D contact method | Resolution typically 5-25 kb; misses sub-loop fine structure |
| Capture Hi-C / PCHi-C (CHiCAGO) | Promoter Capture Hi-C | Asymptotic CHiCAGO score | High-resolution promoter-anchored | Wet-lab cost; promoter capture only |
| EpiMap (Boix 2021) reference | None (pre-computed lookup) | Bulk-derived enhancer-gene predictions in 833 epigenomes | Fast, comprehensive | Cell-type-agnostic for tissues outside the reference set |
| GeneHancer / FANTOM5 (legacy) | None (pre-computed lookup) | Pre-computed; varied methods per database | Comprehensive lookup; widely cited | Older; less reliable than ABC for cell-type-specific |

Methodology evolves; verify against current Engreitz lab releases (ABC), ENCODE 4 publications (ENCODE-rE2G), and Mumbach 2017 (HiChIP) before locking pipelines.

## ABC Mathematics

For each candidate (enhancer E, gene G) pair within the cis window (default 5 Mb):

```
ABC(E -> G) = Activity_E * Contact_E,G / sum_{all e in window}(Activity_e * Contact_e,G)
```

- **Activity_E** = ATAC reads at E * H3K27ac reads at E (geometric mean of normalized signals; reflects "enhancer strength")
- **Contact_E,G** = Hi-C/Micro-C contact frequency from E to G's TSS (after distance-correction)
- **Window** = +/- 5 Mb cis (default; ENCODE-rE2G uses 1 Mb)

Threshold typical: ABC >= 0.02 for high-confidence; >= 0.01 for exploratory.

When Hi-C is unavailable, ABC uses an "average contact" derived from 5 ENCODE Hi-C cell types as proxy. Performance degrades but it's not zero (Fulco 2019).

## ENCODE-rE2G Differences from ABC

ENCODE-rE2G (2024) is a reformulation:

- **Logistic regression** trained on CRISPRi-FlowFISH ground truth (~10 cell types)
- **Features:** ABC score components + 3D contact + distance + activity ratios
- **Multiple feature configurations:** "abc-features", "no-hic-features" for cells without 3D data
- **Output:** Per-pair probability of regulatory connection
- **Pre-trained models** for ENCODE cell lines; logistic params vary by cell type

ENCODE-rE2G generally outperforms ABC at CRISPRi recall, especially at modest distances (50-500 kb). For ENCODE cell types, prefer ENCODE-rE2G; for novel cell types, ABC remains the default.

## Per-Tool Failure Modes

### ABC -- Wrong cell-type-matched Hi-C

**Trigger:** Using K562 Hi-C contact when actual cell type is GM12878.

**Mechanism:** Contact frequencies differ across cell types at compartment and TAD boundaries; using mismatched Hi-C produces wrong ABC scores.

**Symptom:** ABC predictions concentrate at known K562-specific loci even when ATAC data is from GM12878.

**Fix:** Use cell-type-matched Hi-C or Micro-C. If unavailable, ABC's "average HiC" (5-cell-type pooled) is the documented fallback with acknowledged degradation. Document the proxy in methods.

### ABC -- H3K27ac normalization

**Trigger:** H3K27ac ChIP-seq with different sequencing depth than ATAC.

**Mechanism:** ABC's "Activity" is the geometric mean of accessibility and H3K27ac signals; both must be normalized to the same scale.

**Symptom:** Activity scores skewed; some peaks have very high activity from H3K27ac alone, others from ATAC alone.

**Fix:** Normalize both signals to reads-per-million in peaks (RPM-IP) before combining. Use ABC's `--qnorm` flag with a quantile-normalization reference file (e.g. `--qnorm src/EnhancersQNormRef.K562.txt` from the ABC repo).

### ENCODE-rE2G -- Cell type not in pre-trained set

**Trigger:** Running pre-trained model on a primary cell type not in CRISPRi training.

**Mechanism:** Logistic regression coefficients learned from ENCODE cell types may not transfer to primary tissues.

**Fix:** Use the closest ENCODE cell type (myeloid lineage -> K562; lymphoid -> GM12878; hepatic -> HepG2). Document the proxy. For high-stakes use, custom retraining requires CRISPRi-FlowFISH data.

### Cicero -- No Hi-C concordance benchmark

**Trigger:** Reporting Cicero connections as enhancer-gene calls without external validation.

**Mechanism:** Cicero is statistical co-accessibility; correlation with Hi-C 3D contacts is ~30-50%. Many strong Cicero connections are NOT Hi-C-validated.

**Fix:** When Hi-C is available, cross-validate; report both. When only ATAC, use Cicero with the explicit caveat that connections are co-accessibility hypotheses, not contact predictions.

### HiChIP -- Loop calling threshold

**Trigger:** Default FitHiChIP at FDR < 0.05.

**Mechanism:** HiChIP loops are abundant (10k-100k per dataset); FDR alone produces a long tail of weak loops.

**Fix:** Threshold at FDR < 0.05 AND number of contacts per loop >= 5; or use the top N most significant where N = expected number of loops based on cell type.

### EpiMap / GeneHancer -- Cell-type-agnostic limitation

**Trigger:** Using EpiMap or GeneHancer pre-computed pairs for a specific cell type.

**Mechanism:** These references aggregate across many tissues / experiments; cell-type-specific connections are diluted.

**Fix:** Use as a baseline / sanity check, not as the primary call. ABC or ENCODE-rE2G in the actual cell type is preferred.

## Decision Tree by Available Data

| Available data | Recommended method |
|---------------|--------------------|
| ATAC + H3K27ac + matched Hi-C/Micro-C | ABC or ENCODE-rE2G (with cell-type-matched contact) |
| ATAC + H3K27ac, no Hi-C | ABC with average HiC fallback; or ENCODE-rE2G `no-hic` model |
| ATAC only, no H3K27ac | Cicero (atac-seq/co-accessibility); ABC with synthetic activity |
| ATAC + H3K27ac HiChIP | FitHiChIP loops + ABC; intersect for high confidence |
| Multiome (ATAC + RNA same cell) | LinkPeaks (Signac) for direct correlation; SCENIC+ for TF networks |
| ENCODE cell type | Pre-computed ENCODE-rE2G predictions (download) |
| Tissue with limited public data | ABC + acknowledge proxy; do not rely on EpiMap |
| Multi-cell-type scATAC | scBasset (atac-seq/deep-learning-atac) for sequence-based per-cell |
| Want experimental validation | CRISPRi-FlowFISH design; use predictions as targeted hypotheses |

## ABC Standard Pipeline

**Goal:** Compute per-(enhancer, gene) ABC scores combining ATAC accessibility, H3K27ac activity, and Hi-C contact.

**Approach:** Build RPGC-normalized ATAC and H3K27ac bigWigs, define non-promoter candidate enhancers from ATAC peaks, run ABC neighborhoods to compute per-candidate activity, then run ABC predict against a Hi-C contact matrix and threshold the per-pair ABC score.

```bash
# 1. Generate ATAC and H3K27ac signal tracks (BAM -> bigWig)
bamCoverage --bam atac.bam --outFileName atac.bw --binSize 50 --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398 --numberOfProcessors 8

# 2. Define enhancer candidates (typically MACS narrowPeak from ATAC)
# Filter to non-promoter regions
bedtools intersect -v -a atac_peaks.narrowPeak -b promoter_regions.bed > candidate_enhancers.bed

# 3. Run ABC neighborhoods (compute Activity per candidate)
# Script path: legacy ABC = src/run.neighborhoods.py; Snakemake-based modern = workflow/scripts/run.neighborhoods.py
python /path/ABC-Enhancer-Gene-Prediction/workflow/scripts/run.neighborhoods.py \
    --candidate_enhancer_regions candidate_enhancers.bed \
    --genes refseq_protein_coding.bed \
    --H3K27ac h3k27ac.bw \
    --DHS atac.bw \
    --chrom_sizes hg38.chrom.sizes \
    --ubiquitously_expressed_genes Genes.ubiquitously_expressed.txt \
    --cellType MyCellType \
    --outdir abc_out/

# 4. Run ABC predictions (Activity * Contact)
python /path/ABC-Enhancer-Gene-Prediction/workflow/scripts/predict.py \
    --enhancers abc_out/EnhancerList.txt \
    --genes abc_out/GeneList.txt \
    --HiCdir hic_data/ \
    --hic_resolution 5000 \
    --score_column ABC.Score \
    --threshold 0.02 \
    --cellType MyCellType \
    --outdir abc_out/Predictions/

# Output: EnhancerPredictionsAllPutative.txt with per-pair ABC scores
```

ABC.Score >= 0.02 is the standard threshold; >= 0.04 for stringent reporting.

## ENCODE-rE2G

```bash
# Snakemake-based; clone the ENCODE-rE2G repo
git clone https://github.com/EngreitzLab/ENCODE_rE2G
cd ENCODE_rE2G

# Edit config.yaml with cell type, ATAC, H3K27ac paths
# Run with appropriate model (cell-type-matched logistic regression weights)
snakemake --use-conda --cores 16 \
    --config cell_type=K562 \
             atac_bw=atac.bw \
             h3k27ac_bw=h3k27ac.bw \
             hic_directory=hic_data/

# Output: scored_predictions.tsv.gz with per-pair scores and binarized predictions
```

Pre-trained models are at https://github.com/EngreitzLab/ENCODE_rE2G/tree/main/models. Choose by tissue similarity if exact cell type not present.

## CRISPRi-FlowFISH Validation Framework

CRISPRi-FlowFISH (Fulco 2019) is the experimental gold-standard:
1. Design sgRNAs tiling each candidate enhancer
2. Transduce CRISPRi-expressing cells; FACS by gene expression (FlowFISH for endogenous; reporter for ectopic)
3. Sequence sgRNAs in low- vs high-expression bins; compute log2 enrichment per sgRNA
4. Significance: meta-test across sgRNAs in same enhancer

A 2-fold expression decrease (p < 0.05) confirms the enhancer regulates the gene.

For predictions to be publication-grade, ENCODE 4 expects:
- **Test set sensitivity / specificity** against published CRISPRi-FlowFISH catalogs (Fulco 2019: K562; Gasperini 2019; Schraivogel 2020 multi-cell-type)
- **Effect-size correlation** between predicted score and observed expression effect
- **Distance bias check** (predictors over-rank close-distance pairs)

## Reconciling Methods

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ABC and ENCODE-rE2G disagree | Different feature weighting; different training distributions | Both valid; report intersection as high-confidence |
| ABC strong, Cicero weak | Co-accessibility sparse for that cell type | Trust ABC if Hi-C is matched |
| HiChIP loop with no ABC prediction | Loop is below ABC threshold; or peak set too narrow | Lower threshold or expand candidate enhancers |
| ENCODE-rE2G high probability, no CRISPRi support | Could be context-dependent biology or false positive | Prioritize for follow-up; not a publishable claim alone |
| EpiMap pair not in ABC | Pre-computed reference is cell-type-aggregated | Use ABC for cell-type-specific |

**Operational rule for high-confidence reporting:** Predictions used for therapeutic target nomination must be (a) above ABC >= 0.02 OR ENCODE-rE2G >= 0.5, AND (b) consistent across two methods (ABC + ENCODE-rE2G or ABC + HiChIP), AND (c) validated experimentally (CRISPRi-FlowFISH preferred). Single-method high-score predictions are exploratory hypotheses.

## Combining Multiple Predictions

**Goal:** Build a high-confidence enhancer-gene set by intersecting ABC, ENCODE-rE2G, and HiChIP evidence.

**Approach:** Load each method's output, merge ABC and ENCODE-rE2G on enhancer-gene pair above per-method thresholds, then flag pairs with HiChIP loop support for triple-method evidence.

```python
import pandas as pd
abc = pd.read_csv('abc_predictions.tsv', sep='\t')
re2g = pd.read_csv('encode_re2g.tsv.gz', sep='\t')
hichip = pd.read_csv('fithichip_loops.bedpe', sep='\t', header=None,
                     names=['chr1','s1','e1','chr2','s2','e2','name','score'])

# High-confidence intersection
high_conf = abc[abc['ABC.Score'] >= 0.02].merge(
    re2g[re2g['re2g_score'] >= 0.5],
    on=['enhancer_id', 'gene'])

# Add HiChIP support flag
hichip_anchors = ...    # extract enhancer/gene pairs from HiChIP loops
high_conf['hichip_support'] = high_conf['enhancer_id'].isin(hichip_anchors)
```

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| ABC predictions concentrate at TSSs | Did not exclude promoter regions from candidates | Pre-filter `bedtools intersect -v` against promoters |
| Activity scores all very small | H3K27ac or ATAC bigWig in wrong scale | Use RPGC normalization |
| ENCODE-rE2G model not converging | Pre-trained model loaded for wrong cell type | Match training cell type via `cell_type` config |
| Cicero connections used as enhancer-gene calls | Method confusion (co-accessibility vs contact) | Switch to ABC if Hi-C available; or document as co-accessibility hypothesis |
| Hi-C resolution too coarse | Default 25 kb resolution masks fine ABC structure | Use 5 kb or 10 kb if Micro-C available |
| FitHiChIP many loops, low specificity | Default FDR alone | Add contact count threshold; or use ENCODE-rE2G HiChIP-trained model |
| GeneHancer / FANTOM5 used as primary call | Cell-type-agnostic limitation | Use as baseline only |

## References

- Fulco CP et al 2019 Nat Genet 51:1664 (ABC; CRISPRi-FlowFISH validation)
- Nasser J et al 2021 Nature 593:238 (ABC genome-wide application)
- ENCODE Project Consortium 2024 (ENCODE-rE2G)
- Mumbach MR et al 2017 Nat Genet 49:1602 (HiChIP H3K27ac)
- Bhattacharyya S et al 2019 Bioinformatics 35:5012 (FitHiChIP)
- Boix CA et al 2021 Nature 590:300 (EpiMap reference)
- Gasperini M et al 2019 Cell 176:377 (CRISPRi at scale)
- Schraivogel D et al 2020 Nat Methods 17:629 (multi-cell-type CRISPRi)
- Pliner HA et al 2018 Mol Cell 71:858 (Cicero co-accessibility)

## Related Skills

- atac-seq/co-accessibility - Cicero (ATAC-only enhancer-promoter inference)
- atac-seq/atac-peak-calling - Generate enhancer candidates
- atac-seq/consensus-peakset - Fixed-width enhancer regions
- atac-seq/deep-learning-atac - chromBPNet variant effect at predicted enhancers
- atac-seq/single-cell-atac - Per-cell-type scATAC inputs
- hi-c-analysis/loop-calling - Hi-C / Micro-C contact prediction
- hi-c-analysis/contact-pairs - Hi-C / Micro-C input
- chip-seq/peak-calling - H3K27ac peaks
- gene-regulatory-networks/scenic-regulons - Downstream TF -> target inference
