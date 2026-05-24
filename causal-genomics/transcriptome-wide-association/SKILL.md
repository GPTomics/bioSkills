---
name: bio-causal-genomics-transcriptome-wide-association
description: Performs gene-level association from GWAS summary statistics via genetically predicted tissue expression using FUSION, PrediXcan, S-PrediXcan, S-MultiXcan, UTMOST, MOSTWAS, kTWAS, EpiXcan, TIGAR-V2, and probabilistic fine-mapping with FOCUS and MA-FOCUS. Use when running TWAS from GWAS sumstats, prioritising candidate causal genes from a GWAS lead locus, picking single-tissue vs cross-tissue models, identifying LD-induced TWAS false positives, choosing ancestry-matched prediction weights, fine-mapping co-regulated TWAS hits, or triangulating TWAS with cis-eQTL Mendelian randomization and colocalization to nominate a causal gene.
tool_type: mixed
primary_tool: FUSION
---

## Version Compatibility

Reference examples tested with: FUSION (head of `gusevlab/fusion_twas`, scripts dated 2023+), MetaXcan / S-PrediXcan / S-MultiXcan 0.7.5+ (`hakyimlab/MetaXcan`), PrediXcan model files from PredictDB (GTEx v8 elastic-net + MASHR), UTMOST (head of `Joker-Jerome/UTMOST`), pyfocus 0.8+ (`bogdanlab/focus`), MA-FOCUS (head of `mancusolab/ma-focus`), TIGAR-V2 (head of `yanglab-emory/TIGAR`), PLINK 1.9 + PLINK 2.0, R 4.3+, Python 3.9-3.11.

Before using code patterns, verify installed versions match. If versions differ:
- R: `Rscript --version`; for FUSION scripts inspect `--help` flags directly in the source
- Python: `pip show metaxcan pyfocus` then `SPrediXcan.py --help`, `SMulTiXcan.py --help`, `focus finemap --help`
- CLI: `plink2 --version`; FUSION ships as R scripts not a binary

If a script throws an error about an argument that has moved (e.g. `--gwas_file` vs `--gwas-file`) or a model database schema change, introspect the installed script with `--help` and adapt rather than retrying. PredictDB model file paths change with GTEx version; pin the version explicitly in scripts.

# Transcriptome-Wide Association

**"Find genes whose predicted tissue expression is associated with my GWAS trait"** -> Train SNP -> expression prediction models on a reference eQTL panel, apply the per-gene SNP weights to GWAS summary statistics or genotypes, and produce a gene-level Z-score equivalent to a weighted sum of SNP Z-scores. The output is a gene-by-tissue association, but TWAS is NOT direct evidence of causal mediation: an LD-tagged eQTL signal produces the same statistical association as a truly causal one, and the dominant failure modes are LD-induced false positives at gene-dense loci, tissue mis-specification, and ancestry mismatch between GWAS and prediction weights.

- CLI (sumstat TWAS, R): `FUSION.assoc_test.R --sumstats g.sumstats --weights weights.pos --weights_dir wgt/ --ref_ld_chr 1KG/EUR. --chr 22 --out chr22.dat`
- CLI (S-PrediXcan, Python): `SPrediXcan.py --model_db_path gtex_v8.db --covariance gtex_v8.cov --gwas_file g.txt --output_file out.csv`
- CLI (S-MultiXcan joint): `SMulTiXcan.py --models_folder mashr_models/ --gwas_folder gwas/ --metaxcan_folder spredixcan_per_tissue/ --output joint.csv`
- CLI (UTMOST cross-tissue): joint test across tissues via UTMOST's per-tissue GBJ / GBJ2 step
- CLI (FOCUS fine-mapping): `focus finemap chr22.dat 1KG_EUR focus.db --p-threshold 5e-8 --out chr22.focus`
- CLI (MA-FOCUS multi-ancestry): `focus finemap` with colon-separated per-ancestry sumstats / LD / weights and ancestry codes in `--locations`

TWAS, cis-eQTL MR, and coloc operate on overlapping evidence: TWAS asks "is the gene's predicted expression associated with the trait?"; cis-eQTL MR asks "does the eQTL effect on expression mediate the trait effect under IV assumptions?"; coloc asks "do the GWAS and eQTL share a causal variant?". Strong causal claims require triangulation, not single-method significance.

## Algorithmic Taxonomy

| Tool | Model | Input | Output | Strength | Fails when |
|------|-------|-------|--------|----------|------------|
| FUSION (Gusev 2016 Nat Genet 48:245) | Weighted sum of SNP Z-scores; per-gene multi-model (BLUP, lasso, elnet, top1, bslmm) selected by cross-validation R^2 | GWAS sumstats + pre-computed weights (.pos + per-gene RData) | TWAS Z, p; conditional joint analysis | Mature, ENCODE-style, pre-trained weights for many tissues (GTEx, CMC, METSIM, YFS, NTR) | LD-induced false positives at gene-dense loci; needs ancestry-matched LD reference; weights are heritability-thresholded so low-h2 genes drop |
| PrediXcan (Gamazon 2015 Nat Genet 47:1091) | Elastic-net SNP -> expression prediction; individual-level genotypes | PLINK genotypes + GWAS phenotype + GTEx prediction DB | Per-gene association test with full regression machinery | Most flexible (allows covariates, interactions, binary outcomes); same TWAS interpretation | Requires individual-level data; biobank-scale compute |
| S-PrediXcan (Barbeira 2018 Nat Commun 9:1825) | Summary-statistic equivalent of PrediXcan; Z-score weighted sum analogous to FUSION | GWAS sumstats + PredictDB model + covariance | Per-gene Z, p | Public PredictDB models (GTEx v8 elastic-net + MASHR-EUR); minimal compute | Pre-trained weights ancestry-specific (EUR primarily); covariance file must match model DB |
| S-MultiXcan (Barbeira 2019 PLoS Genet 15:e1007889) | Joint multi-tissue test combining per-tissue S-PrediXcan via PCA-regularised regression | Folder of per-tissue S-PrediXcan outputs + model folder | Single joint p per gene + per-tissue significance | Boosts power when causal tissue is unknown; standard for transcriptome-wide screens | Joint test cannot pinpoint causal tissue; correlated tissues produce ill-conditioned regression |
| UTMOST (Hu 2019 Nat Genet 51:568) | Cross-tissue elastic-net (group lasso) for joint tissue prediction + GBJ test | Per-tissue eQTL data + GWAS sumstats | Cross-tissue joint statistic | Often better-powered than S-MultiXcan at cross-tissue genes | Computationally heavier than MetaXcan; tissue weights are less interpretable |
| MOSTWAS (Bhattacharya 2021 PLoS Genet 17:e1009398) | Mediator-aware TWAS adding distal trans-mediating SNPs to cis-only models | GWAS sumstats + MOSTWAS weights | Per-gene Z, p (TWAS + distal mediator extension) | Recovers signal at genes with non-cis genetic regulation | Trans-mediation models need large reference panel; weights less broadly available |
| kTWAS (Cao 2021 Brief Bioinform 22:bbaa270) | Kernel-based TWAS using SKAT-style aggregation | GWAS sumstats or genotypes + per-gene SNP set | Per-gene p | Robust to non-linear and rare-variant contributions | Loses the eQTL-weighting interpretability; less standardised |
| EpiXcan (Zhang 2019 Nat Commun 10:3300) | Adds epigenome-derived per-SNP prior to weight training | eQTL + epigenome + GWAS sumstats | Per-gene Z, p (epigenome-informed) | Higher prediction R^2 in epigenome-rich tissues | Requires matched epigenome data for the prediction tissue |
| TIGAR-V2 (Tang 2021 NAR 49:e64) | Dirichlet process regression (non-parametric Bayes) for SNP -> expression | Reference eQTL + GWAS sumstats | Per-gene Z, p | Captures non-elastic-net effect structures; more flexible weight learning | Slower training; benefits depend on locus genetics |
| FOCUS (Mancuso 2019 Nat Genet 51:675) | Probabilistic gene-level fine-mapping over TWAS Z-scores using gene-by-gene predicted-expression correlation as analog of LD | FUSION/S-PrediXcan TWAS Z + ancestry-matched LD reference | Per-gene PIP + credible gene set | Resolves co-regulated gene clusters into a probabilistic causal gene; standard add-on after TWAS | Requires the same prediction-weight panel that produced TWAS Z; PIPs depend on prior |
| MA-FOCUS (Lu 2022 AJHG 109:1388-1404) | Multi-ancestry FOCUS; joint gene fine-mapping across ancestries with shared causal-gene assumption | Per-ancestry TWAS sumstats + per-ancestry weights + per-ancestry LD | Cross-ancestry PIP | Smaller credible gene sets when AFR/EAS contribute non-EUR LD information | Trans-ethnic gene-effect heterogeneity violated; weights must be ancestry-matched |
| JEPEG / JEPEG-Mix (Lee 2015 / 2016) | Gene-based test combining eQTL and functional weights | GWAS sumstats + JEPEG annotation database | Per-gene p | Lightweight gene-burden alternative to TWAS | Less granular than full PrediXcan/FUSION machinery; minimally updated |

Methodology evolves; the FUSION-vs-PrediXcan landscape has been stable but probabilistic fine-mapping (FOCUS, MA-FOCUS) and integrative methods (MOSTWAS, EpiXcan, OmicsXcan) continue to advance. Verify against the current PredictDB release notes (predictdb.org) and the latest FUSION weight panels (gusevlab.org/projects/fusion) before locking on a tissue or model.

S-PrediXcan and FUSION are mathematically near-identical: both compute a weighted sum of GWAS SNP Z-scores using per-gene SNP-expression weights and an LD-aware variance correction. The practical differences reduce to (a) the weight panel (FUSION elastic-net vs PredictDB MASHR), (b) the LD reference, and (c) the per-gene heritability threshold; the supplement of Barbeira 2018 works through the algebra. Method choice should therefore be driven by panel availability and ancestry-match, not by the underlying algorithm.

### PredictDB Model Choice and GTEx Versioning

| Release | Cohort | Status (2026) | When to use |
|---------|--------|----------------|-------------|
| GTEx v8 (838 donors, 49 tissues, 2020) | EUR-dominant (~85%) | Current PredictDB standard | Default; pre-trained MASHR + elastic-net DBs available |
| GTEx v9 (2023) | Expanded harmonization | Not migrated into PredictDB | Do not use until PredictDB rebuilds |
| GTEx v10 (2024 AnVIL release) | Re-aligned to GRCh38 v44 | Limited harmonization; not PredictDB-default | Wait for community-validated weight panels |

**Operational rule:** Use GTEx v8 unless there is an explicit biological reason to deviate (tissue not in v8, ancestry-specific panel preferred). In methods, pin exactly: "GTEx v8 MASHR-EUR, PredictDB release 2022-01".

### MASHR vs Elastic-Net Models

PredictDB ships two cross-validated model families per tissue (Barbeira 2020 Genome Biol 21:36):

| Model | Construction | Per-gene SNP count | When to use |
|-------|--------------|---------------------|-------------|
| MASHR | Cross-tissue posterior mean from DAP-G fine-mapped SNPs | ~10x sparser | Primary discovery; higher per-gene R^2 in most genes; standard for S-MultiXcan |
| Elastic-net | Per-tissue lasso/ridge mix (alpha = 0.5) | Denser | Tissues where MASHR's cross-tissue prior is mis-specified (ovary, testis, isolated-organ traits) |

**Operational rule:** Never mix MASHR and elastic-net within a single S-MultiXcan run; the inter-tissue covariance and condition number assumptions break. Pick one family and apply consistently across all tissues.

## Decision Tree by Experimental Scenario

| Scenario | Recommended workflow | Why |
|----------|---------------------|-----|
| GWAS summary stats only, EUR, single hypothesis tissue (e.g. liver for LDL) | S-PrediXcan with GTEx v8 MASHR-EUR weights, or FUSION with GTEx liver | Standard pre-trained pipeline; minimal compute |
| GWAS summary stats only, tissue unknown a priori | S-MultiXcan (standard GTEx v8) OR UTMOST (custom panel) -- see S-MultiXcan vs UTMOST table | Joint multi-tissue inflates power; tissue prioritisation requires LDSC-SEG separately |
| Multiple TWAS hits at one locus (gene-dense region) | Run TWAS then FOCUS for probabilistic fine-mapping | LD ties co-regulated genes; FOCUS PIP distinguishes likely causal gene |
| Multi-ancestry GWAS (EUR + EAS + AFR) | Per-ancestry S-PrediXcan with matched weights, then MA-FOCUS to combine | Single-ancestry weights miscalibrated in other ancestries; joint fine-mapping shrinks credible set |
| Individual-level genotypes available (UKB) | PrediXcan (full regression) | Allows covariates, interactions, binary outcomes natively |
| Low-N tissue (GTEx N < 100) | Substitute eQTLGen (whole blood, N ~ 31k) OR skip the tissue | Per-gene CV R^2 unstable below N ~ 100; weights overfit |
| Drug-target prioritisation (TWAS as causal-gene evidence) | TWAS + cis-eQTL MR + coloc + FOCUS triangulation | TWAS alone is associational; triangulation strengthens causal claim |
| Trans-acting / mediator-aware analysis | MOSTWAS | Adds distal trans-mediating SNPs to cis-only models |
| Rare-variant or non-linear gene effects | kTWAS or TIGAR-V2 | Kernel / non-parametric flexibility |
| HLA region (chr6:25-35 Mb hg38) | Exclude or use HLA-specific tools | Long-range LD breaks every gene-by-gene method; standard TWAS PIPs not interpretable |
| Splicing-mediated trait (e.g. neuropsych for sQTL) | sTWAS (sQTL-weighted TWAS) using GTEx splicing models | Splicing mediates many GWAS effects; cis-sQTL panels available in PredictDB |
| Cell-type-specific trait | sc-eQTL-based TWAS (e.g. OneK1K, Yazar 2022 Science 376:eabf3041) | Bulk-tissue TWAS averages over cell types; single-cell eQTL recovers cell-type specificity |

### Tissue Selection Protocol

Tissue choice drives TWAS power and false-positive rate; selecting tissues by inspecting TWAS hit count is circular. Run all three of the following on the GWAS sumstats (independent of any TWAS run) and pick the primary TWAS tissue from the intersection:

1. **Stratified LDSC tissue prioritization** (Finucane 2018 Nat Genet 50:621): `ldsc.py --h2-cts <sumstats> --ref-ld-chr-cts <annot> --w-ld-chr <weights>` against 200+ tissue-specific gene expression annotations
2. **CELLEX** (Timshel 2020 eLife 9:e55851): single-cell tissue / cell-type prioritization on the same GWAS
3. **MAGMA gene-property analysis** (de Leeuw 2016 PLoS Comput Biol 12:e1004219): cheaper substitute when LDSC unavailable

**Operational rule:** Primary TWAS tissue = the tissue with FDR-significant enrichment in at least two of the three methods. Run secondary tissues in S-MultiXcan for cross-tissue replication. Bonferroni for tissue selection alone: 0.05 / 200 annotations = 2.5e-4.

### S-MultiXcan vs UTMOST

| Method | Use case | Rationale |
|--------|----------|-----------|
| S-MultiXcan (Barbeira 2019) | Standard GTEx v8 analysis | Pre-computed MASHR weights; lower compute barrier; PCA-regularised inter-tissue regression |
| UTMOST (Hu 2019 Nat Genet 51:568) | Custom eQTL panel with cross-tissue retraining | Higher power at genes with shared cross-tissue eQTL architecture; group-lasso enforces sparsity across tissues |

Benchmarks: Hu 2019, Barbeira 2019. Choose by panel availability first; the methods recover overlapping but non-identical gene sets.

### Single-Cell and Cell-Type-Resolved TWAS

Bulk-tissue TWAS averages over cell composition; sc-eQTL TWAS recovers cell-type-specific regulation but at lower per-cell-type power.

| Panel | Reference | Cells / tissue |
|-------|-----------|----------------|
| OneK1K | Yazar 2022 Science 376:eabf3041 | PBMC, ~982 donors, 14 cell types |
| HipSci iPSC-eQTL | Kilpinen 2017 Nature 546:370 | iPSC, ~317 donors |
| BLUEPRINT | Chen 2016 Cell 167:1398 | Monocytes, neutrophils, T cells |

**Tooling state (2026):** No fully pre-built scPrediXcan equivalent to MASHR; sc-eQTL weights are panel-specific. Train custom PredictDB or use TIGAR-V2's Bayesian DPR on the sc-eQTL matrix.

**Operational rule:** Run standard bulk-tissue TWAS first; run sc-eQTL TWAS in the prioritized cell type as a secondary analysis; require concordance between bulk and sc results before nominating a cell-type-specific gene. Upstream sc preprocessing: cross-reference single-cell/preprocessing.

## Per-Tool Failure Modes

### LD-induced TWAS false positives (most common pitfall)

**Trigger:** Two or more genes at the same locus have correlated cis-eQTLs (shared causal eQTL SNP or LD-linked eQTL SNPs).

**Mechanism:** TWAS Z-scores are linear combinations of SNP Z-scores weighted by per-gene SNP effects. When two genes share many high-weight SNPs (e.g. nearby genes regulated by the same enhancer or LD-tagged independent eQTLs), their TWAS Z-scores are positively correlated. A single causal GWAS variant therefore produces significant Z at multiple co-regulated genes (Wainberg 2019 Nat Genet 51:592; Mancuso 2019 Nat Genet 51:675).

**Symptom:** A GWAS lead locus shows 3-10 genes all passing genome-wide TWAS significance (p < 2.3e-6 ~ 0.05/22k); per-gene LocusZoom-style plots look near-identical; conditional analysis (FUSION `--coloc_P` or `--joint`) reveals only 1-2 independent gene signals; the genes lie within 1 Mb of each other.

**Fix:** Always run FOCUS after TWAS to obtain per-gene PIPs. Report only genes with PIP >= 0.8 as candidate causal; report co-significant genes with PIP < 0.5 as LD-tagged. Cross-check with cis-eQTL coloc (PP.H4 >= 0.7) for the candidate causal gene. FUSION's `FUSION.post_process.R --joint` performs conditional analysis as a lighter-weight alternative.

### Tissue mis-specification

**Trigger:** Running TWAS in a tissue that does not host the causal regulatory effect (e.g. whole blood for a psychiatric trait; pancreas for an LDL trait).

**Mechanism:** A gene's cis-eQTL effect size varies across tissues; a wrong-tissue model has weaker per-gene prediction R^2 and lower power. Conversely, eQTL effects in the wrong tissue can still tag the GWAS signal via LD and produce spurious associations not present in the causal tissue.

**Symptom:** Strong TWAS signal in a tissue biologically irrelevant to the trait; null in the expected tissue; tissue-prioritisation methods (LDSC-SEG, Finucane 2018 Nat Genet 50:621; CELL-TYPE-SPECIFIC LDSC, Calderon 2017) disagree with the TWAS tissue.

**Fix:** Run S-MultiXcan to combine tissues if causal tissue is unknown. For prioritisation, use LDSC-SEG / CELLEX / EWCE on the GWAS sumstats independently of TWAS, and report TWAS in the prioritised tissues. Never report a single-tissue TWAS hit as causal without independent tissue evidence (single-cell eQTL, chromatin accessibility in matched cell type).

### Ancestry mismatch in prediction weights

**Trigger:** Running TWAS on a non-EUR GWAS using GTEx (~ 85% EUR) weights, or vice versa.

**Mechanism:** Cis-eQTL effect sizes and LD structure are ancestry-specific; prediction weights trained in one ancestry transfer with reduced R^2 and biased Z-scores in another. Power is lost preferentially at loci where the causal eQTL is not shared across ancestries (Patel 2022 Genome Med 14:31).

**Symptom:** Genome-wide TWAS hit count much lower than expected given GWAS power; non-EUR-specific GWAS loci fail to produce TWAS hits; per-gene prediction R^2 substantially reduced.

**Fix:** Use ancestry-matched prediction panels where available: MESA multi-ethnic eQTL (Mogil 2018 PLoS Genet), eQTLGen-Asian, AFGR (Africa) when published, or MAGE (Taliun-style multi-ancestry eQTL). Move to MA-FOCUS for cross-ancestry joint fine-mapping. Document the ancestry assumption explicitly in methods.

### Low-N tissue weights are unstable

**Trigger:** Using a GTEx tissue with N < 100 donors (e.g. several brain sub-regions, kidney cortex in v7).

**Mechanism:** Per-gene elastic-net weights are cross-validated with the available donors. Below ~ 100 donors, the cross-validation R^2 has high variance and the heritability filter (FUSION requires hsq_p < 0.01) drops many genes. Surviving weights overfit, inflating per-gene Z under the null.

**Symptom:** Tissue produces unusually high TWAS hit count or unusually high genomic inflation; per-gene CV R^2 distribution is bimodal with a long heavy tail.

**Fix:** Skip GTEx tissues with N < 100 unless biologically essential. Substitute eQTLGen for whole blood (N ~ 31k, Vosa 2021 Nat Genet 53:1300) where blood is acceptable. For brain, use PsychENCODE (N ~ 1300, Wang 2018 Science 362:eaat8464) or BrainSeq (N ~ 350+) when available; verify the matching prediction-weight panel exists.

### HLA region

**Trigger:** Any gene within chr6:25-35 Mb (hg38; extended MHC) reported by TWAS.

**Mechanism:** Long-range LD (r2 > 0.5 over many Mb) and extreme structural variation mean per-gene prediction weights at HLA capture haplotype rather than gene-specific regulation. Standard TWAS gene-level inference is biologically meaningless here.

**Fix:** Exclude chr6:25-35 Mb from genome-wide TWAS summaries by default. For HLA-driven traits (autoimmune, infection, transplantation), impute classical HLA alleles using one of:

| Tool | Reference | Notes |
|------|-----------|-------|
| HIBAG | Zheng 2014 Pharmacogenomics J 14:192 | R package; pre-trained per-ancestry classifiers |
| SNP2HLA | Jia 2013 PLoS One 8:e64683 | Beagle-based imputation; supports T1DGC reference |
| HLA-TAPAS | Luo 2021 Nat Genet 53:1504 | Current standard; multi-ancestry reference; recommended for new analyses |

Then test classical alleles plus amino-acid residues (Raychaudhuri 2012 Nat Genet 44:291 set the gold standard for residue-level association in MHC). Do NOT run SNP-level TWAS inside the MHC.

### Correlated-expression confounding (co-regulated genes)

**Trigger:** Two or more genes are functionally co-regulated by a single TF or enhancer, producing nearly-identical predicted-expression vectors.

**Mechanism:** Even with separate per-gene cis-eQTL prediction, downstream co-regulation makes predicted expression highly correlated; TWAS cannot distinguish which gene mediates the trait.

**Symptom:** FOCUS credible gene set contains multiple genes with PIP roughly equal (e.g. three genes at 0.3 each); functional follow-up (MPRA, CRISPRi screens) is needed to break the tie.

**Fix:** Acknowledge the limit of statistical resolution; report the full credible gene set and prioritise on orthogonal evidence (CRISPRi/CRISPRa effect size in matched cell type, e.g. Open Targets-style, or MPRA at allelic series; protein-level pQTL coloc if available).

## TWAS - MR - Coloc Triangulation

A TWAS-significant gene is associational, not causal. The strongest defensible claim that a gene mediates a GWAS effect comes from triangulating three orthogonal lines of evidence:

| Line | What it tests | Threshold |
|------|---------------|-----------|
| TWAS | Predicted-expression association with trait | S-MultiXcan joint p < 2.3e-6, OR S-PrediXcan per-tissue p < 4.6e-8 (49-tissue Bonferroni), OR per-tissue FDR < 0.05 |
| cis-eQTL MR | Causal effect of expression on trait under IV assumptions (cross-reference causal-genomics/mendelian-randomization) | Wald-ratio or IVW p < 0.05/n_genes; instrument F > 10 |
| Colocalization | Shared causal variant between GWAS and eQTL (cross-reference causal-genomics/colocalization-analysis) | coloc.abf or coloc.susie PP.H4 >= 0.7 |
| FOCUS | Probabilistic per-gene PIP under TWAS fine-mapping | PIP >= 0.8 |

**Operational rule:** Report a gene as a "strong candidate causal gene" only when 3 of the 4 are concordant (TWAS hit + coloc PP.H4 >= 0.7 + FOCUS PIP >= 0.8, with cis-MR as a supporting fourth). 2-of-4 concordance is "suggestive"; 1-of-4 is "associational only". The combination is more conservative than any single method but matches the standards used in modern GWAS-to-target pipelines (Open Targets Genetics Mountjoy 2021 Nat Genet 53:1527; FinnGen R10 release notes).

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|--------------------|
| S-MultiXcan joint p (gene-wide) | < 2.3e-6 (0.05 / 22k genes) | One p per gene across tissues; standard genome-wide TWAS significance |
| S-PrediXcan per-tissue (49 tissues) | < 4.6e-8 (0.05 / (49 x 22k)) | Cross-tissue x cross-gene Bonferroni |
| Per-tissue FDR (alternative) | BH q < 0.05 within each tissue | Less conservative; report cross-tissue replication of FDR-significant hits |
| FUSION cross-validation R^2 | >= 0.01 (per gene) | FUSION default heuristic; below this, gene is not heritable in tissue and weights drop |
| FUSION heritability p | < 0.01 (hsq_p) | FUSION default; filters out non-heritable expression |
| Coloc PP.H4 for triangulation | >= 0.7 | Open Targets / common practice; >= 0.8 for high-confidence |
| FOCUS gene PIP (causal) | >= 0.8 | Mancuso 2019 convention; PIP >= 0.5 is suggestive |
| cis-eQTL MR instrument F | >= 10 | Standard MR convention to avoid weak-instrument bias |
| Tissue eQTL N (weight stability) | >= 100 donors | Below this, per-gene elastic-net weights unstable |
| cis-window radius | +/- 500 kb of gene TSS/TES (FUSION) or +/- 1 Mb (S-PrediXcan) | Captures most cis-eQTL signal; window choice rarely changes top hits |
| LD reference panel size | >= 500 individuals matched ancestry | 1000 Genomes superpopulation reference is standard |
| MA-FOCUS minimum ancestry N | >= 2 ancestries with significant TWAS | Joint inference requires non-trivial heterogeneity |
| S-MultiXcan condition number | < 30 (correlation matrix) | PCA regularisation kicks in above this; near-collinear tissues collapsed |

## FUSION Pipeline

**Goal:** Run sumstat TWAS using FUSION and conditional joint analysis to identify independently associated genes.

**Approach:** Format GWAS sumstats to FUSION expected columns (SNP A1 A2 Z); run `FUSION.assoc_test.R` per chromosome with pre-computed weights and ancestry-matched LD; post-process with `FUSION.post_process.R --joint` to identify independent genes; flag conditional-significant genes for follow-up.

```bash
# Pre-computed FUSION weights live at http://gusevlab.org/projects/fusion/
# Example: GTEx v8 Whole_Blood; download .pos summary + per-gene RData files into wgt_dir/

# GWAS sumstats expected columns: SNP A1 A2 Z (Z-score on standardised scale)
# Use TwoSampleMR or a custom munger to harmonise alleles upstream

for chr in {1..22}; do
    Rscript FUSION.assoc_test.R \
        --sumstats gwas.sumstats \
        --weights gtex_whole_blood.pos \
        --weights_dir gtex_whole_blood_wgt/ \
        --ref_ld_chr 1000G_EUR_LD/EUR. \
        --chr ${chr} \
        --out twas_chr${chr}.dat
done
cat twas_chr*.dat > twas_all.dat

# Conditional joint analysis at each significant locus
Rscript FUSION.post_process.R \
    --sumstats gwas.sumstats \
    --input twas_all.dat \
    --out twas_joint.dat \
    --ref_ld_chr 1000G_EUR_LD/EUR. \
    --chr 22 \
    --plot --locus_win 100000
# twas_joint.dat reports per-gene conditional Z; genes with joint Z > 4 are independent
```

The `--locus_win 100000` parameter defines the conditioning window; 100 kb is conservative for non-HLA loci. FUSION's `--coloc_P` flag runs single-SNP coloc internally but is less robust than running coloc separately on the per-gene top eQTL.

## S-PrediXcan + S-MultiXcan Pipeline

**Goal:** Run TWAS across all GTEx tissues using pre-trained PredictDB models, then combine via S-MultiXcan for a joint multi-tissue test.

**Approach:** For each tissue, run S-PrediXcan with the matched model DB and covariance file; collect per-tissue outputs into a folder; run S-MultiXcan with the same model folder and GWAS to produce a joint multi-tissue Z and per-tissue significance.

```bash
# PredictDB models: predictdb.org
# GTEx v8 MASHR-EUR is the standard EUR panel
# Each tissue has a .db (model) and .txt.gz (covariance) file pair

mkdir -p spredixcan_out
for tissue in Whole_Blood Liver Brain_Frontal_Cortex_BA9 Adipose_Subcutaneous; do
    python SPrediXcan.py \
        --model_db_path mashr_models/mashr_${tissue}.db \
        --covariance mashr_models/mashr_${tissue}.txt.gz \
        --gwas_file gwas.txt \
        --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 \
        --beta_column BETA --pvalue_column P \
        --output_file spredixcan_out/${tissue}.csv
done

# Joint multi-tissue
python SMulTiXcan.py \
    --models_folder mashr_models/ \
    --models_name_pattern "mashr_(.*)\.db" \
    --snp_covariance gtex_v8_expression_mashr_snp_covariance.txt.gz \
    --metaxcan_folder spredixcan_out/ \
    --metaxcan_filter "(.*)\.csv" \
    --metaxcan_file_name_parse_pattern "(.*)\.csv" \
    --gwas_file gwas.txt \
    --snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 \
    --beta_column BETA --pvalue_column P \
    --output joint_multitissue.csv
```

S-MultiXcan applies PCA regularisation on the inter-tissue correlation matrix; `--regularization 0.1` (default) is the conditioning ridge. Tissues that are nearly collinear with another (e.g. multiple brain sub-regions) are absorbed into shared components and do not contribute independent power.

## FOCUS Probabilistic Fine-Mapping

**Goal:** Resolve a locus with multiple co-significant TWAS genes into a probabilistic causal-gene credible set.

**Approach:** Build (or download) a FOCUS gene-prediction database matching the TWAS weight panel; run `focus finemap` with the TWAS sumstats and ancestry-matched LD reference; report PIPs and credible sets.

```bash
# Install: pip install pyfocus
# Pre-built FOCUS DBs for FUSION/PrediXcan weights live at github.com/bogdanlab/focus

# Convert FUSION TWAS output to FOCUS sumstat format if needed
# FOCUS expects: CHR SNP BP A1 A2 Z P from the underlying GWAS sumstats (NOT TWAS Z)

focus finemap \
    gwas.sumstats \
    1000G_EUR_chr \
    focus_gtex_v8_whole_blood.db \
    --p-threshold 5e-8 \
    --tissue Whole_Blood \
    --out gwas_focus_whole_blood
# Output: per-gene PIP, credible-set membership flag, and locus-level group probability
```

For a custom prediction-weight panel without a pre-built FOCUS DB, construct one from FUSION weights:

```bash
focus import-fusion --gwas gwas.sumstats --pos custom_panel.pos --out custom_focus.db
```

FOCUS PIPs depend on the per-locus prior probability that any gene is causal. Report a sensitivity scan over `--prior_prob`:

| Setting | Use case |
|---------|----------|
| `--prior_prob 1e-3` (default) | Standard genome-wide TWAS; matches Mancuso 2019 |
| `--prior_prob 1e-4` (conservative) | High-prior gene-dense locus where most genes are not causal |
| `--prior_prob 1e-2` (liberal) | Pre-prioritised candidate region where one gene is expected |

Cite the Mancuso 2019 supplement for the sensitivity-scan protocol. MA-FOCUS extends with per-ancestry weights and is invoked via the same `focus finemap` CLI -- multi-ancestry mode is signaled by colon-separated per-ancestry sumstats, LD references, and weight DBs, plus paired ancestry codes in `--locations`:

```bash
focus finemap \
    eur.sumstats.tsv.gz:eas.sumstats.tsv.gz:afr.sumstats.tsv.gz \
    1000G_EUR_chr:1000G_EAS_chr:1000G_AFR_chr \
    focus_eur.db:focus_eas.db:focus_afr.db \
    --chr 22 --locations 38:EUR-EAS-AFR \
    --out gwas_ma_focus
```

MA-FOCUS assumes a shared causal gene across ancestries; gene-specific heterogeneity (e.g. an ancestry-specific eQTL) violates the assumption and produces inflated H0/heterogeneous-group probability.

## When Standard Pipeline is Insufficient

The FUSION / S-PrediXcan / S-MultiXcan + FOCUS pipeline is the default. Escalate to a specialised method only when the standard pipeline misses an expected hit (a strong GWAS signal with no TWAS gene, or a known causal gene without recovery).

| Method | Triggering scenario | Yield |
|--------|--------------------|-------|
| MOSTWAS (Bhattacharya 2021 PLoS Genet 17:e1009398) | Trans-mediator architecture suspected (immune, neuropsych traits) | ~15% additional hits via distal mediator terms |
| EpiXcan (Zhang 2019 Nat Commun 10:3300) | Paired epigenome data available (Roadmap, EpiMap, ENCODE cell-matched DNase/H3K27ac) | Higher prediction R^2 in epigenome-rich tissues |
| TIGAR-V2 (Tang 2021 NAR 49:e64) | Training a custom Bayesian DPR panel (no pre-trained PredictDB for the tissue) | Captures non-elastic-net effect structures |
| kTWAS (Cao 2021 Brief Bioinform 22:bbaa270) | Rare-variant or population-specific contexts; cis-eQTL panels under-powered | Kernel aggregation robust to non-linear and rare effects |

**Operational rule:** Default to FUSION / S-PrediXcan / S-MultiXcan + FOCUS first; document an expected-but-missing hit before escalating.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FUSION significant + S-PrediXcan null | Different weight panels (FUSION vs PredictDB) or different LD references | Re-run with matched panel; check FUSION weight CV-R^2 vs PredictDB MASHR posterior mean |
| S-MultiXcan p << min per-tissue p | Joint test borrowing strength across tissues | Genuine; report the joint p and the top contributing tissue |
| Many co-significant genes at one locus, no FOCUS PIP > 0.5 | LD-tied co-regulated genes, true causal gene not in panel | Expand panel (add brain or tissue-specific weights); functional follow-up needed |
| FOCUS PIP > 0.8 + coloc PP.H4 < 0.5 | Sparse eQTL signal (single SNP drives prediction) + locus has another co-localising signal | Investigate the single top-eQTL SNP; check for fine-mapped credible-set overlap |
| TWAS hit + cis-eQTL MR null | TWAS hit is LD-tagged, not mediated by expression | Trust the cis-eQTL MR result; flag the gene as TWAS-positive but non-causal |
| MA-FOCUS PIP much lower than per-ancestry FOCUS | Cross-ancestry heterogeneity; gene effect is not shared | Report per-ancestry separately; do not pool |
| TWAS significant in wrong tissue | LD-induced via tissue-shared eQTL | Verify with LDSC-SEG tissue prioritisation; treat top-tissue TWAS hit as the trustworthy one |

**Operational rule:** No single-method TWAS result is sufficient evidence of causal gene identity. The minimum reporting standard is (a) TWAS significance threshold met with multiple-testing correction; (b) FOCUS PIP >= 0.8 OR coloc PP.H4 >= 0.7; (c) explicit acknowledgement of tissue choice and ancestry of prediction weights. Triangulation with cis-eQTL MR or independent CRISPR/MPRA validation lifts a finding from "candidate" to "supported".

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error: weight file <gene>.RDat not found` (FUSION) | weights.pos points to relative paths; `--weights_dir` not set correctly | Verify `--weights_dir` matches the .pos `WGT` column directory |
| `KeyError: 'SNP'` or `'A1'` (S-PrediXcan) | GWAS sumstat columns not mapped via `--snp_column` etc. | Pass `--snp_column SNP --effect_allele_column A1 --non_effect_allele_column A2 --beta_column BETA --pvalue_column P` explicitly |
| All TWAS Z's are 0 or NA | Allele coding mismatch between GWAS and weight panel | Harmonise (effect-allele flip); check `--keep_non_rsid` if using non-rsID SNPs |
| Genomic inflation lambda >> 1.05 | Wrong LD reference OR population structure not controlled in upstream GWAS | Fix at the GWAS stage; do not adjust TWAS lambda post-hoc |
| FOCUS reports PIP = 1 for one gene at every locus | Only one gene in panel at locus; degenerate posterior | Expand panel coverage or report locus as panel-limited |
| S-MultiXcan condition number warning | Tissues near-collinear (multiple brain regions) | Increase regularisation (`--regularization 0.5`) or restrict to tissue subset |
| FOCUS DB not matching FUSION weights | Custom weight panel without FOCUS DB | Build FOCUS DB from weights using `focus import-fusion` |
| MA-FOCUS H0 probability dominates | Cross-ancestry heterogeneity at the locus | Run per-ancestry FOCUS separately; do not force joint |
| S-PrediXcan output has effect sizes much larger than expected | `sdY` proxy mis-specified; standardised vs unstandardised mismatch | Confirm GWAS Z and beta scale; rerun with `--additional_output` |

## Required Reporting for Publication

A defensible TWAS report includes every item below in methods or supplement:

- GWAS sumstat source, effective N (Neff), and ancestry composition
- Weight panel and version (e.g. "GTEx v8 MASHR-EUR, PredictDB release 2022-01")
- LD reference panel and version (e.g. "1000 Genomes Phase 3 EUR" or "UK Biobank array")
- Per-tissue list (or explicit "all 49 GTEx v8 tissues")
- Multiple-testing correction strategy (S-MultiXcan joint, per-tissue Bonferroni, or per-tissue FDR)
- FOCUS PIP threshold and prior-probability sensitivity scan
- Coloc PP.H4 threshold (cross-reference causal-genomics/colocalization-analysis)
- cis-eQTL MR estimate, instrument F-statistic, and TwoSampleMR package version (cross-reference causal-genomics/mendelian-randomization)
- Triangulation rule (e.g. "3-of-4 concordance across TWAS, FOCUS, coloc, cis-MR")
- HLA exclusion confirmed (chr6:25-35 Mb dropped from genome-wide summaries)

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "LD-induced false positive at gene-dense locus?" | FOCUS PIP reported per gene; only PIP >= 0.8 carried forward as candidate causal |
| "Tissue mis-specification?" | Tissue Selection Protocol applied (sLDSC + CELLEX + MAGMA); cross-tissue S-MultiXcan reported for replication |
| "Why GTEx and not eQTLGen?" | Justified by trait biology: blood-relevant traits use eQTLGen (N ~ 31k); tissue-specific traits use GTEx v8 MASHR |
| "MHC?" | chr6:25-35 Mb excluded; HLA-TAPAS run separately for HLA-relevant traits |
| "Triangulation?" | TWAS + coloc + cis-MR + FOCUS run; 3-of-4 concordance required for the strong-candidate label |
| "Ancestry transfer?" | Ancestry-matched weights used where available; MA-FOCUS applied for multi-ancestry GWAS |
| "Prior sensitivity in FOCUS?" | `--prior_prob` scanned at 1e-2, 1e-3, 1e-4; PIPs reported across the scan |

## Tool Install Notes

- **FUSION**: gusevlab.org/projects/fusion or `git clone https://github.com/gusevlab/fusion_twas`. R scripts; needs plink (PLINK 1.9), Rscript, and the GBJ R package for omnibus. Pre-computed weights for GTEx v7/v8, CMC, YFS, METSIM, NTR, MESA at the same site.
- **MetaXcan / S-PrediXcan / S-MultiXcan**: `pip install metaxcan` OR `git clone https://github.com/hakyimlab/MetaXcan`. Python scripts in `software/`. Compatible with Python 3.9-3.11.
- **PredictDB models**: predictdb.org. GTEx v8 elastic-net (single-tissue) and MASHR (cross-tissue posterior) databases for EUR; multi-ethnic panels emerging.
- **UTMOST**: `git clone https://github.com/Joker-Jerome/UTMOST`. Python. Needs precomputed cross-tissue weights or training pipeline.
- **FOCUS**: `pip install pyfocus`. CLI `focus`. Pre-built DBs for GTEx v7/v8 panels at github.com/bogdanlab/focus.
- **MA-FOCUS**: `git clone https://github.com/mancusolab/ma-focus && cd ma-focus && pip install .` (no PyPI release; install from source). Same CLI as single-ancestry FOCUS -- `focus finemap` -- with colon-separated per-ancestry sumstats / LD / weight DBs and paired ancestry codes in `--locations`.
- **TIGAR-V2**: `git clone https://github.com/yanglab-emory/TIGAR`. Python + R hybrid; ships with example data.
- **MOSTWAS**: `git clone https://github.com/bhattacharya-a-bt/MOSTWAS`. R package; install via `devtools::install_github`.
- **EpiXcan**: Bitbucket roussos-lab/epixcan_translation. Workflow-style pipeline.

## References

- Gamazon ER, Wheeler HE, Shah KP, Mozaffari SV, Aquino-Michaels K et al 2015 Nat Genet 47:1091 (PrediXcan)
- Gusev A, Ko A, Shi H, Bhatia G, Chung W et al 2016 Nat Genet 48:245 (FUSION TWAS)
- Barbeira AN, Dickinson SP, Bonazzola R, Zheng J, Wheeler HE et al 2018 Nat Commun 9:1825 (S-PrediXcan)
- Barbeira AN, Pividori M, Zheng J, Wheeler HE, Nicolae DL et al 2019 PLoS Genet 15:e1007889 (S-MultiXcan)
- Hu Y, Li M, Lu Q, Weng H, Wang J et al 2019 Nat Genet 51:568 (UTMOST cross-tissue)
- Mancuso N, Freund MK, Johnson R, Shi H, Kichaev G et al 2019 Nat Genet 51:675 (FOCUS gene-level fine-mapping)
- Lu Z, Wang X, Carr M, Kim A, Gazal S et al 2022 AJHG 109:1388-1404 (MA-FOCUS multi-ancestry)
- Wainberg M, Sinnott-Armstrong N, Mancuso N, Barbeira AN, Knowles DA et al 2019 Nat Genet 51:592 (TWAS limitations and LD-induced false positives)
- Bhattacharya A, Li Y, Love MI 2021 PLoS Genet 17:e1009398 (MOSTWAS distal mediation)
- Cao C, Kwok D, Edie S, Li Q, Ding B et al 2021 Brief Bioinform 22:bbaa270 (kTWAS)
- Zhang W, Voloudakis G, Rajagopal VM, Readhead B, Dudley JT et al 2019 Nat Commun 10:3300 (EpiXcan)
- Tang S, Buchman AS, De Jager PL, Bennett DA, Epstein MP, Yang J 2021 Nucleic Acids Res 49:e64 (TIGAR-V2)
- Vosa U, Claringbould A, Westra HJ, Bonder MJ, Deelen P et al 2021 Nat Genet 53:1300 (eQTLGen reference)
- GTEx Consortium 2020 Science 369:1318 (GTEx v8 multi-tissue eQTL)
- Mountjoy E, Schmidt EM, Carmona M, Schwartzentruber J, Peat G et al 2021 Nat Genet 53:1527 (Open Targets Genetics)
- Mogil LS, Andaleon A, Badalamenti A, Dickinson SP, Guo X et al 2018 PLoS Genet 14:e1007586 (MESA multi-ethnic eQTL)
- Patel RA, Musharoff SA, Spence JP, Pimentel H, Tcheandjieu C et al 2022 Genome Med 14:31 (TWAS ancestry transfer)
- Yazar S, Alquicira-Hernandez J, Wing K, Senabouth A, Gordon MG et al 2022 Science 376:eabf3041 (OneK1K sc-eQTL)

## Related Skills

- causal-genomics/fine-mapping - Variant-level credible sets feeding FOCUS gene-level fine-mapping
- causal-genomics/colocalization-analysis - Coloc PP.H4 triangulation with TWAS hits
- causal-genomics/mendelian-randomization - cis-eQTL MR triangulation; drug-target prioritisation
- causal-genomics/effector-gene-prioritization - Downstream gene mapping from TWAS candidate sets
- causal-genomics/proteome-mr-drug-target - Drug-target triangulation using protein-level evidence
- causal-genomics/pleiotropy-detection - Distinguishing horizontal pleiotropy from mediated TWAS signal
- causal-genomics/mediation-analysis - Downstream gene-mediated trait effects given TWAS hits
- population-genetics/association-testing - Upstream GWAS summary statistic generation
- population-genetics/linkage-disequilibrium - LD reference panel construction for FUSION / FOCUS
- differential-expression/deseq2-basics - Generating eQTL count data for custom prediction-weight training
- single-cell/preprocessing - Cell-type-resolved eQTL panels for sc-TWAS
- workflows/gwas-pipeline - End-to-end GWAS pipeline producing TWAS input
- variant-calling/variant-annotation - Functional annotation of TWAS / FOCUS top variants
