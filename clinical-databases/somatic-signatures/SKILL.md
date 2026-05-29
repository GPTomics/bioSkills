---
name: bio-clinical-databases-somatic-signatures
description: Extracts and assigns COSMIC v3.4 mutational signatures (84 SBS / 11 DBS / 18 ID / 24 CN / 16 SV) from somatic VCFs using SigProfilerSuite, MutationalPatterns, MuSiCal mvNMF, SigNet, or HRDetect. Use when characterizing DNA-damage etiology (BRCA1/2 HRD, MMR-D, POLE, APOBEC3A, UV, tobacco, aflatoxin, 5-FU/SBS17b, platinum, colibactin SBS88), routing PARP inhibitor decisions, or auditing de novo extraction vs refit choice for cohort size.
tool_type: mixed
primary_tool: SigProfilerAssignment
---

## Version Compatibility

Reference examples tested with: SigProfilerMatrixGenerator 1.2+ (Bergstrom 2019), SigProfilerExtractor 1.1.24+ (Islam 2022), SigProfilerAssignment 0.1+ (Diaz-Gay 2023), MutationalPatterns 3.12+ (Manders 2022), MuSiCal 0.7+ (Liu 2024), SigNet (Serrano 2023, bioRxiv), HRDetect (Davies 2017 / Degasperi 2022 implementations), pandas 2.2+, R 4.3+. COSMIC v3.4 (September 2024) is the current reference catalog: 84 SBS, 11 DBS, 18 ID, 24 CN, 16 SV signatures.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. COSMIC signature naming evolves: SBS40 was split to SBS40a/b/c in v3.4 (Senkin 2024); SBS17 split to SBS17a/b (5-FU); SBS10 split to SBS10a-d (POLE/POLD1).

# Somatic Mutational Signatures; Etiology, Extraction, Clinical Use

**'Extract mutational signatures from this tumor cohort and identify HRD/MMR/APOBEC processes'** -> Generate 96-context (or DBS/ID/CN/SV) matrix from VCF; choose de novo extraction (NMF) vs refit-to-COSMIC by cohort size; map dominant signatures to etiology; flag clinical actionability.

- Python (recommended): SigProfilerMatrixGenerator -> SigProfilerExtractor (de novo) or SigProfilerAssignment (refit)
- R alternative: `MutationalPatterns::fit_to_signatures()` (strict refit) or `extract_signatures()` (NMF de novo)
- Python (mvNMF for non-uniqueness): MuSiCal (Liu 2024 *Nat Genet*)
- Python (deep learning low-mutation count): SigNet (Serrano 2023)
- R (HRD-specific): HRDetect (Davies 2017 *Nat Med*); the 5-feature BRCA-deficiency classifier

## COSMIC v3.4 Catalog: Evolution and Composition

| Class | Count | Encoding |
|-------|-------|----------|
| **SBS** (Single Base Substitutions) | 84 | 96 trinucleotide contexts (6 substitution types x 16 trinucleotides) |
| **DBS** (Doublet Base Substitutions) | 11 | 78 strand-agnostic doublet classes (Bergstrom 2019) |
| **ID** (Insertion/Deletion) | 18 | 83 categories (indel length x repeat context x microhomology) |
| **CN** (Copy Number) | 24 | 48 channels (total CN x heterozygosity x segment length; Drews 2022 *Nature*) |
| **SV** (Structural Variants) | 16 | 32 channels (cluster x length x type) |

**Recent splits to know:**
- **SBS40 -> SBS40a / SBS40b / SBS40c** (Senkin 2024 *Nature*): pan-cancer active (40a); RCC-specific (40b/c).
- **SBS17 -> SBS17a (T>C uncertain) / SBS17b (T>G in CTT, 5-FU)** (Christensen 2019 *Nat Commun*; Pich 2019 *Nat Genet*).
- **SBS7 -> SBS7a / 7b / 7c / 7d** (UV photoproduct chemistry; Alexandrov 2020).
- **SBS10 -> SBS10a (POLE P286R) / 10b (POLE V411L) / 10c / 10d (POLD1)** (Mertz 2020 *Mol Cell*).

## Etiology Table (Postdoc-grade)

| Signature | Etiology | Clinical implication | Notes |
|-----------|----------|---------------------|-------|
| **SBS1** | Spontaneous 5mC deamination at CpG | Age-correlated; mitotic-rate biomarker | Clock-like |
| **SBS5** | UNKNOWN, clock-like; age-correlated | -- | Reviewer-accepted: "unknown, clock-like"; NOT polymerase fidelity errors |
| **SBS2 / SBS13** | APOBEC (APOBEC3A dominant per Petljak 2022) | Often co-occur; kataegis; ICI response signal | A3A vs A3B via YTCA vs RTCA tetranucleotide ratio |
| **SBS3** | HRD (BRCA1/2 deficient flat profile) | **PARP inhibitor eligibility** | HRDetect 98.7% sensitivity (Davies 2017) |
| **ID6** | HRD microhomology-mediated deletions | PARP eligibility | Pairs with SBS3 |
| **CN17 (HRD-CN1)** | HRD chromosomal instability | PARP eligibility; also BRCA1 promoter hypermethylation | Drews 2022 |
| **SBS6 / 14 / 15 / 20 / 21 / 26 / 44 + ID1 / 2** | MMR-D | ICI eligibility | Lynch typically 6/15/26/44; sporadic MLH1-hyperMet typically 21/26 |
| **SBS14 + SBS20** | POLE+MMR or POLD1+MMR double defect | Ultra-hypermutator; ICI excellent response | >500 mut/Mb |
| **SBS10a / 10b** | POLE-exo P286R / V411L | Hypermutator; ICI excellent response | 100-300 mut/Mb pure POLE |
| **SBS10c / 10d** | POLD1 | -- | Less common |
| **SBS28** | POLE indirect | Often co-extracted with SBS10 | -- |
| **SBS4 + DBS2** | Tobacco smoking; benzo[a]pyrene-G adducts | Lung cancer | C>A bias |
| **SBS7a/b/c/d + DBS1** | UV (CPD vs 6-4 photoproduct chemistry) | Melanoma | CC>TT dipyrimidine, CC>AA |
| **SBS24** | Aflatoxin | HCC (geographic) | C>A at CpC; Schulze 2015 *Nat Genet* |
| **SBS22** | Aristolochic acid | UTC, HCC | T>A at CpTpG; Hoang 2013 *Sci Transl Med* |
| **SBS17b** | 5-Fluorouracil | Therapy-induced | T>G in CTT context |
| **SBS31 / 35 / 86 / 87** | Platinum chemotherapy | Therapy-induced; second cancers | Cisplatin / carboplatin / oxaliplatin |
| **SBS11** | Temozolomide | Glioma post-TMZ | C>T at unmethylated CpC/CpT |
| **SBS88 + ID18** | Colibactin (pks+ E. coli) | CRC etiology; NTHL1-syndrome backgrounds | Pleguezuelos-Manzano 2020 *Nature* |
| **SBS30** | NTHL1 BER deficiency | Lynch-like; cancer predisposition | High cosine to FFPE artifact |
| **SBS-FFPE-artifact** | Formalin-induced C>T (NOT SBS33 as commonly cited) | Sequencing artifact | ~0.90 cosine to SBS30 (formalin-induced C>T characterization in mutational-signatures literature; specific paper attribution removed pending verification) |

**CRITICAL CORRECTION:** The widely-cited "SBS33 = FFPE artifact" is wrong. Modern literature attributes FFPE artifact to a signature resembling SBS30 (NTHL1-BER-deficiency profile); after enzymatic uracil repair the artifact instead resembles SBS1.

## Tool Taxonomy

| Tool | Approach | Class coverage | When to use | Fails when |
|------|----------|----------------|-------------|-----------|
| **SigProfilerSuite** (Alexandrov lab; Bergstrom 2019 *BMC Genomics*; Islam 2022 *Cell Genomics*; Diaz-Gay 2023 *Bioinformatics*) | Matrix gen -> NMF de novo / forward-backward refit | SBS / DBS / ID / CN / SV | Field standard; CPIC-equivalent for signatures | Heavy compute for de novo (100 NMF replicates) |
| **MutationalPatterns** (Manders 2022 *BMC Genomics*) | R-based; strict refit + NMF de novo | SBS / DBS / ID; lesion segregation | R workflows; reproducible refit | Lacks SV signatures |
| **MuSiCal** (Liu 2024 *Nat Genet*) | mvNMF (minimum-volume NMF) addressing NMF non-uniqueness | SBS / DBS / ID | Mid-size cohorts; novel signatures suspected | Less benchmarking at very large scale |
| **SigNet** (Serrano 2023, bioRxiv) | ANN-based signature attribution | SBS | Low mutation counts (best in 13-tool benchmark) | New tool; reproducibility data still maturing |
| **YAPSA** (Hubschmann 2020) | Linear combination decomposition | SBS | Comparison runs | Less widely used |
| **MutSignatures** (Fantini 2020) | Probabilistic refits | SBS | -- | -- |
| **deconstructSigs** (Rosenthal 2016) | NNLS (unregularized) | SBS | **DEPRECATED; never use** | NNLS overfits onto reference set; superseded by SigProfilerAssignment |
| **mSigAct** | Signatures from RNA-seq | SBS | RNA-seq only contexts | Limited resolution |
| **Helmsman** (Carlson 2018) | Fast matrix construction | SBS / DBS / ID | Preprocessing step only | Not for extraction/refit |
| **HRDetect** (Davies 2017 *Nat Med*) | Lasso logistic on 6 features (SBS3 / ID6 / RS3 / RS5 / HRD-LOH) | HRD-specific | BRCA1/2 deficiency classifier | Breast/ovarian-trained; cross-cancer needs revalidation |
| **MutationTimer** (Gerstung 2020 *Nature*) | Mutation timing relative to CN states | SBS | PCAWG-style evolution | Requires Battenberg/ASCAT CN; >=30x coverage |

**The deprecation:** deconstructSigs is the most-cited signature tool in publications but is **operationally deprecated**. NNLS without regularization overfits onto the ~70-signature reference; reviewers flag manuscripts using it without SigProfilerAssignment sensitivity. Replace with SigProfilerAssignment or MutationalPatterns strict refit.

## De Novo vs Refit: The Field's Most-Contested Choice

**Degasperi 2022** *Science* (12,222 WGS, UK 100k Genomes) argued refitting underestimates novel signatures because variance is forced onto existing references. They identified 40 additional SBS and 18 DBS signatures by full de novo extraction.

**Operational rule:**

| Mutation count per sample | Cohort size | Approach |
|---------------------------|-------------|----------|
| > 200 (SBS96) | N >= 50 (or 100 for DBS/ID/CN) | **De novo extraction** (SigProfilerExtractor, MuSiCal); validate via split-sample CV + bootstrap stability |
| > 200 | N < 50 | Refit (SigProfilerAssignment) |
| 50-200 | Any | Refit only; flag low confidence |
| < 50 | Any | **Do not attempt single-sample signature analysis** |

**SigProfilerExtractor stability gates:**
- `nmf_replicates = 100` (default 100, do not reduce)
- minimum stability >= 0.2 per signature
- minimum average stability >= 0.8 across signatures
- combined stability == 1.0 for selected rank

Manuscripts reporting extraction without these stability values are unreviewable.

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Single tumor WGS, > 200 mutations | SigProfilerAssignment refit | Single-sample de novo is unstable |
| Cohort >= 50 WGS, novel etiology suspected | SigProfilerExtractor de novo + cross-validate | Capture potentially novel signatures |
| Cohort >= 50 WGS, established cancer type | SigProfilerAssignment refit | Field consensus; fast |
| Mid-size cohort with novel signatures | MuSiCal mvNMF | Handles NMF non-uniqueness |
| Low mutation count (<100/sample) | SigNet (Serrano 2023) | Best low-count performance per Pancotti et al 2025 *Brief Bioinform* benchmark |
| BRCA1/2 deficiency screen | HRDetect (Davies 2017) | 6-feature lasso classifier; 98.7% sensitivity |
| Tumor evolution / mutation timing | MutationTimer (Gerstung 2020) | Requires Battenberg CN; PCAWG-validated |
| FFPE samples | SigProfilerAssignment with explicit FFPE-artifact handling | SBS30-like artifact; matched fresh-frozen controls ideal |
| WES (not WGS) | SigProfilerMatrixGenerator with `exome=True` | Trinucleotide-context correction for capture bias |
| RNA-seq only | mSigAct | Limited resolution; supplement with DNA-seq if available |
| HRD CN signatures | SigProfilerExtractor CN mode + CN17 (HRD-CN1) | Drews 2022 framework |
| Cross-cancer signature comparison | SigProfilerSuite with strand bias on | Aristolochic-acid SBS22 shows strong transcribed-strand bias |

## Standard Workflow: SigProfilerSuite

**Goal:** Generate 96-context mutation matrix, extract de novo signatures with stability validation, and decompose to COSMIC v3.4 reference.

**Approach:** Three-step pipeline with explicit version pinning and stability gates.

```python
# Step 1: Install reference genome (one-time)
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38')

# Step 2: Generate matrix
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
matrices = matGen.SigProfilerMatrixGeneratorFunc(
    project='cohort_2026',
    genome='GRCh38',
    vcfFiles='/path/to/vcf_directory',
    plot=True,
    exome=False,  # True if WES; corrects trinucleotide capture bias
    bed_file=None,  # Restrict to BED region if panel
    chrom_based=False,
    tsb_stat=True  # Transcribed-strand statistics
)
```

```python
# Step 3a (cohort >= 50): de novo extraction with stability gates
from SigProfilerExtractor import sigpro as sig
sig.sigProfilerExtractor(
    input_type='matrix',
    input_data='cohort_2026/output/SBS/cohort_2026.SBS96.all',
    output='extraction_output',
    reference_genome='GRCh38',
    opportunity_genome='GRCh38',
    minimum_signatures=1,
    maximum_signatures=12,
    nmf_replicates=100,            # Required for stability
    cpu=-1,
    seeds='random',
    matrix_normalization='gmm',
    resample=True,
    batch_size=1,
    refit_denovo_signatures=True,
    cosmic_version=3.4              # Match to current COSMIC release
)
```

```python
# Step 3b (single sample or cohort < 50): refit to COSMIC
from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(
    samples='cohort_2026/output/SBS/cohort_2026.SBS96.all',
    output='assignment_output',
    input_type='matrix',
    genome_build='GRCh38',
    cosmic_version=3.4,
    signature_database='SBS_GRCh38_GRCh38',  # Verify against the SigProfilerAssignment release; the bundled
                                              # COSMIC signature-database identifiers change between versions.
    nnls_add_penalty=0.05,         # Forward-add gate
    nnls_remove_penalty=0.01,      # Backward-remove gate
    initial_remove_penalty=0.05,
    refit_denovo_signatures=False,
    make_plots=True,
    sample_reconstruction_plots=True
)
```

## MutationalPatterns Strict Refit (R Alternative)

**Goal:** Same as SigProfilerAssignment but in R; suited for Bioconductor pipelines.

**Approach:** Cosine-based stopping reduces overfitting vs deconstructSigs.

```r
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load VCFs as GRanges
vcf_files <- list.files('vcf_dir', pattern = '\\.vcf$', full.names = TRUE)
sample_names <- gsub('\\.vcf$', '', basename(vcf_files))
vcfs <- read_vcfs_as_granges(vcf_files, sample_names,
                              ref_genome = 'BSgenome.Hsapiens.UCSC.hg38')

# Generate 96-context matrix
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = 'BSgenome.Hsapiens.UCSC.hg38')

# Fit to COSMIC v3.4 with strict refit (cosine-stopping; avoids deconstructSigs overfit)
signatures <- get_known_signatures(muttype = 'snv', source = 'COSMIC_v3.4',
                                    sig_type = 'reference', genome = 'GRCh38')
strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)

# Plot relative + absolute contributions
plot_contribution(strict_refit$fit_res$contribution, signatures, mode = 'relative')
```

## HRDetect for BRCA1/2 Deficiency

**Goal:** Classify tumors as HRD vs HR-proficient using the 6-feature Davies 2017 lasso.

**Approach:** Compute SBS3, ID6 (microhomology-mediated deletions), RS3 (rearrangement signature 3), RS5, HRD-LOH score; apply lasso classifier.

```r
# Davies 2017 HRDetect framework
# Features: SBS3, ID6/ID8, RS3, RS5, HRD-LOH
# Output: probability of HRD; threshold 0.7 = HRD-positive
library(signature.tools.lib)

hrdetect <- HRDetect_pipeline(
    SNV_vcf_files = snv_vcfs,
    Indels_vcf_files = indel_vcfs,
    SV_bedpe_files = sv_bedpes,
    CNV_tab_files = cnv_tables,
    genome.v = 'hg38',
    nparallel = 8
)

# hrdetect$hrdetect_output has BRCA_prob per sample
# >= 0.7 = HRD-positive; consider PARP inhibitor
```

## Per-Operation Failure Modes

**1. Single-sample de novo extraction**
- Trigger: Run SigProfilerExtractor on a cohort of 1.
- Mechanism: NMF requires multiple samples to find stable rank; single-sample 96-context spectrum has unstable signature decomposition.
- Symptom: Tool runs but signatures are noisy and inconsistent across replicates.
- Fix: Use refit (SigProfilerAssignment) for cohorts < 50; never de novo on single samples.

**2. Sub-100-mutation sample analyzed individually**
- Trigger: Calculate signatures for a tumor with <100 mutations.
- Mechanism: 96-context SBS spectrum needs 200-500 mutations for stable estimation; sub-100 produces signal-to-noise dominated by stochastic context distribution.
- Symptom: Random or implausible signature contributions.
- Fix: Aggregate samples in a meta-tumor for cohort analysis; for single-sample at low count consider SigNet which is optimized for low counts.

**3. FFPE artifact misclassified as SBS30 / SBS33**
- Trigger: Pipeline reports SBS33 (or SBS30) as biologically meaningful.
- Mechanism: FFPE-induced C>T deamination produces a profile resembling SBS30 (~0.90 cosine); after enzymatic uracil repair resembles SBS1. Pre-2022 literature incorrectly cited SBS33.
- Symptom: False NTHL1-BER-deficiency or "unknown SBS33" reports in cohorts using FFPE samples without matched controls.
- Fix: Run matched fresh-frozen controls in cohort; flag FFPE samples for separate analysis; use enzymatic-uracil pretreatment; expect SBS30-like artifact, not SBS33.

**4. WES + signature analysis without trinucleotide correction**
- Trigger: Run SigProfilerExtractor on WES VCFs with `exome=False`.
- Mechanism: WES capture has biased trinucleotide composition vs whole genome.
- Symptom: Apparent signature differences from WGS-derived signatures are artifactual.
- Fix: Set `exome=True`; this triggers trinucleotide-context correction.

**5. Refit chosen for cohort with novel etiology**
- Trigger: Tropical-region cohort with putative novel mutagen exposure; refit to COSMIC.
- Mechanism: Refit constrains variance onto existing catalog; novel signatures appear as residual or are decomposed onto closest-cosine known signatures.
- Symptom: Apparent absence of novel etiology despite biological hypothesis.
- Fix: For cohorts >= 50 run de novo extraction with cross-validation; Senkin 2024 kidney cancer cohort exemplifies the gain.

**6. APOBEC SBS2 vs SBS13 conflation; A3A vs A3B**
- Trigger: Report "APOBEC activity" without subtype.
- Mechanism: Petljak 2022 *Nature* established APOBEC3A as dominant active deaminase; the SBS2/SBS13 ratio reflects REV1-dependent translesion synthesis.
- Symptom: Lose mechanistic insight; potential mis-attribution.
- Fix: Distinguish A3A (YTCA 5' tetranucleotide preference) vs A3B (RTCA); cite Petljak 2022.

**7. SBS5 attributed to "polymerase fidelity errors"**
- Trigger: Manuscript claims SBS5 = replication errors.
- Mechanism: SBS5 etiology is contested; clock-like, age-correlated, modulated by ERCC2/TC-NER but not established as polymerase errors.
- Symptom: Reviewers reject.
- Fix: Report as "unknown, clock-like" until field consensus.

**8. Cross-version comparison without re-extraction**
- Trigger: Compare SBS40 contributions from v3.2 (single signature) vs v3.4 (split into 40a/b/c).
- Mechanism: Signature splits/merges occur between versions; cross-version exposures are not directly comparable.
- Symptom: Apparent "loss" or "gain" of activity due to renaming.
- Fix: Re-run with current COSMIC version; document version in methods.

**9. Strand bias ignored**
- Trigger: SBS22 (aristolochic acid) reported without transcribed-strand bias.
- Mechanism: Aristolochic-acid mutagenesis strongly biased toward transcribed strand; SigProfiler handles via `tsb_stat=True`, MutationalPatterns supports, deconstructSigs ignores.
- Symptom: Mis-attribution; loss of mechanistic evidence.
- Fix: Use SigProfiler or MutationalPatterns; enable strand-bias output.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| SigProfiler de novo vs refit different signatures | De novo finds novel + refit forces onto COSMIC | If cohort >= 50, prefer de novo; document |
| MutationalPatterns vs SigProfilerAssignment different exposures | NNLS vs forward-backward selection differences | Compare cosine similarity to mut_mat; pick better reconstruction |
| HRDetect calls HR-deficient + BRCA wildtype | BRCA1 promoter hypermethylation; PALB2 / FBXW7 / CDK12 alterations | Confirm with HRD-LOH score; assay BRCA1 methylation |
| Cosine to SBS3 high but ID6 absent | Single-feature HRD signal insufficient | Use HRDetect 6-feature classifier, not SBS3 alone |
| APOBEC signature present + low TMB | Cohort has APOBEC but not hypermutator | Both can coexist; YTCA/RTCA discriminates A3A vs A3B |
| FFPE samples produce SBS30 / SBS33-like signal | Almost always artifact | Run matched FF controls; use enzymatic uracil pretreatment |
| Cohort signature contributions implausible | Sub-100-mutation samples included | Stratify by mutation count; report >=200 separately |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| SBS96 stable extraction | >=200 mutations per sample | Alexandrov 2020 |
| De novo extraction cohort | N >= 50 (SBS); N >= 100 (DBS/ID/CN) | Field consensus |
| nmf_replicates | 100 (default; do not reduce) | SigProfilerExtractor |
| Stability gate | minimum stability >= 0.2; average >= 0.8 | SigProfilerExtractor defaults |
| Cosine similarity for "same signature" | > 0.85 (some use 0.90) | Convention |
| SBS-FFPE-artifact cosine to SBS30 | ~0.90 | formalin-induced C>T characterization (mutational-signatures literature; specific primary citation pending verification) |
| HRDetect threshold | BRCA_prob >= 0.7 = HRD-positive | Davies 2017 |
| POLE-exo + MMR mutation count | >500 mut/Mb (ultra-hypermutator) | Alexandrov 2020 |
| Pure POLE-exo mutation count | 100-300 mut/Mb | Alexandrov 2020 |
| MMR-D typical mutation count | 30-50 mut/Mb | Salem 2018 *Cancer Discov* |
| COSMIC version | 3.4 (September 2024) | COSMIC database |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| Tool reports SBS33 in FFPE cohort | Mis-attribution of FFPE artifact | Confirm by examining trinucleotide pattern; FFPE artifact resembles SBS30 in modern catalog |
| Single tumor signature attribution unstable | Sub-200-mutation sample analyzed alone | Aggregate; use SigNet for low-count |
| Refit ignores novel etiology | Forced onto COSMIC reference | Run de novo on cohort if N >= 50 |
| HRDetect false negative | Missing one of 6 features (ID6, RS3, RS5, HRD-LOH) | Confirm all features computed; assay BRCA1 methylation |
| Strand bias not detected | Tool/setting ignores transcribed-strand | Use SigProfilerSuite with `tsb_stat=True` or MutationalPatterns |
| WES vs WGS signatures differ | Capture-bias trinucleotide composition | Set `exome=True` in SigProfilerMatrixGenerator |
| Platinum-treated tumor: SBS31 vs SBS35 confusion | Both attributed to platinum; SBS35 closer to direct Drost lab signature | Report both; cosine to direct |
| Aristolochic-acid signature in non-exposure context | Bias from highly-expressed transcribed-strand artifacts | Check geographic + clinical history |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why not deconstructSigs?" | Operationally deprecated; NNLS overfits. Use SigProfilerAssignment or MutationalPatterns strict refit. |
| "Single tumor signatures meaningless?" | Below 200 mutations: yes. We aggregate cohorts and run refit for low-mutation samples. |
| "The de novo NMF rank choice?" | nmf_replicates=100 with stability gates; minimum >=0.2, average >=0.8, combined =1.0; SigProfilerExtractor defaults. |
| "FFPE samples bias signatures" | Modern attribution: FFPE artifact resembles SBS30 (not SBS33); we run matched FF controls or use enzymatic uracil pretreatment. |
| "SBS5 etiology; 'unknown' is unsatisfying" | Tomasetti-Vogelstein clock model; Druck 2026 FHIT + TC-NER; field has not converged; we report "unknown, clock-like". |
| "Why no APOBEC subtype distinction?" | Reported via YTCA vs RTCA tetranucleotide ratio per Petljak 2022; not all tools surface this; we used SigProfilerTopography. |
| "Cohort cross-comparison with old paper" | Re-extracted with COSMIC v3.4; signature splits (SBS40 -> 40a/b/c, SBS17 -> 17a/b) make pre-2024 exposures non-comparable. |
| "HRDetect cross-cancer validation" | Original Davies 2017 trained on breast cancer; we revalidated in our cohort with cross-cancer HRD-LOH score. |

## References

- Alexandrov LB et al. 2020. The repertoire of mutational signatures in human cancer. *Nature* 578:94. (PCAWG)
- Tate JG et al. 2019. COSMIC: the Catalogue Of Somatic Mutations In Cancer. *Nucleic Acids Res* 47:D941. (COSMIC v89)
- Senkin S et al. 2024. Geographic variation of mutagenic exposures in kidney cancer genomes. *Nature* 629:910. (SBS40a/b/c split)
- Drews RM et al. 2022. A pan-cancer compendium of chromosomal instability. *Nature* 606:976. (CN signatures)
- Petljak M et al. 2022. Mechanisms of APOBEC3 mutagenesis in human cancer cells. *Nature* 607:799. (A3A dominance)
- Bergstrom EN et al. 2019. SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. *BMC Genomics* 20:685.
- Islam SMA et al. 2022. Uncovering novel mutational signatures by de novo extraction with SigProfilerExtractor. *Cell Genomics* 2:100179.
- Diaz-Gay M et al. 2023. SigProfilerAssignment: an R + Python package for assigning known mutational signatures to individual samples. *Bioinformatics* 39:btad756.
- Manders F et al. 2022. MutationalPatterns: the one-stop shop for the analysis of mutational processes. *BMC Genomics* 23:134.
- Liu Q et al. 2024. Identifying mutational signatures with minimum volume non-negative matrix factorization. *Nat Genet*.
- Davies H et al. 2017. HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures. *Nat Med* 23:517.
- Degasperi A et al. 2022. Substitution mutational signatures in whole-genome-sequenced cancers in the UK population. *Science* 376:abl9283.
- Christensen S et al. 2019. 5-Fluorouracil treatment induces characteristic T>G mutations in human cancer. *Nat Commun* 10:4571. (SBS17b)
- Pich O et al. 2019. The mutational footprints of cancer therapies. *Nat Genet* 51:1732.
- Hayward NK et al. 2017. Whole-genome landscapes of major melanoma subtypes. *Nature* 545:175. (UV signatures)
- Schulze K et al. 2015. Exome sequencing of hepatocellular carcinomas. *Nat Genet* 47:505. (Aflatoxin SBS24)
- Hoang ML et al. 2013. Mutational signature of aristolochic acid exposure as revealed by whole-exome sequencing. *Sci Transl Med* 5:197ra102.
- Pleguezuelos-Manzano C et al. 2020. Mutational signature in colorectal cancer caused by genotoxic pks+ E. coli. *Nature* 580:269. (Colibactin SBS88)
- Gerstung M et al. 2020. The evolutionary history of 2,658 cancers. *Nature* 578:122. (MutationTimer)
- Mertz TM et al. 2020. POLE mutation spectra are shaped by the mutant allele identity. *Mol Cell* 78:1116.
- (FFPE-induced C>T mutational artifact: the earlier "Guyard 2022 Nat Commun" attribution could not be verified -- consult current FFPE-artifact literature for a confirmed primary citation.)
- Pancotti C et al. 2025. Deep learning methods for the recovery of mutational signatures. *Brief Bioinform* 26:bbaf545.
- COSMIC Signatures: `https://cancer.sanger.ac.uk/signatures/`

## Related Skills

- clinical-databases/tumor-mutational-burden - TMB and ICI biomarker
- clinical-databases/msi-detection - MSI as MMR-D biomarker (paired with SBS6/15/26/44)
- variant-calling/variant-calling - Somatic VCF input
- variant-calling/variant-calling - Mutect2 / Strelka2 somatic upstream
- data-visualization/heatmaps-clustering - Signature contribution visualization
