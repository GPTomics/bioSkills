---
name: bio-workflows-neoantigen-pipeline
description: Orchestrates neoantigen discovery from somatic variants to ranked vaccine candidates, chaining HLA typing (OptiType/arcasHLA + LOHHLA), VEP annotation (Wildtype+Frameshift plugins) + expression/readcount annotation, proximal-variant phasing, pVACseq MHC-I/II binding, CCF/clonality, and immunogenicity/quality ranking. Use when recognizing that binding is single-digit PPV and the critical steps are downstream (full-resolution HLA + LOH gating, proximal-variant phasing, clonality from purity+CN not raw VAF, expression), sequencing normalize+annotate -> phase -> HLA -> binding -> quality in the defensible order, dropping candidates on LOH-lost alleles, supplying --phased-proximal-variants-vcf so the mutant peptide is real, or ranking WITHIN patient rather than a fixed IC50 threshold. Hands mechanism to the immunoinformatics component skills; not a re-teach of any single step.
tool_type: mixed
primary_tool: pVACtools
goal_approach_exempt: true
workflow: true
depends_on:
  - clinical-databases/hla-typing
  - immunoinformatics/mhc-binding-prediction
  - immunoinformatics/mhc-class-ii-prediction
  - immunoinformatics/neoantigen-prediction
  - immunoinformatics/immunogenicity-scoring
  - immunoinformatics/epitope-prediction
qc_checkpoints:
  - after_hla: "HLA types resolved to 4-digit, coverage adequate"
  - after_binding: "Predictions for ALL alleles (LOH-lost alleles dropped); ranked within patient, not hard-thresholded across patients"
  - after_neoantigen: "Expressed (RNA-confirmed); clonality via CCF from purity+CN, not raw VAF"
  - after_scoring: "Top candidates are a tier-1 hypothesis list for MS + T-cell validation"
---

## Version Compatibility

Reference examples tested with: Ensembl VEP 111+, pVACtools 4.1+ (Frameshift plugin REPLACED the legacy Downstream in 2.0+), MHCflurry 2.1+, NetMHCpan 4.1, OptiType 1.3+ / arcasHLA, LOHHLA, WhatsHap 2.0+ (phasing), matplotlib 3.8+, numpy 1.26+, pandas 2.2+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Neoantigen Pipeline

**"Predict neoantigens from my tumor sequencing data"** -> Orchestrate HLA typing (OptiType), somatic variant calling, pVACtools neoantigen prediction, MHC binding scoring, and immunogenicity-based candidate ranking for personalized cancer immunotherapy.

Complete workflow from somatic variants to ranked neoantigen vaccine candidates for personalized cancer immunotherapy.

## Key Judgment -- binding is the easy part; PPV lives downstream

A binding-only pipeline has single-digit-percent positive predictive value (TESLA; Wells 2020 Cell 183:818). The critical steps are downstream of binding: correct full-resolution HLA typing (wrong allele = confident garbage), HLA loss-of-heterozygosity (run LOHHLA and DROP candidates on a lost allele — it invalidates predictions silently), proximal-variant phasing (supply `--phased-proximal-variants-vcf` or the mutant peptide is wrong), cancer cell fraction for clonality (clonal beats subclonal; use purity + copy number, not raw VAF), expression, and quality features (agretopicity, foreignness). Treat the ranked output as a tier-1 hypothesis list for immunopeptidomics MS and functional T-cell validation, not a final answer. Add MHC class II (CD4) neoantigens for vaccine help (see immunoinformatics/mhc-class-ii-prediction). Note on DAI below: agretopicity is most often the WT/MT binding ratio; whichever form is used, an anchor-position mutation inflates it without changing the TCR-facing surface, and a barely-presented WT makes it unstable — pair it with anchor evaluation.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| HLA typing at full 4-digit resolution (class I + II) | A wrong allele is confident garbage; every binding prediction inherits it; reconcile DNA vs RNA calls |
| Variant source + somatic caller (matched-normal preferred) | Tumor-only calling leaks germline; indels/frameshifts are disproportionately valuable; expression must be RNA-confirmed and annotated INTO the VCF |
| Proximal-variant phasing | Without it the mutant peptide is one the tumor never makes; germline SNPs in cis are especially treacherous |
| HLA-LOH gate | Candidates on a lost allele are silently invalid (~17% pan-cancer, 30%+ HNSCC/NSCLC/cervical) |

## The canonical order and why

Somatic PASS calls -> normalize + VEP-annotate (Wildtype + Frameshift plugins) -> annotate expression + DNA/RNA readcounts INTO the VCF -> PHASE proximal variants -> HLA typing + LOHHLA -> MHC binding -> clonality (CCF from purity+CN) -> quality features -> tier/rank -> pVACview review.

- **Order-trap 1 - normalize + annotate with the RIGHT plugins BEFORE pVACseq.** pVACseq needs the Wildtype plugin (matched WT peptide -> agretopicity) and the Frameshift plugin (novel ORF); Frameshift REPLACED the legacy Downstream in pVACtools 2.0+. Normalize before annotate.
- **Order-trap 2 - PHASE proximal variants BEFORE translating the mutant peptide.** THE review-sinker: editing variants independently yields a peptide the patient never makes. Merge somatic+germline, phase (WhatsHap/GATK), supply `--phased-proximal-variants-vcf`.
- **Order-trap 3 - HLA typing (+ LOHHLA) BEFORE binding.** Binding is per-allele; a wrong or lost allele makes every downstream prediction garbage. Drop LOH-lost alleles before ranking.
- **Order-trap 4 - CCF/clonality from purity+copy-number BEFORE calling something subclonal.** Low purity makes clonal look subclonal; correct VAF to cancer-cell fraction (copy-number/allele-specific-copy-number). Clonal beats subclonal.
- **Order-trap 5 - rank WITHIN patient; do NOT hard-threshold IC50 across patients.** Immunogenicity scores are relative.

## Workflow Overview

```
Somatic VCF (annotated) + Tumor RNA-seq (optional)
        |
        v
[1. HLA Typing] --> arcasHLA / OptiType (if types not provided)
        |
        v
[2. MHC Binding Prediction] --> MHCflurry / NetMHCpan
        |
        v
[3. Neoantigen Calling] --> pVACseq
        |
        v
[4. Immunogenicity Scoring] --> Multi-factor ranking
        |
        v
Ranked Vaccine Candidates (TSV + visualizations)
```

## Prerequisites (Ensembl VEP 111+)

```bash
pip install pvactools mhcflurry vatools

mhcflurry-downloads fetch

conda install -c bioconda ensembl-vep arcas-hla optitype
```

## Primary Path: pVACseq Pipeline

### Step 1: HLA Typing (if not provided)

HLA types are critical for MHC binding prediction. If not already known from clinical testing:

```bash
# From tumor RNA-seq BAM
arcasHLA extract tumor.bam -t 8 -o hla_output/
arcasHLA genotype hla_output/tumor.extracted.1.fq.gz hla_output/tumor.extracted.2.fq.gz \
    -g A,B,C,DRB1,DQB1,DQA1,DPB1,DPA1 -t 8 -o hla_output/   # type the DQA1/DPA1 alpha chains too: NetMHCIIpan needs PAIRED DQ/DP alleles

# Parse results
cat hla_output/tumor.genotype.json
```

```python
import json

with open('hla_output/tumor.genotype.json') as f:
    hla_data = json.load(f)

hla_alleles = []
for gene, alleles in hla_data.items():
    for allele in alleles:
        # arcasHLA emits 3-field alleles (A*01:01:01); pVACseq/IEDB validate 2-field (HLA-A*01:01)
        hla_alleles.append('HLA-' + ':'.join(allele.split(':')[:2]))

# Format for pVACseq: HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,...
hla_string = ','.join(hla_alleles)
print(f'HLA alleles: {hla_string}')
```

### Step 2: VCF Annotation with VEP

pVACseq requires VEP-annotated VCF with specific fields:

```bash
# Annotate somatic VCF
vep --input_file somatic.vcf \
    --output_file somatic.vep.vcf \
    --format vcf --vcf --symbol --terms SO \
    --plugin Frameshift --plugin Wildtype \
    --offline --cache \
    --pick --fork 4

# Add expression data (optional but recommended)
# Positionals: <vcf> <expression_file> {kallisto,stringtie,cufflinks,custom} {gene,transcript}
vcf-expression-annotator somatic.vep.vcf \
    expression.tsv custom gene \
    -s tumor_sample --id-column gene_id --expression-column tpm \
    -o somatic.vep.expression.vcf

# PHASE proximal variants (the review-sinker). Merge somatic + germline, phase with WhatsHap,
# and pass the result to pVACseq via --phased-proximal-variants-vcf so a second variant in the
# same codon-window (esp. a germline SNP in cis) yields the peptide the tumor ACTUALLY makes.
whatshap phase -o phased.vcf.gz --reference reference.fa somatic_plus_germline.vcf.gz tumor.bam
tabix -p vcf phased.vcf.gz
```

### Step 3: Run pVACseq (Ensembl VEP 111+)

```bash
# Basic run with MHC Class I
pvacseq run \
    somatic.vep.vcf \
    tumor_sample \
    "HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*44:02,HLA-C*07:02,HLA-C*05:01" \
    MHCflurry MHCnuggetsI NetMHCpan \
    pvacseq_output/ \
    -e1 8,9,10,11 \
    --iedb-install-directory /path/to/iedb \
    -t 8

# With expression filtering
pvacseq run \
    somatic.vep.expression.vcf \
    tumor_sample \
    "HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*44:02" \
    MHCflurry NetMHCpan \
    pvacseq_output/ \
    -e1 8,9,10,11 \
    --phased-proximal-variants-vcf phased.vcf.gz \
    --tumor-purity 0.7 \
    --tdna-vaf 0.1 \
    --expn-val 1 \
    -t 8
```

Drop candidates on HLA-LOH-lost alleles (run LOHHLA/DASH) BEFORE ranking, and correct clonality to cancer-cell fraction (CCF from purity + copy number, not raw VAF; see copy-number/allele-specific-copy-number) — the raw-VAF filter below is a coarse proxy.

### Step 4: Filter and Rank Candidates

```python
import pandas as pd
import numpy as np

results = pd.read_csv('pvacseq_output/MHC_Class_I/tumor_sample.filtered.tsv', sep='\t')

# Binding affinity filter (IC50 <500nM considered strong binder)
# IC50 <500nM: strong binder; 500-5000nM: weak binder
strong_binders = results[results['Median MT IC50 Score'] < 500].copy()

# Differential agretopicity index (DAI): WT/MT IC50 ratio (== pVACtools Fold Change), matching the
# WT/MT ratio definition. DAI > 1 = MT binds better than WT (mutation created/improved binding); higher = more tumor-specific.
strong_binders['DAI'] = strong_binders['Median WT IC50 Score'] / strong_binders['Median MT IC50 Score']

# Expression filter (if available)
if 'Gene Expression' in strong_binders.columns:
    # TPM >1 ensures detectable expression
    strong_binders = strong_binders[strong_binders['Gene Expression'] > 1]

# VAF filter: prioritize clonal mutations
# VAF >0.1 ensures mutation present in substantial tumor fraction
strong_binders = strong_binders[strong_binders['Tumor DNA VAF'] > 0.1]

# Multi-factor scoring
def immunogenicity_score(row):
    score = 0
    # Strong binding (IC50 <150nM is very strong)
    if row['Median MT IC50 Score'] < 150:
        score += 3
    elif row['Median MT IC50 Score'] < 500:
        score += 2

    # High DAI (tumor-specificity). DAI is the WT/MT IC50 ratio: >1 = MT binds better than WT.
    if row['DAI'] > 10:
        score += 2
    elif row['DAI'] > 2:
        score += 1

    # Clonal mutation (high VAF)
    if row['Tumor DNA VAF'] > 0.3:
        score += 2
    elif row['Tumor DNA VAF'] > 0.15:
        score += 1

    # Expressed (if available)
    if 'Gene Expression' in row.index and row['Gene Expression'] > 10:
        score += 1

    return score

strong_binders['Immunogenicity Score'] = strong_binders.apply(immunogenicity_score, axis=1)

# Rank by composite score
ranked = strong_binders.sort_values('Immunogenicity Score', ascending=False)

# Top candidates for vaccine
top_candidates = ranked.head(20)
top_candidates.to_csv('top_neoantigen_candidates.tsv', sep='\t', index=False)

print(f'Total strong binders: {len(strong_binders)}')
print(f'Top 20 candidates exported')
print(ranked[['Gene Name', 'MT Epitope Seq', 'HLA Allele', 'Median MT IC50 Score', 'DAI', 'Immunogenicity Score']].head(10))
```

### Step 5: MHC Class II Neoantigens (CD4+ T cell help)

```bash
pvacseq run \
    somatic.vep.vcf \
    tumor_sample \
    "DRB1*01:01,DRB1*07:01,DQA1*05:01-DQB1*02:01,DQA1*03:01-DQB1*03:01" \
    MHCnuggetsII NetMHCIIpan \
    pvacseq_class2_output/ \
    -e2 15 \
    --iedb-install-directory /path/to/iedb \
    -t 8
```

## Alternative: Standalone MHCflurry

For quick binding predictions without full pVACseq pipeline:

```python
from mhcflurry import Class1PresentationPredictor

predictor = Class1PresentationPredictor.load()

peptides = ['SIINFEKL', 'GILGFVFTL', 'NLVPMVATV']
alleles = ['HLA-A*02:01', 'HLA-B*07:02']

results = predictor.predict(peptides=peptides, alleles=alleles,
                            include_affinity_percentile=True, verbose=0)
print(results[['peptide', 'best_allele', 'presentation_score', 'affinity', 'affinity_percentile']])
```

## Visualization

```python
import matplotlib.pyplot as plt
import seaborn as sns

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# IC50 distribution
ax1 = axes[0]
ax1.hist(ranked['Median MT IC50 Score'], bins=50, edgecolor='black')
ax1.axvline(500, color='red', linestyle='--', label='500nM threshold')
ax1.set_xlabel('Median MT IC50 (nM)')
ax1.set_ylabel('Count')
ax1.set_title('Binding Affinity Distribution')
ax1.legend()

# DAI vs IC50
ax2 = axes[1]
scatter = ax2.scatter(ranked['Median MT IC50 Score'], ranked['DAI'],
                      c=ranked['Immunogenicity Score'], cmap='viridis', alpha=0.7)
ax2.set_xlabel('MT IC50 (nM)')
ax2.set_ylabel('Differential Agretopicity Index')
ax2.set_title('Tumor Specificity vs Binding')
plt.colorbar(scatter, ax=ax2, label='Immunogenicity Score')

# Top genes
ax3 = axes[2]
gene_counts = ranked['Gene Name'].value_counts().head(15)
gene_counts.plot(kind='barh', ax=ax3)
ax3.set_xlabel('Number of Neoantigens')
ax3.set_title('Top Genes with Neoantigens')

plt.tight_layout()
plt.savefig('neoantigen_summary.pdf')
```

## Parameter Recommendations

| Step | Parameter | Value | Rationale |
|------|-----------|-------|-----------|
| pVACseq | -e1 | 8,9,10,11 | MHC-I binds 8-11mer peptides |
| pVACseq | -e2 | 15 | MHC-II binds 13-25mer, 15 is core |
| Filtering | IC50 | <500nM | Standard strong binder threshold |
| Filtering | VAF | >0.1 | Ensures clonal representation |
| Filtering | Expression | >1 TPM | Detectable transcription |
| Ranking | DAI (WT/MT IC50 ratio) | >2 moderate, >10 strong | MT binds better than WT (>1); higher = more tumor-specific |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Peptides the tumor never makes | Proximal variants edited independently (unphased) | `--phased-proximal-variants-vcf` (WhatsHap/GATK); include germline in cis |
| Frameshift ORFs lost / no agretopicity | Wrong/legacy VEP plugin (Downstream instead of Frameshift; missing Wildtype) | `pvacseq install_vep_plugin`; run `--plugin Wildtype --plugin Frameshift` |
| Confident but invalid predictions | HLA allele wrong or on a LOH-lost haplotype | Full 4-digit typing + LOHHLA drop before ranking |
| `--expn-val`/VAF filters silently pass everything | Expression/readcounts not annotated into the VCF | `vcf-expression-annotator` + `vcf-readcount-annotator` before pVACseq |
| Clonal candidate mis-tiered subclonal | Raw VAF used as clonality on a low-purity tumor | CCF from purity + copy number (copy-number/allele-specific-copy-number) |
| Candidates mis-ranked across patients | Fixed IC50 threshold applied cross-patient | Rank WITHIN patient (immunoinformatics/immunogenicity-scoring) |
| No neoantigens found | Low mutation burden | Lower IC50 threshold to 1000nM; check TMB/MSI first |

## References

- Hundal J, Kiwala S, Feng YY, et al (2020) pVACtools: a computational toolkit to identify and visualize cancer neoantigens. *Cancer Immunology Research* 8:409-420. DOI 10.1158/2326-6066.CIR-19-0401.
- Wells DK, van Buuren MM, Dang KK, et al (2020) Key parameters of tumor epitope immunogenicity revealed through a consortium approach improve neoantigen prediction (TESLA). *Cell* 183:818-834. DOI 10.1016/j.cell.2020.09.015. (single-digit PPV of binding-only.)
- McGranahan N, Rosenthal R, Hiley CT, et al (2017) Allele-specific HLA loss and immune escape in lung cancer evolution. *Cell* 171:1259-1271. DOI 10.1016/j.cell.2017.10.001. (LOHHLA.)
- Wood MA, Nguyen A, Struck AJ, et al (2020) neoepiscope improves neoepitope prediction with multivariant phasing. *Bioinformatics* 36:713-720. DOI 10.1093/bioinformatics/btz653. (phasing matters.)

## Output Files

| File | Description |
|------|-------------|
| `*.filtered.tsv` | pVACseq filtered neoantigens |
| `*.all_epitopes.tsv` | All predicted epitopes |
| `top_neoantigen_candidates.tsv` | Ranked vaccine candidates |
| `neoantigen_summary.pdf` | Visualization figures |

## Related Skills

- immunoinformatics/mhc-binding-prediction - MHCflurry parameters; BA vs EL, %Rank vs nM, abundance bias
- immunoinformatics/mhc-class-ii-prediction - class II (CD4) neoantigens for vaccine help
- immunoinformatics/neoantigen-prediction - pVACtools details; LOHHLA, phasing, clonality
- immunoinformatics/immunogenicity-scoring - rank within patient (don't threshold); fitness-model quality
- immunoinformatics/epitope-prediction - B-cell epitopes
- clinical-databases/hla-typing - HLA typing (T1K is the 2024-2026 all-rounder; OptiType for class I; arcasHLA for RNA-seq); check HLA-LOH via LOHHLA / DASH which abolishes neoantigen presentation in ~17% pan-cancer (~30%+ HNSCC / NSCLC / cervical)
- clinical-databases/tumor-mutational-burden - TMB-H pan-tumor ICI biomarker; check before neoantigen-vaccine candidate selection
- clinical-databases/msi-detection - MSI-H / dMMR pan-tumor ICI biomarker; MSI-H supersedes TMB-H per Sha 2020
- clinical-databases/somatic-signatures - Clonal neoantigen burden (McGranahan 2016 Science) predicts ICI response better than total TMB
- workflows/somatic-variant-pipeline - Upstream somatic calling
