---
name: bio-variant-calling-clinical-interpretation
description: Clinical variant interpretation using ClinVar, ACMG guidelines, and pathogenicity predictors. Prioritize variants for diagnostic and research applications. Use when interpreting clinical significance of variants.
tool_type: mixed
primary_tool: InterVar
---

## Version Compatibility

Reference examples tested with: Entrez Direct 21.0+, bcftools 1.19+, cyvcf2 0.30+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Clinical Variant Interpretation

Prioritize and interpret variants for clinical significance using databases and ACMG/AMP guidelines.

## Interpretation Framework

```
Annotated VCF
    │
    ├── Database Lookup
    │   ├── ClinVar (clinical assertions)
    │   ├── OMIM (disease associations)
    │   └── gnomAD (population frequency)
    │
    ├── Computational Predictions
    │   ├── SIFT, PolyPhen-2
    │   ├── CADD, REVEL
    │   └── SpliceAI
    │
    ├── ACMG Classification
    │   └── Pathogenic → Likely Pathogenic → VUS → Likely Benign → Benign
    │
    └── Prioritized Variant List
```

## ClinVar Annotation

**Goal:** Annotate variants with ClinVar clinical significance and filter by pathogenicity.

**Approach:** Download the ClinVar VCF, add CLNSIG/CLNDN/CLNREVSTAT fields with bcftools annotate, then filter by significance level.

**"Find pathogenic variants in my VCF"** → Cross-reference variants against ClinVar clinical assertions and extract those classified as pathogenic or likely pathogenic.

### Download ClinVar

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

### Annotate with bcftools

```bash
bcftools annotate \
    -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT \
    input.vcf.gz -Oz -o with_clinvar.vcf.gz
```

### Filter Pathogenic Variants

```bash
# Pathogenic or Likely pathogenic
bcftools view -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' \
    with_clinvar.vcf.gz -Oz -o pathogenic.vcf.gz

# Exclude benign
bcftools view -e 'INFO/CLNSIG~"Benign" || INFO/CLNSIG~"Likely_benign"' \
    with_clinvar.vcf.gz -Oz -o not_benign.vcf.gz
```

## ClinVar Significance Levels

| CLNSIG | Meaning | Action |
|--------|---------|--------|
| Pathogenic | Disease-causing | Report |
| Likely_pathogenic | Probably disease-causing | Report with caveat |
| Uncertain_significance | VUS | May report, needs follow-up |
| Likely_benign | Probably not disease-causing | Usually exclude |
| Benign | Not disease-causing | Exclude |
| Conflicting | Multiple interpretations | Manual review |

## ClinVar Review Status

| CLNREVSTAT | Stars | Meaning |
|------------|-------|---------|
| practice_guideline | 4 | Expert panel reviewed |
| reviewed_by_expert_panel | 3 | ClinGen expert reviewed |
| criteria_provided,_multiple_submitters | 2 | Consistent assertions |
| criteria_provided,_single_submitter | 1 | One submitter with criteria |
| no_assertion_criteria | 0 | No criteria provided |

```bash
# Filter for high-confidence assertions (2+ stars)
bcftools view -i 'INFO/CLNREVSTAT~"multiple_submitters" || \
    INFO/CLNREVSTAT~"expert_panel" || \
    INFO/CLNREVSTAT~"practice_guideline"' \
    with_clinvar.vcf.gz -Oz -o high_confidence.vcf.gz
```

## InterVar (ACMG Classification)

**Goal:** Classify variants according to ACMG/AMP guidelines using automated criteria evaluation.

**Approach:** Convert VCF to ANNOVAR format, run InterVar to evaluate 28 ACMG criteria, and output five-tier classification.

Automated ACMG/AMP variant classification.

### Installation

```bash
git clone https://github.com/WGLab/InterVar.git
cd InterVar
# Download databases per documentation
```

### Run InterVar

```bash
python Intervar.py \
    -i input.avinput \
    -o output \
    -b hg38 \
    -d humandb/ \
    --input_type=AVinput
```

### From VCF

```bash
# Convert VCF to ANNOVAR format
convert2annovar.pl -format vcf4 input.vcf > input.avinput

# Run InterVar
python Intervar.py -i input.avinput -o intervar_results -b hg38
```

## ACMG/AMP Criteria

### Pathogenic Criteria

| Code | Type | Description |
|------|------|-------------|
| PVS1 | Very Strong | Null variant in gene where LOF is disease mechanism |
| PS1-4 | Strong | Same AA change, functional studies, etc. |
| PM1-6 | Moderate | Hot spot, absent from controls, etc. |
| PP1-5 | Supporting | Co-segregation, computational evidence |

### Benign Criteria

| Code | Type | Description |
|------|------|-------------|
| BA1 | Stand-alone | AF >5% in gnomAD |
| BS1-4 | Strong | AF greater than expected, functional studies |
| BP1-7 | Supporting | Missense in gene with truncating mechanism |

### Classification Rules

These are guidelines, not absolute rules -- expert review is always required.

**Pathogenic:** (PVS1 AND >=1 PS) OR (>=2 PS) OR (1 PS AND >=3 PM) OR (>=2 PM AND >=2 PP) OR (1 PM AND >=4 PP)

**Likely Pathogenic:** (1 PVS1 AND 1 PM) OR (1 PS AND 1-2 PM) OR (1 PS AND >=2 PP) OR (>=3 PM) OR (2 PM AND >=2 PP) OR (1 PM AND >=4 PP)

**VUS:** Does not meet criteria for pathogenic, likely pathogenic, benign, or likely benign. This is not a classification of uncertainty about pathogenicity per se, but rather insufficient evidence to classify in either direction.

**Likely Benign:** (1 BS AND 1 BP) OR (>=2 BP)

**Benign:** (BA1 alone) OR (>=2 BS)

### Annotation Concordance and PVS1

Different annotation tools (VEP vs. SnpEff) disagree on loss-of-function classification in 33-44% of cases. This directly impacts PVS1 -- the strongest single evidence criterion. For clinical-grade interpretation, annotate with both VEP and SnpEff and flag discrepant LOF calls for manual review. A variant called LOF by only one tool should not receive PVS1 without further investigation.

## Population Frequency Filtering

**Goal:** Restrict to rare variants that could be disease-causing.

**Approach:** Filter by gnomAD allele frequency threshold appropriate for the disease model (dominant vs. recessive). Use filtering allele frequency (FAF) where available, as it accounts for sampling error and is more appropriate than raw AF for clinical filtering.

| Disease Model | AF Threshold | Rationale |
|---------------|-------------|-----------|
| Dominant | < 0.0001 (1/10,000) | Penetrant dominant variants cannot be common |
| Recessive | < 0.01 (1/100) | Carriers can be relatively common |
| BA1 (stand-alone benign) | > 0.05 | Too common to cause rare disease |

Population-specific frequencies must be checked, not just global AF. A variant rare globally may be common in one population (e.g., sickle cell HbS in West African populations). The gnomAD filtering allele frequency (FAF) field adjusts for population substructure and sampling uncertainty and is preferred over raw AF for clinical filtering.

```bash
# Dominant disease model (AF < 0.0001)
bcftools view -i 'INFO/gnomAD_AF<0.0001 || INFO/gnomAD_AF="."' \
    input.vcf.gz -Oz -o dominant_rare.vcf.gz

# Recessive disease model (AF < 0.01)
bcftools view -i 'INFO/gnomAD_AF<0.01 || INFO/gnomAD_AF="."' \
    input.vcf.gz -Oz -o recessive_rare.vcf.gz

# When FAF annotation is available (preferred for clinical use)
bcftools view -i 'INFO/gnomAD_FAF<0.0001 || INFO/gnomAD_FAF="."' \
    input.vcf.gz -Oz -o faf_filtered.vcf.gz
```

## Pathogenicity Score Filtering

**Goal:** Prioritize variants using computational pathogenicity predictors.

**Approach:** Filter by CADD PHRED score (deleteriousness) and REVEL score (missense pathogenicity), alone or in combination with ClinVar. No single predictor is sufficient for clinical classification; ACMG PP3/BP4 criteria allow computational evidence as supporting, not strong, evidence.

### Predictor Interpretation

| Predictor | Threshold | Notes |
|-----------|-----------|-------|
| CADD PHRED | >= 20 (top 1%), >= 30 (top 0.1%) | Low specificity (~12%); use for prioritization, not binary classification |
| REVEL | > 0.5 (balanced), > 0.7 (high specificity) | Best-performing single missense predictor overall |
| AlphaMissense | > 0.564 pathogenic, < 0.340 benign | Protein structure-based (AlphaFold2); not trained on pathogenicity labels, so orthogonal signal to sequence-based predictors |
| SpliceAI | > 0.5 suggests splice-altering | Important for non-coding and synonymous variants near splice sites |
| SIFT | < 0.05 deleterious | Older predictor; largely superseded by REVEL for missense |
| PolyPhen-2 | > 0.85 probably damaging | Missense-only; complements but does not replace REVEL |

CADD scores variants of all types (coding, non-coding, indels) on a unified scale but trades sensitivity for breadth. REVEL is restricted to missense variants but outperforms CADD for that subset. AlphaMissense provides structural context that sequence-based tools miss, making it particularly valuable for variants in poorly conserved regions where REVEL may underperform.

### CADD Scores

```bash
# CADD > 20 (top 1% -- prioritization, not diagnostic)
bcftools view -i 'INFO/CADD_PHRED>20' input.vcf.gz -Oz -o cadd_filtered.vcf.gz

# CADD > 30 (top 0.1% -- higher confidence)
bcftools view -i 'INFO/CADD_PHRED>30' input.vcf.gz -Oz -o highly_deleterious.vcf.gz
```

### REVEL Scores

```bash
# REVEL > 0.5 (balanced sensitivity/specificity)
bcftools view -i 'INFO/REVEL>0.5' input.vcf.gz -Oz -o revel_filtered.vcf.gz

# REVEL > 0.7 (high specificity, fewer false positives)
bcftools view -i 'INFO/REVEL>0.7' input.vcf.gz -Oz -o revel_stringent.vcf.gz
```

### SpliceAI Filtering

```bash
# SpliceAI delta score > 0.5 (likely splice-altering)
bcftools view -i 'INFO/SpliceAI_DS_AG>0.5 || INFO/SpliceAI_DS_AL>0.5 || \
    INFO/SpliceAI_DS_DG>0.5 || INFO/SpliceAI_DS_DL>0.5' \
    input.vcf.gz -Oz -o splice_candidates.vcf.gz
```

### Combined Filtering

```bash
bcftools view -i '(INFO/CADD_PHRED>20 || INFO/REVEL>0.5) && \
    (INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely" || INFO/CLNSIG=".")' \
    input.vcf.gz -Oz -o prioritized.vcf.gz
```

## Python: Clinical Prioritization

**Goal:** Implement a multi-criteria variant classification pipeline in Python.

**Approach:** Combine ClinVar lookups, population frequency, and computational scores (CADD, REVEL, AlphaMissense) into a tiered classification function. This is a prioritization helper for research triage, not a substitute for formal ACMG review. Clinical reporting requires expert curation of each criterion.

```python
from cyvcf2 import VCF, Writer

def classify_variant(variant):
    """Heuristic prioritization tier based on available annotations.
    Not equivalent to formal ACMG classification -- treats computational
    scores as supporting evidence only, consistent with PP3/BP4 criteria.
    """
    clnsig = str(variant.INFO.get('CLNSIG', ''))
    af = variant.INFO.get('gnomAD_AF', 0) or 0
    cadd = variant.INFO.get('CADD_PHRED', 0) or 0
    revel = variant.INFO.get('REVEL', 0) or 0
    alphamissense = variant.INFO.get('AM_pathogenicity', 0) or 0

    if 'Pathogenic' in clnsig and 'Likely' not in clnsig:
        return 'PATHOGENIC'
    if 'Likely_pathogenic' in clnsig:
        return 'LIKELY_PATHOGENIC'

    if 'Benign' in clnsig or af > 0.05:
        return 'BENIGN'

    predictor_hits = sum([cadd > 25, revel > 0.7, alphamissense > 0.564])
    if predictor_hits >= 2 and af < 0.0001:
        return 'LIKELY_PATHOGENIC'
    if predictor_hits >= 1 and af < 0.0001:
        return 'VUS_FAVOR_PATH'
    if predictor_hits >= 1 and af < 0.01:
        return 'VUS_ELEVATED'

    if cadd < 10 and revel < 0.3 and alphamissense < 0.340:
        return 'LIKELY_BENIGN'

    return 'VUS'

vcf = VCF('annotated.vcf.gz')
results = []
report_tiers = {'PATHOGENIC', 'LIKELY_PATHOGENIC', 'VUS_FAVOR_PATH', 'VUS_ELEVATED'}

for variant in vcf:
    classification = classify_variant(variant)
    if classification in report_tiers:
        results.append({
            'chrom': variant.CHROM, 'pos': variant.POS,
            'ref': variant.REF, 'alt': variant.ALT[0],
            'gene': variant.INFO.get('SYMBOL', 'Unknown'),
            'consequence': variant.INFO.get('Consequence', 'Unknown'),
            'classification': classification,
            'clnsig': variant.INFO.get('CLNSIG', '.'),
            'cadd': variant.INFO.get('CADD_PHRED', '.'),
            'revel': variant.INFO.get('REVEL', '.'),
            'alphamissense': variant.INFO.get('AM_pathogenicity', '.'),
            'af': variant.INFO.get('gnomAD_AF', '.')
        })

for r in results:
    print(f'{r["gene"]}\t{r["chrom"]}:{r["pos"]}\t{r["consequence"]}\t{r["classification"]}')
```

## VUS Re-analysis

ClinVar reclassifies approximately 7% of variants annually. ACMG recommends periodic re-analysis of VUS, particularly when new functional studies are published, ClinVar receives new submissions, or new population data becomes available (gnomAD updates). Automated re-analysis pipelines should re-annotate stored VCFs against the latest ClinVar release and flag VUS that have been reclassified since the original report.

```bash
# Re-annotate with updated ClinVar and identify reclassified VUS
bcftools annotate -a clinvar_latest.vcf.gz \
    -c INFO/CLNSIG_NEW:=INFO/CLNSIG \
    original_results.vcf.gz -Oz -o reannotated.vcf.gz

# Find variants where original was VUS but ClinVar now has a classification
bcftools view -i 'INFO/CLNSIG_ORIG~"Uncertain" && \
    (INFO/CLNSIG_NEW~"athogenic" || INFO/CLNSIG_NEW~"enign")' \
    reannotated.vcf.gz -Oz -o reclassified.vcf.gz
```

## Gene Panel Filtering

**Goal:** Restrict analysis to variants within a clinical gene panel.

**Approach:** Filter by BED coordinates or VEP gene symbol annotations to target specific genes.

```bash
# Filter to gene panel
bcftools view -R gene_panel.bed input.vcf.gz -Oz -o panel_variants.vcf.gz

# Or by gene symbol (requires VEP annotation)
bcftools view -i 'INFO/CSQ~"BRCA1" || INFO/CSQ~"BRCA2"' \
    input.vcf.gz -Oz -o brca_variants.vcf.gz
```

## Disease-Specific Resources

| Resource | Content | Use |
|----------|---------|-----|
| ClinVar | Clinical assertions | Primary lookup |
| OMIM | Gene-disease relationships | Gene prioritization |
| HGMD | Published mutations | Literature evidence |
| gnomAD | Population frequencies | Rarity filtering |
| ClinGen | Gene validity/dosage | LOF interpretation |

## Reporting Template

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SYMBOL\t%INFO/Consequence\t\
%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/gnomAD_AF\t%INFO/CADD_PHRED\n' \
    prioritized.vcf.gz > clinical_report.tsv
```

## Complete Workflow

**Goal:** Run an end-to-end clinical variant interpretation pipeline from annotation through reporting.

**Approach:** Chain ClinVar annotation, rare variant filtering, pathogenicity extraction, VUS review, and TSV report generation.

```bash
#!/bin/bash
set -euo pipefail

INPUT=$1
CLINVAR=$2
OUTPUT_PREFIX=$3

echo "=== Add ClinVar annotations ==="
bcftools annotate -a $CLINVAR \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNVC \
    $INPUT -Oz -o ${OUTPUT_PREFIX}_clinvar.vcf.gz

echo "=== Filter rare variants ==="
bcftools view -i 'INFO/gnomAD_AF<0.01 || INFO/gnomAD_AF="."' \
    ${OUTPUT_PREFIX}_clinvar.vcf.gz -Oz -o ${OUTPUT_PREFIX}_rare.vcf.gz

echo "=== Extract pathogenic/likely pathogenic ==="
bcftools view -i 'INFO/CLNSIG~"athogenic"' \
    ${OUTPUT_PREFIX}_rare.vcf.gz -Oz -o ${OUTPUT_PREFIX}_pathogenic.vcf.gz

echo "=== Extract high-impact VUS ==="
bcftools view -i 'INFO/CLNSIG~"Uncertain" && INFO/CADD_PHRED>20' \
    ${OUTPUT_PREFIX}_rare.vcf.gz -Oz -o ${OUTPUT_PREFIX}_vus_review.vcf.gz

echo "=== Generate report ==="
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SYMBOL\t%INFO/Consequence\t\
%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/gnomAD_AF\t%INFO/CADD_PHRED\n' \
    ${OUTPUT_PREFIX}_pathogenic.vcf.gz > ${OUTPUT_PREFIX}_report.tsv

echo "=== Complete ==="
echo "Pathogenic: ${OUTPUT_PREFIX}_pathogenic.vcf.gz"
echo "VUS for review: ${OUTPUT_PREFIX}_vus_review.vcf.gz"
echo "Report: ${OUTPUT_PREFIX}_report.tsv"
```

## Related Skills

- variant-calling/variant-annotation - VEP/SnpEff annotation pipelines for functional consequence prediction
- variant-calling/filtering-best-practices - Quality and artifact filtering prior to clinical interpretation
- variant-calling/vcf-basics - VCF format fundamentals and field extraction
- database-access/entrez-fetch - Programmatic download of ClinVar and OMIM data
- pathway-analysis/go-enrichment - Gene set and pathway analysis for variant gene lists
