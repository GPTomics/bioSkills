# Clinical Interpretation Usage Guide

## Overview

Classify a variant's clinical significance by assembling independent, calibrated evidence under the correct framework for a pinned context. Germline Mendelian variants use ACMG/AMP with the ClinGen SVI refinements (graded PVS1, PM2 downgraded to Supporting, PP5/BP6 retired, calibrated PP3/BP4, a Bayesian points sum); somatic/tumor variants use the AMP/ASCO/CAP actionability tiers and the ClinGen/CGC/VICC oncogenicity system. A ClinVar assertion and a gnomAD frequency are leads to weigh, not answers to adopt, and every classification carries an expiry date.

## Prerequisites

```bash
# ClinVar VCF (match your genome build; GRCh38 shown)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# Tools
conda install -c bioconda bcftools htslib
pip install cyvcf2

# InterVar (automated ACMG, 2015 rules; treat output as a starting point, not final)
git clone https://github.com/WGLab/InterVar.git
```

## Quick Start

Tell your AI agent what you want to do:
- "Classify this germline variant with current ACMG and write the criterion-by-criterion rationale"
- "Decide whether to use the germline or somatic framework for this tumor variant"
- "Check whether there is a ClinGen VCEP specification for this gene before I apply generic ACMG"
- "Is this ClinVar assertion usable as evidence, or just a lead?"
- "Apply the graded PVS1 decision tree to this stop-gained variant"
- "Judge this variant's gnomAD grpmax filtering AF against the disease's maximum credible AF"
- "Reanalyze my prior VUS against the latest ClinVar and gnomAD"

## ACMG points classification

Summed-points thresholds (Tavtigian 2020): Pathogenic >= 10, Likely Pathogenic 6-9, VUS 0-5, Likely Benign -1 to -6, Benign <= -7. Points by strength: Supporting +/-1, Moderate +/-2, Strong +/-4, Very Strong +/-8 (benign criteria subtract). Prefer this over the 2015 verbal combining rules; it accommodates the graded PVS1 and calibrated PP3/BP4 strengths.

## Downloading and annotating ClinVar

```bash
bcftools annotate -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNVC \
    input.vcf.gz -Oz -o with_clinvar.vcf.gz
```

Always carry CLNREVSTAT so a single-submitter (1-star) or conflicting call is never mistaken for evidence. Match the ClinVar build to your VCF build; ClinVar publishes separate GRCh37 and GRCh38 VCFs.

## InterVar (automated ACMG)

InterVar automates the 2015 combining rules. It does NOT apply all of the 2018-2025 ClinGen refinements or any gene-specific VCEP specification, so treat its output as a first pass to be corrected by hand (graded PVS1, PM2_Supporting, calibrated PP3/BP4).

```bash
convert2annovar.pl -format vcf4 input.vcf > input.avinput
python Intervar.py -i input.avinput -o intervar_results -b hg38 -d humandb/ --input_type=AVinput
```

## Gene panel filtering

```bash
# By BED coordinates
bcftools view -R gene_panel.bed input.vcf.gz -Oz -o panel_variants.vcf.gz

# By gene symbol (requires VEP CSQ annotation)
bcftools view -i 'INFO/CSQ~"BRCA1" || INFO/CSQ~"BRCA2"' input.vcf.gz -Oz -o brca_variants.vcf.gz
```

## Reporting and end-to-end workflow

```bash
#!/bin/bash
set -euo pipefail
INPUT=$1; CLINVAR=$2; OUT=$3

# 1. ClinVar leads (keep review status)
bcftools annotate -a "$CLINVAR" -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNVC \
    "$INPUT" -Oz -o "${OUT}_clinvar.vcf.gz"

# 2. Rare by grpmax filtering AF (recessive-model example; keep absent sites)
bcftools view -i 'INFO/fafmax_faf95_max<0.01 || INFO/fafmax_faf95_max="."' \
    "${OUT}_clinvar.vcf.gz" -Oz -o "${OUT}_rare.vcf.gz"

# 3. Pathogenic/likely-pathogenic leads
bcftools view -i 'INFO/CLNSIG~"athogenic"' "${OUT}_rare.vcf.gz" -Oz -o "${OUT}_path.vcf.gz"

# 4. Report (carry CLNREVSTAT and gene/consequence)
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SYMBOL\t%INFO/Consequence\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n' \
    "${OUT}_path.vcf.gz" > "${OUT}_report.tsv"
```

## Key databases

| Database | Purpose | Caveat |
|----------|---------|--------|
| ClinVar | Submitted clinical assertions | Assertions are leads; read stars + submitter + evidence |
| gnomAD | Population/grpmax frequencies | grpmax FAF vs disease-max AF; presence != benign |
| OMIM | Gene-disease relationships | Gene-level, not variant classification |
| ClinGen | Gene validity, dosage, VCEP specs | Check for a VCEP spec before generic ACMG |
| OncoKB / CIViC / COSMIC | Somatic actionability / oncogenicity / recurrence | Somatic only; read the evidence item, not the letter |
| HGMD | Published mutations | Literature lead; curation quality varies |

## Example Prompts

### Germline classification
> "Apply current ACMG (graded PVS1, PM2_Supporting, calibrated PP3/BP4) to this variant and give me the points total and tier"

> "Walk the PVS1 decision tree for this frameshift, including whether it escapes NMD on the MANE Select transcript"

### Evidence quality
> "This variant is Pathogenic 1-star in ClinVar with one submitter -- can I use that as evidence?"

> "Is this REVEL score enough for PP3, and at what strength?"

### Frequency
> "The variant is 0.4% globally but I suspect a founder allele -- check grpmax against the maximum credible AF for this recessive disease"

### Somatic
> "This is a tumor variant -- classify actionability with AMP/ASCO/CAP tiers and oncogenicity with the Horak framework, not ACMG"

### Reanalysis
> "Reanalyze my prior VUS list against the latest ClinVar and gnomAD and flag anything that now crosses a threshold"

## What the Agent Will Do

1. Decide the framework first: germline ACMG vs somatic tiers/oncogenicity, and check for a gene-specific VCEP specification.
2. Annotate ClinVar assertions as leads (carrying CLNREVSTAT) and gnomAD grpmax frequencies.
3. Apply current ACMG criteria: graded PVS1 (mechanism + NMD + MANE transcript), PM2_Supporting, retired PP5/BP6, one calibrated predictor for PP3/BP4.
4. Sum Bayesian points to a tier, or route a tumor variant to Li 2017 tiers and Horak 2022 oncogenicity.
5. Flag VUS for a reanalysis loop and record the pinned context (build, transcript, tool/db versions) on the report.

## Tips

- Pin everything: genome build, MANE Select transcript, annotation tool + version, predictor + version, gnomAD/ClinVar version. A classification is only reproducible relative to these.
- Never use flat 2015 defaults; verify the current ClinGen SVI recommendations and any VCEP spec.
- One predictor, calibrated strength -- SIFT/PolyPhen are decorative, CADD is for non-coding ranking, and stacking correlated tools over-calls pathogenic.
- Frequency is per-disease: grpmax FAF vs the maximum credible AF, not a global 1% line; presence in gnomAD is not benign (recessive, late-onset, clonal hematopoiesis).
- Never apply germline ACMG to a somatic variant; tumor tier is tumor-type-specific.
- A classification expires -- schedule reanalysis as new gnomAD/ClinVar/functional data arrive.

## Related Skills

- variant-calling/variant-annotation - Add VEP/SnpEff/ANNOVAR consequences and MANE Select transcripts feeding PVS1
- variant-calling/variant-normalization - Normalize/left-align variants before ClinVar and HGVS matching
- variant-calling/filtering-best-practices - Filter quality/artifact variants before clinical review
- variant-calling/vcf-basics - VCF field extraction and INFO parsing
- database-access/entrez-fetch - Programmatic ClinVar and OMIM download
