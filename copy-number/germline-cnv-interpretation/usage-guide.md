# Germline CNV Interpretation Usage Guide

## Overview

Classifying a constitutional copy number variant for clinical reporting means applying the 2019 ACMG/ClinGen technical standards: a semiquantitative, points-based rubric that sums evidence into one of five categories (pathogenic, likely pathogenic, VUS, likely benign, benign). There are separate rubrics for copy-number loss and copy-number gain. Tools (ClassifyCNV, AnnotSV) automate the gene-content, dosage-overlap, and population-frequency evidence, but de novo status, segregation, and literature evidence must be scored by the interpreter, which is why unsupervised tool runs systematically return VUS. This skill is for constitutional CNVs only; somatic tumor CNVs use a different framework.

## Prerequisites

```bash
git clone https://github.com/Genotek/ClassifyCNV    # then run its update_clingen.sh
conda install -c bioconda annotsv bedtools
pip install pandas
```

Inputs: constitutional CNVs as a BED (chrom, start, end, type = DEL/DUP) or VCF; current ClinGen dosage databases (update before each batch); for full classification, trio/family data and access to literature/DECIPHER.

## Quick Start

Tell the AI agent what to do:
- "Classify these constitutional CNVs with the ACMG/ClinGen points framework"
- "Score this deletion against the 2019 ACMG/ClinGen copy-number-loss rubric"
- "Identify which of my VUS CNVs would change class if de novo evidence were added"
- "Cross-check ClassifyCNV and AnnotSV classifications and reconcile disagreements"
- "Explain why this tumor CNV should not be scored with the constitutional rubric"

## Example Prompts

### Classification

> "Run ClassifyCNV on this constitutional CNV BED in hg38 and produce a scoresheet, then tell me which sections were scored automatically and which need manual evidence."

> "Score this de novo 16p11.2 deletion against the ACMG/ClinGen loss rubric, including the Section 5 de novo points."

### Triage and review

> "I have 200 CNVs classified VUS by ClassifyCNV. Identify the ones near a tier boundary where segregation or de novo evidence would tip the classification."

> "Cross-check ClassifyCNV against AnnotSV ranks for my callset and reconcile the disagreements against the 2019 standard."

### Scope

> "These CNVs came from a tumor. Explain why the constitutional ACMG/ClinGen framework does not apply and what the correct somatic framework is."

## What the Agent Will Do

1. Confirm the CNVs are constitutional (not somatic)
2. Verify genome build and update ClinGen dosage databases
3. Score the automatable sections (gene content, dosage overlap, frequency) with ClassifyCNV
4. Cross-check with AnnotSV's ACMG-aligned rank
5. Flag VUS CNVs near tier boundaries that need case-specific evidence
6. Add de novo, segregation, and literature points from family/clinical data
7. Produce a documented per-criterion classification

## Tips

- The framework is constitutional-only; never apply it to somatic tumor CNVs.
- Tools score gene content, dosage overlap, and frequency, not de novo status, segregation, or literature; unsupervised runs therefore default to VUS.
- A VUS near a tier boundary is a signal that case-specific evidence is missing, not a final answer.
- Update the ClinGen dosage databases before each batch; record the release date.
- Keep CNV coordinates, the `--GenomeBuild` argument, and all databases on one build.
- Distinguish whole-gene loss from partial-gene overlap; they score under different criteria and a truncating partial deletion can be more damaging.
- Document each scored criterion; the ClinGen web calculator is the reference tally.

## Related Skills

- copy-number/cnv-annotation - Gene, dosage, and database annotation feeding the rubric
- copy-number/gatk-cnv - GATK-gCNV germline CNV calling
- copy-number/cnvkit-analysis - Germline CNV calling from panels and exomes
- clinical-databases/clinvar-lookup - ClinVar CNV records and prior classifications
- clinical-databases/variant-prioritization - Somatic variant tiering
- clinical-databases/gnomad-frequencies - Population frequency for benign evidence
