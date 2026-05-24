---
name: bio-copy-number-germline-cnv-interpretation
description: Classify constitutional (germline) copy number variants for clinical reporting using the 2019 ACMG/ClinGen technical standards points-based framework, with ClassifyCNV and AnnotSV for semi-automated scoring. Covers the separate copy-number-loss and copy-number-gain rubrics, the five-tier classification, ClinGen haploinsufficiency/triplosensitivity and dosage-sensitive regions, de novo and segregation evidence, and population-frequency benign evidence. Use when assigning pathogenic/likely-pathogenic/VUS/likely-benign/benign to a constitutional CNV, scoring a CNV against ACMG/ClinGen criteria, or distinguishing the automatable evidence from the case-specific evidence requiring manual input.
tool_type: mixed
primary_tool: ClassifyCNV
---

## Version Compatibility

Reference examples tested with: ClassifyCNV 1.1+, AnnotSV 3.4+, Python 3.10+ with
pandas 2.2+; bedtools 2.31+.

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `python ClassifyCNV.py --help`, `AnnotSV --version`
- Update the bundled ClinGen/dosage databases — ClassifyCNV ships an `update_clingen.sh`;
  dosage curation changes, and a stale database silently mis-scores.

This skill is for **constitutional/germline** CNVs only. Somatic tumor CNVs use a
different framework (AMP/ASCO/CAP and OncoKB tiers) — do not apply ACMG/ClinGen
constitutional scoring to a tumor.

# Germline CNV Interpretation

**"Is this constitutional CNV pathogenic"** -> Apply the 2019 ACMG/ClinGen technical standards: a semiquantitative, points-based rubric that sums evidence into one of five clinical categories. There are two separate rubrics — one for copy-number **loss**, one for copy-number **gain** — because the evidence for deletion and duplication pathogenicity is different. The total score maps to a five-tier classification.

- CLI: `ClassifyCNV` (automates the observed-evidence sections), `AnnotSV` (ACMG-aligned rank)
- Manual: case-specific evidence (de novo status, segregation, prior literature) is scored by the interpreter, not the tool

## The Points Framework

| Total score | Classification |
|-------------|----------------|
| >= 0.99 | Pathogenic |
| 0.90 to 0.98 | Likely pathogenic |
| -0.89 to 0.89 | Variant of uncertain significance (VUS) |
| -0.90 to -0.98 | Likely benign |
| <= -0.99 | Benign |

Evidence is grouped into sections (the loss and gain rubrics each have five). For
copy-number **loss**: Section 1 — does the CNV contain protein-coding or functionally
important elements; Section 2 — overlap with established haploinsufficient genes/regions
(strong positive) or established benign regions (strong negative); Section 3 — number of
protein-coding genes; Section 4 — detailed case/literature evidence (case-control, prior
probands, phenotype specificity); Section 5 — inheritance (de novo with confirmed
parentage is strong positive; inherited from an unaffected parent is negative). The
**gain** rubric is structured the same way but keyed to triplosensitivity and the
distinct evidence base for duplications.

The decisive postdoc-level point: **a tool can only score the evidence it is given.**
ClassifyCNV and AnnotSV automate Sections 1-3 (gene content, dosage-region overlap,
population frequency) well; Sections 4-5 (de novo status, segregation, literature)
require the interpreter to supply points. An unsupervised tool run therefore
systematically lands CNVs in VUS — the absence of family/literature evidence is not
neutral, it is unscored.

## Classification Workflow

| Step | Source | Automatable |
|------|--------|-------------|
| Gene content, functional elements | RefSeq/GENCODE | Yes (ClassifyCNV/AnnotSV) |
| Established HI/TS gene & region overlap | ClinGen dosage map | Yes |
| Protein-coding gene count | Gene model | Yes |
| Population frequency (benign evidence) | gnomAD-SV, DGV | Yes |
| Case-control / prior probands / phenotype fit | Literature, DECIPHER, internal DB | Partial — interpreter scores |
| De novo status, segregation | Trio/family data | No — interpreter scores |

## Semi-Automated Scoring with ClassifyCNV

**Goal:** Score the automatable ACMG/ClinGen sections for a set of constitutional CNVs.

**Approach:** Provide CNVs as a BED with an explicit DEL/DUP type; ClassifyCNV applies the
2019 rubric against the bundled ClinGen databases and emits a per-CNV scoresheet.

```bash
# Input BED: chrom, start, end, type  (type = DEL or DUP)
python ClassifyCNV.py \
    --infile constitutional_cnvs.bed \
    --GenomeBuild hg38 \
    --precise \
    --outdir classifycnv_out

# Output Scoresheet.txt: per-CNV total score, classification, and per-criterion points.
```

```python
import pandas as pd

def review_classifycnv(scoresheet):
    '''Flag CNVs whose ACMG class likely changes once case-specific evidence is added.'''
    df = pd.read_csv(scoresheet, sep='\t')
    # VUS CNVs near a tier boundary are the ones where de novo / segregation evidence
    # (Sections 4-5, not scored automatically) would tip the classification.
    df['near_boundary'] = df['Total score'].between(0.60, 0.89) | \
                          df['Total score'].between(-0.89, -0.60)
    df['needs_manual_evidence'] = (df['Classification'] == 'VUS') & df['near_boundary']
    return df
```

## Comprehensive Annotation Cross-Check with AnnotSV

```bash
AnnotSV -SVinputFile constitutional_cnvs.vcf -genomeBuild GRCh38 \
    -annotationMode both -outputFile annotsv_out.tsv
# AnnotSV emits an ACMG-aligned rank (1 benign - 5 pathogenic) per SV; use it to
# cross-check ClassifyCNV, not as a standalone clinical classification.
```

## Failure Modes

### Applying constitutional scoring to a somatic CNV

**Trigger:** Running ACMG/ClinGen germline classification on tumor copy number.

**Mechanism:** The 2019 standards are explicitly constitutional; somatic CNV clinical significance uses the AMP/ASCO/CAP tier system and oncology evidence (therapy, prognosis).

**Symptom:** Tumor amplifications classified as "pathogenic germline variants"; clinically meaningless report.

**Fix:** Confirm the CNV is constitutional (present in germline DNA). For tumors, use somatic oncology frameworks — see clinical-databases/variant-prioritization.

### Treating a tool's VUS as a final answer

**Trigger:** Reporting ClassifyCNV/AnnotSV output verbatim without adding case evidence.

**Mechanism:** Tools score gene content, dosage overlap, and frequency, but not de novo status, segregation, or literature; absent that input the score sits in the VUS band.

**Symptom:** Nearly every novel CNV classified VUS; clinically relevant de novo deletions under-called.

**Fix:** Treat tool output as the Section 1-3 baseline. Add Section 4-5 points from trio data, segregation, DECIPHER, and literature before issuing a classification. A VUS near a tier boundary specifically signals missing case evidence.

### Stale ClinGen dosage database

**Trigger:** Using ClassifyCNV/AnnotSV bundled databases without updating.

**Mechanism:** ClinGen dosage curation is ongoing; HI/TS scores and dosage-sensitive regions change. A stale database scores Section 2 wrong.

**Symptom:** A gene with a newly curated HI score 3 is scored as having no dosage evidence; classification too low.

**Fix:** Run the database update script before a classification batch; record the ClinGen release date in the report.

### Genome-build mismatch

**Trigger:** CNV coordinates and the `--GenomeBuild` argument (or annotation databases) on different builds.

**Mechanism:** Coordinates silently shift; the wrong genes and dosage regions are scored.

**Symptom:** Implausible gene content; a known disorder locus scored as gene-poor.

**Fix:** Confirm CNV coordinates, `--GenomeBuild`, and all databases are the same build; verify a landmark CNV.

### Partial-gene overlap scored as whole-gene loss

**Trigger:** Scoring a deletion that removes only part of a haploinsufficient gene as a full-gene loss.

**Mechanism:** The rubric distinguishes whole-gene loss from partial overlap; a deletion of a few exons may create a truncating allele with different (sometimes greater) impact, scored under different criteria.

**Symptom:** Partial-gene CNVs mis-scored; truncating deletions under- or over-weighted.

**Fix:** Record whether the CNV removes the whole gene or part of it, and which exons; apply the rubric's partial-overlap criteria explicitly.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| ClassifyCNV VUS, AnnotSV rank 4 | Different weighting of the same evidence | Re-derive points manually against the 2019 standard |
| Tool says benign, locus is a known disorder | Stale dosage database or build mismatch | Update databases; verify build |
| Two interpreters disagree on a VUS | Section 4-5 evidence weighted differently | Use the ClinGen calculator; document each criterion |
| De novo deletion still VUS | Section 5 points not added | Add confirmed-de-novo points |

**Operational rule:** A clinical CNV classification is final only when (1) the CNV is
confirmed constitutional, (2) databases and builds are current and consistent, (3) the
automatable Sections 1-3 are scored by a tool, and (4) the interpreter has scored
Sections 4-5 from case-specific evidence. Document each criterion and its points; the
ClinGen web calculator is the reference tally.

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| Pathogenic | total score >= 0.99 | Riggs 2020 ACMG/ClinGen technical standards |
| Likely pathogenic | 0.90 to 0.98 | Riggs 2020 |
| VUS | -0.89 to 0.89 | Riggs 2020 |
| Likely benign | -0.90 to -0.98 | Riggs 2020 |
| Benign | <= -0.99 | Riggs 2020 |
| Established dosage sensitivity | ClinGen HI/TS score = 3 | ClinGen: sufficient evidence |
| Common-CNV benign frequency | high population frequency | Section 2/4 benign evidence |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Tumor CNVs classified "pathogenic germline" | Constitutional rubric applied to somatic | Use somatic oncology frameworks |
| Almost everything classified VUS | Sections 4-5 not scored | Add de novo/segregation/literature points |
| Known disorder locus scored benign | Stale dosage DB or build mismatch | Update ClinGen databases; check build |
| Wrong genes scored | Build mismatch | Align coordinates, --GenomeBuild, databases |
| Partial-gene deletion mis-scored | Whole-gene assumption | Apply partial-overlap criteria |
| ClassifyCNV vs AnnotSV disagree | Different evidence weighting | Re-derive against the 2019 standard manually |

## References

- Riggs ER et al 2020. Technical standards for the interpretation and reporting of constitutional copy-number variants: a joint consensus recommendation of ACMG and ClinGen. Genet Med 22:245
- Gurbich TA, Ilinsky VV 2020. ClassifyCNV: a tool for clinical annotation of copy-number variants. Sci Rep 10:20375
- Geoffroy V et al 2018. AnnotSV: an integrated tool for structural variations annotation. Bioinformatics 34:3572
- Rehm HL et al 2015 NEJM 372:2235 (ClinGen launch / framework). Dosage-sensitivity curation methodology is in Riggs ER et al 2012 Genet Med 14:680 (original ClinGen dosage-sensitivity workflow). Current ClinGen Dosage Sensitivity Map: clinicalgenome.org.

## Related Skills

- copy-number/cnv-annotation - Gene, dosage, and database annotation feeding the rubric
- copy-number/gatk-cnv - GATK-gCNV germline CNV calling
- copy-number/cnvkit-analysis - Germline CNV calling from panels/exomes
- clinical-databases/clinvar-lookup - ClinVar CNV records and prior classifications
- clinical-databases/variant-prioritization - Somatic variant tiering (the non-germline path)
- clinical-databases/gnomad-frequencies - Population frequency for benign evidence
