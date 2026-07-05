---
name: bio-variant-calling-clinical-interpretation
description: Classify variant clinical significance with the ACMG/AMP germline framework and its 2018-2025 ClinGen refinements (graded PVS1 decision tree, PM2 downgraded to Supporting, PP5/BP6 retired, calibrated PP3/BP4, Bayesian points), the AMP/ASCO/CAP somatic tiers and ClinGen oncogenicity system, ClinVar star-rating and gnomAD grpmax filtering-AF interpretation. Use when deciding germline-vs-somatic framework, applying current (not flat-2015) ACMG points, checking for a gene-specific VCEP specification, judging whether a ClinVar assertion or gnomAD frequency is usable evidence, calibrating a pathogenicity predictor, evaluating PVS1 on the MANE Select transcript, or building a VUS reanalysis loop. Not for functional annotation itself (see variant-calling/variant-annotation).
tool_type: mixed
primary_tool: bcftools
---

## Version Compatibility

Reference examples tested with: bcftools 1.19+, cyvcf2 0.30+, InterVar 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: interpretation guidance evolves. Flat 2015 ACMG defaults are OUT OF DATE; verify the current ClinGen SVI recommendations and any gene-specific VCEP specification before classifying. Never apply germline ACMG to a somatic variant.

# Clinical Variant Interpretation

**"Classify this variant / write the ACMG rationale"** -> Assemble independent, calibrated evidence lines and combine them under the correct framework for a pinned (gene, transcript, disease) context.
- Germline Mendelian: ACMG/AMP + ClinGen SVI refinements (below).
- Somatic/tumor: AMP/ASCO/CAP tiers + ClinGen/CGC/VICC oncogenicity (never ACMG).

## The governing principle

A variant's clinical significance is NOT a database lookup. It is a Bayesian sum of INDEPENDENT, CALIBRATED evidence relative to a pinned context (genome build, MANE Select transcript, gene disease-mechanism, disease prevalence, framework version). Three traps sink most naive pipelines:

1. **Flat 2015 defaults are obsolete.** A current classifier applies the ClinGen SVI refinements (graded PVS1, PM2 downgraded, PP5/BP6 retired, calibrated PP3/BP4, Bayesian points). Using the raw 2015 combining rules is a known error.
2. **Germline and somatic are different questions with different frameworks.** "Does this cause a Mendelian disorder" (ACMG) vs "is this an actionable/oncogenic tumor variant" (Li tiers, Horak oncogenicity). Applying ACMG to a somatic variant is a category error.
3. **A ClinVar assertion is a LEAD, not evidence.** Concordance between submitters is not independence; 1-star is not usable; re-derive from the underlying data.

## Pick the framework FIRST

| Variant origin | Framework | Question answered | Cite |
|----------------|-----------|-------------------|------|
| Germline (constitutional) | ACMG/AMP + ClinGen SVI | Pathogenic..Benign for a Mendelian disorder | Richards 2015; Abou Tayoun 2018; Tavtigian 2020; Pejaver 2022 |
| Somatic (tumor), actionability | AMP/ASCO/CAP tiers I-IV | Diagnostic/prognostic/therapeutic significance in THIS tumor type | Li 2017 |
| Somatic, oncogenicity | ClinGen/CGC/VICC points | Oncogenic..Benign (is it a driver) | Horak 2022 |

Before applying generic ACMG, **check for a ClinGen Variant Curation Expert Panel (VCEP) specification for the gene** (e.g. hearing loss, RASopathy, cardiomyopathy, ENIGMA BRCA1/2). A VCEP spec reweights and constrains criteria and OVERRIDES generic defaults; a 3-star ClinVar assertion often reflects one.

## ACMG/AMP germline: the 2015 baseline and its mandatory refinements

The 2015 consensus (Richards 2015 *Genet Med* 17:405-424) defines five tiers (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign) and 28 coded criteria at default strengths: PVS1 (very strong), PS1-4 (strong), PM1-6 (moderate), PP1-5 (supporting); BA1 (stand-alone), BS1-4 (strong), BP1-7 (supporting). A director does NOT interpret with raw 2015 anymore. Apply these ClinGen SVI corrections:

| Refinement | What changed | Consequence for the classifier |
|------------|--------------|-------------------------------|
| Graded PVS1 (Abou Tayoun 2018 *Hum Mutat* 39:1517) | PVS1 is a decision tree, not automatic for any null | Emit PVS1 at Very Strong / Strong / Moderate / Supporting per NMD + mechanism (below) |
| PM2 -> Supporting (ClinGen SVI PM2 v1.0, approved Sept 2020) | Absence from gnomAD is WEAK | Apply PM2 at Supporting, never Moderate |
| PP5 / BP6 RETIRED (Biesecker & Harrison 2018 *Genet Med* 20:1687) | An assertion cannot substitute for evidence | Never use PP5/BP6; cite the underlying data instead |
| Calibrated PP3 / BP4 (Pejaver 2022 *AJHG* 109:2163) | Computational evidence is graded, not flat-Supporting | Use ONE calibrated predictor at its calibrated strength (below) |
| Bayesian points (Tavtigian 2018/2020) | Verbal combining rules approximate naive Bayes | Sum points; graded/fractional strengths are coherent |

### Bayesian points system (Tavtigian 2020 *Hum Mutat* 41:1734)

**Goal:** Combine graded evidence into a tier reproducibly instead of matching verbal rule patterns.

**Approach:** Assign each met criterion a point value by strength (benign subtracts), sum, and threshold. This underlies the emerging points-based ACMG/AMP/CAP/ClinGen overhaul, so prefer it over the 2015 verbal table.

| Strength | Points (P side) | OddsPath (Tavtigian 2018, prior ~0.10) |
|----------|-----------------|-----------------------------------------|
| Supporting | +1 | ~2.08 |
| Moderate | +2 | ~4.33 |
| Strong | +4 | ~18.7 |
| Very Strong | +8 | ~350 |

Classification by summed points: Pathogenic >= 10, Likely Pathogenic 6-9, VUS 0-5, Likely Benign -1 to -6, Benign <= -7 (confirm the exact benign cutpoints against Tavtigian 2020 before hard-coding). Benign criteria (BA1/BS/BP) contribute negative points at the same magnitudes.

### PVS1 decision tree and the NMD 50-nt rule

**Goal:** Assign PVS1 the CORRECT strength for a null variant instead of firing it on any "HIGH impact" call.

**Approach:** Route by gene LOF mechanism, then variant type, then NMD prediction and exon location (Abou Tayoun 2018). Evaluate on the MANE Select transcript, not whichever isoform maximizes severity.

- **Gene mechanism gate:** PVS1 applies ONLY where loss of function is the established disease mechanism (haploinsufficiency). For gain-of-function / dominant-negative genes a null may be benign -- PVS1 must not fire.
- **NMD 50-55 nt rule:** a premature termination codon >~50-55 nt upstream of the last exon-exon junction triggers nonsense-mediated decay (true LOF -> full strength). A PTC in the LAST exon, within ~50 nt of the final junction, or in a single-exon gene ESCAPES NMD -- protein is made; downgrade PVS1 (Strong/Moderate/Supporting) by how much functional protein / which domains are lost.
- **Transcript relevance:** confirm the affected exon is in biologically expressed transcripts; a canonical-splice change in a minor non-expressed isoform is not PVS1.
- "HIGH impact stop_gained" from SnpEff/ANNOVAR is NOT PVS1 -- impact buckets know nothing about NMD or mechanism. Evaluate PVS1 on MANE Select; do not use the worst-consequence transcript. See variant-calling/variant-annotation.

## ClinVar: assertions are leads, not evidence

**"Look up this variant in ClinVar"** -> Read WHO submitted, at what review status, on WHAT evidence -- then re-derive, do not adopt the conclusion.

| CLNREVSTAT | Stars | Usable as evidence? |
|------------|-------|---------------------|
| practice_guideline | 4 | Strongest single-DB signal; still verify vs current evidence |
| reviewed_by_expert_panel | 3 | VCEP; strong, often implies a gene specification |
| criteria_provided,_multiple_submitters,_no_conflicts | 2 | Consensus; check submitters shared no common error |
| criteria_provided,_single_submitter | 1 | A LEAD only -- not usable as evidence |
| criteria_provided,_conflicting_classifications | 1 | Conflict is an informative signal, not noise to average |
| no_assertion_criteria_provided | 0 | No weight |

Rules: 1-star / no-criteria is not evidence. Conflicting interpretations flag genuinely hard variants (penetrance, ancestry, mechanism) -- investigate, do not average. Concordance is not independence (two submitters can copy one original error). PP5/BP6 are retired precisely because an assertion cannot be an evidence input.

### Annotate and read ClinVar fields (bcftools / cyvcf2)

**Goal:** Attach ClinVar assertions as LEADS and surface review status alongside significance.

**Approach:** Annotate CLNSIG/CLNDN/CLNREVSTAT from the ClinVar VCF, then always carry CLNREVSTAT so a 1-star call is never mistaken for evidence. Download the build-matched ClinVar VCF first (usage-guide.md).

```bash
bcftools annotate -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT input.vcf.gz -Oz -o with_clinvar.vcf.gz

# Surface P/LP leads WITH their review status (never drop CLNREVSTAT)
bcftools view -i 'INFO/CLNSIG~"athogenic"' with_clinvar.vcf.gz \
  | bcftools query -f '%CHROM:%POS %REF>%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\n'
```

## Population frequency: grpmax filtering-AF, not a global cutoff

**Goal:** Decide BA1/BS1 (or PM2_Supporting) correctly for THIS disease, not with a universal 1% line.

**Approach:** Compare the gnomAD grpmax filtering allele frequency to the maximum credible population AF derived from disease prevalence, heterogeneity, inheritance and penetrance (Whiffin 2017 *Genet Med* 19:1151). A flat cutoff is wrong in both directions.

- **Filtering AF (FAF)** is the LOWER bound of the 95% CI of the grpmax (genetic-ancestry-group max) AF -- gnomAD v4 exposes it as the `fafmax_faf95_max` INFO field (`fafmax_faf95_max_joint` in the joint exome+genome VCF). Using grpmax, not global AF, avoids diluting a variant common in one ancestry across the whole cohort; using the CI lower bound guards against a noisy small-subpopulation estimate.
- **Rule:** if FAF > the disease's maximum credible population AF, apply BA1/BS1. This is per-disease.
- **Presence in gnomAD is NOT benign.** Exceptions a director watches for: recessive carriers are healthy (pathogenic alleles sit at carrier frequency, e.g. CFTR); late-onset / reduced-penetrance alleles appear in adult cohorts (BRCA, Lynch); somatic / clonal-hematopoiesis contamination leaks low-AF calls in DNMT3A/TET2; artifacts in homopolymer/segdup regions -- respect gnomAD PASS/quality flags, not raw AF.
- gnomAD ancestry groups are unevenly sampled, so "absent" is much weaker evidence for an under-represented ancestry; PM2/BS1 strength is implicitly ancestry-dependent. gnomAD v2.1.1 is GRCh37 (Karczewski 2020 *Nature* 581:434); v3/v4 are GRCh38 -- never eyeball "absent" across builds without liftover.

```bash
# Illustrative: filter on a grpmax filtering-AF field, keeping absent sites (annotation-dependent)
bcftools view -i 'INFO/fafmax_faf95_max<0.0001 || INFO/fafmax_faf95_max="."' \
    input.vcf.gz -Oz -o faf_filtered.vcf.gz
```

## Pathogenicity predictors: ONE, calibrated

**Goal:** Convert a computational score into PP3/BP4 at a defensible strength without double-counting.

**Approach:** Pick ONE predictor that reached >= Strong in the ClinGen calibration and apply it at its calibrated threshold (Pejaver 2022). Stacking correlated tools fakes independence and silently over-calls pathogenic.

- PP3 and BP4 are graded (Supporting/Moderate/Strong) and mutually exclusive. For REVEL (Ioannidis 2016 *AJHG* 99:877) the well-reproduced SUPPORTING thresholds are PP3 >= 0.644 and BP4 <= 0.290; higher-strength (Moderate/Strong) cutoffs exist -- read them from the Pejaver 2022 supplement or the current ClinGen SVI table rather than hard-coding.
- Use only ONE tool. REVEL is an ensemble of 13 scores (incl. SIFT, PolyPhen), so "REVEL agrees with PolyPhen" is not corroboration -- PolyPhen is INSIDE REVEL.
- **SIFT and PolyPhen-2 did not reach even Supporting** in the calibration -- a "damaging" call is decorative, not evidence. **Raw CADD did not reach Supporting for PP3** (CADD>=20 calibrated to benign-Moderate -- mild evidence AGAINST missense pathogenicity, the opposite of how it is usually invoked); CADD is for genome-wide/non-coding ranking, not missense PP3.
- **AlphaMissense** (Cheng 2023 *Science* 381:eadg7492) is proteome-wide and not trained on ClinVar labels, but its developer class cutoffs are NOT ACMG strengths -- check the current ClinGen SVI tool list for its calibrated PP3/BP4 thresholds before assigning a strength.
- **Splicing (SpliceAI, Jaganathan 2019 *Cell* 176:535):** delta scores 0-1, developer guidance 0.2 recall / 0.5 recommended / 0.8 precision. A high delta is a PREDICTION; converting it to PS3/PP3 strength needs the ClinGen splicing calibration, and the default scoring window is narrow -- widen it (deep-intronic/pseudoexon variants are otherwise missed). SpliceAI does not report the mis-splicing OUTCOME (exon skip vs intron retention), which determines PVS1 applicability.

### Python: research-triage prioritization (NOT formal ACMG)

**Goal:** Rank candidate variants for review triage using available annotations.

**Approach:** Combine ClinVar leads, grpmax frequency and a single calibrated predictor into a tier. This is a triage helper, not an ACMG classification -- computational scores are supporting only, and stacking here is for RANKING, not evidence.

```python
from cyvcf2 import VCF

def triage_tier(variant):
    # Triage ranking ONLY; not equivalent to ACMG. ClinVar is a lead (carry review status
    # separately), scores are PP3/BP4-supporting, and stacking predictors here just ranks.
    clnsig = str(variant.INFO.get('CLNSIG', ''))
    faf = variant.INFO.get('fafmax_faf95_max', 0) or 0
    revel = variant.INFO.get('REVEL', 0) or 0  # single calibrated predictor

    if 'Pathogenic' in clnsig and 'Likely' not in clnsig:
        return 'PATHOGENIC_LEAD'
    if 'Likely_pathogenic' in clnsig:
        return 'LIKELY_PATHOGENIC_LEAD'
    if 'Benign' in clnsig or faf > 0.05:  # BA1 territory; confirm vs disease-max credible AF
        return 'BENIGN_LEAD'
    if revel >= 0.644 and faf < 0.0001:    # REVEL PP3_Supporting threshold (Pejaver 2022)
        return 'VUS_FAVOR_PATH'
    if revel <= 0.290:                     # REVEL BP4_Supporting threshold
        return 'VUS_FAVOR_BENIGN'
    return 'VUS'

vcf = VCF('annotated.vcf.gz')
report = {'PATHOGENIC_LEAD', 'LIKELY_PATHOGENIC_LEAD', 'VUS_FAVOR_PATH'}
for v in vcf:
    tier = triage_tier(v)
    if tier in report:
        gene = v.INFO.get('SYMBOL', 'NA')
        print(f'{gene}\t{v.CHROM}:{v.POS}\t{tier}\t{v.INFO.get("CLNREVSTAT", ".")}')
```

## Somatic variants: a separate framework

**"Interpret this tumor variant"** -> Ask about actionability and oncogenicity in THIS tumor type, never germline pathogenicity. Tier is tumor-type-specific (BRAF V600E is Tier I in melanoma, lower elsewhere) -- a context-dependence with no germline analog.

**AMP/ASCO/CAP tiers (Li 2017 *J Mol Diagn* 19:4)** -- clinical actionability:
- Tier I: strong significance (FDA-approved therapy for this variant + tumor type, or in guidelines).
- Tier II: potential significance (therapy in another tumor type; trial evidence; multiple studies).
- Tier III: unknown clinical significance (the somatic "VUS").
- Tier IV: benign/likely benign (common, no oncogenic role).

**ClinGen/CGC/VICC oncogenicity (Horak 2022 *Genet Med* 24:986)** -- a SEPARATE points-based axis (Oncogenic..Benign) using cancer-specific codes (hotspot recurrence, functional oncogenic data, tumor frequency). Oncogenicity != actionability: an oncogenic driver may have no drug (Tier III despite oncogenic).

Knowledgebase evidence levels: OncoKB Level 1-4 + R1/R2 (therapeutic), CIViC evidence A-E (read the evidence item, not just the letter), COSMIC recurrence (a hotspot SIGNAL, not clinical actionability). **Tumor-only** assays cannot cleanly separate somatic from germline -- a ~50%/~100% VAF variant may be germline; filter and disclose explicitly, or use paired tumor-normal.

## Classification has an expiry date

A classification is a snapshot relative to the evidence available on its date. Build a reanalysis loop: periodically re-annotate stored VCFs against the latest ClinVar and gnomAD releases and flag VUS whose evidence changed (new functional/segregation data, a new VCEP spec, a frequency that now crosses BA1/BS1). A one-time classification without reanalysis is a latent error.

**Goal:** Re-score stored VUS against a newer ClinVar release and surface those whose assertion has since become definitive.

**Approach:** Re-annotate the prior results with the current ClinVar under a distinct INFO tag, then select records that were Uncertain but now carry a pathogenic/benign assertion.

```bash
# Re-annotate against a newer ClinVar; find VUS that now carry a definitive assertion
bcftools annotate -a clinvar_latest.vcf.gz -c INFO/CLNSIG_NEW:=INFO/CLNSIG \
    prior_results.vcf.gz -Oz -o reannotated.vcf.gz
bcftools view -i 'INFO/CLNSIG~"Uncertain" && (INFO/CLNSIG_NEW~"athogenic" || INFO/CLNSIG_NEW~"enign")' \
    reannotated.vcf.gz -Oz -o reclassified.vcf.gz
```

## Common Errors

| Symptom / mistake | Cause | Fix |
|-------------------|-------|-----|
| PVS1 fired on any stop_gained | Used SnpEff HIGH-impact bucket | Route through the Abou Tayoun tree: mechanism + NMD + MANE transcript |
| PM2 applied at Moderate | Flat 2015 default | PM2_Supporting (ClinGen SVI 2020) |
| Over-called pathogenic | Stacked SIFT+PolyPhen+REVEL | One calibrated predictor at its calibrated strength; the others are inside REVEL |
| Adopted a 1-star ClinVar "Pathogenic" | Treated an assertion as evidence | 1-star is a lead; re-derive; carry CLNREVSTAT |
| Benign called on global AF > 1% | Ignored grpmax + disease context | grpmax FAF vs disease max-credible AF (Whiffin) |
| Common founder allele benignized | Global AF hid an ancestry-specific frequency | Use grpmax; presence in gnomAD != benign |
| ACMG applied to a tumor variant | Wrong framework | Li 2017 tiers + Horak 2022 oncogenicity |
| "Absent in gnomAD" across versions | v2 is GRCh37, v3/v4 GRCh38 | Liftover the variant; check site callability |

## Related Skills

- variant-calling/variant-annotation - VEP/SnpEff/ANNOVAR consequence calls, MANE Select transcripts, tool concordance feeding PVS1
- variant-calling/variant-normalization - left-align/normalize before ClinVar/HGVS matching
- variant-calling/filtering-best-practices - quality/artifact filtering before clinical review
- variant-calling/vcf-basics - VCF field extraction and INFO parsing
- database-access/entrez-fetch - programmatic ClinVar/OMIM download

## References

- Richards S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the ACMG and the AMP. *Genetics in Medicine*. 2015;17(5):405-424.
- Abou Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Human Mutation*. 2018;39(11):1517-1524.
- Tavtigian SV, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genetics in Medicine*. 2018;20(9):1054-1060.
- Tavtigian SV, et al. Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Human Mutation*. 2020;41(10):1734-1737.
- Biesecker LG, Harrison SM. The ACMG/AMP reputable source criteria for the interpretation of sequence variants. *Genetics in Medicine*. 2018;20(12):1687-1688.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *American Journal of Human Genetics*. 2022;109(12):2163-2177.
- Whiffin N, et al. Using high-resolution variant frequencies to empower clinical genome interpretation. *Genetics in Medicine*. 2017;19(10):1151-1158.
- Karczewski KJ, et al. The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*. 2020;581(7809):434-443.
- Ioannidis NM, et al. REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *American Journal of Human Genetics*. 2016;99(4):877-885.
- Cheng J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*. 2023;381(6664):eadg7492.
- Jaganathan K, et al. Predicting splicing from primary sequence with deep learning. *Cell*. 2019;176(3):535-548.
- Morales J, et al. A joint NCBI and EMBL-EBI transcript set for clinical genomics and research (MANE). *Nature*. 2022;604:310-315.
- Li MM, et al. Standards and guidelines for the interpretation and reporting of sequence variants in cancer: a joint consensus recommendation of AMP, ASCO, and CAP. *Journal of Molecular Diagnostics*. 2017;19(1):4-23.
- Horak P, et al. Standards for the classification of pathogenicity of somatic variants in cancer (oncogenicity): joint recommendations of ClinGen, CGC, and VICC. *Genetics in Medicine*. 2022;24(5):986-998.
- Landrum MJ, et al. ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Research*. 2018;46(D1):D1062-D1067.
