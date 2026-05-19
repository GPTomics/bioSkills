---
name: bio-clinical-databases-tumor-mutational-burden
description: Calculates tumor mutational burden from WES/WGS/panel data with Friends of Cancer Research harmonization equations, per-assay calibration (FDA 10/Mb = 7.8 TSO500 = 8.4 OncomineTML), synonymous/indel/germline filtering, hypermutator tiering, blood TMB, and integration with HLA-LOH and neoantigen quality (Luksza 2017 fitness). Use when assessing ICI eligibility under tumor-specific cutoffs (McGrail 2021), comparing tissue vs bTMB, or auditing TMB-H reporting against ESMO 2024 and FDA pembrolizumab pan-tumor 2020.
tool_type: python
primary_tool: cyvcf2
---

## Version Compatibility

Reference examples tested with: cyvcf2 0.30+, VEP 111+ (or snpEff 5.2+), pandas 2.2+, numpy 1.26+, LOHHLA 1.0+ (Marty 2017), DASH 1.0+ (Montesion 2021). v4.1 (May 2024) gnomAD is current for germline subtraction. Friends of Cancer Research TMB harmonization framework (Vega 2021 *Ann Oncol*) and ESMO 2024 (Mosele *Ann Oncol*) define the operational thresholds.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying. TMB calculation requires VCF with VEP / snpEff / Funcotator consequence annotations; the panel size used as denominator MUST match the assay's actual scored region, NOT the panel's total content.

# Tumor Mutational Burden; Calculation, Harmonization, ICI Eligibility

**'Calculate TMB from this somatic VCF and apply ICI eligibility cutoff'** -> Count nonsynonymous coding variants passing VAF/depth/germline filters; divide by assay scored region in Mb; apply assay-calibrated TMB-H cutoff; integrate with MSI / HLA-LOH / neoantigen quality.

- Python: `cyvcf2.VCF()` + VEP/snpEff consequence parsing + panel-size normalization
- CLI: `bcftools view` filtering + custom counting
- HLA-LOH: LOHHLA (Marty 2017 *Cell*) or DASH (Montesion 2021 *Cancer Discov*)
- Neoantigen quality: pVAC-tools, NetMHCpan-4.1, Luksza 2017 fitness model

## Regulatory and Trial Landscape

| Event | Year | Threshold | Notes |
|-------|------|-----------|-------|
| **KEYNOTE-158 + FDA pembrolizumab pan-tumor approval** | 2020 | TMB-H >= 10 mut/Mb | FoundationOne CDx companion diagnostic; 10 cohorts |
| **Friends of Cancer Research TMB harmonization Phase I (Merino 2020)** | 2020 | -- | 11 panels vs WES truth; 3-fold panel-specific differences |
| **Friends of Cancer Research Phase II (Vega 2021)** | 2021 | Calibration equations | 19 platforms; per-assay calibration to WES-aligned TMB-Mb |
| **ESMO 2024 (Mosele *Ann Oncol*)** | 2024 | TMB-H >= 10/Mb retained | NOT endorsed for breast, prostate, glioma |
| **KEYNOTE-189 (NSCLC + pembrolizumab + chemo)** | 2018 | -- | TMB-H did NOT enrich for benefit with chemo backbone |
| **POSEIDON / KEYNOTE-021 / KEYNOTE-407** | 2019-2022 | -- | TMB inconsistent with chemo backbones |
| **B-F1RST + BFAST Cohort C (bTMB)** | 2022 | bTMB >= 16/Mb | BFAST Cohort C FAILED primary endpoint |

## Friends of Cancer Research Harmonization: Cross-Panel Calibration

Merino 2020 *J Immunother Cancer*: in silico panel sampling from TCGA WES truth showed panel-specific TMB can differ 3-fold for identical samples. Vega 2021 *Ann Oncol* derived per-panel calibration equations to translate panel TMB to WES-aligned TMB-Mb.

**Per-panel calibration to FoundationOne 10/Mb sensitivity:**

| Panel | Scored region (Mb) | Equivalent threshold for FDA 10/Mb pan-tumor | Fails when |
|-------|---------------------|----------------------------------------------|-----------|
| **FoundationOne CDx** | 0.8 Mb scored (NOT 1.1 Mb total) | 10 mut/Mb (FDA reference; F1CDx companion) | Using 1.1 Mb panel total inflates TMB ~37%; pipeline excludes synonymous (F1CDx includes them) |
| **MSK-IMPACT v3** | 0.98 Mb | ~10 (full Vega 2021 calibration recommended) | Tumor purity < 30%; non-paired-normal mode |
| **MSK-IMPACT v4** | 1.22 Mb | ~10 | -- |
| **TruSight Oncology 500** | ~1.3 Mb scored (from 1.94 Mb total) | **7.8 mut/Mb** | Pipeline uses 10/Mb instead of Vega 2021 calibrated 7.8 |
| **Oncomine Tumor Mutation Load** | 1.2 Mb | **8.4 mut/Mb** | Pipeline uses 10/Mb instead of Vega 2021 calibrated 8.4 |
| **Caris MI Tumor Seek** | ~1.2 Mb |; (verify Caris docs) | -- |
| **Tempus xT v3** | 0.6 Mb | -- | Below 0.8 Mb minimum reliability threshold |
| **Predicine ATLAS** | ~0.6 Mb | -- | Below 0.8 Mb minimum; high sampling variance |

**TMB =/= TMB across vendors.** Manuscripts that compare TMB across panels without per-assay calibration are unreviewable. Use the Vega 2021 calibration equations or WES re-projection.

## Variant-Counting Subtleties

These choices alter TMB by 5-20%:

| Variable | Convention | Notes |
|----------|-----------|-------|
| **Synonymous variants** | **FoundationOne CDx INCLUDES synonymous** (rationale: reduces sampling noise); MSK-IMPACT and most academic pipelines exclude | The FDA companion diagnostic counts synonymous; frequent misconception |
| **Indels** | FoundationOne includes; some assays exclude frameshift only | 5-15% TMB impact |
| **Germline subtraction** | Paired-normal (gold standard); else gnomAD AF <=0.5% (sometimes 1%) for tumor-only | Population-stratified gnomAD AF for ancestry-diverse cohorts |
| **VAF threshold** | FoundationOne >=5%; >=10% for tumor-only no UMI; down to 2% with paired-normal | Lower VAF risks contamination/artifacts |
| **Hotspots** | COSMIC-confirmed driver hotspots typically EXCLUDED (not random) | Inflates TMB if included |
| **Tumor purity** | FoundationOne >=20%; MSK-IMPACT >=30% | Below floor erodes VAF-based filtering |
| **VEP version** | Pin to assay's annotation version | gnomAD v4 uses VEP 105 |

## Hypermutator Tiering

| Class | Threshold | Common etiology |
|-------|-----------|----------------|
| **TMB-H (FDA pan-cancer)** | >= 10 mut/Mb | Variable; ICI eligible |
| **Hypermutator (research)** | >= 100 mut/Mb | MMR-D, POLE-exo |
| **Ultra-hypermutator** | >= 500 mut/Mb | POLE+MMR concurrent |

MSI-H typically 30-50 mut/Mb; pure POLE-exo P286R 100-300 mut/Mb; POLE-exo + MMR-D exceeds 500. MSI-H and TMB-H overlap substantially in CRC and endometrial (~80% of MSI-H are TMB-H) but only ~16% of TMB-H solid tumors are MSI-H (Salem 2018 *Cancer Discov*).

## The Tumor-Type-Specific Cutoff Debate

**McGrail 2021** *Ann Oncol* is the most damning paper for the universal 10/Mb cutoff. TMB-H predicts ICI response in melanoma, NSCLC, bladder; but FAILS in breast, prostate, glioma. ORR in TMB-H melanoma/NSCLC/bladder was 39.8%; TMB-H breast/prostate/glioma was 15.3%. Mechanistic explanation: TMB only predicts when baseline CD8 T-cell infiltrate is present.

**Sha 2020** *Cell Rep Med*: TMB-H predicts ICI benefit in MSS subset but adds nothing on top of MSI-H (because MSI-H is uniformly hypermutator and uniformly responsive).

**Samstein 2019** *Nat Genet* (MSK-IMPACT 1,662 ICI-treated): cancer-specific TMB cutoffs (top 20% within each tumor type) outperform universal 10/Mb.

**ESMO 2024** retained TMB-H >= 10/Mb pan-tumor but explicitly noted exceptions: **NOT endorsed for breast, prostate, glioma** based on negative real-world data.

## Blood TMB (bTMB): The Negative-Trial Story

**Gandara 2018** *Nat Med*: bTMB on Foundation Medicine FoundationACT panel; POPLAR + OAK retrospective. bTMB >= 16 mut/Mb showed PFS benefit with atezolizumab in NSCLC.

**B-F1RST (Kim 2022)**: bTMB >= 16 prospectively predictive for atezolizumab first-line NSCLC.

**BFAST Cohort C (Dziadziuszko 2022)**: FAILED primary endpoint; atezolizumab vs chemo in bTMB-H NSCLC did not improve investigator-assessed PFS. Dominant confounder: low ctDNA shed fraction produces false-negative bTMB.

**Operational state:** bTMB is research-grade in tissue-naive settings; tissue TMB remains the regulatory standard.

## Neoantigen Quality: Beyond Raw TMB

**Luksza 2017** *Nature*: neoantigen fitness model. Combines "non-selfness" (TCR recognition probability via IEDB similarity) + "selfness" (MHC binding affinity differential vs WT peptide). Pancreatic-cancer validation (Balachandran 2017 *Nature*): long-term survivors had higher-quality neoantigens. Luksza 2022 *Nature*: immunoediting over 10 years.

**McGranahan 2016** *Science*: **clonal neoantigen burden** (mutations present in all tumor cells) predicts ICI response better than total. Subclonal-rich tumors evade despite high TMB.

**HLA-LOH** (Marty 2017 *Cell*, LOHHLA; Montesion 2021 *Cancer Discov*, DASH): HLA-LOH occurs in ~40% of NSCLC and abolishes neoantigen presentation for the lost allele. ~17% pan-cancer; >30% in HNSCC / NSCLC / cervical. Co-occurs with high subclonal burden + APOBEC + immune escape.

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Pan-tumor ICI eligibility (FDA pembrolizumab) | TMB-H >= 10/Mb on FoundationOne CDx | FDA companion diagnostic |
| Non-FoundationOne panel | Apply Vega 2021 calibration to FoundationOne 10/Mb equivalent | TSO500 = 7.8; Oncomine = 8.4 |
| WES TMB | Compute directly; threshold per ESMO 2024 = 10/Mb | WES is reference standard |
| Tissue-naive bTMB | Caution: BFAST Cohort C failed | Research-grade; check ctDNA shed fraction |
| Breast / prostate / glioma | TMB-H NOT endorsed per ESMO 2024 + McGrail 2021 | Tumor-type-specific cutoffs |
| MSI-H + TMB-H concurrence | MSI-H supersedes for ICI biomarker decision | Sha 2020 |
| Hypermutator characterization (>=100/Mb) | Confirm MMR-D or POLE-exo via signatures + IHC | Co-occurrence is common |
| Neoantigen quality (research) | Luksza fitness + HLA-LOH (LOHHLA / DASH) + clonality (McGranahan) | Beyond raw TMB |
| Cross-panel comparison | Vega 2021 calibration equations OR WES re-projection | Direct comparison invalid |

## Standard Workflow

**Goal:** Compute TMB from a VEP-annotated somatic VCF with full filtering.

**Approach:** Parse cyvcf2; apply VAF + depth + germline (gnomAD) filters; count nonsynonymous coding consequences; divide by scored Mb.

```python
from cyvcf2 import VCF
import re

NONSYNONYMOUS_CONSEQUENCES = {
    'missense_variant', 'stop_gained', 'stop_lost', 'start_lost', 'start_retained',
    'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
    'splice_donor_variant', 'splice_acceptor_variant',
    'protein_altering_variant', 'initiator_codon_variant'
}

# Vega 2021-calibrated scored regions (Mb)
PANEL_SCORED_REGION = {
    'FoundationOne_CDx': 0.8,           # Scored region; NOT 1.1 panel total
    'MSK_IMPACT_v3': 0.98,
    'MSK_IMPACT_v4': 1.22,
    'TSO500': 1.3,                       # Scored from 1.94 total
    'Oncomine_TML': 1.2,
    'Caris_MI': 1.2,
    'Tempus_xT_v3': 0.6,                 # Borderline reliability
    'WES': 30.0,
    'WGS': 3000.0
}

# Vega 2021 calibration: equivalent thresholds for FDA 10/Mb FoundationOne sensitivity
ASSAY_TMB_H_CUTOFF = {
    'FoundationOne_CDx': 10.0,
    'TSO500': 7.8,
    'Oncomine_TML': 8.4,
    'MSK_IMPACT_v3': 10.0,  # Approximate; full Vega 2021 calibration recommended
    'MSK_IMPACT_v4': 10.0,
    'WES': 10.0
}


def parse_consequences_from_vep(csq_field, csq_header):
    '''Parse VEP CSQ INFO field; returns list of per-transcript consequence types.'''
    if not csq_field:
        return []
    cons_idx = csq_header.index('Consequence')
    out = []
    for transcript in csq_field.split(','):
        fields = transcript.split('|')
        if len(fields) > cons_idx:
            out.append(fields[cons_idx])
    return out


def is_nonsynonymous(consequences, include_synonymous=False):
    '''Check if variant has nonsynonymous coding consequence.

    FoundationOne CDx convention INCLUDES synonymous (set include_synonymous=True).
    MSK-IMPACT and most academic pipelines exclude.
    '''
    target = set(NONSYNONYMOUS_CONSEQUENCES)
    if include_synonymous:
        target.add('synonymous_variant')
    for cons_str in consequences:
        for cons in cons_str.split('&'):
            if cons in target:
                return True
    return False


def calculate_tmb(vcf_path, scored_region_mb, csq_header,
                   min_vaf=0.05, min_depth=100, max_gnomad_af=0.005,
                   include_synonymous=False, exclude_hotspots=True,
                   hotspot_bed=None):
    '''Calculate TMB with filtering per Vega 2021 harmonization.

    Args:
        scored_region_mb: panel's SCORED region (NOT total panel)
        min_vaf: 0.05 (FoundationOne) to 0.10 (tumor-only no UMI)
        max_gnomad_af: 0.005 (0.5%) typical for tumor-only germline filter
        include_synonymous: True for FoundationOne CDx-compatible; False for MSK-IMPACT
        exclude_hotspots: COSMIC drivers excluded (not random mutations)
    '''
    vcf = VCF(vcf_path)
    cons_idx = csq_header.index('Consequence') if 'Consequence' in csq_header else 1
    nonsyn_count = 0
    total_pass = 0

    for v in vcf:
        if v.FILTER is not None:  # FILTER == None means PASS in cyvcf2
            continue
        depth = v.INFO.get('DP', 0)
        if depth < min_depth:
            continue
        vaf = _get_vaf(v)
        if vaf is None or vaf < min_vaf:
            continue
        gnomad_af = v.INFO.get('gnomAD_AF', 0) or v.INFO.get('AF_popmax', 0) or 0
        if gnomad_af > max_gnomad_af:
            continue
        total_pass += 1

        csq = v.INFO.get('CSQ', '')
        consequences = parse_consequences_from_vep(csq, csq_header)
        if is_nonsynonymous(consequences, include_synonymous=include_synonymous):
            nonsyn_count += 1

    tmb = nonsyn_count / scored_region_mb
    return {
        'tmb': round(tmb, 2),
        'nonsynonymous_count': nonsyn_count,
        'total_passing_filters': total_pass,
        'scored_region_mb': scored_region_mb
    }


def _get_vaf(variant):
    '''Extract VAF from genotype FORMAT fields (Mutect2 AD or AF).'''
    try:
        ad = variant.format('AD')
        if ad is not None and len(ad) > 0:
            ad0 = ad[0]
            total = sum(ad0)
            return ad0[1] / total if total > 0 else None
    except Exception:
        pass
    try:
        af = variant.format('AF')
        if af is not None and len(af) > 0:
            return float(af[0])
    except Exception:
        pass
    return None


def classify_tmb(tmb_value, assay='FoundationOne_CDx'):
    '''Apply ESMO 2024 / FDA pembrolizumab cutoff with Vega 2021 calibration per assay.'''
    cutoff = ASSAY_TMB_H_CUTOFF.get(assay, 10.0)
    if tmb_value >= 500:
        category = 'Ultra-hypermutator (>=500/Mb; POLE+MMR likely)'
    elif tmb_value >= 100:
        category = 'Hypermutator (>=100/Mb; MMR-D or POLE)'
    elif tmb_value >= cutoff:
        category = f'TMB-H (>= {cutoff}/Mb {assay}-calibrated; pan-tumor ICI eligible per FDA 2020)'
    else:
        category = 'TMB-low'
    return category
```

## TMB-MSI Concordance and Reconciliation

**Goal:** When MSI-H is present, TMB-H adds no information (Sha 2020).

```python
def tmb_msi_reconcile(tmb_value, msi_status, tumor_type=None):
    '''Reconcile TMB + MSI for ICI decision.'''
    tmb_high = tmb_value >= 10
    msi_high = msi_status == 'MSI-H'

    if msi_high:
        return ('ICI eligible by MSI-H (FDA 2017 pembrolizumab); TMB-H adds no information '
                '(Sha 2020 Cell Rep Med).')
    if tmb_high and tumor_type in ('breast', 'prostate', 'glioma'):
        return ('TMB-H present but NOT endorsed for this tumor type (ESMO 2024; McGrail 2021). '
                'Tumor-specific cutoffs recommended.')
    if tmb_high:
        return ('TMB-H pan-tumor; ICI eligible (FDA pembrolizumab 2020). '
                'Confirm baseline CD8 infiltrate; check HLA-LOH (Marty 2017 LOHHLA).')
    return 'TMB-low; MSS. Standard-of-care chemo.'
```

## Per-Operation Failure Modes

**1. Using panel total size as denominator (NOT scored region)**
- Trigger: Compute TMB = nonsynonymous count / 1.1 Mb for FoundationOne.
- Mechanism: FoundationOne CDx total panel is 1.1 Mb; SCORED region (counted for TMB denominator) is 0.8 Mb.
- Symptom: TMB underestimated by ~37%.
- Fix: Use 0.8 Mb for FoundationOne CDx scored region per Vega 2021.

**2. Cross-panel comparison without calibration**
- Trigger: TSO500 reports TMB = 9.5; compared to FoundationOne 10/Mb cutoff.
- Mechanism: Vega 2021 showed equivalent threshold is 7.8/Mb on TSO500 (not 10/Mb).
- Symptom: TSO500 TMB-H called positive at incorrect threshold.
- Fix: Apply assay-specific Vega 2021 calibration; TSO500 cutoff = 7.8/Mb.

**3. FoundationOne synonymous mis-handling**
- Trigger: Compare academic pipeline (no synonymous) to FoundationOne reference (synonymous included).
- Mechanism: FoundationOne CDx counts synonymous; MSK-IMPACT and most academic pipelines exclude.
- Symptom: Academic pipeline TMB systematically lower than FoundationOne by ~10-20%.
- Fix: Match the counting convention to the comparison reference; document explicitly.

**4. Tumor-only TMB inflated**
- Trigger: Tumor-only WES with naive germline filter (gnomAD AF > 1%).
- Mechanism: Population-specific common variants leak through if gnomAD AF threshold not stratified by ancestry.
- Symptom: AFR/EAS patient TMB inflated 1.5-3x; misclassified as TMB-H.
- Fix: Stratify gnomAD AF by patient ancestry; use grpmax FAF95; threshold <= 0.5%.

**5. bTMB applied without ctDNA shed check**
- Trigger: Report bTMB low in a metastatic patient.
- Mechanism: Low ctDNA shed fraction produces false-negative bTMB (BFAST Cohort C failure mechanism).
- Symptom: Patient with high tissue TMB labeled bTMB-low; ICI not offered.
- Fix: Check tumor fraction (e.g., ichorCNA, MAF of known driver) before trusting bTMB-low; consider tissue TMB.

**6. TMB-H applied to breast / prostate / glioma**
- Trigger: ICI prescribed for TMB-H breast cancer based on pan-tumor approval.
- Mechanism: McGrail 2021 demonstrated TMB fails to enrich for ICI response in breast, prostate, glioma; ESMO 2024 explicitly excludes these.
- Symptom: ICI offered with low expectation of benefit; patient bears unnecessary toxicity.
- Fix: Apply tumor-type-specific cutoffs (Samstein 2019); document ESMO 2024 caveat in report.

**7. Hotspots inflating TMB**
- Trigger: Include BRAF V600E and KRAS G12C in TMB count.
- Mechanism: Driver hotspots are non-random; including biases TMB upward in driver-mutated samples.
- Symptom: TMB inflated in samples with strong drivers.
- Fix: Exclude COSMIC-confirmed hotspots (provide hotspot BED).

**8. MSI-H -> add TMB-H -> additive ICI confidence**
- Trigger: Report TMB-H as additional support for ICI in MSI-H patient.
- Mechanism: MSI-H is uniformly hypermutator + uniformly ICI-responsive; adding TMB-H is statistical tautology (Sha 2020).
- Symptom: Reviewer flag.
- Fix: Report MSI-H + TMB-H concurrence but explicitly note TMB-H is NOT additive given MSI-H.

**9. Ignoring HLA-LOH**
- Trigger: TMB-H + neoantigen prediction without LOH check.
- Mechanism: HLA-LOH abolishes neoantigen presentation for lost allele in ~17% pan-cancer (>30% HNSCC / NSCLC / cervical).
- Symptom: Apparent neoantigen burden inflated.
- Fix: Run LOHHLA (Marty 2017) or DASH (Montesion 2021); flag HLA-LOH-positive tumors.

## Reconciliation: When Sources Disagree

| Pattern | Likely cause | Action |
|---------|-------------|--------|
| Vendor TMB vs WES TMB differ 2-3x | Panel-specific scored region + counting convention | Apply Vega 2021 calibration |
| FoundationOne vs MSK-IMPACT same sample differ | Synonymous handling differs | Document both; cite Vega 2021 |
| Tissue TMB vs bTMB differ | ctDNA shed fraction low; tumor heterogeneity | Trust tissue; check ctDNA fraction for bTMB confidence |
| TMB-H + MSI-H | Expected concurrence | MSI-H is the primary biomarker; TMB-H not additive |
| TMB-H + clinical PD-L1-negative | Independent biomarkers | Report both; ICI eligibility still per TMB-H pan-tumor |
| Patient with TMB-H but PR rate low | Tumor-type-specific cutoff; HLA-LOH | Apply Samstein 2019 cancer-specific cutoff; check HLA-LOH |
| POLE-exo + MMR-D | Ultra-hypermutator | ICI excellent response expected |

## Quantitative Thresholds and Conventions

| Threshold | Convention | Source |
|-----------|-----------|--------|
| FDA pembrolizumab pan-tumor | TMB-H >= 10 mut/Mb on FoundationOne CDx | FDA 2020 |
| TSO500 equivalent cutoff | 7.8 mut/Mb | Vega 2021 |
| Oncomine TML equivalent cutoff | 8.4 mut/Mb | Vega 2021 |
| Hypermutator | >= 100 mut/Mb | Research convention |
| Ultra-hypermutator | >= 500 mut/Mb | POLE+MMR; ICI excellent |
| MSI-H typical TMB | 30-50 mut/Mb | Salem 2018 |
| MSI-H + TMB-H overlap | ~80% MSI-H are TMB-H in CRC/endometrial; ~16% TMB-H are MSI-H | Salem 2018 |
| Tumor purity floor | FoundationOne >=20%; MSK-IMPACT >=30% | Vendor documentation |
| Min VAF | FoundationOne 5%; tumor-only no UMI 10% | Vendor documentation |
| Tumor-only germline filter | gnomAD AF <=0.5% (sometimes 1%) | Convention |
| Panel size minimum | >= 0.8 Mb workable; >= 1.0 Mb preferred; < 0.5 Mb unreliable | Vega 2021 |
| HLA-LOH frequency | ~17% pan-cancer; >30% HNSCC / NSCLC / cervical | Montesion 2021 |

## Common Errors

| Symptom | Cause | Solution |
|---------|-------|----------|
| TMB much lower than FoundationOne report | Used panel total (1.1) instead of scored (0.8) | Use 0.8 Mb for FoundationOne |
| Academic TMB systematically lower | Excluded synonymous; FoundationOne includes | Match counting convention |
| AFR / EAS tumor-only TMB inflated | gnomAD AF filter EUR-only | Use grpmax FAF95; stratify by patient ancestry |
| bTMB negative but tissue positive | Low ctDNA shed | Use tissue TMB; check fraction |
| TMB-H in breast cancer with poor response | ESMO 2024 / McGrail 2021 exclusion | Use tumor-type-specific cutoff |
| MSI-H + TMB-H reported as additive | Tautology | MSI-H is primary biomarker |
| POLE-exo + low TMB | Tumor sequencing artifact OR low tumor purity | Check VAF distribution; re-call if purity low |
| Variant counting differs across replicates | Random VAF sampling at borderline thresholds | Set explicit VAF floor + replicate-stable filter |

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Why panel-specific cutoffs?" | Vega 2021 demonstrated 3-fold panel variance; FoundationOne 10/Mb = TSO500 7.8 = Oncomine 8.4. Universal 10/Mb is wrong across non-F1 platforms. |
| "TMB-H is supposed to be tumor-agnostic" | FDA pan-tumor approval based on KEYNOTE-158; ESMO 2024 retained but excluded breast/prostate/glioma. McGrail 2021 + Samstein 2019 demonstrate tumor-type-specific cutoffs. |
| "Synonymous variants?" | FoundationOne CDx counts synonymous; academic pipelines exclude. We document the counting convention and apply Vega 2021 calibration. |
| "Why exclude hotspots?" | Driver hotspots are non-random; including biases TMB upward in driver-mutated samples vs cohort comparator. |
| "Tumor-only TMB unreliable" | Acknowledged; we apply stringent gnomAD grpmax FAF95 filtering stratified by patient ancestry; report paired-normal-validated subset separately. |
| "Why HLA-LOH integration?" | Marty 2017 + Montesion 2021 show ~17% pan-cancer (>30% HNSCC / NSCLC / cervical) lose HLA via LOH; apparent neoantigen burden over-estimated without LOH check. |
| "bTMB?" | BFAST Cohort C failed primary endpoint (Dziadziuszko 2022); bTMB is research-grade in tissue-naive only; we use tissue TMB as regulatory standard. |
| "Why ultra-hypermutator distinction?" | POLE+MMR (>=500 mut/Mb) shows superior ICI response per multiple case series; mechanistically distinct from MMR-D alone. |

## References

- Marabelle A et al. 2020. Association of TMB with efficacy of pembrolizumab in advanced solid tumours from the phase 2 KEYNOTE-158 study. *Lancet Oncol* 21:1353.
- Merino DM et al. 2020. Establishing guidelines to harmonize tumor mutational burden (TMB). *J Immunother Cancer* 8:e000147. (FoC Phase I)
- Vega DM et al. 2021. Aligning TMB analytical validation efforts: TMB harmonization Phase II project. *Ann Oncol* 32:1626.
- Mosele MF et al. 2024. Recommendations for the use of NGS for patients with metastatic cancers. *Ann Oncol* 35:588. (ESMO 2024)
- McGrail DJ et al. 2021. High TMB predicts response to immunotherapy independent of clinical or pathologic features. *Ann Oncol* 32:661.
- Sha D et al. 2020. TMB as a predictive biomarker in solid tumors. *Cell Rep Med* 1:100043.
- Samstein RM et al. 2019. TMB and survival after immunotherapy across cancer types. *Nat Genet* 51:202.
- Chalmers ZR et al. 2017. Analysis of 100,000 human cancer genomes reveals the landscape of TMB. *Genome Med* 9:34.
- Yarchoan M et al. 2017. Tumor mutational burden and response rate to PD-1 inhibition. *NEJM* 377:2500.
- Gandara DR et al. 2018. Blood-based TMB as a predictor of response to atezolizumab in NSCLC. *Nat Med* 24:1441.
- Dziadziuszko R et al. 2022. Atezolizumab vs platinum-based chemotherapy in blood-based TMB-positive NSCLC: BFAST Cohort C. *Nat Med* 28:2541.
- Salem ME et al. 2018. Landscape of TMB in different cancers: implications for ICI response. *Cancer Discov* 8:1136.
- Luksza M et al. 2017. A neoantigen fitness model predicts tumour response to checkpoint blockade immunotherapy. *Nature* 551:517.
- Luksza M et al. 2022. Neoantigen quality predicts immunoediting in survivors of pancreatic cancer. *Nature* 606:389.
- McGranahan N et al. 2016. Clonal neoantigens elicit T cell immunoreactivity and sensitivity to immune checkpoint blockade. *Science* 351:1463.
- Marty R et al. 2017. MHC-I genotype restricts the oncogenic mutational landscape. *Cell* 171:1272. (LOHHLA)
- Montesion M et al. 2021. Somatic HLA class I loss is a widespread mechanism of immune evasion which refines the use of TMB as a biomarker. *Cancer Discov* 11:282. (DASH)
- Friends of Cancer Research TMB harmonization resources: `https://friendsofcancerresearch.org/tmb/`

## Related Skills

- clinical-databases/somatic-signatures - Mutational signatures including HRD (PARP) and MMR-D (ICI)
- clinical-databases/msi-detection - MSI-H is the related ICI biomarker
- clinical-databases/hla-typing - HLA typing for neoantigen prediction and LOH
- variant-calling/variant-calling - Mutect2 / Strelka2 somatic upstream
- variant-calling/clinical-interpretation - ACMG / AMP cancer framework
