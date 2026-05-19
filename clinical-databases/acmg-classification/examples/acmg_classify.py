'''ACMG/AMP classification with Tavtigian point system + Pejaver 2022 calibration.

Reference: requests 2.31+, pandas 2.2+ | Verify GeneBe/InterVar API if version differs.
Implements: Tavtigian 2018/2020 point system; Pejaver 2022 REVEL PP3/BP4; Walker 2023 SpliceAI;
Whiffin 2017 BS1; Abou Tayoun 2018 PVS1 (decision-tree sketch); Brnich 2020 PS3/BS3 OddsPath.
'''
import requests
import pandas as pd


# Pejaver 2022 AJHG -- calibrated thresholds for REVEL.
# Pattern: lower bound <= score < upper bound triggers the indicated code.
REVEL_THRESHOLDS = [
    (float('-inf'), 0.003, 'BP4_VeryStrong'),
    (0.003, 0.016, 'BP4_Strong'),
    (0.016, 0.290, 'BP4_Moderate'),
    (0.290, 0.644, 'BP4_Supporting'),
    (0.644, 0.773, 'PP3_Supporting'),
    (0.773, 0.932, 'PP3_Moderate'),
    (0.932, float('inf'), 'PP3_Strong'),
]

BAYESDEL_NOAF_THRESHOLDS = [
    (float('-inf'), -0.36, 'BP4_Strong'),
    (-0.36, -0.18, 'BP4_Moderate'),
    (-0.18, -0.08, 'BP4_Supporting'),
    (-0.08, 0.13, None),
    (0.13, 0.27, 'PP3_Supporting'),
    (0.27, 0.50, 'PP3_Moderate'),
    (0.50, float('inf'), 'PP3_Strong'),
]

# Tavtigian 2020 Hum Mutat -- point assignments per strength.
STRENGTH_POINTS = {
    'PVS1_VeryStrong': 8, 'PVS1_Strong': 4, 'PVS1_Moderate': 2, 'PVS1_Supporting': 1,
    'PS1': 4, 'PS2': 4, 'PS3': 4, 'PS3_Moderate': 2, 'PS3_Supporting': 1, 'PS4': 4,
    'PS4_Moderate': 2, 'PS4_Supporting': 1,
    'PM1': 2, 'PM2_Supporting': 1, 'PM2': 1,  # PM2 was downgraded to Supporting in SVI 2020
    'PM3': 2, 'PM3_Strong': 4, 'PM3_VeryStrong': 8, 'PM4': 2, 'PM5': 2, 'PM6': 2,
    'PP1': 1, 'PP1_Moderate': 2, 'PP1_Strong': 4, 'PP2': 1,
    'PP3_Supporting': 1, 'PP3_Moderate': 2, 'PP3_Strong': 4, 'PP4': 1, 'PP5': 1,
    # Benign codes (Tavtigian negative-signed)
    'BA1': -100,                # Standalone benign
    'BS1': -4, 'BS2': -4, 'BS3': -4, 'BS3_Moderate': -2, 'BS3_Supporting': -1, 'BS4': -4,
    'BP1': -1, 'BP2': -1, 'BP3': -1,
    'BP4_Supporting': -1, 'BP4_Moderate': -2, 'BP4_Strong': -4, 'BP4_VeryStrong': -8,
    'BP5': -1, 'BP6': -1, 'BP7': -1
}


def pejaver_calibrate(score, thresholds):
    '''Map a score to ACMG strength per Pejaver 2022 calibration. ONE predictor only.'''
    if score is None:
        return None
    for lo, hi, code in thresholds:
        if lo <= score < hi:
            return code
    return None


def revel_pp3_bp4(revel_score):
    '''REVEL score -> PP3/BP4 strength per Pejaver 2022.'''
    return pejaver_calibrate(revel_score, REVEL_THRESHOLDS)


def bayesdel_pp3_bp4(bayesdel_noaf_score):
    '''BayesDel-noAF score -> PP3/BP4 strength per Pejaver 2022.'''
    return pejaver_calibrate(bayesdel_noaf_score, BAYESDEL_NOAF_THRESHOLDS)


def spliceai_walker2023(ds_max, has_corroborating_evidence=False):
    '''Walker 2023 SVI Splicing Subgroup framework.

    SpliceAI DS_max >= 0.5 + corroborating evidence -> PP3_Strong.
    SpliceAI >= 0.20 -> minimum threshold for ANY splicing PP3 (PP3_Supporting default).
    SpliceAI < 0.1 -> BP4_Moderate.
    '''
    if ds_max is None:
        return None
    if ds_max >= 0.5 and has_corroborating_evidence:
        return 'PP3_Strong'
    if ds_max >= 0.5:
        return 'PP3_Moderate'  # Without corroboration, downgrade
    if ds_max >= 0.20:
        return 'PP3_Supporting'
    if ds_max < 0.1:
        return 'BP4_Moderate'
    return None


def alphamissense_supporting_only(am_score):
    '''AlphaMissense is supporting-only as of May 2026 (ClinGen not endorsed PP3 calibration).

    Cheng 2023 developer threshold 0.564 is NOT the Pejaver-style PP3 calibration.
    Schmidt 2025 calibration is preliminary; use REVEL/BayesDel for PP3.
    '''
    if am_score is None:
        return None
    if am_score >= 0.7:
        return 'PP3_Supporting'
    if am_score <= 0.2:
        return 'BP4_Supporting'
    return None


def pvs1_decision_tree(variant_type, is_in_5prime_or_pre_last_exon_junction,
                        coding_pct_removed, in_critical_region,
                        is_alt_isoform_in_disease_tissue=True):
    '''Abou Tayoun 2018 PVS1 decision tree (simplified -- check VCEP-specific CSpec for refinements).

    Args:
        variant_type: 'nonsense' | 'frameshift' | 'splice_donor' | 'splice_acceptor' |
                       'initiation' | 'single_exon_del' | 'multi_exon_del'
        is_in_5prime_or_pre_last_exon_junction: True if NMD-predicted
        coding_pct_removed: fraction of coding sequence removed (e.g., 0.15 for 15%)
        in_critical_region: True if removes critical functional domain
        is_alt_isoform_in_disease_tissue: True if alt isoform is disease-relevant

    Returns: 'PVS1_VeryStrong' / 'PVS1_Strong' / 'PVS1_Moderate' / 'PVS1_Supporting' / None
    '''
    if variant_type in ('nonsense', 'frameshift'):
        if is_in_5prime_or_pre_last_exon_junction:
            return 'PVS1_VeryStrong'  # NMD-triggered
        if in_critical_region or coding_pct_removed > 0.10:
            return 'PVS1_Strong'
        return 'PVS1_Moderate'
    if variant_type in ('splice_donor', 'splice_acceptor'):
        if is_in_5prime_or_pre_last_exon_junction:
            return 'PVS1_VeryStrong'
        return 'PVS1_Strong'
    if variant_type == 'initiation':
        return 'PVS1_Moderate'   # Per Abou Tayoun 2018
    if variant_type in ('single_exon_del', 'multi_exon_del'):
        if in_critical_region:
            return 'PVS1_VeryStrong'
        return 'PVS1_Strong'
    return None


def whiffin_max_credible_af(prevalence, max_allelic_contribution=1.0,
                              max_genetic_contribution=1.0, penetrance=1.0):
    '''Whiffin 2017 max-credible-AF for BS1 application.

    For dominant inheritance. For autosomal recessive, transform appropriately.
    '''
    return (prevalence * max_genetic_contribution * max_allelic_contribution) / (penetrance * 2)


def bs1_ba1(grpmax_faf95, max_credible_af, ba1_threshold=0.05):
    '''Apply ClinGen SVI BS1/BA1/PM2_Supporting from gnomAD grpmax FAF95.

    BA1 default 5% in non-bottleneck group (VCEP override common; Hearing Loss = 0.5% AR).
    BS1 = gene-specific Whiffin max-credible-AF.
    PM2_Supporting = absent or ultra-rare (downgraded from PM2_Moderate in SVI 2020).
    '''
    if grpmax_faf95 is None or grpmax_faf95 == 0.0:
        return 'PM2_Supporting'
    if grpmax_faf95 > ba1_threshold:
        return 'BA1'
    if grpmax_faf95 > max_credible_af:
        return 'BS1'
    return None


def ps3_oddspath(odds_path):
    '''Brnich 2020 PS3/BS3 OddsPath calibration.

    OddsPath > 18.7 -> Very Strong; 4.3-18.7 -> Strong; 2.1-4.3 -> Moderate; 1.2-2.1 -> Supporting.
    Mirror values for benignity (BS3).
    '''
    if odds_path is None:
        return None
    if odds_path > 18.7:
        return 'PS3'  # Strong (Tavtigian 4 points); upgrade to VeryStrong via custom
    if odds_path > 4.3:
        return 'PS3'  # Strong (default)
    if odds_path > 2.1:
        return 'PS3_Moderate'
    if odds_path > 1.2:
        return 'PS3_Supporting'
    # Benign side
    if odds_path < (1 / 18.7):
        return 'BS3'
    if odds_path < (1 / 4.3):
        return 'BS3'
    if odds_path < (1 / 2.1):
        return 'BS3_Moderate'
    if odds_path < (1 / 1.2):
        return 'BS3_Supporting'
    return None


def tavtigian_classify(criteria_assigned):
    '''Sum Tavtigian points and classify P / LP / VUS / LB / B.'''
    points = sum(STRENGTH_POINTS.get(c, 0) for c in criteria_assigned)
    if any(c == 'BA1' for c in criteria_assigned):
        return {'classification': 'Benign', 'points': points,
                'criteria': criteria_assigned, 'rationale': 'BA1 standalone'}
    if points >= 10:
        category = 'Pathogenic'
    elif points >= 6:
        category = 'Likely Pathogenic'
    elif points >= 0:
        category = 'VUS'
    elif points >= -6:
        category = 'Likely Benign'
    else:
        category = 'Benign'
    return {'classification': category, 'points': points, 'criteria': criteria_assigned}


def classify_with_subsumption(criteria_assigned):
    '''Apply Abou Tayoun PVS1 subsumption: PVS1 + PP3 -> only PVS1 (PP3 subsumed).

    Also: PVS1 + PM4 -> only PVS1.
    '''
    has_pvs1 = any(c.startswith('PVS1') for c in criteria_assigned)
    if has_pvs1:
        filtered = [c for c in criteria_assigned
                    if not c.startswith('PP3') and c != 'PM4']
    else:
        filtered = criteria_assigned
    return tavtigian_classify(filtered)


def genebe_api(hgvs):
    '''Query GeneBe automated ACMG classifier (Stawinski 2024).

    Open-source; Tavtigian-point-system-based; comparable to VarSome (commercial).
    '''
    r = requests.get('https://api.genebe.net/cloud/api-public/v1/variant',
                     params={'variant': hgvs, 'genome': 'hg38'},
                     timeout=30)
    r.raise_for_status()
    return r.json()


# Cancer somatic: Li 2017 AMP/ASCO/CAP Tier I-IV
def cancer_amp_tier(biomarker_in_oncokb_level, in_nccn_guidelines=False,
                     same_tumor_type=False, evidence_quality='preclinical'):
    '''Apply AMP/ASCO/CAP 2017 Tier I-IV.

    Tier I-A: FDA drug for same tumor + biomarker -> on-label
    Tier I-B: professional guideline (NCCN/ESMO) for same tumor
    Tier II-C: FDA drug in different tumor (off-label / basket)
    Tier II-D: preclinical / investigational
    Tier III: VUS-somatic; Tier IV: benign-somatic
    '''
    if biomarker_in_oncokb_level == 1 and same_tumor_type:
        return 'Tier I-A: FDA drug + same tumor + biomarker (on-label)'
    if in_nccn_guidelines and same_tumor_type:
        return 'Tier I-B: professional guideline (NCCN/ESMO)'
    if biomarker_in_oncokb_level <= 3 and not same_tumor_type:
        return 'Tier II-C: FDA drug in different tumor (off-label / basket)'
    if evidence_quality in ('preclinical', 'case_report'):
        return 'Tier II-D: investigational / preclinical'
    return 'Tier III/IV: VUS or benign-somatic'


if __name__ == '__main__':
    # Example: BRCA1 missense
    revel_score = 0.95
    am_score = 0.85
    spliceai_max = 0.05
    grpmax_faf95 = 0.0
    odds_path_ps3 = 8.0

    criteria = []
    revel_call = revel_pp3_bp4(revel_score)
    if revel_call:
        criteria.append(revel_call)
    spliceai_call = spliceai_walker2023(spliceai_max)
    if spliceai_call:
        criteria.append(spliceai_call)
    ps3_call = ps3_oddspath(odds_path_ps3)
    if ps3_call:
        criteria.append(ps3_call)
    bs1_call = bs1_ba1(grpmax_faf95, max_credible_af=1e-6, ba1_threshold=0.05)
    if bs1_call:
        criteria.append(bs1_call)

    result = classify_with_subsumption(criteria)
    print(f'Criteria: {criteria}')
    print(f'Tavtigian points: {result["points"]} -> {result["classification"]}')

    # Subsumption test: PVS1 + PP3 -> only PVS1
    test_criteria = ['PVS1_VeryStrong', 'PP3_Strong', 'PM2_Supporting']
    result2 = classify_with_subsumption(test_criteria)
    print(f'\nSubsumption example (PVS1 + PP3 -> PVS1 only):')
    print(f'  Input: {test_criteria}')
    print(f'  Final criteria: {[c for c in test_criteria if not c.startswith("PP3")]}')
    print(f'  Classification: {result2["classification"]} (points: {result2["points"]})')
