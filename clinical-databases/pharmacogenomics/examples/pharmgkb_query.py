'''Pharmacogenomic workflow: PharmGKB + CYP2D6 activity score (Caudle 2020) + DPYD AS 2024.

Reference: requests 2.31+, pandas 2.2+ | Verify CPIC guideline versions if differs.
Activity values follow Caudle 2020 *Clin Transl Sci* (CYP2D6 *10 reset to 0.25).
DPYD follows CPIC 2024 (Lam et al. *Clin Pharmacol Ther*) activity-score framework.
'''
import requests
import time
import pandas as pd

PHARMGKB = 'https://api.pharmgkb.org/v1'

# Caudle 2020 activity values; *10 reset from 0.5 (pre-2020) to 0.25 (2020+).
# East-Asian populations carrying *10 were reclassified from NM to IM.
CYP2D6_ACTIVITY = {
    '*1': 1.0, '*2': 1.0, '*35': 1.0,
    '*3': 0.0, '*4': 0.0, '*5': 0.0, '*6': 0.0, '*7': 0.0, '*8': 0.0,
    '*11': 0.0, '*12': 0.0, '*14': 0.0, '*15': 0.0, '*19': 0.0, '*20': 0.0,
    '*36': 0.0, '*40': 0.0, '*42': 0.0, '*68': 0.0, '*13': 0.0,
    '*9': 0.5, '*41': 0.5, '*17': 0.5, '*29': 0.5,
    '*10': 0.25,
}

# CYP2C19
CYP2C19_ACTIVITY = {
    '*1': 1.0, '*17': 1.5,
    '*2': 0.0, '*3': 0.0, '*4': 0.0, '*5': 0.0, '*6': 0.0, '*7': 0.0, '*8': 0.0
}

# CYP2C9 (Daneshjou 2014 emphasizes including *5/*6/*8/*11 for African-ancestry warfarin)
CYP2C9_ACTIVITY = {
    '*1': 1.0,
    '*2': 0.5, '*3': 0.0,
    '*5': 0.0, '*6': 0.0, '*8': 0.5, '*11': 0.5  # AFR-common; the COAG failure variants
}

# DPYD CPIC 2024 activity values (Lam et al. Clin Pharmacol Ther)
DPYD_2024_ACTIVITY = {
    'c.1905+1G>A': 0.0,  # DPYD*2A; splice donor
    'c.1679T>G': 0.0,    # DPYD*13; p.I560S
    'c.2846A>T': 0.5,    # p.D949V
    'HapB3': 0.5,        # c.1129-5923C>G + c.1236G>A linked
    # c.85T>C (DPYD*9A) is NOT in CPIC 2024 actionable set; evidence does not support clinical decrement
}


def _allele_activity(allele_str, table):
    '''Look up per-allele activity. Handle xN copy-number suffix correctly.

    Key footgun: *4xN is clinically silent. *4 has activity 0; 0 * N = 0.
    Only functional alleles (*1, *2, *35) become UM when amplified.
    '''
    if 'x' in allele_str:
        base, n = allele_str.split('x')
        copies = int(n) if n.isdigit() else 2  # 'N' usually >=2; assume 2 unless quantified
        return table.get(base, 1.0) * copies
    return table.get(allele_str, 1.0)


def cyp2d6_phenotype(diplotype):
    '''CYP2D6 diplotype -> activity score + Caudle 2020 phenotype bin.

    Diplotype string format: '*1/*4' or '*4xN/*10' or '*2xN/*17'
    '''
    left, right = diplotype.split('/')
    left_score = _allele_activity(left, CYP2D6_ACTIVITY)
    right_score = _allele_activity(right, CYP2D6_ACTIVITY)
    total = left_score + right_score
    if total == 0:
        phenotype = 'Poor Metabolizer'
    elif total < 1.25:
        phenotype = 'Intermediate Metabolizer'
    elif total <= 2.25:
        phenotype = 'Normal Metabolizer'
    else:
        phenotype = 'Ultrarapid Metabolizer'
    return {'diplotype': diplotype, 'left_activity': left_score,
            'right_activity': right_score, 'activity_score': total, 'phenotype': phenotype}


def cyp2c19_phenotype(diplotype):
    '''CYP2C19 diplotype -> phenotype. Note RM phenotype bin between NM and UM.'''
    left, right = diplotype.split('/')
    total = CYP2C19_ACTIVITY.get(left, 1.0) + CYP2C19_ACTIVITY.get(right, 1.0)
    if total == 0:
        phenotype = 'Poor Metabolizer'
    elif total < 1.5:
        phenotype = 'Intermediate Metabolizer'
    elif total <= 2.0:
        phenotype = 'Normal Metabolizer'
    elif total <= 2.5:
        phenotype = 'Rapid Metabolizer'
    else:
        phenotype = 'Ultrarapid Metabolizer'
    return {'diplotype': diplotype, 'activity_score': total, 'phenotype': phenotype}


def dpyd_activity(variants):
    '''CPIC 2024 DPYD gene activity score.

    Sum the two lowest activities across the two alleles.
    AS 2.0: full dose; AS 1.5: 50% start + TDM; AS 1.0: 50% start + TDM; AS 0: avoid.
    '''
    activities = sorted([DPYD_2024_ACTIVITY.get(v, 1.0) for v in variants])
    gene_as = sum(activities[:2])
    if gene_as >= 1.99:
        dose = 'Full dose'
    elif gene_as >= 1.0:
        dose = '50% starting dose + therapeutic drug monitoring'
    else:
        dose = 'Avoid fluoropyrimidines'
    return {'variants': variants, 'gene_activity_score': gene_as, 'dosing': dose}


def pharmgkb_clinical_annotations(gene_symbol):
    '''PharmGKB clinical annotations for a gene; PharmGKB level 1A/1B/2A/2B/3/4.'''
    r = requests.get(f'{PHARMGKB}/data/clinicalAnnotation',
                     params={'view': 'base', 'location.genes.symbol': gene_symbol},
                     timeout=30)
    r.raise_for_status()
    return r.json().get('data', [])


def cpic_guidelines(gene_symbol):
    '''CPIC guidelines for a gene via PharmGKB API.'''
    r = requests.get(f'{PHARMGKB}/data/guideline',
                     params={'view': 'base', 'relatedGenes.symbol': gene_symbol, 'source': 'CPIC'},
                     timeout=30)
    r.raise_for_status()
    return r.json().get('data', [])


def hla_pgx_screen(hla_alleles_4field):
    '''Screen 4-field HLA alleles against established PGx contraindication table.

    Critical: requires 4-field resolution. B*57:01 (abacavir risk) vs B*57:03 (no risk).
    HLA-B*35:02 (minocycline DILI) vs B*35:01 (TMP-SMX DILI) -- different drugs.
    '''
    pgx_alleles = {
        'B*57:01': {'drug': 'Abacavir', 'reaction': 'HSS', 'population': 'All (5-8% NFE)',
                    'or': 100, 'cite': 'Mallal 2008 NEJM (PREDICT-1)'},
        'B*15:02': {'drug': 'Carbamazepine/oxcarbazepine', 'reaction': 'SJS/TEN',
                    'population': 'Han Chinese, Thai, Malay, Indian', 'or': 2500,
                    'cite': 'Chung 2004 Nature; FDA black-box 2007'},
        'A*31:01': {'drug': 'Carbamazepine', 'reaction': 'DRESS/SJS', 'population': 'EUR, Japanese',
                    'or': 12, 'cite': 'McCormack 2011 NEJM'},
        'B*58:01': {'drug': 'Allopurinol', 'reaction': 'SJS/TEN', 'population': 'Han Chinese, Korean, Thai',
                    'or': 580, 'cite': 'Hung 2005 PNAS'},
        'B*13:01': {'drug': 'Dapsone', 'reaction': 'DDS', 'population': 'Han Chinese, SE Asian',
                    'cite': 'Zhang 2013 NEJM'},
        'B*35:02': {'drug': 'Minocycline', 'reaction': 'DILI', 'population': 'All',
                    'cite': 'Urban 2017; NOT *35:01'},
        'B*35:01': {'drug': 'TMP-SMX', 'reaction': 'DILI', 'population': 'Mixed',
                    'cite': 'Li 2021 Hepatology'},
        'B*14:01': {'drug': 'TMP-SMX', 'reaction': 'DILI', 'population': 'African',
                    'cite': 'Li 2021'},
        'B*15:13': {'drug': 'Phenytoin', 'reaction': 'SJS', 'population': 'Malaysian',
                    'cite': 'Chang 2017'},
    }
    matches = []
    for allele in hla_alleles_4field:
        normalized = allele.replace('HLA-', '')
        if normalized in pgx_alleles:
            matches.append({'allele': allele, **pgx_alleles[normalized]})
    return matches


if __name__ == '__main__':
    print('=== CYP2D6 phenotype examples (Caudle 2020) ===')
    examples = ['*1/*1', '*1/*4', '*4/*4', '*1/*10', '*4xN/*10', '*1xN/*4', '*2xN/*2']
    for dip in examples:
        result = cyp2d6_phenotype(dip)
        print(f"  {dip:20s} AS={result['activity_score']:.2f}  {result['phenotype']}")

    print('\n=== DPYD CPIC 2024 ===')
    case1 = dpyd_activity(['c.1905+1G>A'])  # Heterozygous *2A
    case2 = dpyd_activity(['c.2846A>T', 'c.2846A>T'])  # Homozygous p.D949V
    print(f"  Het *2A: AS={case1['gene_activity_score']}, dosing: {case1['dosing']}")
    print(f"  Hom D949V: AS={case2['gene_activity_score']}, dosing: {case2['dosing']}")

    print('\n=== HLA PGx screen example ===')
    hits = hla_pgx_screen(['HLA-B*57:01', 'HLA-A*02:01', 'HLA-B*15:02'])
    for h in hits:
        print(f"  {h['allele']}: {h['drug']} -- {h['reaction']} ({h.get('or', '?')}x OR; {h['cite']})")

    print('\n=== PharmGKB clinical annotations for CYP2D6 (first 5) ===')
    annotations = pharmgkb_clinical_annotations('CYP2D6')
    for ann in annotations[:5]:
        drug = ann['chemicals'][0]['name'] if ann.get('chemicals') else '?'
        level = ann.get('levelOfEvidence', '?')
        print(f"  CYP2D6 + {drug}: Level {level}")
