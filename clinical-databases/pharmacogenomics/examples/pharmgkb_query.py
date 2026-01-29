'''Query PharmGKB for pharmacogenomic annotations and interpret star alleles'''

import requests
import time

# --- PharmGKB API ---
# Free access, no API key required
# Rate limit: Be respectful, add delays for batch queries

def get_pharmgkb_annotations(gene_symbol):
    '''Get PharmGKB clinical annotations for a gene'''
    url = 'https://api.pharmgkb.org/v1/data/clinicalAnnotation'
    params = {'view': 'base', 'location.genes.symbol': gene_symbol}
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json().get('data', [])
    return []

def get_cpic_guidelines(gene_symbol):
    '''Get CPIC dosing guidelines for a gene'''
    url = 'https://api.pharmgkb.org/v1/data/guideline'
    params = {'view': 'base', 'relatedGenes.symbol': gene_symbol, 'source': 'CPIC'}
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json().get('data', [])
    return []

# --- Star Allele Activity Scores ---
# Based on CPIC Clinical Pharmacogenetics Implementation Consortium guidelines
# https://cpicpgx.org/

# CYP2D6 activity values
# Score 0 = no function, 0.5 = decreased, 1.0 = normal
CYP2D6_ACTIVITY = {
    '*1': 1.0,   # Normal (reference)
    '*2': 1.0,   # Normal
    '*3': 0.0,   # No function - frameshift
    '*4': 0.0,   # No function - splicing defect (most common loss-of-function in Europeans)
    '*5': 0.0,   # Gene deletion
    '*6': 0.0,   # No function - frameshift
    '*7': 0.0,   # No function
    '*8': 0.0,   # No function
    '*9': 0.5,   # Decreased function
    '*10': 0.25, # Decreased function (most common in East Asians)
    '*14': 0.0,  # No function
    '*17': 0.5,  # Decreased function
    '*29': 0.5,  # Decreased function
    '*41': 0.5,  # Decreased function
}

# CYP2C19 activity values
CYP2C19_ACTIVITY = {
    '*1': 1.0,   # Normal (reference)
    '*2': 0.0,   # No function (c.681G>A, most common loss-of-function)
    '*3': 0.0,   # No function
    '*4': 0.0,   # No function
    '*17': 1.5,  # Increased function
}

# CYP2C9 activity values
CYP2C9_ACTIVITY = {
    '*1': 1.0,   # Normal (reference)
    '*2': 0.5,   # Decreased function
    '*3': 0.0,   # No function (or very reduced)
}


def calculate_activity_score(allele1, allele2, gene='CYP2D6'):
    '''Calculate activity score from diplotype

    Args:
        allele1, allele2: Star alleles (e.g., '*1', '*4')
        gene: Gene name for activity table lookup
    '''
    activity_tables = {
        'CYP2D6': CYP2D6_ACTIVITY,
        'CYP2C19': CYP2C19_ACTIVITY,
        'CYP2C9': CYP2C9_ACTIVITY,
    }

    table = activity_tables.get(gene, {})
    # Default to normal function if allele not in table
    score1 = table.get(allele1, 1.0)
    score2 = table.get(allele2, 1.0)
    return score1 + score2


def get_metabolizer_phenotype(activity_score, gene='CYP2D6'):
    '''Convert activity score to metabolizer phenotype

    Thresholds based on CPIC guidelines for CYP2D6:
    - PM (Poor Metabolizer): 0
    - IM (Intermediate Metabolizer): >0 to <=1.25
    - NM (Normal Metabolizer): >1.25 to <=2.25
    - UM (Ultrarapid Metabolizer): >2.25 (requires gene duplication)
    '''
    if gene == 'CYP2D6':
        if activity_score == 0:
            return 'Poor Metabolizer (PM)'
        elif activity_score <= 1.25:
            return 'Intermediate Metabolizer (IM)'
        elif activity_score <= 2.25:
            return 'Normal Metabolizer (NM)'
        else:
            return 'Ultrarapid Metabolizer (UM)'

    elif gene == 'CYP2C19':
        if activity_score == 0:
            return 'Poor Metabolizer (PM)'
        elif activity_score < 1.5:
            return 'Intermediate Metabolizer (IM)'
        elif activity_score <= 2.0:
            return 'Normal Metabolizer (NM)'
        elif activity_score <= 2.5:
            return 'Rapid Metabolizer (RM)'
        else:
            return 'Ultrarapid Metabolizer (UM)'

    return 'Unknown'


# Example: Interpret CYP2D6 diplotype
print('=== CYP2D6 Metabolizer Status ===')
diplotypes = [
    ('*1', '*1'),   # Normal/Normal
    ('*1', '*4'),   # Normal/No function
    ('*4', '*4'),   # No function/No function (PM)
    ('*1', '*10'),  # Normal/Decreased
    ('*1', '*41'),  # Normal/Decreased
]

for a1, a2 in diplotypes:
    score = calculate_activity_score(a1, a2, 'CYP2D6')
    phenotype = get_metabolizer_phenotype(score, 'CYP2D6')
    print(f'  {a1}/{a2}: Score={score:.2f} -> {phenotype}')

print()

# Query PharmGKB for CYP2D6 drug interactions
print('=== CYP2D6 Drug Annotations (top 5) ===')
annotations = get_pharmgkb_annotations('CYP2D6')
for ann in annotations[:5]:
    gene = ann['location']['genes'][0]['symbol'] if ann.get('location', {}).get('genes') else 'N/A'
    drug = ann['chemicals'][0]['name'] if ann.get('chemicals') else 'N/A'
    level = ann.get('levelOfEvidence', 'N/A')
    print(f'  {gene} + {drug}: Level {level}')

time.sleep(0.5)  # Rate limit courtesy

# Get CPIC guidelines
print()
print('=== CPIC Guidelines for CYP2D6 ===')
guidelines = get_cpic_guidelines('CYP2D6')
for g in guidelines[:5]:
    name = g.get('name', 'N/A')
    drug = g['chemicals'][0]['name'] if g.get('chemicals') else 'N/A'
    print(f'  {name}')
