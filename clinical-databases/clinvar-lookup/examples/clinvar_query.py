'''Query ClinVar for variant pathogenicity via REST API and local VCF'''

import requests
import time

def search_clinvar(term, retmax=100):
    '''Search ClinVar and return variant IDs'''
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {'db': 'clinvar', 'term': term, 'retmax': retmax, 'retmode': 'json'}
    response = requests.get(url, params=params)
    data = response.json()
    return data.get('esearchresult', {}).get('idlist', [])

def get_clinvar_summary(ids):
    '''Get summary for ClinVar variation IDs'''
    if not ids:
        return []
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    params = {'db': 'clinvar', 'id': ','.join(ids), 'retmode': 'json'}
    response = requests.get(url, params=params)
    data = response.json()
    results = data.get('result', {})
    return [results[vid] for vid in ids if vid in results]

# Search for pathogenic BRCA1 variants
# Limit to 10 for demonstration; typical searches may return hundreds
gene = 'BRCA1'
ids = search_clinvar(f'{gene}[gene] AND pathogenic[clinical_significance]', retmax=10)
print(f'Found {len(ids)} pathogenic {gene} variants (showing first 10)')

# NCBI rate limit: 3 requests/second without API key
time.sleep(0.5)

summaries = get_clinvar_summary(ids)
for s in summaries[:5]:
    title = s.get('title', 'N/A')
    clinsig = s.get('clinical_significance', {}).get('description', 'N/A')
    review = s.get('clinical_significance', {}).get('review_status', 'N/A')
    print(f'  {title}: {clinsig} ({review})')

# Parse clinical significance into categories
def categorize_significance(clinsig):
    '''Categorize ClinVar clinical significance

    Categories:
    - pathogenic: Disease-causing (Pathogenic, Likely pathogenic)
    - benign: Not disease-causing (Benign, Likely benign)
    - vus: Uncertain significance
    - conflicting: Multiple labs disagree
    '''
    clinsig_lower = clinsig.lower()
    if 'pathogenic' in clinsig_lower and 'likely' not in clinsig_lower:
        return 'pathogenic'
    elif 'likely pathogenic' in clinsig_lower:
        return 'likely_pathogenic'
    elif 'benign' in clinsig_lower and 'likely' not in clinsig_lower:
        return 'benign'
    elif 'likely benign' in clinsig_lower:
        return 'likely_benign'
    elif 'conflicting' in clinsig_lower:
        return 'conflicting'
    else:
        return 'vus'

# Review status to star rating
# Higher stars = more confidence in classification
REVIEW_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, single submitter': 1,
    'no assertion criteria provided': 0,
    'no assertion provided': 0
}

print('\nClinVar review status guide:')
for status, stars in REVIEW_STARS.items():
    print(f'  {"*" * stars if stars else "-"} {status}')
