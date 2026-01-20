'''Metabolite annotation from m/z values'''
import pandas as pd
import requests

HMDB_API = 'https://hmdb.ca/metabolites/search'

def query_hmdb(mass, ppm=10):
    '''Query HMDB by exact mass'''
    tolerance = mass * ppm / 1e6
    params = {
        'query[mass_start]': mass - tolerance,
        'query[mass_end]': mass + tolerance,
        'format': 'json'
    }

    # Note: Actual HMDB API may differ; this is illustrative
    # Consider using a local database for production
    try:
        response = requests.get(HMDB_API, params=params, timeout=10)
        if response.ok:
            return response.json()
    except Exception:
        pass
    return []

def annotate_feature(mz, adduct='[M+H]+', ppm=10):
    '''Annotate a single feature'''
    adduct_masses = {
        '[M+H]+': 1.007276,
        '[M+Na]+': 22.989218,
        '[M-H]-': -1.007276,
    }

    neutral_mass = mz - adduct_masses.get(adduct, 0)
    matches = query_hmdb(neutral_mass, ppm)

    if matches:
        best = matches[0]
        return {
            'compound_name': best.get('name'),
            'compound_id': best.get('hmdb_id'),
            'formula': best.get('formula'),
            'confidence': 4  # Mass match only
        }
    return {'compound_name': None, 'compound_id': None, 'formula': None, 'confidence': 5}


# Load feature table
features = pd.read_csv('feature_table.csv')

# Annotate
print('Annotating features...')
annotations = [annotate_feature(mz) for mz in features['mz']]
annotation_df = pd.DataFrame(annotations)

# Combine
result = pd.concat([features, annotation_df], axis=1)

# Summary
print(f'\nAnnotation summary:')
print(result['confidence'].value_counts().sort_index())

# Save
result.to_csv('annotated_features.csv', index=False)
print('\nSaved to annotated_features.csv')
