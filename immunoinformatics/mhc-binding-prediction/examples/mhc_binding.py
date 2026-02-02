'''MHC binding prediction with MHCflurry'''

import pandas as pd


def predict_mhc_binding(peptides, alleles):
    '''Predict MHC class I binding using MHCflurry

    MHCflurry outputs:
    - mhcflurry_affinity: Predicted IC50 in nM
    - mhcflurry_affinity_percentile: Percentile rank among random peptides
    - mhcflurry_presentation_score: Combined binding + processing

    Interpretation:
    - IC50 < 50 nM: Strong binder
    - IC50 50-500 nM: Moderate binder
    - Percentile < 0.5%: Strong binder
    - Percentile 0.5-2%: Moderate binder
    '''
    try:
        from mhcflurry import Class1PresentationPredictor
        predictor = Class1PresentationPredictor.load()

        results = []
        for peptide in peptides:
            for allele in alleles:
                pred = predictor.predict(peptides=[peptide], alleles=[allele])
                results.append({
                    'peptide': peptide,
                    'allele': allele,
                    'affinity_nM': pred['mhcflurry_affinity'].values[0],
                    'percentile': pred['mhcflurry_affinity_percentile'].values[0],
                    'presentation': pred['mhcflurry_presentation_score'].values[0]
                })

        return pd.DataFrame(results)
    except ImportError:
        print('MHCflurry not installed. Install with: pip install mhcflurry')
        return None


def classify_binding(affinity_nm):
    '''Classify binding strength by IC50

    Standard thresholds used in literature:
    - Strong: <50 nM (high-affinity binder)
    - Moderate: 50-500 nM (potential epitope)
    - Weak: 500-5000 nM (unlikely functional)
    - Non-binder: >5000 nM
    '''
    if affinity_nm < 50:
        return 'Strong'
    elif affinity_nm < 500:
        return 'Moderate'
    elif affinity_nm < 5000:
        return 'Weak'
    else:
        return 'Non-binder'


def scan_protein(protein_seq, alleles, lengths=[8, 9, 10, 11], percentile_threshold=2.0):
    '''Scan protein for potential MHC epitopes

    MHC class I typically binds 8-11 amino acid peptides.
    9-mers are most common.

    Args:
        percentile_threshold: Include peptides with rank below this
                             2.0% is standard cutoff for potential binders
    '''
    try:
        from mhcflurry import Class1PresentationPredictor
        predictor = Class1PresentationPredictor.load()

        epitopes = []
        for length in lengths:
            for i in range(len(protein_seq) - length + 1):
                peptide = protein_seq[i:i + length]

                for allele in alleles:
                    pred = predictor.predict(peptides=[peptide], alleles=[allele])
                    percentile = pred['mhcflurry_affinity_percentile'].values[0]

                    if percentile < percentile_threshold:
                        epitopes.append({
                            'peptide': peptide,
                            'position': i + 1,
                            'length': length,
                            'allele': allele,
                            'affinity_nM': pred['mhcflurry_affinity'].values[0],
                            'percentile': percentile,
                            'binding': classify_binding(pred['mhcflurry_affinity'].values[0])
                        })

        return pd.DataFrame(epitopes).sort_values('percentile')
    except ImportError:
        print('MHCflurry not installed')
        return None


# Common HLA alleles for population-level analysis
COMMON_ALLELES = {
    'HLA-A': ['HLA-A*02:01', 'HLA-A*01:01', 'HLA-A*03:01', 'HLA-A*24:02', 'HLA-A*11:01'],
    'HLA-B': ['HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*44:02', 'HLA-B*15:01', 'HLA-B*35:01'],
}


if __name__ == '__main__':
    print('MHC Binding Prediction Example')
    print('=' * 50)

    # Test peptides (known epitopes)
    peptides = [
        'SIINFEKL',   # Ovalbumin, strong HLA-A*02:01 binder (mouse)
        'GILGFVFTL',  # Influenza M1, HLA-A*02:01 epitope
        'NLVPMVATV',  # CMV pp65, HLA-A*02:01 epitope
    ]
    alleles = ['HLA-A*02:01', 'HLA-A*03:01']

    print('\nPredicting binding for known epitopes...')
    results = predict_mhc_binding(peptides, alleles)

    if results is not None:
        results['binding'] = results['affinity_nM'].apply(classify_binding)
        print(results.to_string(index=False))

        # Protein scan example
        print('\n' + '=' * 50)
        print('Scanning short protein for epitopes...')
        test_protein = 'MSIINFEKLAAAGIAGFVFTLVSSAYNLVPMVATVQTLNF'

        epitopes = scan_protein(test_protein, ['HLA-A*02:01'], lengths=[9])
        if epitopes is not None and len(epitopes) > 0:
            print(f'\nFound {len(epitopes)} potential epitopes:')
            print(epitopes[['peptide', 'position', 'affinity_nM', 'binding']].head(10).to_string(index=False))
    else:
        # Simulated results if MHCflurry not available
        print('\nSimulated results (MHCflurry not installed):')
        for pep in peptides:
            for allele in alleles:
                print(f'  {pep} + {allele}: ~100-500 nM (moderate binder)')
