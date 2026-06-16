'''Expand MaxQuant Phospho (STY)Sites multiplicity and filter to class I localization.

Self-contained: writes a tiny synthetic site table to a tempdir, processes it, cleans up.
The point is the multiplicity expansion (Intensity___1/___2/___3, THREE underscores) and
the protein-adjustment caveat - quantifying on the collapsed base Intensity mixes
phospho-states and can fake dephosphorylation.
'''
# Reference: pandas 2.2+, numpy 1.26+ | Verify API if version differs
import os
import tempfile
import pandas as pd
import numpy as np

CLASS_I_PROB = 0.75   # Olsen 2006 class-I convention; comparability standard, not a calibrated FLR

synthetic = pd.DataFrame({
    'Gene names': ['AKT1;AKT1b', 'MAPK1', np.nan, 'GSK3B'],
    'Protein': ['P31749', 'P28482', 'Q9Y6K9', 'P49841'],
    'Amino acid': ['S', 'T', 'Y', 'S'],
    'Position': [473, 185, 99, 9],
    'Localization prob': [0.99, 0.82, 0.61, 0.95],   # last-but-one is class II, filtered out
    'Sequence window': ['_' * 31, '_' * 31, '_' * 31, '_' * 31],
    'Reverse': ['', '', '', ''],
    'Potential contaminant': ['', '', '', ''],
    'Intensity___1': [1.0e8, 5.0e7, 2.0e7, 8.0e7],
    'Intensity___2': [3.0e7, 0.0, 0.0, 2.0e7],
    'Intensity___3': [0.0, 0.0, 0.0, 0.0],
})

with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, 'Phospho (STY)Sites.txt')
    synthetic.to_csv(path, sep='\t', index=False)

    phospho = pd.read_csv(path, sep='\t', low_memory=False)
    print(f'Total sites: {len(phospho)}')

    contaminant_col = 'Potential contaminant' if 'Potential contaminant' in phospho.columns else 'Contaminant'
    phospho = phospho[(phospho['Reverse'] != '+') & (phospho[contaminant_col] != '+')]

    phospho['loc_class'] = pd.cut(phospho['Localization prob'], bins=[0, 0.25, 0.5, 0.75, 1.0], labels=['IV', 'III', 'II', 'I'], include_lowest=True)
    print('\nLocalization class distribution:')
    print(phospho['loc_class'].value_counts())

    confident = phospho[phospho['Localization prob'] >= CLASS_I_PROB].copy()
    print(f'\nClass I sites (prob >= {CLASS_I_PROB}): {len(confident)}')

    gene = confident['Gene names'].where(confident['Gene names'].notna(), confident['Protein'])
    confident['site_id'] = gene.str.split(';').str[0] + '_' + confident['Amino acid'] + confident['Position'].astype(int).astype(str)

    mult_cols = [c for c in confident.columns if c.startswith('Intensity') and '___' in c and c.split('___')[-1] in {'1', '2', '3'}]
    long = confident.melt(id_vars=['site_id', 'Amino acid', 'Position', 'Localization prob'], value_vars=mult_cols, var_name='run_multiplicity', value_name='intensity')
    long['multiplicity'] = long['run_multiplicity'].str.split('___').str[-1]
    long = long[long['intensity'] > 0].copy()
    long['log2_intensity'] = np.log2(long['intensity'])

    print('\nMultiplicity-resolved site observations (collapsing these would fake dephosphorylation):')
    print(long[['site_id', 'multiplicity', 'log2_intensity']].to_string(index=False))

    print('\nNext step (not shown): feed the PTM site table AND a paired global proteome to')
    print('MSstatsPTM groupComparisonPTM and call only ADJUSTED.Model hits regulated.')
