'''Prepare a subject-level analysis dataset from CDISC SDTM domain files'''
# Reference: pyreadstat 1.2+, pandas 2.1+ | Verify API if version differs

import pandas as pd
import numpy as np

dm = pd.read_csv('DM.csv')
ae = pd.read_csv('AE.csv')

print(f'DM: {len(dm)} subjects, columns: {list(dm.columns)}')
print(f'AE: {len(ae)} events across {ae["USUBJID"].nunique()} subjects')

assert dm['USUBJID'].is_unique, 'Duplicate USUBJIDs in DM'

severity_map = {'MILD': 1, 'MODERATE': 2, 'SEVERE': 3, 'LIFE THREATENING': 4, 'FATAL': 5}
ae['AESEV_NUM'] = ae['AESEV'].map(severity_map).fillna(0)

ae_summary = ae.groupby('USUBJID').agg(
    ae_count=('AETERM', 'count'),
    had_serious=('AESER', lambda x: (x == 'Y').any()),
    max_severity=('AESEV_NUM', 'max')
).reset_index()

analysis = dm[['USUBJID', 'ARM', 'ARMCD', 'AGE', 'SEX', 'RACE']].copy()
analysis = analysis.merge(ae_summary, on='USUBJID', how='left')

analysis['ae_count'] = analysis['ae_count'].fillna(0).astype(int)
analysis['had_serious'] = analysis['had_serious'].fillna(False)
analysis['max_severity'] = analysis['max_severity'].fillna(0).astype(int)

analysis['had_any_ae'] = (analysis['ae_count'] > 0).astype(int)

print(f'\nAnalysis dataset: {len(analysis)} subjects, {len(analysis.columns)} columns')
print(f'\nSubjects by treatment arm:')
print(analysis['ARM'].value_counts())

print(f'\nAdverse event summary by arm:')
arm_summary = analysis.groupby('ARM').agg(
    n_subjects=('USUBJID', 'count'),
    pct_any_ae=('had_any_ae', 'mean'),
    pct_serious=('had_serious', 'mean'),
    mean_ae_count=('ae_count', 'mean')
)
arm_summary['pct_any_ae'] = (arm_summary['pct_any_ae'] * 100).round(1)
arm_summary['pct_serious'] = (arm_summary['pct_serious'] * 100).round(1)
arm_summary['mean_ae_count'] = arm_summary['mean_ae_count'].round(1)
print(arm_summary)

print(f'\nAge distribution:')
print(analysis.groupby('ARM')['AGE'].describe().round(1))

analysis.to_csv('analysis_dataset.csv', index=False)
print(f'\nSaved analysis_dataset.csv ({len(analysis)} rows, {len(analysis.columns)} columns)')
