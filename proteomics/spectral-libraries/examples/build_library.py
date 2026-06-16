'''Calibrate predicted iRT, merge libraries on the full transition key, and QC.

Operates on a small in-memory library so it runs without network or input files.
Koina, MS2PIP, and DeepLC prediction calls are shown in comments for reference;
swap the in-memory table for real predictions in production.
'''
# Reference: pandas 2.2+, numpy 1.26+, scipy 1.12+ | Verify API if version differs
import numpy as np
import pandas as pd
from scipy import stats

# --- Predicted-library generation (reference; needs network/local model) ---
# from koinapy import Koina  # verify constructor signature with help(Koina)
# inputs = pd.DataFrame({'peptide_sequences': ['LGGNEQVTR'], 'precursor_charges': [2], 'collision_energies': [30]})
# fragments = Koina('Prosit_2019_intensity', 'koina.wilhelmlab.org:443').predict(inputs)
# import ms2pip; result = ms2pip.predict_batch(psms, model='HCD')  # v4 module-level API, ProcessingResult objects
# from deeplc import DeepLC; dlc = DeepLC(); dlc.calibrate_preds(seq_df=cal_df); rt = dlc.make_preds(seq_df=pep_df)

R2_MIN = 0.95  # below this the iRT-to-RT fit is untrustworthy and extraction windows misplace
FRAGMENTS_PER_PRECURSOR = 6  # confident peak-group scoring without inviting interference
TRANSITION_KEY = ['ModifiedSequence', 'PrecursorCharge', 'FragmentType', 'FragmentSeriesNumber', 'FragmentCharge']

IRT_PEPTIDES = {'LGGNEQVTR': -24.92, 'GAGSSEPVTGLDAK': 0.00, 'VEATFGVDESNAK': 12.39,
                'YILAGVENSK': 19.79, 'TPVISGGPYEYR': 28.71, 'TPVITGAPYEYR': 33.38,
                'DGLDAASYYAPVR': 42.26, 'ADVTPADFSEWSK': 54.62, 'GTFIIDPGGVIR': 70.52,
                'GTFIIDPAAVIR': 87.23, 'LFLQFGAQGSPFLK': 100.00}

def calibrate_irt(anchor_irt, observed_rt):
    slope, intercept, r, _, _ = stats.linregress(anchor_irt, observed_rt)
    if r ** 2 < R2_MIN:
        raise ValueError(f'iRT fit R^2={r**2:.3f} < {R2_MIN}; gradient may be nonlinear, use LOWESS')
    return slope, intercept, r ** 2

def merge_libraries(libs):
    combined = pd.concat(libs, ignore_index=True)
    combined['precursor_total'] = combined.groupby(['ModifiedSequence', 'PrecursorCharge'])['LibraryIntensity'].transform('sum')
    combined = combined.sort_values('precursor_total', ascending=False)
    return combined.drop_duplicates(subset=TRANSITION_KEY).drop(columns='precursor_total')

def library_stats(lib):
    n_prec = lib.groupby(['ModifiedSequence', 'PrecursorCharge']).ngroups
    return {'precursors': n_prec, 'proteins': lib['ProteinId'].nunique(),
            'transitions_per_precursor': round(len(lib) / n_prec, 1)}

anchors = pd.DataFrame({'peptide': list(IRT_PEPTIDES), 'irt': list(IRT_PEPTIDES.values())})
anchors['observed_rt'] = 5.0 + 0.18 * anchors['irt'] + np.random.default_rng(0).normal(0, 0.1, len(anchors))
slope, intercept, r2 = calibrate_irt(anchors['irt'], anchors['observed_rt'])
print(f'iRT calibration: RT = {slope:.3f} * iRT + {intercept:.3f}, R^2 = {r2:.3f}')

rng = np.random.default_rng(1)
def make_lib(seq, charge, protein, n):
    return pd.DataFrame({'ModifiedSequence': seq, 'PrecursorCharge': charge, 'ProteinId': protein,
                         'FragmentType': ['y'] * n, 'FragmentSeriesNumber': np.arange(1, n + 1),
                         'FragmentCharge': 1, 'LibraryIntensity': rng.uniform(0.1, 1.0, n)})

lib1 = make_lib('LGGNEQVTR', 2, 'P1', FRAGMENTS_PER_PRECURSOR)
lib2 = pd.concat([make_lib('LGGNEQVTR', 3, 'P1', FRAGMENTS_PER_PRECURSOR),
                  make_lib('VEATFGVDESNAK', 2, 'P2', FRAGMENTS_PER_PRECURSOR)], ignore_index=True)

merged = merge_libraries([lib1, lib2])
print('Merged library stats:', library_stats(merged))
print(f'Charge-2 and charge-3 of LGGNEQVTR both retained: '
      f'{set(merged.loc[merged.ModifiedSequence == "LGGNEQVTR", "PrecursorCharge"]) == {2, 3}}')
