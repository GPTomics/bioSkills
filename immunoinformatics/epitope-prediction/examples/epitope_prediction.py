'''Epitope prediction with the right confidence per epitope class.

The load-bearing decision is which problem is being solved. T-cell epitope
prediction is mature (it reduces to MHC presentation); B-cell linear prediction is
unreliable (~AUC 0.6) because ~90% of natural epitopes are conformational, so the
defensible B-cell path is structure-based DiscoTope-3.0 on an AlphaFold model,
gated by pLDDT. This script routes a request to the appropriate method, classifies
BepiPred-3.0 output at its real default threshold (0.1512, not 0.5), and filters
DiscoTope calls by structure confidence.
'''
# Reference: BepiPred-3.0, pandas 2.2+ | Verify API if version differs

import pandas as pd

BEPIPRED3_DEFAULT_THRESHOLD = 0.1512  # Clifford 2022; balances sens/spec - NOT 0.5


def route_epitope_request(target, has_structure, native_response):
    '''Pick a defensible method. native_response=True means predicting the antibody
    response against folded antigen (the conformational case linear models cannot do).'''
    if target == 't_cell_cd8':
        return 'NetMHCpan-4.1 EL / MHCflurry (see mhc-binding-prediction); skip NetChop by default'
    if target == 't_cell_cd4':
        return 'NetMHCIIpan-4.3 (see mhc-class-ii-prediction); less reliable than class I'
    if has_structure or native_response:
        return 'DiscoTope-3.0 on the (folded) structure; gate by pLDDT'
    return 'BepiPred-3.0 linear - only for peptide/denatured targets; misses ~90% native epitopes'


def call_bepipred_regions(per_residue, threshold=BEPIPRED3_DEFAULT_THRESHOLD, min_len=5):
    '''Collapse per-residue BepiPred-3.0 probabilities into contiguous epitope regions.
    per_residue: DataFrame with columns Position, Residue, Score.'''
    df = per_residue.copy()
    df['is_epitope'] = df['Score'] >= threshold
    regions, run = [], []
    for _, r in df.iterrows():
        if r['is_epitope']:
            run.append(r)
        elif run:
            if len(run) >= min_len:
                regions.append({'start': run[0]['Position'], 'end': run[-1]['Position'],
                                'sequence': ''.join(x['Residue'] for x in run),
                                'mean_score': sum(x['Score'] for x in run) / len(run)})
            run = []
    return pd.DataFrame(regions)


def gate_discotope_by_plddt(df, plddt_col='pLDDT', score_col='DiscoTope-3.0 score', min_plddt=70):
    '''DiscoTope-3.0 AUC-PR is only ~0.22 (many false positives) and accuracy drops
    in low-pLDDT loops; keep calls only in confidently-folded regions.'''
    return df[df[plddt_col] >= min_plddt].sort_values(score_col, ascending=False)


if __name__ == '__main__':
    route_epitope_request('t_cell_cd8', has_structure=False, native_response=False)
    route_epitope_request('b_cell', has_structure=False, native_response=True)

    demo = pd.DataFrame({
        'Position': range(1, 13),
        'Residue': list('MKTAYIAKQRQI'),
        'Score': [0.05, 0.20, 0.31, 0.18, 0.40, 0.22, 0.08, 0.05, 0.16, 0.19, 0.25, 0.30],
    })
    call_bepipred_regions(demo)
