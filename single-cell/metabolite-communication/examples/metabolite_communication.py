'''Analyze metabolite-mediated cell-cell communication with MEBOCOST'''
# Reference: mebocost 1.0+, scanpy 1.10+ | Verify API if version differs

import scanpy as sc
import pandas as pd
import numpy as np


def prepare_data_for_mebocost(adata, cell_type_col='cell_type', min_cells=50):
    '''Prepare AnnData for MEBOCOST.

    Requirements: log-normalized expression, gene SYMBOLS (not Ensembl IDs),
    and enough cells per type. min_cells=50 keeps per-group statistics stable;
    fewer cells give unreliable permutation FDRs.
    '''
    if adata.X.max() > 50:
        print('Data appears not log-normalized; normalizing')
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    counts = adata.obs[cell_type_col].value_counts()
    valid = counts[counts >= min_cells].index.tolist()
    if len(valid) < len(counts):
        print(f'Removing rare cell types (<{min_cells} cells): {set(counts.index) - set(valid)}')
    adata = adata[adata.obs[cell_type_col].isin(valid)].copy()
    print(f'Final dataset: {adata.n_obs} cells, {len(valid)} cell types')
    return adata


def run_mebocost_analysis(adata, cell_type_col='cell_type', species='human',
                          config_path='./mebocost.conf', n_shuffle=1000):
    '''Run MEBOCOST metabolite communication inference.

    Flow: Sender -> enzyme (synthesis) -> Metabolite -> Sensor -> Receiver.
    config_path points to mebocost.conf listing the metabolite-enzyme-sensor database.
    cutoff_prop=0.15 is the dropout floor; n_shuffle builds the permutation FDR.
    '''
    from mebocost import mebocost

    mebo = mebocost.create_obj(adata=adata, group_col=cell_type_col, condition_col=None,
                               met_est='mebocost', config_path=config_path, species=species,
                               cutoff_exp='auto', cutoff_met='auto', cutoff_prop=0.15,
                               sensor_type='All', thread=8)
    commu_res = mebo.infer_commu(n_shuffle=n_shuffle, seed=12345, Return=True,
                                 min_cell_number=10, pval_method='permutation_test_fdr',
                                 pval_cutoff=0.05, thread=None)
    return mebo, commu_res


def summarize_results(commu_res, fdr_threshold=0.05):
    '''Summarize significant communications. Filter on permutation_test_fdr, not raw p-value.

    Result columns are capitalized: Sender, Receiver, Metabolite_Name, Sensor,
    Annotation (Transporter/Enzyme), Commu_Score, Norm_Commu_Score, permutation_test_fdr.
    '''
    sig = commu_res[commu_res['permutation_test_fdr'] < fdr_threshold].copy()
    print(f'\nMetabolite communication (FDR < {fdr_threshold})')
    print(f'Total tested: {len(commu_res)}  Significant: {len(sig)}')
    if len(sig) == 0:
        return sig

    print(f"Unique metabolites: {sig['Metabolite_Name'].nunique()}")
    print('\nTop metabolites:')
    for met, n in sig['Metabolite_Name'].value_counts().head(10).items():
        print(f'  {met}: {n}')

    sig['pair'] = sig['Sender'] + ' -> ' + sig['Receiver']
    print('\nTop sender -> receiver pairs:')
    for pair, n in sig['pair'].value_counts().head(5).items():
        print(f'  {pair}: {n}')

    transporter = sig[sig['Annotation'] == 'Transporter']
    print(f'\nTransporter-based (bidirectional, lower-confidence) calls: {len(transporter)}')
    return sig


def analyze_specific_metabolite(commu_res, metabolite, fdr_threshold=0.05):
    '''Detail communications for one metabolite. Output is a hypothesis: the metabolite
    is inferred from enzyme expression, never measured.'''
    sig = commu_res[(commu_res['Metabolite_Name'] == metabolite) &
                    (commu_res['permutation_test_fdr'] < fdr_threshold)]
    if len(sig) == 0:
        print(f'No significant {metabolite} communication')
        return None
    print(f'\n{metabolite}: machinery consistent with sender -> receiver flow')
    for _, row in sig.iterrows():
        print(f"  {row['Sender']} -> {row['Receiver']} via {row['Sensor']} "
              f"(score {row['Commu_Score']:.3f}, FDR {row['permutation_test_fdr']:.4f})")
    return sig


if __name__ == '__main__':
    print('MEBOCOST metabolite communication workflow')
    print('Run on real log-normalized scRNA-seq with gene symbols and a mebocost.conf:')
    print('  adata = prepare_data_for_mebocost(adata)')
    print('  mebo, commu_res = run_mebocost_analysis(adata, config_path="./mebocost.conf")')
    print('  sig = summarize_results(commu_res)')
    print('Validate every hit with metabolomics / MSI / tracing before claiming production.')
