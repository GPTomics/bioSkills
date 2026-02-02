'''Analyze metabolite-mediated cell-cell communication'''

import scanpy as sc
import pandas as pd
import numpy as np


def prepare_data_for_mebocost(adata, cell_type_col='cell_type', min_cells=50):
    '''Prepare AnnData for MeboCost analysis

    Requirements:
    - Log-normalized expression
    - Cell type annotations with sufficient cells per type
    - Gene symbols (not Ensembl IDs)

    Args:
        min_cells: Minimum cells per cell type (50 typical)
                  Fewer cells = unreliable statistics
    '''
    # Verify log normalization
    if adata.X.max() > 50:
        print('Data appears not log-normalized. Normalizing...')
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Filter rare cell types
    cell_counts = adata.obs[cell_type_col].value_counts()
    valid_types = cell_counts[cell_counts >= min_cells].index.tolist()

    if len(valid_types) < len(cell_counts):
        removed = set(cell_counts.index) - set(valid_types)
        print(f'Removing rare cell types (<{min_cells} cells): {removed}')

    adata = adata[adata.obs[cell_type_col].isin(valid_types)].copy()

    print(f'Final dataset: {adata.n_obs} cells, {len(valid_types)} cell types')
    return adata


def run_mebocost_analysis(adata, cell_type_col='cell_type', species='human', n_perms=1000):
    '''Run MeboCost metabolite communication analysis

    MeboCost infers communication via:
    Sender cell -> Enzyme (secretes) -> Metabolite -> Receptor -> Receiver cell

    Args:
        n_perms: Permutations for p-value calculation
                1000: Standard for publication
                100: Quick exploration
                0: No significance testing (fastest)
    '''
    try:
        import mebocost as mbc
    except ImportError:
        print('MeboCost not installed. Install with: pip install mebocost')
        return None

    mebo = mbc.create_obj(
        adata=adata,
        group_col=cell_type_col,
        species=species
    )

    mebo.infer_commu(n_permutations=n_perms, seed=42)

    return mebo


def summarize_results(mebo, pval_threshold=0.05):
    '''Summarize metabolite communication results

    Significance threshold:
    - 0.05: Standard (5% FDR)
    - 0.01: Stringent
    - 0.1: Exploratory
    '''
    results = mebo.commu_res.copy()
    sig = results[results['pval'] < pval_threshold]

    print(f'\nMetabolite Communication Summary (p < {pval_threshold})')
    print('=' * 50)
    print(f"Total interactions tested: {len(results)}")
    print(f"Significant interactions: {len(sig)} ({len(sig)/len(results):.1%})")

    if len(sig) > 0:
        print(f"\nUnique metabolites involved: {sig['metabolite'].nunique()}")
        print(f"Unique sender cell types: {sig['sender'].nunique()}")
        print(f"Unique receiver cell types: {sig['receiver'].nunique()}")

        print('\nTop 10 metabolites by frequency:')
        for met, count in sig['metabolite'].value_counts().head(10).items():
            print(f'  {met}: {count} interactions')

        print('\nTop sender-receiver pairs:')
        sig['pair'] = sig['sender'] + ' -> ' + sig['receiver']
        for pair, count in sig['pair'].value_counts().head(5).items():
            print(f'  {pair}: {count} interactions')

    return sig


def analyze_specific_metabolite(mebo, metabolite, pval_threshold=0.05):
    '''Detailed analysis for a specific metabolite'''
    sig = mebo.commu_res[
        (mebo.commu_res['metabolite'] == metabolite) &
        (mebo.commu_res['pval'] < pval_threshold)
    ]

    if len(sig) == 0:
        print(f'No significant {metabolite} communication found')
        return None

    print(f'\n{metabolite} Communication Flow')
    print('=' * 50)

    # Group by sender
    for sender in sig['sender'].unique():
        sender_data = sig[sig['sender'] == sender]
        receivers = sender_data['receiver'].unique()
        print(f'\n{sender} secretes {metabolite} to:')
        for _, row in sender_data.iterrows():
            print(f"  -> {row['receiver']} via {row['receptor']} "
                  f"(score: {row['commu_score']:.3f}, p={row['pval']:.4f})")

    return sig


# Simulated example (MeboCost requires real data)
if __name__ == '__main__':
    print('MeboCost Metabolite Communication Analysis')
    print('=' * 50)

    # Simulate a small dataset for demonstration
    # In practice, use real scRNA-seq data
    np.random.seed(42)
    n_cells = 500
    n_genes = 100

    # Create mock AnnData
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['T_cell', 'Macrophage', 'Fibroblast', 'Epithelial'],
                                      size=n_cells)
    })
    var = pd.DataFrame(index=[f'Gene{i}' for i in range(n_genes)])

    adata = sc.AnnData(X=X.astype(float), obs=obs, var=var)

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print('\nDataset prepared:')
    print(f'  Cells: {adata.n_obs}')
    print(f'  Genes: {adata.n_vars}')
    print(f'  Cell types: {adata.obs["cell_type"].nunique()}')

    # Note: Full MeboCost analysis requires the package and real gene names
    # This is a demonstration of the workflow structure
    print('\nTo run full analysis:')
    print('  1. Install mebocost: pip install mebocost')
    print('  2. Use real scRNA-seq data with gene symbols')
    print('  3. Call run_mebocost_analysis(adata)')
