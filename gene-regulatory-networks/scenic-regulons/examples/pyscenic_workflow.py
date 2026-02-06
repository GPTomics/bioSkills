'''
Complete pySCENIC three-step pipeline: GRNBoost2, cisTarget, AUCell.

Requires:
- filtered.loom: scRNA-seq expression in loom format
- allTFs_hg38.txt: TF list for human
- hg38 ranking databases (.feather files)
- motifs-v9-nr.hgnc-m0.001-o0.0.tbl: motif annotations
'''

import os
import glob
import pickle
import subprocess
import pandas as pd
import numpy as np
import loompy
from pyscenic.utils import load_tf_names
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.rss import regulon_specificity_scores
from ctxcore.rnkdb import FeatherRankingDatabase

LOOM_FILE = 'filtered.loom'
TF_FILE = 'allTFs_hg38.txt'
MOTIF_ANNOTATIONS = 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
DB_GLOB = '*.genes_vs_motifs.rankings.feather'
NUM_WORKERS = 8


def step1_grn_inference(loom_file, tf_file, output='adj.tsv', num_workers=8):
    '''
    GRN inference with GRNBoost2.
    Uses arboreto_with_multiprocessing.py to avoid dask >= 2.0 compatibility issues.
    '''
    # Locate the bundled script
    import pyscenic
    scenic_dir = os.path.dirname(pyscenic.__file__)
    arboreto_script = os.path.join(scenic_dir, 'arboreto_with_multiprocessing.py')

    if not os.path.exists(arboreto_script):
        arboreto_script = 'arboreto_with_multiprocessing.py'

    subprocess.run([
        'python', arboreto_script,
        loom_file, tf_file,
        '--method', 'grnboost2',
        '--output', output,
        '--num_workers', str(num_workers),
        '--seed', '42'
    ], check=True)
    print(f'Step 1 complete: {output}')
    return output


def step2_cistopic_pruning(adj_file, db_glob, motif_annotations, loom_file):
    '''Prune co-expression modules by cis-regulatory motif enrichment.'''
    db_fnames = glob.glob(db_glob)
    dbs = [FeatherRankingDatabase(fname) for fname in db_fnames]
    print(f'Loaded {len(dbs)} ranking databases')

    adjacencies = pd.read_csv(adj_file, sep='\t')
    print(f'Loaded {len(adjacencies)} TF-target adjacencies')

    df = prune2df(dbs, adjacencies, motif_annotations)
    regulons = df2regulons(df)

    with open('regulons.pkl', 'wb') as f:
        pickle.dump(regulons, f)

    print(f'Step 2 complete: {len(regulons)} regulons')
    for reg in sorted(regulons, key=lambda r: -len(r))[:10]:
        print(f'  {reg.name}: {len(reg)} targets')

    return regulons


def step3_aucell_scoring(loom_file, regulons, auc_threshold=0.05, num_workers=8):
    '''Score regulon activity per cell with AUCell.'''
    ds = loompy.connect(loom_file)
    expr_matrix = pd.DataFrame(ds[:, :], index=ds.ra.Gene, columns=ds.ca.CellID).T
    ds.close()

    # auc_threshold 0.05: consider top 5% of ranked genes per cell
    auc_mtx = aucell(expr_matrix, regulons, auc_threshold=auc_threshold, num_workers=num_workers)
    auc_mtx.to_csv('auc_matrix.csv')

    print(f'Step 3 complete: scored {auc_mtx.shape[1]} regulons across {auc_mtx.shape[0]} cells')
    return auc_mtx


def compute_rss(auc_mtx, cell_type_series):
    '''Regulon specificity scores per cell type.'''
    rss = regulon_specificity_scores(auc_mtx, cell_type_series)
    rss.to_csv('rss_scores.csv')

    for ct in rss.columns:
        top = rss[ct].sort_values(ascending=False).head(5)
        print(f'\n{ct}:')
        for reg, score in top.items():
            print(f'  {reg}: {score:.3f}')

    return rss


if __name__ == '__main__':
    # Step 1: GRN inference
    # adj_file = step1_grn_inference(LOOM_FILE, TF_FILE, num_workers=NUM_WORKERS)

    # Step 2: Regulon pruning
    # regulons = step2_cistopic_pruning('adj.tsv', DB_GLOB, MOTIF_ANNOTATIONS, LOOM_FILE)

    # Step 3: AUCell scoring
    # with open('regulons.pkl', 'rb') as f:
    #     regulons = pickle.load(f)
    # auc_mtx = step3_aucell_scoring(LOOM_FILE, regulons, num_workers=NUM_WORKERS)

    # RSS analysis (requires cell type labels)
    # cell_types = pd.read_csv('cell_types.csv', index_col=0)['cell_type']
    # rss = compute_rss(auc_mtx, cell_types)

    print('Uncomment sections above to run with actual data')
