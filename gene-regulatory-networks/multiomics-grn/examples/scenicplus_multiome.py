'''
SCENIC+ multiome workflow for enhancer-driven GRN inference.

Requires:
- 10x Multiome output from CellRanger ARC (filtered_feature_bc_matrix.h5, atac_fragments.tsv.gz)
- hg38-blacklist.v2.bed: ENCODE blacklist regions
- cisTarget databases and TF list
'''

import subprocess
import scanpy as sc
import pandas as pd
import numpy as np
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models, evaluate_models


def prepare_rna(matrix_h5):
    '''Load and preprocess scRNA-seq from CellRanger ARC.'''
    adata = sc.read_10x_h5(matrix_h5, gex_only=True)
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    # MT% < 20: standard cutoff for most tissues
    adata = adata[adata.obs.pct_counts_mt < 20].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.8)
    sc.tl.umap(adata)

    print(f'RNA: {adata.n_obs} cells, {adata.n_vars} genes')
    return adata


def call_peaks_macs3(fragments_file, output_dir='macs3_output'):
    '''Call peaks from ATAC fragments with MACS3.'''
    subprocess.run([
        'macs3', 'callpeak',
        '-t', fragments_file,
        '-f', 'BEDPE',
        '--nomodel', '--shift', '-75', '--extsize', '150',
        '-g', 'hs',
        '--keep-dup', 'all',
        '-n', 'multiome_peaks',
        '--outdir', output_dir
    ], check=True)
    print(f'Peaks called in {output_dir}/')
    return f'{output_dir}/multiome_peaks_peaks.narrowPeak'


def run_cistopic(fragments_file, peaks_file, blacklist_file, n_cpu=8):
    '''Run cisTopic LDA topic modeling on scATAC data.'''
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments=fragments_file,
        path_to_regions=peaks_file,
        path_to_blacklist=blacklist_file,
        n_cpu=n_cpu
    )

    # Test multiple topic numbers; select by coherence/perplexity
    models = run_cgs_models(cistopic_obj, n_topics=[10, 20, 30, 40, 50],
                            n_cpu=n_cpu, n_iter=300, random_state=42)
    model = evaluate_models(models, return_model=True)
    cistopic_obj.add_LDA_model(model)

    print(f'cisTopic: {model.n_topic} topics selected')
    return cistopic_obj


def run_scenicplus(adata_rna, cistopic_obj, tf_file, save_path='scenicplus_output/'):
    '''Run SCENIC+ eRegulon inference.'''
    from scenicplus.scenicplus_class import SCENICPLUS
    from scenicplus.wrappers.run_scenicplus import run_scenicplus as _run

    scplus_obj = SCENICPLUS(adata_rna, cistopic_obj, menr=None)

    _run(
        scplus_obj,
        variable=['GeneExpressionLevel'],
        species='hsapiens',
        assembly='hg38',
        tf_file=tf_file,
        save_path=save_path,
        biomart_host='http://www.ensembl.org',
        upstream=[1000, 150000],
        downstream=[1000, 150000],
        calculate_TF_eGRN_correlation=True,
        calculate_DEGs_DARs=True,
        export_to_loom_file=True,
        export_to_UCSC_file=True,
        n_cpu=8
    )

    eregulons = scplus_obj.uns['eRegulon_metadata']
    print(f'Found {eregulons["Region_signature_name"].nunique()} eRegulons')

    ereg_summary = eregulons.groupby('TF').agg(
        n_regions=('Region', 'nunique'),
        n_genes=('Gene', 'nunique')
    ).sort_values('n_genes', ascending=False)
    print(ereg_summary.head(20))

    return scplus_obj


if __name__ == '__main__':
    # adata_rna = prepare_rna('filtered_feature_bc_matrix.h5')
    # peaks_file = call_peaks_macs3('atac_fragments.tsv.gz')
    # cistopic_obj = run_cistopic('atac_fragments.tsv.gz', peaks_file, 'hg38-blacklist.v2.bed')
    # scplus_obj = run_scenicplus(adata_rna, cistopic_obj, 'allTFs_hg38.txt')

    print('Uncomment sections above to run with actual 10x multiome data')
