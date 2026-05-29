'''
SCENIC+ multiome workflow for enhancer-driven (eRegulon) GRN inference.

Uses the CURRENT Snakemake-based SCENIC+ pipeline. The legacy manual SCENICPLUS-object
API (run_scenicplus wrapper) is deprecated -- pre-2024 tutorials are stale. The Python
here covers data prep (RNA QC, consensus peak calling); SCENIC+ itself is driven by its
Snakemake workflow. Verify command names against scenicplus.readthedocs.io at run time.

Requires:
- 10x Multiome output from CellRanger ARC (filtered_feature_bc_matrix.h5, atac_fragments.tsv.gz)
- A cell-type annotation (consensus peaks are called from per-cell-type pseudobulk)
- cisTarget databases (region-based) and a TF list
'''
# Reference: scenicplus (Snakemake workflow), pycisTopic 2.0+, scanpy 1.10+, macs3 3.0+ | Verify API if version differs

import glob
import subprocess
import scanpy as sc
import pandas as pd


def prepare_rna(matrix_h5):
    '''Load and preprocess scRNA-seq from CellRanger ARC; annotate cell types before peak calling.'''
    adata = sc.read_10x_h5(matrix_h5, gex_only=True)
    adata.var_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 20].copy()      # MT% < 20: standard tissue cutoff

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
    '''Call peaks from ATAC fragments with MACS3 for the region universe.'''
    subprocess.run([
        'macs3', 'callpeak',
        '-t', fragments_file, '-f', 'BEDPE',
        '--nomodel', '--shift', '-75', '--extsize', '150',
        '-g', 'hs', '--keep-dup', 'all',
        '-n', 'multiome_peaks', '--outdir', output_dir
    ], check=True)
    print(f'Peaks called in {output_dir}/')
    return f'{output_dir}/multiome_peaks_peaks.narrowPeak'


def run_scenicplus_snakemake(out_dir='scenicplus_run', cores=16):
    '''Scaffold and run the SCENIC+ Snakemake pipeline.

    Edit out_dir/config/config.yaml between scaffolding and running to point at the
    scATAC fragments, the scRNA AnnData, the cell-type annotation (for pseudobulk peak
    calling), the cisTarget databases, and the output paths. The pipeline performs topic
    modeling, motif enrichment, region-to-gene and TF-to-gene inference, and eRegulon
    assembly -- the steps the deprecated run_scenicplus wrapper used to bundle.
    '''
    subprocess.run(['scenicplus', 'init_snakemake', '--out_dir', out_dir], check=True)
    print(f'Edit {out_dir}/Snakemake/config/config.yaml, then run the pipeline.')
    # Run from inside the Snakemake dir so config-relative paths resolve.
    subprocess.run(['snakemake', '--cores', str(cores)],
                   cwd=f'{out_dir}/Snakemake', check=True)
    # Output name/dir are set in config.yaml (output_data); the exact filename spelling has
    # varied across versions, so resolve the direct (high-confidence) eRegulon table by glob.
    return glob.glob(f'{out_dir}/**/eRegulon*direct*.tsv', recursive=True)[0]


def summarize_eregulons(eregulon_tsv):
    '''Summarize eRegulons (TF -> region -> gene triplets) by target counts per TF.'''
    eregulons = pd.read_csv(eregulon_tsv, sep='\t')
    summary = (eregulons.groupby('TF')
               .agg(n_regions=('Region', 'nunique'), n_genes=('Gene', 'nunique'))
               .sort_values('n_genes', ascending=False))
    print(summary.head(20))
    return summary


if __name__ == '__main__':
    # adata_rna = prepare_rna('filtered_feature_bc_matrix.h5')   # annotate cell types before peaks
    # peaks_file = call_peaks_macs3('atac_fragments.tsv.gz')
    # ereg_tsv = run_scenicplus_snakemake('scenicplus_run')
    # summarize_eregulons(ereg_tsv)

    print('Uncomment sections above to run with actual 10x multiome data')
