#!/usr/bin/env python3
'''
cfDNA methylation analysis for cancer detection.
'''

import subprocess
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import nnls


def extract_methylation(bam_file, reference, output_prefix, min_depth=5):
    '''Extract methylation from bisulfite BAM using MethylDackel.'''
    subprocess.run([
        'MethylDackel', 'extract',
        reference,
        bam_file,
        '-o', output_prefix,
        '--minDepth', str(min_depth),
        '--mergeContext'
    ], check=True)

    bedgraph = f'{output_prefix}_CpG.bedGraph'
    meth = pd.read_csv(bedgraph, sep='\t', header=None,
                       names=['chrom', 'start', 'end', 'meth_pct', 'meth', 'unmeth'])
    meth['beta'] = meth['meth'] / (meth['meth'] + meth['unmeth'])

    return meth


def find_dmrs(cancer_samples, normal_samples, min_diff=0.2, min_sites=5):
    '''Find differentially methylated regions.'''
    all_cancer = pd.concat(cancer_samples)
    all_normal = pd.concat(normal_samples)

    cancer_beta = all_cancer.groupby(['chrom', 'start', 'end'])['beta'].mean()
    normal_beta = all_normal.groupby(['chrom', 'start', 'end'])['beta'].mean()

    merged = pd.DataFrame({'cancer': cancer_beta, 'normal': normal_beta}).dropna()
    merged['diff'] = merged['cancer'] - merged['normal']

    dmrs = merged[abs(merged['diff']) >= min_diff]
    return dmrs.sort_values('diff', key=abs, ascending=False)


def tissue_deconvolution(sample_meth, reference_atlas):
    '''Deconvolve tissue composition from methylation.'''
    common = sample_meth.index.intersection(reference_atlas.index)
    sample_vec = sample_meth.loc[common].values
    ref_matrix = reference_atlas.loc[common].values

    proportions, _ = nnls(ref_matrix, sample_vec)
    proportions = proportions / proportions.sum()

    return dict(zip(reference_atlas.columns, proportions))


def calculate_cancer_score(sample_meth, cancer_markers, normal_range):
    '''Calculate cancer detection score from marker regions.'''
    scores = []

    for marker in cancer_markers:
        if marker in sample_meth.index:
            beta = sample_meth.loc[marker]
            normal_mean, normal_std = normal_range[marker]
            z = (beta - normal_mean) / normal_std
            scores.append(z)

    return np.mean(scores) if scores else np.nan


if __name__ == '__main__':
    print('cfDNA Methylation Analysis')
    print('=' * 40)
    print('1. extract_methylation() - Extract from bisulfite BAM')
    print('2. find_dmrs() - Find differential methylation')
    print('3. tissue_deconvolution() - Deconvolve tissue origins')
    print('4. calculate_cancer_score() - Score cancer markers')
