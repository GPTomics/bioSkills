#!/usr/bin/env python3
'''cfDNA methylation analysis: MethylDackel extraction, beta values, NNLS tissue
deconvolution, and region-level DMR discovery with explicit BH FDR.'''
# Reference: numpy 1.26+, pandas 2.2+, scipy 1.12+ | Verify API if version differs

import subprocess
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls
from statsmodels.stats.multitest import multipletests


def choose_trimming(reference, bam, prefix):
    'Run MethylDackel mbias to write M-bias plots and print the suggested --OT/--OB string.'
    subprocess.run(['MethylDackel', 'mbias', reference, bam, prefix], check=True)


def extract_methylation(reference, bam, prefix, min_depth=1, trim=None):
    '''Extract per-CpG methylation with MethylDackel and parse the fixed bedGraph.
    bedGraph columns: chrom / start / end / methylation%(int) / count_meth / count_unmeth.
    min_depth kept low because per-CpG depth gates erase low-coverage cfDNA signal.'''
    cmd = ['MethylDackel', 'extract', reference, bam, '-o', prefix,
           '--mergeContext', '--minDepth', str(min_depth)]
    if trim:
        cmd += ['--OT', trim, '--OB', trim]
    subprocess.run(cmd, check=True)

    meth = pd.read_csv(f'{prefix}_CpG.bedGraph', sep='\t', skiprows=1, header=None,
                       names=['chrom', 'start', 'end', 'meth_pct', 'count_meth', 'count_unmeth'])
    meth['beta'] = meth['count_meth'] / (meth['count_meth'] + meth['count_unmeth'])
    return meth


def deconvolve_tissue(sample_beta, atlas):
    '''Atlas-based tissue-of-origin via NNLS with the simplex constraint (w>=0, sum w=1).
    sample_beta: Series indexed by marker. atlas: DataFrame markers x cell_types.'''
    markers = sample_beta.index.intersection(atlas.index)
    w, _ = nnls(atlas.loc[markers].values, sample_beta.loc[markers].values)
    w = w / w.sum()
    return dict(zip(atlas.columns, w))


def region_dmrs(cancer, normal, region_col='region', min_samples=3):
    '''Region-level DMR discovery (correct altitude vs per-CpG t-tests).
    cancer/normal: long DataFrames with [region_col, beta], one row per sample-region.
    FDR via Benjamini-Hochberg specified explicitly (statsmodels default is Holm-Sidak).'''
    rows = []
    for region, c in cancer.groupby(region_col)['beta']:
        n = normal.loc[normal[region_col] == region, 'beta'].dropna()
        c = c.dropna()
        if len(c) < min_samples or len(n) < min_samples:
            continue
        _, p = stats.mannwhitneyu(c, n, alternative='two-sided')
        rows.append((region, c.mean() - n.mean(), p))

    res = pd.DataFrame(rows, columns=['region', 'delta_beta', 'pvalue'])
    if len(res):
        res['fdr'] = multipletests(res['pvalue'], method='fdr_bh')[1]
    return res.sort_values('fdr') if len(res) else res


def haplotype_concordance_score(read_states, background_p=0.05):
    '''Read-level detection primitive: a fragment concordantly methylated across k CpGs in a
    tumor block has background probability ~background_p**k under the WBC null. Returns the
    per-read negative-log10 background probability; large values flag tumor-derived molecules.
    read_states: list of arrays of 0/1 methylation calls per CpG along each fragment.'''
    scores = []
    for states in read_states:
        states = np.asarray(states)
        if states.mean() > 0.9:
            k = len(states)
            scores.append(-np.log10(background_p ** k))
        else:
            scores.append(0.0)
    return np.array(scores)


if __name__ == '__main__':
    atlas = pd.DataFrame(
        {'wbc': [0.85, 0.10, 0.90, 0.05], 'liver': [0.10, 0.88, 0.20, 0.92], 'colon': [0.15, 0.80, 0.10, 0.85]},
        index=['m1', 'm2', 'm3', 'm4'])
    sample = pd.Series([0.80, 0.18, 0.84, 0.13], index=['m1', 'm2', 'm3', 'm4'])
    fractions = deconvolve_tissue(sample, atlas)
    print('Tissue-of-origin fractions:', {k: round(v, 3) for k, v in fractions.items()})

    reads = [[1, 1, 1, 1, 1, 1], [0, 1, 0, 1, 0, 0], [1, 1, 1, 1, 1]]
    print('Read haplotype concordance scores:', haplotype_concordance_score(reads).round(2))
