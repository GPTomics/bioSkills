#!/usr/bin/env python3
'''
Longitudinal ctDNA monitoring for treatment response.
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def analyze_tf_dynamics(patient_data):
    '''Analyze tumor fraction dynamics over treatment.'''
    patient_data = patient_data.sort_values('timepoint').copy()

    baseline_tf = patient_data.iloc[0]['tumor_fraction']
    patient_data['log2_fc'] = np.log2(patient_data['tumor_fraction'] / baseline_tf)

    min_tf = patient_data['tumor_fraction'].min()
    min_idx = patient_data['tumor_fraction'].idxmin()

    metrics = {
        'baseline_tf': baseline_tf,
        'nadir_tf': min_tf,
        'nadir_timepoint': patient_data.loc[min_idx, 'timepoint'],
        'max_reduction': baseline_tf - min_tf,
        'log2_max_reduction': np.log2(baseline_tf / min_tf) if min_tf > 0 else np.inf
    }

    return patient_data, metrics


def define_response(tf_series, baseline, criteria='2log'):
    '''Define molecular response based on reduction criteria.'''
    if criteria == '2log':
        threshold = baseline / 100
    elif criteria == '1log':
        threshold = baseline / 10
    elif criteria == 'undetectable':
        threshold = 0.001

    return tf_series < threshold


def track_mutations(mutation_data):
    '''Track mutations across timepoints.'''
    pivot = mutation_data.pivot_table(
        index='timepoint', columns='mutation', values='vaf', aggfunc='first'
    )
    pivot['mean_vaf'] = pivot.mean(axis=1)

    last = pivot.index.max()
    cleared = [m for m in pivot.columns if m != 'mean_vaf' and
               (pivot.loc[last, m] < 0.001 or pd.isna(pivot.loc[last, m]))]

    return pivot, cleared


def calculate_clearance_kinetics(mutation_data, mutation):
    '''Calculate mutation clearance half-life.'''
    data = mutation_data[mutation_data['mutation'] == mutation].sort_values('timepoint')

    if len(data) < 3:
        return None

    x = data['timepoint'].values
    y = np.log(data['vaf'].values + 1e-6)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    half_life = np.log(2) / abs(slope) if slope < 0 else np.inf

    return {'mutation': mutation, 'half_life': half_life, 'r_squared': r_value**2}


def detect_relapse(tf_series, threshold_increase=2):
    '''Detect molecular relapse from rising ctDNA.'''
    tf_series = tf_series.sort_values('timepoint')

    nadir_idx = tf_series['tumor_fraction'].idxmin()
    nadir_tf = tf_series.loc[nadir_idx, 'tumor_fraction']
    nadir_time = tf_series.loc[nadir_idx, 'timepoint']

    post_nadir = tf_series[tf_series['timepoint'] > nadir_time]

    for idx, row in post_nadir.iterrows():
        if row['tumor_fraction'] > nadir_tf * threshold_increase:
            return {'relapse': True, 'timepoint': row['timepoint'], 'nadir_tf': nadir_tf}

    return {'relapse': False, 'nadir_tf': nadir_tf}


def plot_ctdna_dynamics(patient_data, treatment_lines=None, output_file=None):
    '''Plot ctDNA dynamics over treatment.'''
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.semilogy(patient_data['timepoint'], patient_data['tumor_fraction'],
                'o-', linewidth=2, markersize=8, color='steelblue')

    if treatment_lines:
        for time, label in treatment_lines:
            ax.axvline(x=time, color='gray', linestyle='--', alpha=0.5)
            ax.text(time, ax.get_ylim()[1], label, rotation=90, va='top', fontsize=9)

    ax.axhline(y=0.01, color='red', linestyle=':', alpha=0.5, label='LOD (1%)')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Tumor Fraction')
    ax.set_title('ctDNA Dynamics During Treatment')
    ax.legend()

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    return fig, ax


if __name__ == '__main__':
    print('Longitudinal ctDNA Monitoring')
    print('=' * 40)
    print('1. analyze_tf_dynamics() - Track tumor fraction')
    print('2. track_mutations() - Monitor specific variants')
    print('3. detect_relapse() - Identify molecular relapse')
    print('4. plot_ctdna_dynamics() - Visualize trends')
