#!/usr/bin/env python3
'''Longitudinal ctDNA monitoring: binary MRD trajectory, censoring-aware clearance
kinetics, confirmed-trend relapse calling, and a monitoring report. Below-LoD draws
are treated as left-censored at the per-sample limit of detection, never as zero.'''
# Reference: numpy 1.26+, pandas 2.2+, scipy 1.12+, matplotlib 3.8+ | Verify API if version differs

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def summarize_trajectory(df):
    '''Baseline, nadir, and baseline-referenced log2 change with below-LoD points censored.
    df columns: timepoint, tumor_fraction, per_sample_lod.'''
    df = df.sort_values('timepoint').copy()
    df['censored'] = df['tumor_fraction'] <= df['per_sample_lod']
    df['detected'] = ~df['censored']
    df['tf_for_log'] = np.where(df['censored'], df['per_sample_lod'], df['tumor_fraction'])
    baseline = df.iloc[0]['tf_for_log']
    df['log2_fc_baseline'] = np.log2(df['tf_for_log'] / baseline)
    detected = df[df['detected']]
    nadir = detected['tumor_fraction'].min() if len(detected) else np.nan
    metrics = {'baseline_tf': baseline, 'nadir_tf': nadir, 'n_censored': int(df['censored'].sum())}
    return df, metrics


def clearance_half_life(df, lod):
    '''Half-life from log-linear decay over the uncensored phase up to the nadir.
    df columns: timepoint, vaf (single mutation). Post-nadir rebound and censored points
    are dropped (so a relapse does not flatten the slope); never log(0)ed.'''
    df = df.sort_values('timepoint')
    n_censored = int((df['vaf'] <= lod).sum())
    uncensored = df[df['vaf'] > lod].reset_index(drop=True)
    if len(uncensored) < 3:  # log-linear slope needs >=3 points for a meaningful CI
        return None
    decay = uncensored.iloc[:uncensored['vaf'].idxmin() + 1]
    if len(decay) < 3:
        return None
    fit = stats.linregress(decay['timepoint'].values, np.log(decay['vaf'].values))
    half_life = np.log(2) / -fit.slope if fit.slope < 0 else np.inf
    return {'half_life_days': half_life, 'slope': fit.slope, 'r_squared': fit.rvalue ** 2,
            'n_points': len(decay), 'n_censored_excluded': n_censored}


def call_molecular_relapse(df, rise_factor=2.0, min_consecutive=2):
    '''Relapse requires detection above LoD and >rise_factor*nadir on >=min_consecutive
    consecutive post-nadir draws (a single excursion is not enough; multiple testing).
    The nadir is floored at the per-sample LoD so the rise bar is never anchored on a
    below-LoD (censored) value. df columns: timepoint, tumor_fraction, per_sample_lod.'''
    df = df.sort_values('timepoint').copy()
    df['detected'] = df['tumor_fraction'] > df['per_sample_lod']
    nadir_time = df.loc[df['tumor_fraction'].idxmin(), 'timepoint']
    nadir_tf = max(df['tumor_fraction'].min(), df['per_sample_lod'].min())
    post = df[df['timepoint'] > nadir_time]
    hit = (post['detected'] & (post['tumor_fraction'] > nadir_tf * rise_factor)).astype(int)
    run = int(hit.groupby((hit == 0).cumsum()).cumsum().max()) if len(hit) else 0
    return {'relapse': bool(run >= min_consecutive), 'nadir_tf': nadir_tf, 'confirmed_consecutive': run}


def classify_response(current_tf, baseline_tf, undetectable, fold_cutoff=100.0):
    '''Response by the supplied fold cutoff (default 2-log/100x; non-harmonized convention).
    undetectable is the binary below-LoD flag for the current draw.'''
    if undetectable:
        return 'molecular complete response (undetectable, LoD-conditional)'
    if current_tf <= baseline_tf / fold_cutoff:
        return 'major molecular response'
    if current_tf < baseline_tf:
        return 'partial molecular response'
    return 'stable / progressive'


def generate_monitoring_report(patient_id, tf_df, mut_df):
    '''Aggregate trajectory, clearance, and relapse into one report.
    tf_df: timepoint, tumor_fraction, per_sample_lod. mut_df: timepoint, mutation, vaf, per_sample_lod.'''
    traj, traj_metrics = summarize_trajectory(tf_df)
    last = traj.iloc[-1]
    relapse = call_molecular_relapse(tf_df)
    response = classify_response(last['tumor_fraction'], traj_metrics['baseline_tf'],
                                 bool(last['censored']))
    kinetics = {}
    for mut, sub in mut_df.groupby('mutation'):
        lod = sub['per_sample_lod'].median()
        res = clearance_half_life(sub[['timepoint', 'vaf']], lod)
        if res is not None:
            kinetics[mut] = res
    return {'patient_id': patient_id, 'baseline_tf': traj_metrics['baseline_tf'],
            'nadir_tf': traj_metrics['nadir_tf'], 'n_censored_draws': traj_metrics['n_censored'],
            'current_response': response, 'relapse_call': relapse, 'clearance_kinetics': kinetics}


def plot_trajectory(traj, lod_band=None, output_file=None):
    '''Log-scale trajectory; censored points drawn at the LoD with a distinct marker, never at zero.'''
    fig, ax = plt.subplots(figsize=(10, 6))
    det = traj[traj['detected']]
    cen = traj[traj['censored']]
    ax.semilogy(det['timepoint'], det['tumor_fraction'], 'o-', color='steelblue',
                linewidth=2, markersize=8, label='detected')
    ax.semilogy(cen['timepoint'], cen['per_sample_lod'], 'v', color='gray',
                markersize=9, label='below LoD (censored)')
    if lod_band is not None:
        ax.axhspan(ax.get_ylim()[0], lod_band, color='red', alpha=0.08, label='below-LoD band')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Tumor fraction')
    ax.set_title('ctDNA MRD trajectory')
    ax.legend()
    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    return fig, ax


if __name__ == '__main__':
    tf_df = pd.DataFrame({
        'timepoint': [0, 28, 56, 90, 150, 210, 250],
        'tumor_fraction': [0.040, 0.006, 0.0008, 0.0003, 0.0003, 0.0012, 0.0030],
        'per_sample_lod': [0.0005, 0.0005, 0.0005, 0.0005, 0.0004, 0.0005, 0.0005]})
    mut_df = pd.DataFrame({
        'timepoint': [0, 28, 56, 90] * 2,
        'mutation': ['TP53_R175H'] * 4 + ['KRAS_G12D'] * 4,
        'vaf': [0.020, 0.003, 0.0004, 0.0002, 0.018, 0.0025, 0.0003, 0.0002],
        'per_sample_lod': [0.0005] * 8})

    report = generate_monitoring_report('PT-001', tf_df, mut_df)
    print('Monitoring report for', report['patient_id'])
    print('  baseline TF:', report['baseline_tf'], '| nadir TF:', report['nadir_tf'])
    print('  censored draws:', report['n_censored_draws'])
    print('  response:', report['current_response'])
    print('  relapse:', report['relapse_call'])
    for mut, k in report['clearance_kinetics'].items():
        print(f'  {mut}: t1/2={k["half_life_days"]:.1f}d R2={k["r_squared"]:.3f} '
              f'(n={k["n_points"]}, censored excluded={k["n_censored_excluded"]})')
