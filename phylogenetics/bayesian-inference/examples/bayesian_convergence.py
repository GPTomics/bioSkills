'''Check MCMC convergence from MrBayes .p (parameter) files across two independent runs'''
# Reference: biopython 1.83+, numpy 1.24+, pandas 2.0+ | Verify API if version differs

import sys
import numpy as np
import pandas as pd


def parse_mrbayes_pfile(filepath):
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('['):
                continue
            if line.startswith('Gen'):
                header = line.split('\t')
                continue
            rows.append(line.split('\t'))
    df = pd.DataFrame(rows, columns=header)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col])
    return df


def compute_ess(series):
    n = len(series)
    if n < 10:
        return 0.0
    mean = series.mean()
    centered = series.values - mean
    variance = np.var(centered, ddof=1)
    if variance == 0:
        return float(n)
    max_lag = min(n // 2, 1000)
    autocorr = np.correlate(centered, centered, mode='full')
    autocorr = autocorr[n - 1:] / (variance * n)
    cumulative = 0.0
    for lag in range(1, max_lag):
        if autocorr[lag] < 0.05:
            break
        cumulative += autocorr[lag]
    ess = n / (1.0 + 2.0 * cumulative)
    return max(ess, 1.0)


def compute_psrf(series1, series2):
    n1, n2 = len(series1), len(series2)
    n = min(n1, n2)
    s1, s2 = series1.iloc[:n], series2.iloc[:n]
    mean1, mean2 = s1.mean(), s2.mean()
    grand_mean = (mean1 + mean2) / 2.0
    between_var = n * ((mean1 - grand_mean) ** 2 + (mean2 - grand_mean) ** 2)
    within_var = (s1.var(ddof=1) + s2.var(ddof=1)) / 2.0
    if within_var == 0:
        return 1.0
    var_estimate = ((n - 1) / n) * within_var + (1 / n) * between_var
    psrf = np.sqrt(var_estimate / within_var)
    return psrf


def assess_convergence(pfile_run1, pfile_run2, burnin_fraction=0.25):
    print(f'Parsing {pfile_run1} and {pfile_run2}...\n')
    df1_raw = parse_mrbayes_pfile(pfile_run1)
    df2_raw = parse_mrbayes_pfile(pfile_run2)

    burnin1 = int(len(df1_raw) * burnin_fraction)
    burnin2 = int(len(df2_raw) * burnin_fraction)
    df1 = df1_raw.iloc[burnin1:].reset_index(drop=True)
    df2 = df2_raw.iloc[burnin2:].reset_index(drop=True)
    print(f'Samples after {burnin_fraction:.0%} burn-in: run1={len(df1)}, run2={len(df2)}\n')

    skip_cols = {'Gen'}
    param_cols = [c for c in df1.columns if c not in skip_cols and c in df2.columns]

    results = []
    for col in param_cols:
        ess1 = compute_ess(df1[col])
        ess2 = compute_ess(df2[col])
        psrf = compute_psrf(df1[col], df2[col])
        min_ess = min(ess1, ess2)
        converged = min_ess >= 200 and psrf < 1.05
        results.append({'parameter': col, 'ess_run1': ess1, 'ess_run2': ess2, 'min_ess': min_ess, 'psrf': psrf, 'converged': converged})

    results_df = pd.DataFrame(results)
    print('Parameter Convergence Summary')
    print('=' * 80)
    for _, row in results_df.iterrows():
        status = 'OK' if row['converged'] else 'FAIL'
        print(f"  {row['parameter']:20s}  ESS: {row['ess_run1']:8.1f} / {row['ess_run2']:8.1f}  PSRF: {row['psrf']:.4f}  [{status}]")

    failed = results_df[~results_df['converged']]
    print(f'\n{"=" * 80}')
    if len(failed) == 0:
        print('All parameters converged (ESS >= 200, PSRF < 1.05)')
    else:
        print(f'{len(failed)} parameter(s) NOT converged:')
        for _, row in failed.iterrows():
            reasons = []
            if row['min_ess'] < 200:
                reasons.append(f'low ESS ({row["min_ess"]:.1f})')
            if row['psrf'] >= 1.05:
                reasons.append(f'high PSRF ({row["psrf"]:.4f})')
            print(f"  {row['parameter']}: {', '.join(reasons)}")
        print('\nRecommendation: Run chains longer (increase ngen). Do NOT just increase samplefreq.')

    return results_df


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: python bayesian_convergence.py run1.p run2.p [burnin_fraction]')
        print('  burnin_fraction: proportion to discard (default 0.25)')
        sys.exit(1)
    pfile1 = sys.argv[1]
    pfile2 = sys.argv[2]
    burnin = float(sys.argv[3]) if len(sys.argv) > 3 else 0.25
    assess_convergence(pfile1, pfile2, burnin)
