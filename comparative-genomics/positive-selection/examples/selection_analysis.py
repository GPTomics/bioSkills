'''Positive selection analysis with PAML codeml'''

import os
import subprocess
from scipy import stats


def write_codeml_control(alignment, tree, outfile, model='M8', output_dir='.'):
    '''Generate codeml control file for site models

    Model comparison for positive selection:
    - M7 (null): beta distribution, 0 < omega < 1
    - M8 (alternative): beta + omega > 1 allowed
    - LRT with df=2, chi2 critical value 5.99 (p<0.05)

    More stringent comparison:
    - M8a (null): beta + omega fixed at 1
    - M8 (alternative): beta + omega estimated
    - LRT with df=1, chi2 critical value 3.84 (p<0.05)
    '''
    nssites_map = {
        'M0': 0,   # One ratio
        'M1a': 1,  # Nearly neutral
        'M2a': 2,  # Positive selection
        'M7': 7,   # Beta
        'M8': 8,   # Beta + omega
        'M8a': 8,  # Beta + omega=1 (fix_omega=1)
    }

    fix_omega = 1 if model == 'M8a' else 0
    nssites = nssites_map.get(model, 8)

    ctl = f'''      seqfile = {alignment}
     treefile = {tree}
      outfile = {outfile}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0

        model = 0
      NSsites = {nssites}

        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = {fix_omega}
        omega = 1

       getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 1
       method = 0
'''
    ctl_path = os.path.join(output_dir, f'codeml_{model}.ctl')
    with open(ctl_path, 'w') as f:
        f.write(ctl)

    return ctl_path


def run_codeml(ctl_file, work_dir='.'):
    '''Execute codeml with control file'''
    original_dir = os.getcwd()
    os.chdir(work_dir)

    result = subprocess.run(['codeml', ctl_file], capture_output=True, text=True)

    os.chdir(original_dir)
    return result.returncode == 0


def parse_mlc_file(mlc_file):
    '''Parse codeml output file for results

    Extracts:
    - Log-likelihood (lnL)
    - Omega estimates
    - Positively selected sites (BEB analysis)
    '''
    results = {'lnL': None, 'np': None, 'omega_classes': [], 'selected_sites': []}

    with open(mlc_file) as f:
        lines = f.readlines()

    in_beb = False
    for i, line in enumerate(lines):
        # Log-likelihood line: lnL(ntime: X  np: Y):  Z
        if line.startswith('lnL'):
            parts = line.split()
            for j, p in enumerate(parts):
                if p.startswith('np:'):
                    results['np'] = int(parts[j + 1].rstrip('):'))
                if ')' in p and j + 1 < len(parts):
                    try:
                        results['lnL'] = float(parts[j + 1])
                    except (ValueError, IndexError):
                        pass

        # Omega class proportions (for site models)
        if 'proportion' in line.lower() and 'omega' in lines[i + 1].lower() if i + 1 < len(lines) else False:
            pass  # Complex parsing for omega classes

        # BEB (Bayes Empirical Bayes) results
        if 'Bayes Empirical Bayes' in line:
            in_beb = True
            continue

        if in_beb and line.strip():
            parts = line.split()
            if len(parts) >= 3:
                try:
                    site = int(parts[0])
                    aa = parts[1]
                    prob = float(parts[2])
                    # BEB probability thresholds:
                    # > 0.95: Site under positive selection (*)
                    # > 0.99: Highly significant (**)
                    if prob > 0.5:
                        results['selected_sites'].append({
                            'site': site,
                            'aa': aa,
                            'prob': prob,
                            'sig': '**' if prob > 0.99 else '*' if prob > 0.95 else ''
                        })
                except ValueError:
                    if 'site' not in line.lower():
                        in_beb = False

    return results


def likelihood_ratio_test(lnL_null, lnL_alt, df):
    '''Perform likelihood ratio test for model comparison

    LRT = 2 * (lnL_alt - lnL_null)

    Critical values (chi-square distribution):
    df=1: 3.84 (p=0.05), 6.63 (p=0.01), 10.83 (p=0.001)
    df=2: 5.99 (p=0.05), 9.21 (p=0.01), 13.82 (p=0.001)
    '''
    if lnL_alt < lnL_null:
        print('Warning: Alternative model has lower likelihood than null')
        return {'LRT': 0, 'p_value': 1.0, 'significant': False}

    lrt = 2 * (lnL_alt - lnL_null)
    p_value = 1 - stats.chi2.cdf(lrt, df)

    return {
        'LRT': lrt,
        'df': df,
        'p_value': p_value,
        'significant': p_value < 0.05,
        'highly_significant': p_value < 0.01
    }


def summarize_selection_test(m7_results, m8_results):
    '''Summarize M8 vs M7 selection test'''
    print('Positive Selection Test: M8 vs M7')
    print('=' * 50)

    print(f"\nM7 (null - beta): lnL = {m7_results['lnL']:.2f}")
    print(f"M8 (alt - beta+omega): lnL = {m8_results['lnL']:.2f}")

    # M8 vs M7: df = 2 (omega class + proportion)
    lrt = likelihood_ratio_test(m7_results['lnL'], m8_results['lnL'], df=2)

    print(f"\nLikelihood Ratio Test:")
    print(f"  LRT statistic: {lrt['LRT']:.4f}")
    print(f"  Degrees of freedom: {lrt['df']}")
    print(f"  P-value: {lrt['p_value']:.6f}")

    if lrt['significant']:
        print('\n*** Positive selection DETECTED ***')
    else:
        print('\nNo significant evidence for positive selection')

    # Report selected sites
    if m8_results['selected_sites']:
        sig_sites = [s for s in m8_results['selected_sites'] if s['sig']]
        print(f"\nPositively selected sites (BEB P > 0.95): {len(sig_sites)}")
        for site in sig_sites:
            print(f"  Site {site['site']}: {site['aa']} (P = {site['prob']:.3f}) {site['sig']}")

    return lrt


if __name__ == '__main__':
    print('PAML Codeml Selection Analysis')
    print('=' * 50)

    # Example results (simulated)
    m7_example = {'lnL': -1234.56, 'np': 13, 'selected_sites': []}
    m8_example = {
        'lnL': -1228.12,
        'np': 15,
        'selected_sites': [
            {'site': 45, 'aa': 'K', 'prob': 0.982, 'sig': '*'},
            {'site': 89, 'aa': 'R', 'prob': 0.996, 'sig': '**'},
            {'site': 123, 'aa': 'D', 'prob': 0.871, 'sig': ''},
            {'site': 156, 'aa': 'N', 'prob': 0.965, 'sig': '*'},
        ]
    }

    summarize_selection_test(m7_example, m8_example)

    print('\n\nTo run on real data:')
    print('1. Prepare codon alignment (in-frame, PHYLIP format)')
    print('2. Generate phylogenetic tree')
    print('3. ctl_m7 = write_codeml_control(aln, tree, "m7.mlc", "M7")')
    print('4. ctl_m8 = write_codeml_control(aln, tree, "m8.mlc", "M8")')
    print('5. run_codeml(ctl_m7); run_codeml(ctl_m8)')
    print('6. m7 = parse_mlc_file("m7.mlc")')
    print('7. m8 = parse_mlc_file("m8.mlc")')
    print('8. summarize_selection_test(m7, m8)')
