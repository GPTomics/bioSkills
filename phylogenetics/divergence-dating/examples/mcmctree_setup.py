'''Generate MCMCTree control files and calibrated tree files for divergence time estimation'''
# Reference: PAML 4.10+ | Verify API if version differs

import os
from pathlib import Path


def write_control_file(outpath, seqfile, treefile, outfile='mcmctree.out',
                       ndata=1, usedata=2, clock=2, model=4, alpha=0.5,
                       ncatg=4, cleandata=0, bdparas='1 1 0.1',
                       rgene_gamma='2 20 1', sigma2_gamma='1 10 1',
                       burnin=50000, sampfreq=50, nsample=20000, seed=-1):
    lines = [
        f'seed = {seed}',
        f'seqfile = {seqfile}',
        f'treefile = {treefile}',
        f'outfile = {outfile}',
        '',
        f'ndata = {ndata}',
        f'usedata = {usedata}',
        f'clock = {clock}',
        f'RootAge = <1.0',
        '',
        f'model = {model}',
        f'alpha = {alpha}',
        f'ncatG = {ncatg}',
        f'cleandata = {cleandata}',
        '',
        f'BDparas = {bdparas}',
        f'kappa_gamma = 6 2',
        f'alpha_gamma = 1 1',
        f'rgene_gamma = {rgene_gamma}',
        f'sigma2_gamma = {sigma2_gamma}',
        '',
        f'print = 1',
        f'burnin = {burnin}',
        f'sampfreq = {sampfreq}',
        f'nsample = {nsample}',
    ]
    Path(outpath).write_text('\n'.join(lines) + '\n')
    return outpath


def format_calibration(cal_type, *args):
    '''Format a single MCMCTree calibration string.

    cal_type: one of 'B' (bounds), 'L' (lower/minimum), 'U' (upper/maximum)
    args: numeric parameters matching MCMCTree notation
      B(tL, tU, pL, pU) - soft lower and upper bounds
      L(tL, p, c, pL)   - minimum bound with Cauchy tail
      U(tU, pU)          - maximum bound
    '''
    params = ', '.join(str(a) for a in args)
    return f"'{cal_type}({params})'"


def build_calibrated_tree(newick, calibrations):
    '''Insert calibration annotations into a Newick tree string.

    calibrations: dict mapping node labels to calibration strings,
      e.g. {'human_chimp': "B(0.06, 0.08, 0.025, 0.025)"}
    '''
    result = newick
    for label, cal_str in calibrations.items():
        cal_notation = f"'{cal_str}'"
        result = result.replace(label, cal_notation)
    return result


def generate_approx_likelihood_configs(seqfile, treefile, outdir='.'):
    '''Generate both step 1 (in.BV) and step 2 (MCMC) control files.

    Step 1: usedata=3 generates the gradient/Hessian file (in.BV).
    Step 2: usedata=2 runs MCMC with approximate likelihood.
    '''
    os.makedirs(outdir, exist_ok=True)

    step1_path = os.path.join(outdir, 'mcmctree_step1.ctl')
    write_control_file(step1_path, seqfile, treefile, usedata=3,
                       outfile='mcmctree_step1.out')

    step2_path = os.path.join(outdir, 'mcmctree_step2.ctl')
    write_control_file(step2_path, seqfile, treefile, usedata=2,
                       outfile='mcmctree_step2.out')

    return step1_path, step2_path


if __name__ == '__main__':
    ctl_path = write_control_file(
        'mcmctree.ctl',
        seqfile='alignment.phy',
        treefile='calibrated_tree.nwk',
        clock=2,
        model=7,
        burnin=100000,
        sampfreq=100,
        nsample=10000,
    )
    print(f'Control file written to: {ctl_path}')

    bounds_cal = format_calibration('B', 0.06, 0.08, 0.025, 0.025)
    lower_cal = format_calibration('L', 0.12, 0.05, 1.0, 0.025)
    upper_cal = format_calibration('U', 1.0, 0.025)
    print(f'Bounds calibration: {bounds_cal}')
    print(f'Lower bound calibration: {lower_cal}')
    print(f'Upper bound calibration: {upper_cal}')

    tree = '((((human, chimp) human_chimp, gorilla) ape_root, mouse), rat);'
    calibrations = {
        'human_chimp': 'B(0.06, 0.08, 0.025, 0.025)',
        'ape_root': 'L(0.12, 0.05, 1.0, 0.025)',
    }
    calibrated = build_calibrated_tree(tree, calibrations)
    print(f'Calibrated tree: {calibrated}')

    step1, step2 = generate_approx_likelihood_configs(
        'alignment.phy', 'calibrated_tree.nwk', outdir='mcmctree_run'
    )
    print(f'Step 1 (in.BV generation): {step1}')
    print(f'Step 2 (MCMC sampling): {step2}')
