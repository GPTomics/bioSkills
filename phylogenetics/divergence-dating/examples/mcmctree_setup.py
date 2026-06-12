'''Generate MCMCTree control files and calibrated tree strings for divergence-time estimation.

A divergence date is mostly a product of the calibration prior, so this helper builds
both the prior-only (usedata=0) and approximate-likelihood (usedata=2) control files and
encodes every fossil as a soft-bounded MINIMUM, never a point. Outputs go to a temp dir.
'''
# Reference: PAML/MCMCTree 4.10+ | Verify API if version differs

import os
import tempfile
from pathlib import Path


def write_control_file(outpath, seqfile, treefile, outfile='mcmctree.out',
                       ndata=1, usedata=2, clock=2, model=4, alpha=0.5,
                       ncatg=4, cleandata=0, bdparas='1 1 0.1',
                       rgene_gamma='2 20 1', sigma2_gamma='1 10 1',
                       burnin=50000, sampfreq=50, nsample=20000, seed=-1):
    '''Write an MCMCTree .ctl file.

    usedata: 0=prior only (effective prior), 1=exact likelihood, 2=approximate, 3=generate in.BV.
    clock: 1=strict, 2=independent rates, 3=autocorrelated rates.
    model: 0=JC69, 4=HKY85, 7=REV (GTR). RootAge is mandatory or MCMCTree refuses to run.
    '''
    lines = [
        f'seed = {seed}',
        f'seqfile = {seqfile}',
        f'treefile = {treefile}',
        f'outfile = {outfile}',
        '',
        f'ndata = {ndata}',
        f'usedata = {usedata}',
        f'clock = {clock}',
        'RootAge = <1.0',
        '',
        f'model = {model}',
        f'alpha = {alpha}',
        f'ncatG = {ncatg}',
        f'cleandata = {cleandata}',
        '',
        f'BDparas = {bdparas}',
        'kappa_gamma = 6 2',
        'alpha_gamma = 1 1',
        f'rgene_gamma = {rgene_gamma}',
        f'sigma2_gamma = {sigma2_gamma}',
        '',
        'print = 1',
        f'burnin = {burnin}',
        f'sampfreq = {sampfreq}',
        f'nsample = {nsample}',
    ]
    Path(outpath).write_text('\n'.join(lines) + '\n')
    return outpath


def format_calibration(cal_type, *args):
    '''Format one MCMCTree calibration string in B()/L()/U() notation.

    B(tL, tU, pL, pU) soft lower and upper bounds; L(tL, p, c, pL) minimum with Cauchy tail;
    U(tU, pU) maximum only. Times are in units of 100 Myr by convention (0.6 = 60 Ma).
    The B()/L()/U() syntax avoids the >/< parsing bug. pL=pU=0.025 is the canonical soft tail.
    '''
    params = ', '.join(str(a) for a in args)
    return f'{cal_type}({params})'


def build_calibrated_tree(newick, calibrations):
    '''Replace internal-node labels with quoted MCMCTree calibration strings.

    calibrations maps node labels to B()/L()/U() strings, each a soft-bounded minimum.
    '''
    result = newick
    for label, cal_str in calibrations.items():
        result = result.replace(label, f"'{cal_str}'")
    return result


def generate_prior_and_posterior_configs(seqfile, treefile, outdir):
    '''Write the mandatory prior-only run (usedata=0) and the approximate-likelihood run.

    Comparing the prior-only effective prior to the posterior on each calibrated node is the
    single most important interpretation step; a posterior that matches the prior is uninformed.
    The in.BV gradient/Hessian (usedata=3) feeds the usedata=2 approximate-likelihood MCMC.
    '''
    os.makedirs(outdir, exist_ok=True)
    prior_path = write_control_file(os.path.join(outdir, 'mcmctree_prior.ctl'),
                                    seqfile, treefile, usedata=0, outfile='prior.out')
    bv_path = write_control_file(os.path.join(outdir, 'mcmctree_step1_bv.ctl'),
                                 seqfile, treefile, usedata=3, outfile='step1.out')
    post_path = write_control_file(os.path.join(outdir, 'mcmctree_step2_post.ctl'),
                                   seqfile, treefile, usedata=2, outfile='posterior.out')
    return prior_path, bv_path, post_path


if __name__ == '__main__':
    outdir = tempfile.mkdtemp(prefix='mcmctree_demo_')

    bounds_cal = format_calibration('B', 0.06, 0.08, 0.025, 0.025)
    lower_cal = format_calibration('L', 0.12, 0.05, 1.0, 0.025)
    upper_cal = format_calibration('U', 1.0, 0.025)
    [bounds_cal, lower_cal, upper_cal]

    tree = '((((human, chimp) human_chimp, gorilla) ape_root, mouse), rat);'
    calibrations = {'human_chimp': 'B(0.06, 0.08, 0.025, 0.025)', 'ape_root': 'L(0.12, 0.05, 1.0, 0.025)'}
    calibrated = build_calibrated_tree(tree, calibrations)

    prior, bv, post = generate_prior_and_posterior_configs('alignment.phy', 'calibrated_tree.nwk', outdir)
    print(f'calibrated tree: {calibrated}')
    print(f'prior-only ctl (usedata=0): {prior}')
    print(f'in.BV ctl (usedata=3): {bv}')
    print(f'posterior ctl (usedata=2): {post}')
