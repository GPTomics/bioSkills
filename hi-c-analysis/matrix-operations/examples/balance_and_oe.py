'''Balance a Hi-C cooler, compute cis expected (P(s)), build an O/E matrix, and take the P(s) log-derivative.'''
# Reference: cooler 0.10+, cooltools 0.7+, bioframe 0.7+ | Verify API if version differs
import cooler
import cooltools
import bioframe
import numpy as np

COOL_URI = 'matrix.mcool::/resolutions/10000'
REGION = 'chr1'
MAD_MAX = 5        # drop bins >5 MAD below median log-marginal; near-empty bins otherwise get exploding weights
IGNORE_DIAGS = 2   # diag 0 = self-ligation/dangling, diag 1 = undigested/religated: ligation chemistry, not 3D contact
SMOOTH_SIGMA = 0.1 # cooltools default: Gaussian std in log10(distance) units for P(s) smoothing


def balance(uri):
    clr = cooler.Cooler(uri)
    bias, stats = cooler.balance_cooler(clr, cis_only=True, mad_max=MAD_MAX, ignore_diags=IGNORE_DIAGS, store=True)
    print('converged:', stats['converged'], '| scale:', stats['scale'], '| divisive:', stats['divisive_weights'])
    return cooler.Cooler(uri)


def cis_expected(clr):
    view_df = bioframe.make_viewframe(clr.chromsizes)
    return cooltools.expected_cis(clr, view_df=view_df, smooth=True, aggregate_smoothed=True, ignore_diags=IGNORE_DIAGS)


def oe_matrix(clr, region, cvd):
    obs = clr.matrix(balance=True).fetch(region)
    chrom = region.split(':')[0] if isinstance(region, str) else region[0]
    exp = cvd[cvd['region1'] == chrom].set_index('dist')['balanced.avg']
    exp_by_dist = exp.reindex(range(obs.shape[0])).to_numpy()
    idx = np.arange(obs.shape[0])
    expected = exp_by_dist[np.abs(np.subtract.outer(idx, idx))]
    return obs / expected


def ps_log_derivative(cvd):
    agg = cvd[(cvd['region1'] == cvd['region2']) & (cvd['dist'] > 0)].drop_duplicates('dist_bp')
    return agg['dist_bp'].to_numpy(), np.gradient(np.log(agg['balanced.avg.smoothed.agg']), np.log(agg['dist_bp']))


clr = balance(COOL_URI)
cvd = cis_expected(clr)
print('expected columns:', list(cvd.columns))

oe = oe_matrix(clr, REGION, cvd)
log_oe = np.log2(oe)
print(f'log2(O/E) {REGION}: {np.nanmin(log_oe):.2f} to {np.nanmax(log_oe):.2f}')

dist_bp, slope = ps_log_derivative(cvd)
mid = (dist_bp > 1e5) & (dist_bp < 1e6)   # reference window where slope ~ -1 (crumpled globule)
print(f'P(s) mean log-slope over 0.1-1 Mb: {np.nanmean(slope[mid]):.2f}')
