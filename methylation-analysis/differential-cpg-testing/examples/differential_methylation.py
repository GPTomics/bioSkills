# Reference: scipy 1.13+, statsmodels 0.14+, pandas 2.2+, numpy 1.26+ | Verify API if version differs
# Per-CpG differential methylation QUICK-LOOK for ARRAY / CONTINUOUS beta matrices.
# This is NOT the recommended test for bisulfite sequencing COUNTS: a continuous test on
# beta = M/Cov discards coverage (an 8-read site weighs the same as an 800-read site).
# For sequencing counts use a beta-binomial count model (DSS / methylKit overdispersion='MN').
# The discipline here mirrors the count path: test on M-values, report effect on beta (delta-beta),
# gate hits on the intersection of FDR and |delta-beta|.

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

MIN_DELTA_BETA = 0.2   # effect-size floor; pairs with FDR to avoid significance-without-effect at large n
FDR = 0.05             # conventional FDR threshold for differential methylation
M_OFFSET = 1e-3        # boundary-safe logit offset so beta in {0,1} does not blow up the M-value


def beta_to_m(beta):
    return np.log2((beta + M_OFFSET) / (1 - beta + M_OFFSET))


def test_dmc(beta_case, beta_ctrl):
    m_case, m_ctrl = beta_to_m(beta_case), beta_to_m(beta_ctrl)
    # equal_var=False -> Welch (scipy default is Student's); test on M-values, not raw beta
    _, pvalues = ttest_ind(m_case, m_ctrl, axis=1, equal_var=False, nan_policy='omit')
    # method='fdr_bh' -> Benjamini-Hochberg (statsmodels default is 'hs', Holm-Sidak)
    _, padj, _, _ = multipletests(pvalues, alpha=FDR, method='fdr_bh')
    delta_beta = beta_case.mean(axis=1) - beta_ctrl.mean(axis=1)   # effect reported on the beta scale
    hit = (padj < FDR) & (np.abs(delta_beta) >= MIN_DELTA_BETA)
    return pd.DataFrame({'cpg_id': beta_case.index, 'mean_case_beta': beta_case.mean(axis=1).values,
                         'mean_ctrl_beta': beta_ctrl.mean(axis=1).values, 'delta_beta': delta_beta.values,
                         'pvalue': pvalues, 'padj': padj, 'hit': hit.values})


def simulate(n_cpg=2000, n_per_group=6, n_true=100, seed=0):
    rng = np.random.default_rng(seed)
    base = rng.uniform(0.05, 0.95, n_cpg)
    case = np.clip(base[:, None] + rng.normal(0, 0.03, (n_cpg, n_per_group)), 0.001, 0.999)
    ctrl = np.clip(base[:, None] + rng.normal(0, 0.03, (n_cpg, n_per_group)), 0.001, 0.999)
    shift = rng.choice([0.3, -0.3], n_true)
    case[:n_true] = np.clip(case[:n_true] + shift[:, None], 0.001, 0.999)
    idx = pd.Index([f'cpg_{i}' for i in range(n_cpg)], name='cpg_id')
    cols_case = [f'case_{i}' for i in range(n_per_group)]
    cols_ctrl = [f'ctrl_{i}' for i in range(n_per_group)]
    return pd.DataFrame(case, index=idx, columns=cols_case), pd.DataFrame(ctrl, index=idx, columns=cols_ctrl), n_true


if __name__ == '__main__':
    beta_case, beta_ctrl, n_true = simulate()
    res = test_dmc(beta_case, beta_ctrl)
    n_hit = int(res['hit'].sum())
    recovered = int(res['hit'].iloc[:n_true].sum())
    print(f'CpGs tested: {len(res)}; hits (FDR<{FDR} AND |delta_beta|>={MIN_DELTA_BETA}): {n_hit}')
    print(f'true DMCs spiked: {n_true}; recovered among hits: {recovered}')
