"""
Power and sample size calculations for clinical trials.

Reference: statsmodels 0.14+ | Verify API if version differs
Reference: scipy 1.12+ | Verify API if version differs

Covers superiority, non-inferiority, equivalence, and time-to-event designs.
Includes the critical δ-vs-MCID distinction and FDA 2016 NI double discount.

NOTE: For confirmatory regulatory work, use R packages gsDesign, rpact, simtrial.
"""

import math
import numpy as np
from statsmodels.stats.power import tt_ind_solve_power


# ------------------------------------------------------------------
# 1. Two-sample t-test sample size (continuous)
# ------------------------------------------------------------------
def two_sample_t_sample_size(mean_diff, sd, alpha=0.05, power=0.80, dropout=0.0):
    """Two-sided two-sample t-test n per arm. Inflated for dropout."""
    effect_size = abs(mean_diff) / sd  # Cohen's d
    n = tt_ind_solve_power(
        effect_size=effect_size, alpha=alpha, power=power, alternative='two-sided'
    )
    n_per_arm = math.ceil(n / (1 - dropout))
    return n_per_arm


# ------------------------------------------------------------------
# 2. Two-sample proportions (binary) -- continuity correction debate
# ------------------------------------------------------------------
def two_proportions_sample_size(p1, p2, alpha=0.05, power=0.80,
                                  alternative='two-sided', continuity_correction=False):
    """
    Sample size per arm for two-sample proportion test.

    continuity_correction=True applies Fleiss's correction; debated.
    For n<100/arm, simulation-based SS is preferred over either formula
    because the asymptotic approximation breaks down (D'Agostino 1988).
    """
    z_alpha = 1.96 if alternative == 'two-sided' else 1.645
    z_beta = 0.84  # power=0.80
    pbar = (p1 + p2) / 2

    n = ((z_alpha * math.sqrt(2 * pbar * (1 - pbar)) +
          z_beta * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) / abs(p1 - p2)) ** 2

    if continuity_correction:
        correction = (1 + math.sqrt(1 + 4 / (n * abs(p1 - p2)))) ** 2 / 4
        n = n * correction

    return math.ceil(n)


# ------------------------------------------------------------------
# 3. Schoenfeld 1981 -- events for two-sample log-rank under PH
# ------------------------------------------------------------------
def schoenfeld_events(hr, alpha=0.05, power=0.80, allocation_ratio=1.0, one_sided=False):
    """
    Schoenfeld 1981 events formula.
    WARNING: assumes proportional hazards. Under non-PH (immuno-oncology),
    under-estimates required events by 20-50%; use Lakatos or simulation.
    """
    z_alpha = 1.96 if not one_sided else 1.645
    z_beta = 0.84
    p = 1 / (1 + allocation_ratio)
    events = (z_alpha + z_beta) ** 2 / (p * (1 - p) * math.log(hr) ** 2)
    return math.ceil(events)


def schoenfeld_n_total(hr, expected_event_prob, dropout=0.15, **kwargs):
    """Convert events to n per arm assuming uniform event probability."""
    events = schoenfeld_events(hr, **kwargs)
    n_per_arm_no_dropout = math.ceil(events / (2 * expected_event_prob))
    n_per_arm = math.ceil(n_per_arm_no_dropout / (1 - dropout))
    return events, n_per_arm


# ------------------------------------------------------------------
# 4. Non-inferiority margin: FDA 2016 double discount
# ------------------------------------------------------------------
def ni_margin_double_discount(historical_effect_point, historical_effect_lower_ci,
                                retention=0.50, mcid=None):
    """
    FDA 2016 double discount: M2 = retention × lower CI bound of historical M1.

    Always verify M2 < MCID (margin must be smaller than clinically meaningful diff).
    """
    M1 = historical_effect_lower_ci  # lower bound of historical effect CI
    M2 = retention * M1

    print(f'Historical effect point estimate: {historical_effect_point}')
    print(f'Historical effect lower CI bound (M1): {M1}')
    print(f'Retention fraction: {retention} (FDA 2016: 50% conservative)')
    print(f'NI margin M2: {M2}')

    if mcid is not None:
        print(f'MCID: {mcid}')
        if M2 > mcid:
            print(f'WARNING: M2 ({M2}) > MCID ({mcid}). Margin is too liberal; '
                  f'patients within "non-inferior" CI may be clinically inferior.')
        else:
            print(f'OK: M2 ({M2}) <= MCID ({mcid}).')

    return M2


def ni_sample_size_binary(p_control, p_test, margin, alpha=0.025, power=0.80):
    """NI sample size per arm for binary endpoint (fixed-margin approach)."""
    z_alpha = 1.96 if alpha == 0.025 else 1.645
    z_beta = 0.84
    variance_term = p_control * (1 - p_control) + p_test * (1 - p_test)
    n = (z_alpha + z_beta) ** 2 * variance_term / margin ** 2
    return math.ceil(n)


# ------------------------------------------------------------------
# 5. TOST equivalence sample size
# ------------------------------------------------------------------
def tost_sample_size(mu_diff, sigma, margin, alpha=0.05, power=0.80):
    """
    TOST equivalence: two one-sided tests at alpha each.
    Margins ±margin; assumed true |mean_diff| < margin.
    """
    z_alpha = 1.645  # one-sided
    z_beta = 0.84
    effect = abs(margin) - abs(mu_diff)
    if effect <= 0:
        raise ValueError(f'|mu_diff| ({abs(mu_diff)}) must be < margin ({abs(margin)})')
    n = 2 * sigma ** 2 * (z_alpha + z_beta) ** 2 / effect ** 2
    return math.ceil(n)


# ------------------------------------------------------------------
# 6. MCID vs δ distinction
# ------------------------------------------------------------------
def delta_mcid_check(delta, mcid):
    """Warn if δ is too close to MCID for adequate power margin."""
    print(f'Delta (effect to detect): {delta}')
    print(f'MCID (clinically meaningful difference): {mcid}')
    ratio = abs(delta) / abs(mcid)
    print(f'Ratio delta/MCID: {ratio:.2f}')
    if ratio < 1.0:
        print('ERROR: delta < MCID. Trial is underpowered for clinically meaningful effect.')
    elif ratio < 1.5:
        print('WARNING: delta < 1.5 × MCID. Postdoc rule: delta >= 1.5 × MCID for safety margin.')
    else:
        print('OK: delta >= 1.5 × MCID.')


# ------------------------------------------------------------------
# 7. Stratified randomisation efficiency gain (Senn 2013)
# ------------------------------------------------------------------
def stratified_efficiency_gain(r_squared):
    """
    Efficiency gain from stratification on a prognostic factor with R² against outcome.
    Effective n = n / (1 - r²)
    """
    if r_squared >= 1 or r_squared < 0:
        raise ValueError('r_squared must be in [0, 1)')
    gain = 1 / (1 - r_squared)
    print(f'R² of stratification factor against outcome: {r_squared}')
    print(f'Efficiency gain: {gain:.3f} (i.e., effective sample size is {gain:.2f}x)')
    return gain


# ------------------------------------------------------------------
# 8. Demonstration scenarios
# ------------------------------------------------------------------
if __name__ == '__main__':
    print('=' * 60)
    print('Scenario 1: Continuous endpoint superiority')
    print('=' * 60)
    n_cont = two_sample_t_sample_size(mean_diff=5, sd=12, alpha=0.05, power=0.80, dropout=0.15)
    print(f'n per arm (inflated for 15% dropout): {n_cont}')

    print('\n' + '=' * 60)
    print('Scenario 2: Binary endpoint (continuity correction debate)')
    print('=' * 60)
    n_uncorr = two_proportions_sample_size(p1=0.40, p2=0.30, continuity_correction=False)
    n_corr = two_proportions_sample_size(p1=0.40, p2=0.30, continuity_correction=True)
    print(f'n per arm uncorrected: {n_uncorr}')
    print(f'n per arm with continuity correction: {n_corr}')
    print(f'(Correction wastes ~{100*(n_corr-n_uncorr)/n_uncorr:.0f}% per D\'Agostino 1988)')

    print('\n' + '=' * 60)
    print('Scenario 3: Oncology OS Schoenfeld vs Lakatos (PH plausible)')
    print('=' * 60)
    events, n_per_arm = schoenfeld_n_total(hr=0.65, expected_event_prob=0.40,
                                            alpha=0.025, power=0.90, one_sided=True)
    print(f'Schoenfeld events needed: {events}')
    print(f'n per arm (40% event rate, 15% dropout): {n_per_arm}')
    print('NOTE: under non-PH (immuno-oncology), use Lakatos or simulation; '
          'Schoenfeld under-estimates by 20-50%')

    print('\n' + '=' * 60)
    print('Scenario 4: Non-inferiority margin selection (FDA 2016 double discount)')
    print('=' * 60)
    M2 = ni_margin_double_discount(
        historical_effect_point=0.25, historical_effect_lower_ci=0.20,
        retention=0.50, mcid=0.15
    )
    n_ni = ni_sample_size_binary(p_control=0.85, p_test=0.85, margin=M2,
                                   alpha=0.025, power=0.80)
    print(f'NI sample size per arm: {n_ni}')

    print('\n' + '=' * 60)
    print('Scenario 5: TOST equivalence')
    print('=' * 60)
    n_tost = tost_sample_size(mu_diff=0.0, sigma=10, margin=5, alpha=0.05, power=0.80)
    print(f'TOST n per arm: {n_tost}')

    print('\n' + '=' * 60)
    print('Scenario 6: Delta vs MCID check')
    print('=' * 60)
    delta_mcid_check(delta=5, mcid=3)
    delta_mcid_check(delta=2, mcid=3)

    print('\n' + '=' * 60)
    print('Scenario 7: Stratification efficiency gain (Senn 2013)')
    print('=' * 60)
    stratified_efficiency_gain(r_squared=0.25)
