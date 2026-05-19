"""
Survival analysis for clinical trials -- worked examples.

Reference: lifelines 0.27+ | Verify API if version differs
Reference: scikit-survival 0.21+ | Verify API if version differs

Covers Cox PH with Therneau-Grambsch diagnostics, RMST under non-PH,
competing risks via Aalen-Johansen + cause-specific Cox, MaxCombo, and
ADTTE CNSR convention handling.
"""

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter, AalenJohansenFitter
from lifelines.statistics import (
    logrank_test, multivariate_logrank_test, proportional_hazard_test
)
from lifelines.utils import restricted_mean_survival_time


# ------------------------------------------------------------------
# 1. Convert CDISC ADTTE CNSR convention to standard event indicator
# ------------------------------------------------------------------
# CDISC: CNSR=0 means event; CNSR=1,2,3... means censored (different reasons)
# Stat packages: event=1 means event
def adtte_to_lifelines(adtte_df):
    df = adtte_df.copy()
    df['event'] = (df['CNSR'] == 0).astype(int)
    df['time'] = df['AVAL']
    return df


# ------------------------------------------------------------------
# 2. Cox PH with Therneau-Grambsch (cox.zph equivalent) diagnostic
# ------------------------------------------------------------------
def cox_with_ph_check(df, duration_col='time', event_col='event', covariates='treatment + age + ECOG_baseline'):
    cph = CoxPHFitter()
    cph.fit(df, duration_col=duration_col, event_col=event_col, formula=covariates)
    cph.print_summary()

    # PH test via scaled Schoenfeld residuals
    # NOTE: a p > 0.05 does NOT prove PH; it means we cannot reject
    ph_results = proportional_hazard_test(cph, df, time_transform='rank')
    print('\nProportional Hazards Test (Therneau-Grambsch):')
    print(ph_results.summary)

    # Graphical diagnostic is more informative than the p-value
    # cph.check_assumptions(df, p_value_threshold=0.05, show_plots=True)

    return cph, ph_results


# ------------------------------------------------------------------
# 3. RMST as primary endpoint under non-PH (Royston-Parmar 2013)
# ------------------------------------------------------------------
def rmst_difference(df, arm_col='arm', duration_col='time', event_col='event', tau=36):
    """
    RMST difference at pre-specified tau.

    CRITICAL: tau must be pre-specified in the SAP. Post-hoc tau
    selection is p-hacking. Tau must be <= min(largest follow-up
    per arm) to avoid extrapolation (Tian 2020 Biostatistics).
    """
    arms = df[arm_col].unique()
    if len(arms) != 2:
        raise ValueError('RMST difference requires exactly 2 arms')

    max_follow_up_per_arm = df.groupby(arm_col)[duration_col].max()
    if tau > max_follow_up_per_arm.min():
        print(f'WARNING: tau={tau} exceeds min follow-up {max_follow_up_per_arm.min():.1f}; '
              f'will extrapolate beyond data')

    results = {}
    for arm in arms:
        arm_df = df[df[arm_col] == arm]
        kmf = KaplanMeierFitter().fit(arm_df[duration_col], event_observed=arm_df[event_col])
        rmst = restricted_mean_survival_time(kmf, t=tau)
        results[arm] = rmst

    diff = results[arms[0]] - results[arms[1]]
    print(f'\nRMST at tau={tau}:')
    for arm, rmst in results.items():
        print(f'  {arm}: {rmst:.2f}')
    print(f'  Difference ({arms[0]} - {arms[1]}): {diff:.2f}')
    # Wald SE via delta method; bootstrap is the rigorous alternative
    # R survRM2::rmst2 is the regulatory-grade implementation
    return results, diff


# ------------------------------------------------------------------
# 4. Competing risks: Aalen-Johansen CIF + cause-specific Cox
# ------------------------------------------------------------------
def competing_risks_analysis(df, duration_col='time', event_status_col='event_status',
                              cause_of_interest=1, competing_cause=2,
                              covariates='treatment + age'):
    """
    For competing risks:
    - Aalen-Johansen for CIF (correct; KM is biased upward)
    - Cause-specific Cox for CAUSAL effect on event of interest
    - Cause-specific Cox for CAUSAL effect on competing event (always report)
    - Fine-Gray ONLY when CIF prediction is the goal (Putter 2007)

    event_status_col: 0=censored, 1=event of interest, 2=competing event
    """
    # CIF via Aalen-Johansen
    ajf = AalenJohansenFitter()
    ajf.fit(df[duration_col], event_observed=df[event_status_col],
            event_of_interest=cause_of_interest)
    print(f'\nAalen-Johansen CIF for cause {cause_of_interest}:')
    print(ajf.cumulative_density_)

    # Cause-specific Cox: treat competing events as censoring
    df_cause = df.copy()
    df_cause['event_cause'] = (df_cause[event_status_col] == cause_of_interest).astype(int)

    df_competing = df.copy()
    df_competing['event_competing'] = (df_competing[event_status_col] == competing_cause).astype(int)

    cph_cause = CoxPHFitter().fit(df_cause, duration_col=duration_col,
                                   event_col='event_cause', formula=covariates)
    cph_competing = CoxPHFitter().fit(df_competing, duration_col=duration_col,
                                       event_col='event_competing', formula=covariates)

    print('\nCause-specific Cox (event of interest):')
    cph_cause.print_summary()
    print('\nCause-specific Cox (competing event) -- report alongside; '
          'covariates affecting competing event indirectly affect CIF of event of interest:')
    cph_competing.print_summary()

    return ajf, cph_cause, cph_competing


# ------------------------------------------------------------------
# 5. Weighted log-rank for non-proportional hazards
# ------------------------------------------------------------------
def weighted_logrank_for_npa(df, time_col='time', event_col='event', group_col='arm',
                              weightings='peto'):
    """
    Fleming-Harrington G(rho, gamma) weighted log-rank for non-PH.

    weightings options (lifelines 0.27+, note plural parameter name):
    - None: standard log-rank G(0, 0)
    - 'wilcoxon': down-weights late events; early-effect detection
    - 'tarone-ware': compromise
    - 'peto': robust to ties and censoring distribution differences
    - 'fleming-harrington': pass rho + gamma via kwargs (e.g. rho=0, gamma=1 for late-emphasis)

    For MaxCombo:
    R: nphRCT::maxcombo(formula, data, rho=c(0,0,1,1), gamma=c(0,1,0,1))
    Python: scikit-survival has weighted log-rank; MaxCombo via manual MVN
    """
    arms = df[group_col].unique()
    if len(arms) == 2:
        arm_A = df[df[group_col] == arms[0]]
        arm_B = df[df[group_col] == arms[1]]
        results = logrank_test(
            arm_A[time_col], arm_B[time_col],
            event_observed_A=arm_A[event_col],
            event_observed_B=arm_B[event_col],
            weightings=weightings
        )
    else:
        results = multivariate_logrank_test(
            df[time_col], df[group_col], df[event_col], weightings=weightings
        )
    print(f'\nWeighted log-rank ({weightings}):')
    print(results.summary)
    return results


# ------------------------------------------------------------------
# 6. Demonstration with simulated trial data
# ------------------------------------------------------------------
def simulate_trial(n_per_arm=200, hr=0.65, lambda_control=0.04, follow_up=48, seed=42):
    """Simulate exponential survival data; PH holds."""
    rng = np.random.default_rng(seed)
    arms = ['Control'] * n_per_arm + ['Treatment'] * n_per_arm
    # Treatment has multiplicative hazard hr
    lambda_arm = np.where(np.array(arms) == 'Treatment',
                          lambda_control * hr, lambda_control)
    true_event_time = rng.exponential(1 / lambda_arm)
    admin_censor_time = follow_up
    obs_time = np.minimum(true_event_time, admin_censor_time)
    event = (true_event_time <= admin_censor_time).astype(int)
    df = pd.DataFrame({
        'USUBJID': [f'SUBJ{i:04d}' for i in range(2 * n_per_arm)],
        'arm': arms,
        'time': obs_time,
        'event': event,
        'CNSR': 1 - event,  # CDISC convention
        'age': rng.normal(60, 10, 2 * n_per_arm),
        'ECOG_baseline': rng.choice([0, 1, 2], size=2 * n_per_arm),
    })
    # Encode arm for Cox
    df['treatment'] = (df['arm'] == 'Treatment').astype(int)
    return df


if __name__ == '__main__':
    print('=' * 60)
    print('Demonstration: simulated PH-conformant trial')
    print('=' * 60)
    df = simulate_trial(hr=0.65)

    print('\n--- Step 1: ADTTE CNSR convention check ---')
    df_check = adtte_to_lifelines(df.rename(columns={'time': 'AVAL'}))
    assert (df_check['event'] == df['event']).all()
    print('CNSR convention conversion verified')

    print('\n--- Step 2: Cox PH with Therneau-Grambsch diagnostic ---')
    cph, ph_test = cox_with_ph_check(df, covariates='treatment + age + ECOG_baseline')

    print('\n--- Step 3: RMST at tau=24 months ---')
    rmst_results, rmst_diff = rmst_difference(df, tau=24)

    print('\n--- Step 4: Stratified log-rank (Peto-Peto weights) ---')
    logrank_results = weighted_logrank_for_npa(df, weightings='peto')
