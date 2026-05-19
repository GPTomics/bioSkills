"""
Missing data sensitivity analyses for clinical trials -- Python worked examples.

Reference: scikit-learn 1.4+ | Verify API if version differs
Reference: statsmodels 0.14+ | Verify API if version differs

NOTE: For confirmatory regulatory work, use R packages mmrm and rbmi.
Python is suitable for exploratory analyses and pipeline development.

Covers:
1. Multiple imputation with Rubin's rules pooling
2. Tipping-point delta-adjustment (Permutt 2016)
3. Simple J2R-style imputation pattern
4. Pre-specified SAP fallback hierarchy
"""

import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import statsmodels.formula.api as smf


# ------------------------------------------------------------------
# 1. Examine missingness pattern first -- the DS domain analog
# ------------------------------------------------------------------
def characterise_dropout(df, arm_col='arm', dropout_col='dropout_reason'):
    """
    Tabulate dropout reasons by arm. If pattern differs by arm and
    reasons are clinically informative, MAR is suspect.
    """
    pattern = pd.crosstab(df[dropout_col], df[arm_col],
                          margins=True, margins_name='Total')
    print('Dropout reasons by arm:')
    print(pattern)

    # Test if dropout rates differ by arm
    overall_dropout = df.groupby(arm_col)[dropout_col].apply(
        lambda x: (~x.isin(['COMPLETED', 'Completed'])).mean()
    )
    print(f'\nOverall dropout rate by arm:\n{overall_dropout}')

    # 10pp threshold: NRC 2010 implicit; dropout asymmetry of this magnitude
    # is the canonical signal that MAR is implausible (Mallinckrodt 2008;
    # Aducanumab AdCom 2021 precedent for differential-dropout concerns).
    if abs(overall_dropout.iloc[0] - overall_dropout.iloc[1]) > 0.10:
        print('\nWARNING: dropout differs by arm > 10pp; MAR is suspect.')
        print('Consider treatment-policy estimand with reference-based MI as primary.')

    return pattern


# ------------------------------------------------------------------
# 2. Multiple imputation with Rubin's rules (MAR primary)
# ------------------------------------------------------------------
def multiple_imputation_logit(df, formula, numeric_cols, n_imputations=20):
    """
    MI with sklearn IterativeImputer + Rubin's rules pooling.

    CRITICAL: sample_posterior=True only works with BayesianRidge.
    For confirmatory work, use R rbmi or mice.
    """
    # Verify estimator
    imputer = IterativeImputer(max_iter=10, random_state=0, sample_posterior=True)
    # imputer.estimator is BayesianRidge by default -- do NOT change

    results = []
    for i in range(n_imputations):
        imputer.set_params(random_state=i)
        imputed_numeric = pd.DataFrame(
            imputer.fit_transform(df[numeric_cols]),
            columns=numeric_cols,
            index=df.index
        )
        # Reattach non-imputed categorical columns
        imputed = imputed_numeric.copy()
        for col in df.columns:
            if col not in numeric_cols:
                imputed[col] = df[col].values

        model = smf.logit(formula, data=imputed).fit(disp=0)
        results.append({
            'coef': model.params.iloc[1],
            'se': model.bse.iloc[1]
        })

    # Rubin's rules pooling
    pooled_coef = np.mean([r['coef'] for r in results])
    within_var = np.mean([r['se'] ** 2 for r in results])
    between_var = np.var([r['coef'] for r in results], ddof=1)
    total_var = within_var + (1 + 1 / n_imputations) * between_var
    pooled_se = np.sqrt(total_var)

    pooled_or = np.exp(pooled_coef)
    z_score = pooled_coef / pooled_se
    # FMI: fraction of missing info
    fmi = (1 + 1 / n_imputations) * between_var / total_var
    # Adequate m: m >= 100 * FMI (von Hippel 2020)
    adequate_m = int(np.ceil(100 * fmi))

    print(f'Pooled coefficient: {pooled_coef:.4f}')
    print(f'Pooled SE: {pooled_se:.4f}')
    print(f'Pooled OR: {pooled_or:.3f}')
    print(f'95% CI: ({np.exp(pooled_coef - 1.96 * pooled_se):.3f}, '
          f'{np.exp(pooled_coef + 1.96 * pooled_se):.3f})')
    print(f'FMI: {fmi:.3f}; adequate m >= {adequate_m}; used m = {n_imputations}')
    if n_imputations < adequate_m:
        print('WARNING: m insufficient; consider increasing per von Hippel 2020')

    # NOTE: Rubin's rules variance shown here is the Cro 2019 information-anchored variance.
    # For reference-based MI (J2R/CR/CIR), Bartlett 2021 + Wolbers 2022 recommend
    # also reporting the frequentist variance via CMI+jackknife (R rbmi::analyse with
    # method_condmean). Active regulatory debate; report both for safety. The
    # frequentist variance is not feasibly computed via sklearn; use R rbmi for
    # confirmatory reference-based MI work.

    return pooled_coef, pooled_se, pooled_or, fmi


# ------------------------------------------------------------------
# 3. Tipping-point analysis (Permutt 2016)
# ------------------------------------------------------------------
def tipping_point_analysis(df, formula, outcome_col, arm_col,
                            arm_to_adjust='Active', delta_range=None,
                            n_imputations=20):
    """
    Apply delta only to active-arm imputed values; find minimum delta
    that flips p > 0.05.

    delta_range: list of delta values to scan (in units of residual SD)
    """
    if delta_range is None:
        delta_range = np.arange(0, 21, 1)

    # Fit MAR model to get residual SD
    complete_df = df.dropna(subset=[outcome_col])
    base_model = smf.logit(formula, data=complete_df).fit(disp=0)
    # For binary outcomes, "residual SD" is conceptually different;
    # for continuous outcomes use the actual residual SD

    p_values = []
    for delta in delta_range:
        # Impute under MAR, then add delta to active-arm imputed values only
        df_adj = df.copy()
        missing_mask = df_adj[outcome_col].isna()
        active_mask = df_adj[arm_col] == arm_to_adjust

        # Simple imputation for demo; real implementation uses rbmi delta_template
        imputer = IterativeImputer(random_state=0, sample_posterior=True)
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        imputed_numeric = pd.DataFrame(
            imputer.fit_transform(df_adj[numeric_cols]),
            columns=numeric_cols, index=df_adj.index
        )
        # Apply delta to active-arm imputed (missing) outcomes
        imputed_numeric.loc[missing_mask & active_mask, outcome_col] += delta

        for col in df_adj.columns:
            if col not in numeric_cols:
                imputed_numeric[col] = df_adj[col].values

        model = smf.logit(formula, data=imputed_numeric).fit(disp=0)
        # Treatment coefficient p-value
        p = model.pvalues.iloc[1]
        p_values.append(p)

    # Report
    print(f'\nTipping-point analysis: delta applied to {arm_to_adjust} arm')
    print(f'{"delta":>8} {"p-value":>10}')
    for d, p in zip(delta_range, p_values):
        print(f'{d:>8.1f} {p:>10.4f}')

    # Find tipping point
    significant = [p < 0.05 for p in p_values]
    if all(significant):
        print(f'\nResult ROBUST: significant across all delta up to {delta_range[-1]}')
        tipping_delta = float('inf')
    elif not any(significant):
        print(f'\nResult NOT ROBUST: not significant even at delta=0')
        tipping_delta = 0
    else:
        # Find smallest delta where p exceeds 0.05
        for i, sig in enumerate(significant):
            if not sig:
                tipping_delta = delta_range[i]
                break
        print(f'\nTipping delta: {tipping_delta}')
        print(f'Judge clinical plausibility: is delta={tipping_delta} larger than MCID?')

    return list(zip(delta_range, p_values)), tipping_delta


# ------------------------------------------------------------------
# 4. SAP fallback hierarchy demonstration
# ------------------------------------------------------------------
def sap_fallback_hierarchy_demo():
    """
    The pre-specified MMRM convergence fallback per industry SAP:
    1. UN with KR
    2. UN with Satterthwaite
    3. Heterogeneous Toeplitz
    4. AR(1) with heterogeneous variances
    5. CS with heterogeneous variances

    Document deviation if invoked. Python implementation here is
    schematic; use R `mmrm` for production.
    """
    print('MMRM SAP fallback hierarchy:')
    hierarchy = [
        ('UN', 'Kenward-Roger-Linear', 'preferred; matches SAS PROC MIXED'),
        ('UN', 'Satterthwaite', 'fallback if KR fails'),
        ('Toeplitz heterogeneous', 'Kenward-Roger', 'k+1 parameters; structured'),
        ('AR(1) heterogeneous', 'Kenward-Roger', 'imposes stationarity'),
        ('CS heterogeneous', 'Kenward-Roger', 'last resort; assumes equal correlation'),
    ]
    for i, (cov, df_method, note) in enumerate(hierarchy, 1):
        print(f'  {i}. {cov} covariance, {df_method} DF -- {note}')
    print('\nDocument deviation in CSR if invoked.')


# ------------------------------------------------------------------
# 5. Simulation: trial data with differential dropout
# ------------------------------------------------------------------
def simulate_trial_with_dropout(n_per_arm=150, dropout_active=0.30, dropout_placebo=0.15, seed=42):
    rng = np.random.default_rng(seed)
    arms = ['Active'] * n_per_arm + ['Placebo'] * n_per_arm
    age = rng.normal(60, 10, 2 * n_per_arm)
    baseline = rng.normal(50, 10, 2 * n_per_arm)

    # True effect: active reduces outcome by 5 units
    true_change = np.where(np.array(arms) == 'Active',
                            rng.normal(-8, 8, 2 * n_per_arm),
                            rng.normal(-3, 8, 2 * n_per_arm))

    # Dropout: active higher rate, predominantly due to AEs (informative)
    dropout_prob = np.where(np.array(arms) == 'Active', dropout_active, dropout_placebo)
    dropped = rng.random(2 * n_per_arm) < dropout_prob

    # Dropout reason: differential
    dropout_reason = np.full(2 * n_per_arm, 'COMPLETED', dtype=object)
    for i in range(2 * n_per_arm):
        if dropped[i]:
            if arms[i] == 'Active':
                dropout_reason[i] = rng.choice(['AE', 'WITHDREW', 'LFU'], p=[0.7, 0.2, 0.1])
            else:
                dropout_reason[i] = rng.choice(['AE', 'WITHDREW', 'LFU'], p=[0.2, 0.4, 0.4])

    outcome_obs = np.where(dropped, np.nan, true_change)

    df = pd.DataFrame({
        'USUBJID': [f'SUBJ{i:04d}' for i in range(2 * n_per_arm)],
        'arm': arms,
        'arm_numeric': (np.array(arms) == 'Active').astype(int),
        'age': age,
        'baseline': baseline,
        'change': outcome_obs,
        'dropout_reason': dropout_reason,
        # Outcome dichotomised for logit demo
        'outcome': np.where(np.isnan(outcome_obs), np.nan, (outcome_obs < -5).astype(float)),
    })
    return df


if __name__ == '__main__':
    print('=' * 60)
    print('Demonstration: differential-dropout trial')
    print('=' * 60)
    df = simulate_trial_with_dropout()

    print('\n--- Step 1: Characterise dropout pattern ---')
    characterise_dropout(df)

    print('\n--- Step 2: MI primary (assuming MAR despite WARNING) ---')
    pooled = multiple_imputation_logit(
        df, formula='outcome ~ arm_numeric + age + baseline',
        numeric_cols=['age', 'baseline', 'outcome'],
        n_imputations=20
    )

    print('\n--- Step 3: Tipping-point analysis ---')
    tp_results, tip_delta = tipping_point_analysis(
        df, formula='outcome ~ arm_numeric + age + baseline',
        outcome_col='outcome', arm_col='arm',
        arm_to_adjust='Active', delta_range=np.arange(0, 0.51, 0.05),
        n_imputations=10
    )

    print('\n--- Step 4: SAP fallback hierarchy reference ---')
    sap_fallback_hierarchy_demo()
