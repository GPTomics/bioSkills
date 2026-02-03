'''Cox proportional hazards regression for survival modeling'''

import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt

clinical = pd.read_csv('clinical.csv')
expr = pd.read_csv('expression.csv', index_col=0)

# Prepare Cox data with clinical + omics features
cox_df = pd.DataFrame({
    'time': clinical['survival_time'],
    'event': clinical['event'],
    'age': clinical['age'],
    'stage': clinical['stage'],
    'BRCA1': expr.loc['BRCA1'],
    'TP53': expr.loc['TP53']
})
cox_df = cox_df.dropna()
print(f'Cox data: {len(cox_df)} patients, {cox_df["event"].sum()} events')

# penalizer=0.1: L2 regularization, use higher values for more features
cph = CoxPHFitter(penalizer=0.1)
cph.fit(cox_df, duration_col='time', event_col='event')

print('\n=== Cox PH Model Summary ===')
cph.print_summary()

print(f'\nConcordance index: {cph.concordance_index_:.3f}')
# C-index: 0.5=random, 0.7+=useful, 0.8+=good

hr_df = cph.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']]
hr_df.columns = ['HR', 'HR_lower', 'HR_upper', 'p_value']
hr_df.to_csv('cox_hazard_ratios.csv')
print('\nSaved hazard ratios to cox_hazard_ratios.csv')

# Check proportional hazards assumption (uncomment for full check)
# cph.check_assumptions(cox_df, show_plots=True)

# Risk stratification using Cox partial hazard
risk_scores = cph.predict_partial_hazard(cox_df)
cox_df['risk_score'] = risk_scores
cox_df['risk_group'] = (risk_scores > risk_scores.median()).map({True: 'high', False: 'low'})

# KM plot for risk groups
fig, ax = plt.subplots(figsize=(8, 6))
for group, color in [('high', 'red'), ('low', 'blue')]:
    mask = cox_df['risk_group'] == group
    kmf = KaplanMeierFitter()
    kmf.fit(cox_df.loc[mask, 'time'], event_observed=cox_df.loc[mask, 'event'], label=f'{group} risk')
    kmf.plot_survival_function(ax=ax, color=color)

high = cox_df[cox_df['risk_group'] == 'high']
low = cox_df[cox_df['risk_group'] == 'low']
lr = logrank_test(high['time'], low['time'], event_observed_A=high['event'], event_observed_B=low['event'])

ax.set_xlabel('Time (months)')
ax.set_ylabel('Survival probability')
ax.set_title(f'Cox risk stratification (C-index={cph.concordance_index_:.3f})\nLog-rank p={lr.p_value:.4e}')
plt.tight_layout()
plt.savefig('cox_risk_stratification.png', dpi=150)
print('Saved cox_risk_stratification.png')
