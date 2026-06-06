'''Under competing risks, 1 - Kaplan-Meier overestimates incidence; use the CIF.

Runs end-to-end on synthetic data with two competing event types. Treating the
competing event as censoring and reporting 1 - KM inflates the estimated incidence
of the event of interest relative to the Aalen-Johansen cumulative incidence
function (CIF). The gap grows with the competing-event rate.
'''
# Reference: lifelines 0.30+, numpy 1.26+ | Verify API if version differs

import numpy as np
from lifelines import KaplanMeierFitter, AalenJohansenFitter

rng = np.random.default_rng(0)
n = 3000
t1 = rng.exponential(4.0, n)                  # time to event of interest (cause 1)
t2 = rng.exponential(1.5, n)                  # time to competing event (cause 2), more frequent
admin = rng.uniform(0, 8, n)                  # administrative censoring
time = np.minimum.reduce([t1, t2, admin])
cause = np.where((t1 <= t2) & (t1 <= admin), 1, np.where((t2 < t1) & (t2 <= admin), 2, 0))  # 0 = censored

horizon = 3.0

# WRONG: treat the competing event (cause 2) as censoring, report 1 - KM for cause 1.
km = KaplanMeierFitter().fit(time, event_observed=(cause == 1))
km_incidence = 1 - float(km.predict(horizon))

# RIGHT: Aalen-Johansen CIF for cause 1, treating cause 2 as a competing risk.
ajf = AalenJohansenFitter().fit(time, cause, event_of_interest=1)
cif = float(ajf.predict(horizon))

print(f'Competing event (cause 2) is {np.mean(cause == 2):.0%} of subjects')
print(f'1 - KM incidence of cause 1 at t={horizon} (WRONG): {km_incidence:.3f}')
print(f'Aalen-Johansen CIF of cause 1 at t={horizon} (RIGHT): {cif:.3f}')
print(f'1 - KM overestimates incidence by {km_incidence - cif:.3f} (it ignores competing risk)')
