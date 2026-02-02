'''Flux balance analysis with COBRApy'''

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba

# BiGG models available at http://bigg.ucsd.edu
# Common models:
# - 'textbook': E. coli core (95 reactions) - good for learning
# - 'iJO1366': Full E. coli K-12 (2583 reactions)
# - 'iMM904': S. cerevisiae (1577 reactions)
# - 'Recon3D': Human metabolism (13543 reactions)

model = cobra.io.load_model('textbook')

print('=== Basic FBA ===')
solution = model.optimize()

# Growth rate interpretation:
# >0.8 h^-1: Fast growth expected (rich media conditions)
# 0.3-0.8 h^-1: Moderate growth
# <0.3 h^-1: Slow growth or stress conditions
# 0: No growth (lethal condition or missing nutrients)
print(f'Growth rate: {solution.objective_value:.4f} h^-1')
print(f'Solver status: {solution.status}')


print('\n=== Key Fluxes ===')
key_reactions = ['GLNS', 'GLUDy', 'PFK', 'PYK', 'ATPS4r']
for rxn_id in key_reactions:
    if rxn_id in model.reactions:
        flux = solution.fluxes[rxn_id]
        print(f'{rxn_id}: {flux:.4f} mmol/gDW/h')


print('\n=== Flux Variability Analysis ===')
# FVA at 90% optimal growth reveals flux flexibility
# fraction_of_optimum=0.9: allows solutions within 10% of max growth
# This identifies essential reactions vs flexible pathways
fva = flux_variability_analysis(model, fraction_of_optimum=0.9,
                                 reaction_list=key_reactions)
print(fva[['minimum', 'maximum']])


print('\n=== Compare Carbon Sources ===')
carbon_sources = {
    'EX_glc__D_e': 'Glucose',
    'EX_ac_e': 'Acetate',
    'EX_succ_e': 'Succinate'
}

for ex_id, name in carbon_sources.items():
    with model:  # Context manager resets changes
        # Close all carbon sources
        for rxn in model.exchanges:
            if any(met.formula and 'C' in met.formula for met in rxn.metabolites):
                rxn.lower_bound = 0

        # Open specific carbon source
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -10  # mmol/gDW/h uptake

            sol = model.optimize()
            print(f'{name}: {sol.objective_value:.4f} h^-1')


print('\n=== Parsimonious FBA ===')
# pFBA minimizes total flux while achieving optimal growth
# Produces more biologically realistic flux distributions
pfba_sol = pfba(model)
fba_total = sum(abs(solution.fluxes))
pfba_total = sum(abs(pfba_sol.fluxes))
print(f'Standard FBA total flux: {fba_total:.1f} mmol/gDW/h')
print(f'pFBA total flux: {pfba_total:.1f} mmol/gDW/h')
print(f'Reduction: {(1 - pfba_total/fba_total)*100:.1f}%')
