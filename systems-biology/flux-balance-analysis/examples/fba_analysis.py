'''FBA, FVA, pFBA, and flux sampling with COBRApy, framed around alternate optima.'''
# Reference: cobrapy 0.29+ | Verify API if version differs
# The __main__ guard is required: FVA and sampling spawn worker processes, which re-import
# this module; without the guard, spawn-based platforms (macOS, Windows) recurse.

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba, loopless_solution
from cobra.sampling import sample

GLUCOSE_UPTAKE = 10   # mmol/gDW/h; standard E. coli aerobic glucose uptake (iJO1366)
FLUX_TOL = 1e-6       # fluxes below this are treated as zero (solver feasibility tolerance)


def main():
    model = cobra.io.load_model('textbook')   # E. coli core (e_coli_core)

    print('=== Basic FBA ===')
    solution = model.optimize()
    # A growth rate is interpretable only against the stated medium and biomass reaction;
    # compare to a measured doubling time (mu = ln2 / t_double), not a fixed fast/slow scale.
    print(f'Growth: {solution.objective_value:.4f} /h  status: {solution.status}')

    print('\n=== Alternate optima: FVA says which single fluxes are trustworthy ===')
    key_reactions = ['PFK', 'PYK', 'ATPS4r', 'PGI', 'FUM']
    fva = flux_variability_analysis(model, reaction_list=key_reactions, fraction_of_optimum=1.0)
    fva['range'] = fva['maximum'] - fva['minimum']
    fva['determined'] = fva['range'].abs() < FLUX_TOL
    print(fva[['minimum', 'maximum', 'determined']].to_string())
    print('A wide range means the plain-FBA value for that reaction was one arbitrary vertex.')

    print('\n=== Compare carbon sources on a defined minimal medium ===')
    minimal = {'EX_o2_e', 'EX_h2o_e', 'EX_h_e', 'EX_nh4_e', 'EX_pi_e', 'EX_so4_e', 'EX_k_e', 'EX_co2_e'}
    for ex_id, name in [('EX_glc__D_e', 'glucose'), ('EX_ac_e', 'acetate'), ('EX_succ_e', 'succinate')]:
        with model:
            for rxn in model.exchanges:
                rxn.lower_bound = 0
            for m in minimal:
                if m in model.reactions:
                    model.reactions.get_by_id(m).lower_bound = -1000
            if ex_id in model.reactions:
                model.reactions.get_by_id(ex_id).lower_bound = -GLUCOSE_UPTAKE
            print(f'{name}: growth = {model.slim_optimize():.4f} /h')

    print('\n=== pFBA and loopless: two principled single representatives ===')
    pfba_sol = pfba(model)
    loopless = loopless_solution(model)
    print(f'FBA total flux     : {solution.fluxes.abs().sum():.1f}')
    print(f'pFBA total flux    : {pfba_sol.fluxes.abs().sum():.1f} (minimal enzyme-cost proxy)')
    print(f'loopless total flux: {loopless.fluxes.abs().sum():.1f} (internal cycles removed)')

    print('\n=== Flux sampling: the distribution, not one optimum ===')
    samples = sample(model, n=500, method='optgp', thinning=100, seed=1)
    print(samples[['PFK', 'PYK', 'ATPS4r']].describe().loc[['mean', '50%', 'std']].to_string())
    print('Check that multiple seeds agree (mixing) before trusting these distributions.')


if __name__ == '__main__':
    main()
