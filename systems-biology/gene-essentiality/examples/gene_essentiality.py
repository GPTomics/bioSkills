'''In-silico gene essentiality: GPR-aware knockouts, cutoff sweep, MOMA, synthetic lethality.'''
# Reference: cobrapy 0.29+ | Verify API if version differs
# The __main__ guard is required: deletion functions spawn workers that re-import this module.

import cobra
from cobra.flux_analysis import single_gene_deletion, double_gene_deletion, moma
from cobra.util.solver import linear_reaction_coefficients

CUTOFFS = (0.01, 0.02, 0.05, 0.10)   # essential = KO growth below this fraction of WT; sweep it
SL_GENE_CAP = 40                     # cap the O(n^2) double-deletion sweep


def essential_calls(model):
    '''Single-gene screen; returns the results DataFrame with a gene column and WT growth.'''
    wt_growth = model.slim_optimize()
    results = single_gene_deletion(model)
    results['gene'] = results['ids'].apply(lambda s: list(s)[0])   # ids elements are gene-id strings
    results['relative'] = results['growth'] / wt_growth
    return results, wt_growth


def confidence_split(results, cutoffs=CUTOFFS):
    '''Core essentials (essential at every cutoff) vs boundary (call depends on the cutoff).'''
    calls = {c: set(results.loc[results['relative'] < c, 'gene']) for c in cutoffs}
    core = set.intersection(*calls.values())
    boundary = set.union(*calls.values()) - core
    return core, boundary


def main():
    model = cobra.io.load_model('textbook')
    results, wt = essential_calls(model)
    print(f'Wild-type growth: {wt:.4f} /h')

    core, boundary = confidence_split(results)
    print(f'High-confidence essential (all cutoffs): {len(core)}')
    print(f'Low-confidence, cutoff-dependent        : {len(boundary)}  {sorted(boundary)}')

    print('\n=== FBA vs MOMA on a fresh knockout (immediate mutant) ===')
    biomass = list(linear_reaction_coefficients(model))[0]   # the objective (biomass) reaction
    wt_sol = model.optimize()
    gene = 'b2276'
    with model:
        model.genes.get_by_id(gene).knock_out()
        fba_growth = model.slim_optimize()
        # moma().objective_value is the QP adjustment objective, NOT growth; read the biomass flux.
        moma_growth = moma(model, solution=wt_sol, linear=True).fluxes[biomass.id]
    print(f'{gene} KO: FBA re-optimized = {fba_growth:.4f}, MOMA immediate = {moma_growth:.4f} /h')

    print('\n=== Synthetic-lethal pairs among viable genes ===')
    viable = list(results.loc[results['relative'] > CUTOFFS[0], 'gene'])[:SL_GENE_CAP]
    dbl = double_gene_deletion(model, gene_list1=viable, gene_list2=viable)
    dbl['n'] = dbl['ids'].apply(len)
    sl = dbl[(dbl['n'] == 2) & (dbl['growth'] / wt < CUTOFFS[0])]
    print(f'Synthetic-lethal pairs: {len(sl)}')
    for s in sl['ids'].head(5):
        print('  ' + ' + '.join(sorted(s)))


if __name__ == '__main__':
    main()
