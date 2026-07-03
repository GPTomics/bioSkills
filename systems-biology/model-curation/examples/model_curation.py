'''Curation checks COBRApy can run without memote: energy-generating cycles, mass/charge balance,
dead ends. memote (consistency + annotation scoring) is shown as the CLI it wraps.
'''
# Reference: cobrapy 0.29+, memote 0.17+ | Verify API if version differs

import numpy as np
import cobra
from cobra.util.solver import linear_reaction_coefficients

BAL_TOL = 1e-6   # element/charge sums below this are treated as balanced

# Moiety-conserving dissipation reactions per energy currency (charged -> discharged). Maximizing
# each with ALL uptake closed must give ~0; a positive flux is an EGC charging that currency for
# free. EGCs are NOT ATP-only (Fritzemeier 2017), so sweep every currency present in the model.
# Genome-scale models should also test GTP/CTP/UTP and quinones (q8h2)/FMNH2. Proton-motive-force
# cycles are subtler (they need both-direction closure) -- memote's dedicated EGC test handles them.
DISSIPATION = {
    'atp': {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'pi_c': 1, 'h_c': 1},
    'nadh': {'nadh_c': -1, 'nad_c': 1, 'h_c': 1},
    'nadph': {'nadph_c': -1, 'nadp_c': 1, 'h_c': 1},
    'fadh2': {'fadh2_c': -1, 'fad_c': 1, 'h_c': 2},
}


def energy_generating_cycles(model, dissipations=DISSIPATION):
    '''Max free-charging flux per energy currency with ALL uptake closed and the NGAM floor removed.
    ~0/infeasible per currency is correct; a positive finite value is an erroneous energy-generating cycle.'''
    results = {}
    with model:
        for ex in model.exchanges:
            ex.lower_bound = 0
        if 'ATPM' in model.reactions:
            model.reactions.get_by_id('ATPM').lower_bound = 0   # NGAM floor would otherwise force infeasibility
        for name, stoich in dissipations.items():
            if any(m not in model.metabolites for m in stoich):
                continue                                        # currency not in this model
            with model:
                rxn = cobra.Reaction(f'EGC_{name}')
                rxn.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
                rxn.bounds = (0, 1000)
                model.add_reactions([rxn])
                model.objective = rxn
                results[name] = model.slim_optimize()
    return results


def imbalance(reaction):
    '''Element imbalance dict and charge imbalance for a reaction.'''
    mass = {}
    charge = 0
    for met, coef in reaction.metabolites.items():
        if met.formula:
            for element, n in met.elements.items():
                mass[element] = mass.get(element, 0) + coef * n
        if met.charge is not None:
            charge += coef * met.charge
    return {e: v for e, v in mass.items() if abs(v) > BAL_TOL}, charge


def unbalanced_reactions(model):
    # exchange/demand/sink and the biomass pseudo-reaction are intentionally unbalanced -> skip them
    pseudo = set(model.boundary) | set(linear_reaction_coefficients(model))
    out = []
    for rxn in model.reactions:
        if rxn in pseudo:
            continue
        mass, charge = imbalance(rxn)
        if mass or abs(charge) > BAL_TOL:
            out.append((rxn.id, mass, charge))
    return out


def dead_ends(model):
    out = []
    for met in model.metabolites:
        produced = any(r.metabolites[met] > 0 for r in met.reactions)
        consumed = any(r.metabolites[met] < 0 for r in met.reactions)
        if not (produced and consumed):
            out.append(met.id)
    return out


def main():
    # memote (consistency + annotation score) is a CLI step:
    #   memote report snapshot model.xml --filename report.html
    # The score measures well-formedness, NOT biological correctness; read which tests fail.
    model = cobra.io.load_model('textbook')

    print('=== Energy-generating cycle test per currency (the defect a high score hides) ===')
    for currency, flux in energy_generating_cycles(model).items():
        cycle = np.isfinite(flux) and flux > 1e-3   # infeasible/nan or ~0 => no cycle (correct)
        print(f'  {currency:7s}: {flux:.4f}  ->  {"ENERGY-GENERATING CYCLE" if cycle else "OK"}')

    print('\n=== Mass/charge balance ===')
    ub = unbalanced_reactions(model)
    print(f'Unbalanced non-exchange reactions: {len(ub)}')
    for rid, mass, charge in ub[:3]:
        print(f'  {rid}: mass {mass} charge {charge}')

    print('\n=== Dead-end metabolites ===')
    de = dead_ends(model)
    print(f'Dead-end metabolites: {len(de)}  (network gaps or wrong stoichiometry)')

    print('\nReminder: these consistency checks do NOT confirm the model predicts biology;')
    print('validate against measured growth/essentiality on the matched medium separately.')


if __name__ == '__main__':
    main()
