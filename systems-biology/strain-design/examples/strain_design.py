'''Growth-coupled strain design with StrainDesign (OptKnock) and verification via the production envelope.

A growth-coupled design makes product secretion obligatory for growth, so selection in the bioreactor
maintains production. OptKnock is optimistic (assumes the cell cooperates at its growth optimum);
verify the coupling with the production envelope and prefer RobustKnock when a design must be robust.
'''
# Reference: StrainDesign 1.15+, cobrapy 0.29+ | Verify API if version differs
# The __main__ guard is required: solver/parallel steps may re-import this module.

import numpy as np
import cobra
import straindesign as sd

BIOMASS = 'Biomass_Ecoli_core'   # the textbook model's actual biomass reaction id
PRODUCT = 'EX_succ_e'            # target: anaerobic succinate secretion (a classic growth-coupling case)
MIN_GROWTH = 0.1                 # keep the designed strain viable (anaerobic growth is lower)
MAX_KO = 4                       # intervention budget (MILP cost); more knockouts = harder search
KO_MARKER = -1.0                 # StrainDesign encodes a knockout as reaction_id -> -1.0


def min_product_at_max_growth(model, frac=0.99):
    '''Pin growth near its max, then minimize product flux. 0 => the cell can grow making none
    (not growth-coupled); nonzero => product secretion is obligatory at high growth (coupled).'''
    with model:
        mu = model.slim_optimize()
        if not np.isfinite(mu):
            return float('nan')
        model.reactions.get_by_id(BIOMASS).lower_bound = frac * mu
        model.objective = PRODUCT
        model.objective_direction = 'min'
        return model.slim_optimize()


def main():
    model = cobra.io.load_model('textbook')
    model.reactions.EX_o2_e.lower_bound = 0   # anaerobic: fermentation makes succinate coupling natural

    optknock = sd.SDModule(
        model, sd.OPTKNOCK,
        inner_objective=BIOMASS,               # the cell maximizes growth
        outer_objective=PRODUCT,               # the engineer maximizes succinate
        constraints=[f'{BIOMASS} >= {MIN_GROWTH}'],
    )

    print('=== OptKnock search (MILP; bounded to keep it fast) ===')
    solutions = sd.compute_strain_designs(
        model, sd_modules=[optknock],
        max_cost=MAX_KO, max_solutions=3, solver='glpk', time_limit=120,
    )
    designs = [[rid for rid, mark in d.items() if mark == KO_MARKER] for d in (solutions.reaction_sd or [])]
    designs = [d for d in designs if d]
    print(f'status: {getattr(solutions, "status", "?")}  designs found: {len(designs)}')
    for i, kos in enumerate(designs[:3], 1):
        print(f'  design {i}: knock out {kos}')

    print('\n=== Verify growth-coupling with the production envelope ===')
    print(f'wild-type minimum {PRODUCT} at max growth: {min_product_at_max_growth(model):.3f} (0 = free to make none)')
    if designs:
        with model:
            for rid in designs[0]:
                model.reactions.get_by_id(rid).knock_out()
            coupled = min_product_at_max_growth(model)
        print(f'design 1 minimum {PRODUCT} at max growth: {coupled:.3f} (nonzero = growth-coupled)')


if __name__ == '__main__':
    main()
