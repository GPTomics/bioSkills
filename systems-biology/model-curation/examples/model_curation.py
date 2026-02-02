'''Validate and curate metabolic models'''

import cobra


def find_dead_end_metabolites(model):
    '''Find metabolites that cannot be produced or consumed

    Dead-end metabolites indicate:
    - Missing reactions in the network
    - Incomplete pathways
    - Incorrectly specified compartments
    '''
    dead_ends = []
    for met in model.metabolites:
        producing = [r for r in met.reactions if r.get_coefficient(met) > 0]
        consuming = [r for r in met.reactions if r.get_coefficient(met) < 0]

        if not producing:
            dead_ends.append({'metabolite': met.id, 'issue': 'no producers'})
        elif not consuming:
            dead_ends.append({'metabolite': met.id, 'issue': 'no consumers'})

    return dead_ends


def check_reaction_balance(reaction):
    '''Check mass and charge balance for a reaction

    Common causes of imbalance:
    - Missing protons (H+)
    - Wrong stoichiometry
    - Missing water molecules
    '''
    mass_balance = {}
    charge_balance = 0

    for met, coef in reaction.metabolites.items():
        if met.formula:
            for element, count in met.elements.items():
                mass_balance[element] = mass_balance.get(element, 0) + coef * count
        if met.charge is not None:
            charge_balance += coef * met.charge

    mass_imbalance = {k: v for k, v in mass_balance.items() if abs(v) > 1e-6}

    return {
        'mass_balanced': len(mass_imbalance) == 0,
        'charge_balanced': abs(charge_balance) < 1e-6,
        'mass_imbalance': mass_imbalance,
        'charge_imbalance': charge_balance
    }


def find_unbalanced_reactions(model):
    '''Find all unbalanced reactions in the model'''
    unbalanced = []
    for rxn in model.reactions:
        if rxn in model.exchanges:
            continue  # Skip exchange reactions
        result = check_reaction_balance(rxn)
        if not result['mass_balanced'] or not result['charge_balanced']:
            unbalanced.append({
                'reaction': rxn.id,
                'equation': rxn.reaction,
                **result
            })
    return unbalanced


def identify_orphan_reactions(model):
    '''Find reactions without gene associations

    Acceptable orphans:
    - Exchange reactions (boundary)
    - Spontaneous reactions (non-enzymatic)
    - Transport (may be unannotated)
    '''
    orphans = {'exchange': [], 'transport': [], 'enzymatic': []}

    for rxn in model.reactions:
        if not rxn.genes:
            if rxn in model.exchanges:
                orphans['exchange'].append(rxn.id)
            elif 'transport' in rxn.name.lower() or any(c in rxn.id.lower() for c in ['_t', 'abc', 'pts']):
                orphans['transport'].append(rxn.id)
            else:
                orphans['enzymatic'].append(rxn.id)

    return orphans


def model_summary(model):
    '''Generate comprehensive model summary'''
    # Basic counts
    summary = {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes),
        'exchanges': len(model.exchanges),
    }

    # Dead-ends
    dead_ends = find_dead_end_metabolites(model)
    summary['dead_end_metabolites'] = len(dead_ends)

    # Orphans
    orphans = identify_orphan_reactions(model)
    summary['orphan_enzymatic'] = len(orphans['enzymatic'])

    # Growth test
    try:
        sol = model.optimize()
        summary['can_grow'] = sol.status == 'optimal' and sol.objective_value > 0.001
        summary['growth_rate'] = sol.objective_value if sol.status == 'optimal' else 0
    except:
        summary['can_grow'] = False
        summary['growth_rate'] = 0

    return summary


if __name__ == '__main__':
    model = cobra.io.load_model('textbook')

    print('Model Curation Report')
    print('=' * 50)

    # Summary statistics
    summary = model_summary(model)
    print('\nModel Statistics:')
    for key, value in summary.items():
        if isinstance(value, float):
            print(f'  {key}: {value:.4f}')
        else:
            print(f'  {key}: {value}')

    # Dead-end metabolites
    print('\nDead-End Metabolites:')
    dead_ends = find_dead_end_metabolites(model)
    if dead_ends:
        for de in dead_ends[:5]:
            print(f"  {de['metabolite']}: {de['issue']}")
        if len(dead_ends) > 5:
            print(f'  ... and {len(dead_ends) - 5} more')
    else:
        print('  None found')

    # Unbalanced reactions
    print('\nUnbalanced Reactions:')
    unbalanced = find_unbalanced_reactions(model)
    if unbalanced:
        for ub in unbalanced[:3]:
            print(f"  {ub['reaction']}: {ub['mass_imbalance']}")
        if len(unbalanced) > 3:
            print(f'  ... and {len(unbalanced) - 3} more')
    else:
        print('  All reactions balanced')

    # Quality assessment
    print('\nQuality Assessment:')
    if summary['can_grow']:
        print('  [PASS] Model produces biomass')
    else:
        print('  [FAIL] Model cannot grow')

    if summary['dead_end_metabolites'] == 0:
        print('  [PASS] No dead-end metabolites')
    else:
        print(f"  [WARN] {summary['dead_end_metabolites']} dead-end metabolites")

    orphan_frac = summary['orphan_enzymatic'] / summary['reactions']
    if orphan_frac < 0.1:
        print(f'  [PASS] Orphan reactions: {orphan_frac:.1%}')
    else:
        print(f'  [WARN] High orphan fraction: {orphan_frac:.1%}')
