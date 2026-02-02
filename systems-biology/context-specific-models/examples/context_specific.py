'''Build tissue and condition-specific metabolic models'''

import cobra
import numpy as np


def gimme(model, expression_data, threshold_percentile=25, min_growth_fraction=0.1):
    '''Gene Inactivity Moderated by Metabolism and Expression (GIMME)

    Creates context-specific model by penalizing flux through lowly-expressed
    reactions while maintaining minimum growth.

    Args:
        expression_data: dict mapping gene_id -> normalized expression (0-1)
        threshold_percentile: Percentile below which genes considered inactive
                             25 = bottom quartile inactive (typical)
        min_growth_fraction: Minimum growth as fraction of wild-type
                            0.1 = at least 10% of max growth required

    This is the simplest context-specific algorithm. Use when:
    - Quick analysis needed
    - Expression data quality is moderate
    - Interpretability is important
    '''
    # Get wild-type growth
    wt_growth = model.optimize().objective_value

    # Calculate expression threshold
    values = list(expression_data.values())
    cutoff = np.percentile(values, threshold_percentile)

    # Create context model
    context_model = model.copy()

    # Set minimum growth
    for rxn in context_model.reactions:
        if rxn.objective_coefficient > 0:  # Biomass reaction
            rxn.lower_bound = min_growth_fraction * wt_growth

    # Constrain lowly-expressed reactions
    constrained = 0
    for rxn in context_model.reactions:
        if not rxn.genes:
            continue

        # Get expression for reaction (max of gene expressions)
        gene_expr = [expression_data.get(g.id, 0.5) for g in rxn.genes]
        rxn_expr = max(gene_expr)

        if rxn_expr < cutoff:
            # Constrain to minimal flux
            rxn.upper_bound = min(rxn.upper_bound, 0.1)
            rxn.lower_bound = max(rxn.lower_bound, -0.1)
            constrained += 1

    print(f'Constrained {constrained} reactions based on low expression')
    return context_model


def simulate_expression_data(model, active_pathways=None):
    '''Generate simulated expression data for testing

    Args:
        active_pathways: List of pathway keywords that should be highly expressed
    '''
    if active_pathways is None:
        active_pathways = ['glycolysis', 'tca', 'oxidative']

    expression = {}
    for gene in model.genes:
        # Check if gene is in active pathway
        rxns = list(gene.reactions)
        in_active = any(any(pw in rxn.name.lower() for pw in active_pathways) for rxn in rxns)

        if in_active:
            expression[gene.id] = np.random.uniform(0.7, 1.0)
        else:
            expression[gene.id] = np.random.uniform(0.0, 0.5)

    return expression


def compare_models(original, context):
    '''Compare flux distributions between models'''
    orig_sol = original.optimize()
    context_sol = context.optimize()

    comparison = {
        'original_growth': orig_sol.objective_value,
        'context_growth': context_sol.objective_value,
        'growth_ratio': context_sol.objective_value / orig_sol.objective_value
    }

    # Find reactions with changed flux
    changed = []
    for rxn_id in original.reactions.list_attr('id'):
        orig_flux = orig_sol.fluxes[rxn_id]
        context_flux = context_sol.fluxes.get(rxn_id, 0)
        if abs(orig_flux - context_flux) > 0.1:
            changed.append({
                'reaction': rxn_id,
                'original': orig_flux,
                'context': context_flux,
                'change': context_flux - orig_flux
            })

    comparison['changed_reactions'] = len(changed)
    comparison['top_changes'] = sorted(changed, key=lambda x: abs(x['change']), reverse=True)[:10]

    return comparison


if __name__ == '__main__':
    model = cobra.io.load_model('textbook')

    print('Context-Specific Model Building')
    print('=' * 50)

    # Simulate expression data
    print('\nGenerating simulated expression data...')
    expression = simulate_expression_data(model, active_pathways=['glycolysis', 'pentose'])

    # Create context model
    print('\nApplying GIMME algorithm...')
    context_model = gimme(model, expression, threshold_percentile=25)

    # Compare models
    print('\nComparing models:')
    comparison = compare_models(model, context_model)
    print(f"  Original growth: {comparison['original_growth']:.4f}")
    print(f"  Context growth: {comparison['context_growth']:.4f}")
    print(f"  Growth ratio: {comparison['growth_ratio']:.2%}")
    print(f"  Changed reactions: {comparison['changed_reactions']}")

    # Show top flux changes
    print('\nTop flux changes:')
    for change in comparison['top_changes'][:5]:
        print(f"  {change['reaction']}: {change['original']:.2f} -> {change['context']:.2f}")
