'''Context-specific extraction: GPR mapping, an illustrative GIMME-style stub, and the threshold sweep.

For real extraction use troppo (GIMME/iMAT/tINIT/FASTCORE) or corda (CORDA) in Python, or the
MATLAB COBRA Toolbox / RAVEN. The stub here only illustrates the objective-protected pruning shape.
'''
# Reference: cobrapy 0.29+, numpy 1.26+ | Verify API if version differs

import cobra
import numpy as np

SEED = 1
LOW_QUANTILES = (0.10, 0.25, 0.50)   # threshold strategy dominates the model; always sweep it
GROWTH_FLOOR = 0.1                    # GIMME-family protects the objective at >= this fraction of WT


def reaction_activity(rxn, gene_expr, default=0.5):
    '''min over AND (complex, limiting subunit), max over OR (isozyme). Simplified OR aggregation.'''
    if not rxn.genes:
        return default
    return max(gene_expr.get(g.id, default) for g in rxn.genes)


def gimme_style_stub(model, gene_expr, low_quantile, growth_floor=GROWTH_FLOOR):
    '''Illustrative objective-protected pruning. NOT a faithful GIMME; do not report as GIMME output.'''
    cutoff = np.quantile(list(gene_expr.values()), low_quantile)
    low = {g for g, v in gene_expr.items() if v < cutoff}
    ctx = model.copy()
    biomass = str(model.objective.expression).split('*')[1].split()[0]
    ctx.reactions.get_by_id(biomass).lower_bound = growth_floor * model.slim_optimize()
    pruned = 0
    for rxn in ctx.reactions:
        genes = {g.id for g in rxn.genes}
        if genes and genes <= low:
            rxn.bounds = (max(rxn.lower_bound, -1.0), min(rxn.upper_bound, 1.0))
            pruned += 1
    return ctx, pruned


def simulate_expression(model, active=('glycolysis', 'pentose')):
    rng = np.random.default_rng(SEED)
    expr = {}
    for gene in model.genes:
        in_active = any(any(a in r.name.lower() for a in active) for r in gene.reactions)
        expr[gene.id] = rng.uniform(0.7, 1.0) if in_active else rng.uniform(0.0, 0.5)
    return expr


def main():
    model = cobra.io.load_model('textbook')
    expr = simulate_expression(model)

    print('=== GPR-mapped reaction activity (min-AND / max-OR) ===')
    for rid in ['PFK', 'CS', 'PGI']:
        r = model.reactions.get_by_id(rid)
        print(f'  {rid}: activity {reaction_activity(r, expr):.2f}  GPR: {r.gene_reaction_rule or "(none)"}')

    print('\n=== Threshold sweep: the same data, three thresholds, different models ===')
    kept_sets = {}
    for q in LOW_QUANTILES:
        ctx, pruned = gimme_style_stub(model, expr, low_quantile=q)
        kept_sets[q] = {r.id for r in ctx.reactions if ctx.reactions.get_by_id(r.id).bounds != (-1.0, 1.0)}
        growth = ctx.slim_optimize()
        state = f'{growth:.4f} /h' if np.isfinite(growth) else 'INFEASIBLE (floor unmet after pruning)'
        print(f'  quantile {q:.2f}: pruned {pruned} reactions, context growth {state}')

    stable = set.intersection(*kept_sets.values())
    threshold_dependent = set.union(*kept_sets.values()) - stable
    print(f'\nReactions stable across thresholds: {len(stable)}')
    print(f'Threshold-dependent (report as hypotheses): {len(threshold_dependent)}')


if __name__ == '__main__':
    main()
