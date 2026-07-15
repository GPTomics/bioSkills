# Reference: AiZynthFinder 4.4+, RDKit 2024.09+ | Verify API if version differs
# Batch retrosynthesis screening for generative-design feasibility

import pandas as pd
from rdkit import Chem


def setup_finder(config_path, expansion_policy, stock, filter_policy=None):
    '''Initialize AiZynthFinder and explicitly select configured resources.'''
    from aizynthfinder.aizynthfinder import AiZynthFinder
    finder = AiZynthFinder(configfile=config_path)
    finder.expansion_policy.select(expansion_policy)
    finder.stock.select(stock)
    if filter_policy is not None:
        finder.filter_policy.select(filter_policy)
    return finder


def run_retrosynthesis_single(finder, smi, iteration_limit=100, time_limit=120):
    '''Plan retrosynthesis for a single target; return route summary list.'''
    finder.target_smiles = smi
    finder.config.search.iteration_limit = iteration_limit
    finder.config.search.time_limit = time_limit
    try:
        finder.tree_search()
        finder.build_routes()
    except Exception as e:
        return {'error': str(e), 'smiles': smi, 'routes': []}

    routes = []
    for tree, score in zip(finder.routes.reaction_trees, finder.routes.scores):
        leaves = list(tree.leafs())
        routes.append({
            'n_steps': len(list(tree.reactions())),
            'score': score,
            'is_solved': tree.is_solved,
            'n_leaves': len(leaves),
            'in_stock_leaves': sum(1 for n in leaves if tree.in_stock(n)),
            'leaf_smiles': [n.smiles for n in leaves],
        })
    return {'smiles': smi, 'routes': routes}


def classify_feasibility(routes, short_route_max_steps=None):
    '''Classify solved state, optionally labeling a user-defined short route.'''
    if not routes:
        return 'unsolved'
    solved = [route for route in routes if route['is_solved']]
    if not solved:
        return 'unsolved'
    if (
        short_route_max_steps is not None
        and min(route['n_steps'] for route in solved) <= short_route_max_steps
    ):
        return 'short_solved'
    return 'solved'


def batch_screen(
    target_smiles_list,
    config_path,
    expansion_policy,
    stock,
    filter_policy=None,
    iteration_limit=100,
    short_route_max_steps=None,
):
    '''Batch retrosynthesis feasibility classification.'''
    finder = setup_finder(
        config_path,
        expansion_policy=expansion_policy,
        stock=stock,
        filter_policy=filter_policy,
    )
    rows = []
    for smi in target_smiles_list:
        if Chem.MolFromSmiles(smi) is None:
            rows.append({
                'smiles': smi,
                'feasibility': 'parse_failure',
                'best_steps': None,
                'n_routes': 0,
                'error': 'invalid SMILES',
            })
            continue
        result = run_retrosynthesis_single(finder, smi, iteration_limit=iteration_limit)
        if result.get('error'):
            rows.append({
                'smiles': smi,
                'feasibility': 'execution_error',
                'best_steps': None,
                'n_routes': 0,
                'error': result['error'],
            })
            continue
        feasibility = classify_feasibility(
            result['routes'],
            short_route_max_steps=short_route_max_steps,
        )
        rows.append({
            'smiles': smi,
            'feasibility': feasibility,
            'best_steps': min(
                (r['n_steps'] for r in result['routes'] if r['is_solved']),
                default=None,
            ),
            'n_routes': len(result['routes']),
            'error': None,
        })
    return pd.DataFrame(rows)


if __name__ == '__main__':
    # Demonstration only; requires config.yaml + template + stock files
    print('Example: configure and select resources per https://molecularai.github.io/aizynthfinder/python_interface.html')
