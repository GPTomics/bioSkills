# Reference: AiZynthFinder 4.4+, RDKit 2024.09+ | Verify API if version differs
# Batch retrosynthesis screening for generative-design feasibility

import pandas as pd
from rdkit import Chem


def setup_finder(config_path):
    '''Initialize AiZynthFinder from a YAML config file.'''
    from aizynthfinder.aizynthfinder import AiZynthFinder
    finder = AiZynthFinder(configfile=config_path)
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
    for r in finder.routes:
        leaves = [n for n in r.leafs()]
        routes.append({
            'depth': r.depth,
            'score': float(r.score) if hasattr(r, 'score') else None,
            'n_leaves': len(leaves),
            'in_stock_leaves': sum(1 for n in leaves if n.in_stock),
            'leaf_smiles': [n.smiles for n in leaves],
        })
    return {'smiles': smi, 'routes': routes}


def classify_feasibility(routes):
    '''Bucket compound as 'easy', 'feasible', or 'unsolved'.'''
    if not routes:
        return 'unsolved'
    best = max(routes, key=lambda r: r['score'] if r['score'] is not None else -1)
    # Easy: short route AND all leaves in stock
    if best['depth'] <= 3 and best['in_stock_leaves'] == best['n_leaves']:
        return 'easy'
    # Feasible: route exists with at least half of leaves in stock
    if best['n_leaves'] > 0 and best['in_stock_leaves'] >= best['n_leaves'] / 2:
        return 'feasible'
    return 'unsolved'


def batch_screen(target_smiles_list, config_path, iteration_limit=100):
    '''Batch retrosynthesis feasibility classification.'''
    finder = setup_finder(config_path)
    rows = []
    for smi in target_smiles_list:
        if Chem.MolFromSmiles(smi) is None:
            rows.append({'smiles': smi, 'feasibility': 'parse_failure', 'best_depth': None, 'n_routes': 0})
            continue
        result = run_retrosynthesis_single(finder, smi, iteration_limit=iteration_limit)
        feasibility = classify_feasibility(result['routes'])
        rows.append({
            'smiles': smi,
            'feasibility': feasibility,
            'best_depth': min((r['depth'] for r in result['routes']), default=None),
            'n_routes': len(result['routes']),
        })
    return pd.DataFrame(rows)


if __name__ == '__main__':
    # Demonstration only; requires config.yaml + template + stock files
    print('Example: configure AiZynthFinder per https://molecularai.github.io/aizynthfinder/howto.html')
