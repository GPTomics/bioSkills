'''In silico gene knockout analysis'''

import cobra
from cobra.flux_analysis import single_gene_deletion
import pandas as pd


def analyze_single_knockouts(model, growth_threshold=0.01):
    '''Perform single gene deletion analysis

    Args:
        model: COBRApy model
        growth_threshold: Below this = essential (default 0.01)
                         This allows for numerical solver tolerance

    Returns:
        DataFrame with knockout results and classifications
    '''
    wt_growth = model.optimize().objective_value

    results = single_gene_deletion(model)
    results['relative_growth'] = results['growth'] / wt_growth
    results['gene_id'] = results['ids'].apply(lambda x: list(x)[0] if x else None)

    # Classification thresholds:
    # Essential: <1% of WT (effectively lethal)
    # Reduced: 1-50% of WT (significant growth defect)
    # Non-essential: >50% of WT
    results['classification'] = 'non-essential'
    results.loc[results['relative_growth'] < 0.5, 'classification'] = 'reduced'
    results.loc[results['relative_growth'] < growth_threshold, 'classification'] = 'essential'

    return results


def find_synthetic_lethal(model, gene_subset=None, growth_threshold=0.01):
    '''Find synthetic lethal gene pairs

    Synthetic lethality: Two genes where:
    - Single KO of either gene is viable
    - Double KO is lethal

    Warning: O(n^2) complexity. For 100 genes = 5000 pairs.
    Full E. coli model (~1500 genes) is impractical.
    '''
    from cobra.flux_analysis import double_gene_deletion

    if gene_subset is None:
        gene_subset = [g.id for g in model.genes[:30]]  # Limit for speed

    # Get viable single knockouts
    single = single_gene_deletion(model, gene_list=gene_subset)
    viable = single[single['growth'] > growth_threshold]
    viable_genes = [list(ids)[0] for ids in viable['ids']]

    print(f'Testing {len(viable_genes)} viable genes ({len(viable_genes)**2//2} pairs)')

    # Double deletions
    double = double_gene_deletion(model, gene_list1=viable_genes, gene_list2=viable_genes)

    # Find synthetic lethals
    sl_pairs = []
    for _, row in double.iterrows():
        genes = list(row['ids'])
        if len(genes) == 2 and row['growth'] < growth_threshold:
            sl_pairs.append({
                'gene1': genes[0],
                'gene2': genes[1],
                'growth': row['growth']
            })

    return sl_pairs


def compare_essentiality_conditions(model):
    '''Compare essential genes under different conditions

    Example: Aerobic vs anaerobic growth
    '''
    results = {}

    # Aerobic
    with model:
        model.reactions.EX_o2_e.lower_bound = -20  # Allow oxygen
        single = single_gene_deletion(model)
        results['aerobic'] = set(list(ids)[0] for ids in single[single['growth'] < 0.01]['ids'])

    # Anaerobic
    with model:
        model.reactions.EX_o2_e.lower_bound = 0  # No oxygen
        single = single_gene_deletion(model)
        results['anaerobic'] = set(list(ids)[0] for ids in single[single['growth'] < 0.01]['ids'])

    # Analysis
    core = results['aerobic'] & results['anaerobic']
    aerobic_only = results['aerobic'] - results['anaerobic']
    anaerobic_only = results['anaerobic'] - results['aerobic']

    return {
        'core_essential': core,
        'aerobic_specific': aerobic_only,
        'anaerobic_specific': anaerobic_only
    }


if __name__ == '__main__':
    model = cobra.io.load_model('textbook')

    print('Gene Essentiality Analysis')
    print('=' * 50)

    # Wild-type growth
    wt = model.optimize()
    print(f'Wild-type growth rate: {wt.objective_value:.4f} h^-1')

    # Single knockouts
    print('\nRunning single gene deletions...')
    results = analyze_single_knockouts(model)

    # Summary
    counts = results['classification'].value_counts()
    print(f'\nGene Classification:')
    for cls, count in counts.items():
        print(f'  {cls}: {count} genes')

    # List essential genes
    essential = results[results['classification'] == 'essential']
    print(f'\nEssential genes ({len(essential)}):')
    for _, row in essential.head(10).iterrows():
        print(f"  {row['gene_id']}: growth = {row['growth']:.4f}")

    # Synthetic lethality (small subset for speed)
    print('\nRunning synthetic lethality analysis...')
    gene_subset = [g.id for g in model.genes[:20]]
    sl_pairs = find_synthetic_lethal(model, gene_subset)
    print(f'Found {len(sl_pairs)} synthetic lethal pairs')

    if sl_pairs:
        print('\nTop synthetic lethal pairs:')
        for pair in sl_pairs[:5]:
            print(f"  {pair['gene1']} + {pair['gene2']}")
