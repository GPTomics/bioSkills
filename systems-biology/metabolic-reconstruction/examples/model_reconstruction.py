'''Build metabolic models from genome sequences'''

import subprocess
import cobra
from pathlib import Path


def reconstruct_with_carveme(fasta_path, output_path, media='M9', grampos=False, verbose=False):
    '''Run CarveMe reconstruction

    Args:
        fasta_path: Path to protein FASTA file
        output_path: Output model file (.xml or .json)
        media: Gap-filling media (M9, LB, M9[glc], etc.)
        grampos: True for Gram-positive bacteria

    Model size expectations by organism type:
    - Bacteria: 1000-2500 reactions
    - Archaea: 800-1500 reactions
    - Fungi: 1500-3000 reactions
    '''
    cmd = ['carve', str(fasta_path), '-o', str(output_path)]

    if media:
        cmd.extend(['--gapfill', media])
    if grampos:
        cmd.append('--grampos')
    if verbose:
        cmd.append('-v')

    subprocess.run(cmd, check=True)
    print(f'Model saved to: {output_path}')
    return output_path


def assess_model_quality(model):
    '''Assess quality of reconstructed model

    Key metrics:
    - Gene-reaction ratio: 1.5-2.5 typical for bacteria
    - Orphan reactions: Should be <10% of total
    - Growth: Must produce biomass on target media
    '''
    metrics = {
        'reactions': len(model.reactions),
        'metabolites': len(model.metabolites),
        'genes': len(model.genes),
        'exchanges': len(model.exchanges),
    }

    # Gene-reaction ratio (typical: 1.5-2.5 for bacteria)
    # Too high may indicate many generic reactions
    # Too low may indicate missing annotations
    if model.genes:
        metrics['gene_reaction_ratio'] = len(model.reactions) / len(model.genes)
    else:
        metrics['gene_reaction_ratio'] = None

    # Orphan reactions (no gene association)
    # Should be <10% of total reactions
    orphans = [r for r in model.reactions if not r.genes]
    metrics['orphan_reactions'] = len(orphans)
    metrics['orphan_fraction'] = len(orphans) / len(model.reactions)

    # Test growth capability
    try:
        sol = model.optimize()
        metrics['can_grow'] = sol.status == 'optimal' and sol.objective_value > 0.001
        metrics['growth_rate'] = sol.objective_value
    except:
        metrics['can_grow'] = False
        metrics['growth_rate'] = 0

    return metrics


def compare_models(model1, model2):
    '''Compare two reconstructed models

    Useful for comparing strains or reconstruction methods.
    '''
    rxns1 = {r.id for r in model1.reactions}
    rxns2 = {r.id for r in model2.reactions}

    return {
        'model1_only': len(rxns1 - rxns2),
        'model2_only': len(rxns2 - rxns1),
        'shared': len(rxns1 & rxns2),
        'jaccard_similarity': len(rxns1 & rxns2) / len(rxns1 | rxns2)
    }


if __name__ == '__main__':
    # Example: Load and assess a model
    # This uses the textbook model as an example
    model = cobra.io.load_model('textbook')

    print('Model Quality Assessment')
    print('=' * 40)
    metrics = assess_model_quality(model)
    for key, value in metrics.items():
        if isinstance(value, float):
            print(f'{key}: {value:.3f}')
        else:
            print(f'{key}: {value}')

    # Quality interpretation
    print('\nQuality Interpretation:')
    if metrics['can_grow']:
        print('  [PASS] Model can produce biomass')
    else:
        print('  [FAIL] Model cannot grow - needs gap-filling')

    if metrics['orphan_fraction'] < 0.1:
        print(f"  [PASS] Orphan reactions: {metrics['orphan_fraction']:.1%} (<10%)")
    else:
        print(f"  [WARN] High orphan fraction: {metrics['orphan_fraction']:.1%}")

    if metrics['gene_reaction_ratio']:
        if 1.5 <= metrics['gene_reaction_ratio'] <= 2.5:
            print(f"  [PASS] Gene-reaction ratio: {metrics['gene_reaction_ratio']:.2f}")
        else:
            print(f"  [WARN] Unusual gene-reaction ratio: {metrics['gene_reaction_ratio']:.2f}")
