#!/usr/bin/env python3
'''
Metabolic modeling pipeline: genome to flux predictions.
Covers reconstruction, curation, FBA, and gene essentiality.
Requires: cobra, carveme, memote, pandas, numpy, matplotlib
'''
import subprocess
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import cobra
    from cobra.flux_analysis import (flux_variability_analysis,
                                      single_gene_deletion,
                                      gapfill)
except ImportError:
    print('Install cobra: pip install cobra')
    raise

# Configuration
PROTEIN_FASTA = 'genome.faa'
OUTPUT_DIR = Path('metabolic_modeling_results')
GRAM_TYPE = 'neg'  # 'neg' or 'pos'

# Thresholds
MIN_GROWTH = 0.01          # Minimum viable growth rate (h^-1)
MEMOTE_THRESHOLD = 0.50    # Minimum acceptable memote score
ESSENTIALITY_THRESHOLD = 0.10  # Growth ratio below this = essential
FVA_FRACTION = 0.90        # Allow 90% of optimal growth for FVA


def run_carveme(protein_fasta, output_path, gram_type='neg', gapfill_media='M9'):
    '''Run CarveMe for automated model reconstruction.'''
    cmd = ['carve', protein_fasta, '-o', str(output_path)]
    if gram_type:
        cmd.append(f'--gram-{gram_type}')
    if gapfill_media:
        cmd.extend(['--gapfill', gapfill_media])

    print(f'Running CarveMe: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    return output_path


def run_memote(model_path, report_path):
    '''Run memote QC on model.'''
    cmd = ['memote', 'run', '--filename', str(report_path), str(model_path)]
    subprocess.run(cmd, check=True)

    # Also generate snapshot for programmatic access
    snapshot_path = report_path.with_suffix('.json')
    cmd_snapshot = ['memote', 'report', 'snapshot', '--filename', str(snapshot_path), str(model_path)]
    subprocess.run(cmd_snapshot, check=True)

    return snapshot_path


def diagnose_model(model):
    '''Diagnose common model issues.'''
    issues = {'dead_ends': [], 'blocked': [], 'no_gpr': []}

    # Dead-end metabolites
    for met in model.metabolites:
        producing = sum(1 for r in met.reactions if r.get_coefficient(met.id) > 0)
        consuming = sum(1 for r in met.reactions if r.get_coefficient(met.id) < 0)
        if producing > 0 and consuming == 0:
            issues['dead_ends'].append((met.id, 'not consumed'))
        elif producing == 0 and consuming > 0:
            issues['dead_ends'].append((met.id, 'not produced'))

    # Blocked reactions (zero flux at optimum)
    try:
        fva = flux_variability_analysis(model, fraction_of_optimum=0)
        blocked = fva[(fva['minimum'] == 0) & (fva['maximum'] == 0)]
        issues['blocked'] = list(blocked.index)
    except Exception:
        pass

    # Reactions without GPR
    for rxn in model.reactions:
        if not rxn.gene_reaction_rule and not rxn.id.startswith('EX_'):
            issues['no_gpr'].append(rxn.id)

    return issues


def run_fba(model, objective=None):
    '''Run FBA and return solution with key metrics.'''
    if objective:
        model.objective = objective

    solution = model.optimize()

    results = {
        'growth_rate': solution.objective_value,
        'status': solution.status,
        'fluxes': solution.fluxes,
        'active_reactions': sum(abs(solution.fluxes) > 1e-6)
    }

    # Key exchange fluxes
    exchange_fluxes = {r.id: solution.fluxes[r.id] for r in model.exchanges
                       if abs(solution.fluxes[r.id]) > 0.01}
    results['exchanges'] = exchange_fluxes

    return results


def run_fva(model, reaction_list=None, fraction=0.9):
    '''Run FVA to determine flux ranges.'''
    fva = flux_variability_analysis(model, reaction_list=reaction_list,
                                     fraction_of_optimum=fraction)
    fva['range'] = fva['maximum'] - fva['minimum']
    fva['midpoint'] = (fva['maximum'] + fva['minimum']) / 2
    return fva


def predict_essentiality(model, growth_threshold=0.1):
    '''Predict gene essentiality via single knockouts.'''
    wt_growth = model.optimize().objective_value

    ko_results = single_gene_deletion(model)
    ko_results['growth_ratio'] = ko_results['growth'] / wt_growth
    ko_results['essential'] = ko_results['growth_ratio'] < growth_threshold

    essential_genes = [list(g)[0].id for g, row in ko_results.iterrows() if row['essential']]

    return {
        'results': ko_results,
        'essential_genes': essential_genes,
        'n_essential': len(essential_genes),
        'n_total': len(model.genes)
    }


def visualize_results(fba_results, fva_results, essentiality_results, output_dir):
    '''Generate summary visualizations.'''
    output_dir = Path(output_dir)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Exchange fluxes
    ax1 = axes[0, 0]
    exchanges = pd.Series(fba_results['exchanges'])
    if len(exchanges) > 0:
        colors = ['#e74c3c' if v < 0 else '#2ecc71' for v in exchanges]
        exchanges.plot(kind='barh', ax=ax1, color=colors)
        ax1.axvline(0, color='black', linewidth=0.5)
        ax1.set_xlabel('Flux (mmol/gDW/h)')
        ax1.set_title('Exchange Fluxes (red=uptake, green=secretion)')

    # 2. FVA ranges for key reactions
    ax2 = axes[0, 1]
    # Select variable reactions
    variable = fva_results[fva_results['range'] > 0.1].head(20)
    if len(variable) > 0:
        ax2.barh(range(len(variable)), variable['range'],
                left=variable['minimum'], alpha=0.7, color='steelblue')
        ax2.set_yticks(range(len(variable)))
        ax2.set_yticklabels(variable.index, fontsize=8)
        ax2.set_xlabel('Flux range (mmol/gDW/h)')
        ax2.set_title('Flux Variability (Top 20 Variable Reactions)')

    # 3. Essentiality distribution
    ax3 = axes[1, 0]
    ko_data = essentiality_results['results']
    ax3.hist(ko_data['growth_ratio'], bins=50, edgecolor='black', alpha=0.7)
    ax3.axvline(ESSENTIALITY_THRESHOLD, color='red', linestyle='--',
                label=f'Threshold ({ESSENTIALITY_THRESHOLD})')
    ax3.set_xlabel('Growth Ratio (KO / WT)')
    ax3.set_ylabel('Number of Genes')
    ax3.set_title('Gene Essentiality Distribution')
    ax3.legend()

    # 4. Summary stats
    ax4 = axes[1, 1]
    ax4.axis('off')
    summary_text = f'''
    Model Summary
    =============

    Growth rate: {fba_results['growth_rate']:.4f} h^-1
    Active reactions: {fba_results['active_reactions']}

    Essential genes: {essentiality_results['n_essential']} / {essentiality_results['n_total']}
    Essentiality rate: {essentiality_results['n_essential']/essentiality_results['n_total']*100:.1f}%

    Highly variable reactions: {sum(fva_results['range'] > 1)}
    Fixed reactions: {sum(fva_results['range'] < 1e-6)}
    '''
    ax4.text(0.1, 0.5, summary_text, transform=ax4.transAxes, fontsize=12,
             verticalalignment='center', fontfamily='monospace')

    plt.tight_layout()
    plt.savefig(output_dir / 'model_analysis_summary.pdf', dpi=150)
    plt.savefig(output_dir / 'model_analysis_summary.png', dpi=150)
    print(f'Saved visualizations to {output_dir}')


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print('=== Metabolic Modeling Pipeline ===')

    # Step 1: Reconstruction (skip if model exists)
    model_path = OUTPUT_DIR / 'model_draft.xml'
    if not model_path.exists():
        print('\n=== Step 1: Model Reconstruction ===')
        run_carveme(PROTEIN_FASTA, model_path, gram_type=GRAM_TYPE)
    else:
        print(f'\nUsing existing model: {model_path}')

    # Load model
    print('\n=== Loading Model ===')
    model = cobra.io.read_sbml_model(str(model_path))
    print(f'Model: {model.id}')
    print(f'Reactions: {len(model.reactions)}')
    print(f'Metabolites: {len(model.metabolites)}')
    print(f'Genes: {len(model.genes)}')

    # Step 2: Model validation
    print('\n=== Step 2: Model Validation ===')
    report_path = OUTPUT_DIR / 'memote_report.html'

    # Quick diagnostic
    issues = diagnose_model(model)
    print(f"Dead-end metabolites: {len(issues['dead_ends'])}")
    print(f"Blocked reactions: {len(issues['blocked'])}")
    print(f"Reactions without GPR: {len(issues['no_gpr'])}")

    # Quick growth test
    solution = model.optimize()
    print(f'Growth rate: {solution.objective_value:.4f} h^-1')

    if solution.objective_value < MIN_GROWTH:
        print(f'WARNING: Growth below threshold ({MIN_GROWTH}). Model may need gap-filling.')

    # Step 3: FBA Analysis
    print('\n=== Step 3: FBA Analysis ===')
    fba_results = run_fba(model)
    print(f"Growth rate: {fba_results['growth_rate']:.4f} h^-1")
    print(f"Active reactions: {fba_results['active_reactions']}")
    print('\nTop exchanges:')
    for rxn_id, flux in sorted(fba_results['exchanges'].items(), key=lambda x: x[1])[:10]:
        direction = 'uptake' if flux < 0 else 'secretion'
        print(f'  {rxn_id}: {flux:.2f} ({direction})')

    # Save fluxes
    fba_results['fluxes'].to_csv(OUTPUT_DIR / 'fba_fluxes.tsv', sep='\t')

    # Step 4: FVA
    print('\n=== Step 4: Flux Variability Analysis ===')
    fva_results = run_fva(model, fraction=FVA_FRACTION)
    print(f'Rigid reactions (fixed flux): {sum(fva_results["range"] < 1e-6)}')
    print(f'Flexible reactions (range > 1): {sum(fva_results["range"] > 1)}')
    fva_results.to_csv(OUTPUT_DIR / 'fva_results.tsv', sep='\t')

    # Step 5: Gene Essentiality
    print('\n=== Step 5: Gene Essentiality Prediction ===')
    essentiality_results = predict_essentiality(model, growth_threshold=ESSENTIALITY_THRESHOLD)
    print(f"Essential genes: {essentiality_results['n_essential']} / {essentiality_results['n_total']}")
    print(f"Essentiality rate: {essentiality_results['n_essential']/essentiality_results['n_total']*100:.1f}%")

    # Save essential genes
    with open(OUTPUT_DIR / 'essential_genes.txt', 'w') as f:
        f.write('\n'.join(essentiality_results['essential_genes']))

    essentiality_results['results'].to_csv(OUTPUT_DIR / 'gene_essentiality.tsv', sep='\t')

    # Visualization
    print('\n=== Generating Visualizations ===')
    visualize_results(fba_results, fva_results, essentiality_results, OUTPUT_DIR)

    # Save curated model
    cobra.io.write_sbml_model(model, str(OUTPUT_DIR / 'model_analyzed.xml'))

    print('\n=== Pipeline Complete ===')
    print(f'Results saved to: {OUTPUT_DIR}/')
    print('\nKey outputs:')
    print(f'  Model: {OUTPUT_DIR}/model_analyzed.xml')
    print(f'  FBA fluxes: {OUTPUT_DIR}/fba_fluxes.tsv')
    print(f'  FVA results: {OUTPUT_DIR}/fva_results.tsv')
    print(f'  Essential genes: {OUTPUT_DIR}/essential_genes.txt')
    print(f'  Summary figures: {OUTPUT_DIR}/model_analysis_summary.pdf')


if __name__ == '__main__':
    main()
