'''Score and rank epitopes for immunogenicity'''

import pandas as pd
import numpy as np


def calculate_immunogenicity_score(peptide_data):
    '''Calculate composite immunogenicity score

    Factors and typical weights:
    - Binding affinity (25%): Lower IC50 = better
    - Agretopicity (20%): MT/WT ratio for neoantigens
    - Processing (10%): Proteasomal cleavage likelihood
    - Expression (15%): Gene expression level
    - Clonality (15%): VAF for tumor neoantigens
    - Foreignness (15%): Dissimilarity to self

    Score range: 0-1 (higher = more immunogenic)
    '''
    scores = {}

    # Binding: transform IC50 to 0-1 (lower = better)
    ic50 = peptide_data.get('ic50_nM', 500)
    scores['binding'] = max(0, 1 - ic50 / 5000)

    # Agretopicity: ratio capped at 10
    agretopicity = peptide_data.get('agretopicity', 1.0)
    scores['agretopicity'] = min(agretopicity / 10, 1)

    # Processing: direct score 0-1
    scores['processing'] = peptide_data.get('processing_score', 0.5)

    # Expression: log-transform, cap at reasonable max
    expression = peptide_data.get('expression_tpm', 10)
    scores['expression'] = min(np.log10(expression + 1) / 3, 1)

    # Clonality: direct VAF
    scores['clonality'] = peptide_data.get('vaf', 0.5)

    # Foreignness: 1 - self_similarity
    self_sim = peptide_data.get('self_similarity', 0.5)
    scores['foreignness'] = 1 - self_sim

    # Weighted sum
    weights = {
        'binding': 0.25,
        'agretopicity': 0.20,
        'processing': 0.10,
        'expression': 0.15,
        'clonality': 0.15,
        'foreignness': 0.15
    }

    total = sum(scores[k] * weights.get(k, 0) for k in scores)

    return total, scores


def rank_epitopes(epitope_df, top_n=20):
    '''Rank epitopes by immunogenicity score'''
    results = []

    for _, row in epitope_df.iterrows():
        total, factors = calculate_immunogenicity_score(row.to_dict())
        result = row.to_dict()
        result['immunogenicity_score'] = total
        result.update({f'{k}_score': v for k, v in factors.items()})
        results.append(result)

    ranked = pd.DataFrame(results)
    ranked = ranked.sort_values('immunogenicity_score', ascending=False)

    # Assign confidence tiers
    n = len(ranked)
    ranked['tier'] = 'low'
    ranked.loc[ranked.index[:max(1, int(n * 0.20))], 'tier'] = 'medium'
    ranked.loc[ranked.index[:max(1, int(n * 0.05))], 'tier'] = 'high'

    return ranked.head(top_n)


def check_anchor_residues(peptide, allele='HLA-A*02:01'):
    '''Check MHC anchor residue preferences

    HLA-A*02:01 preferences:
    - Position 2: L, M, I, V, A (hydrophobic)
    - Position C-terminus: V, L, I (hydrophobic)

    Good anchors improve binding stability
    '''
    if len(peptide) < 2:
        return {'anchor_score': 0}

    pos2_preferred = set('LMIVA')
    c_term_preferred = set('VLI')

    pos2_good = peptide[1] in pos2_preferred
    c_term_good = peptide[-1] in c_term_preferred

    return {
        'pos2_residue': peptide[1],
        'pos2_preferred': pos2_good,
        'c_term_residue': peptide[-1],
        'c_term_preferred': c_term_good,
        'anchor_score': int(pos2_good) + int(c_term_good)
    }


def summarize_ranking(ranked_df):
    '''Summarize epitope ranking results'''
    print('Immunogenicity Ranking Summary')
    print('=' * 50)

    print(f"\nTotal candidates: {len(ranked_df)}")
    tier_counts = ranked_df['tier'].value_counts()
    for tier in ['high', 'medium', 'low']:
        if tier in tier_counts:
            print(f"  {tier.capitalize()} confidence: {tier_counts[tier]}")

    print(f"\nScore range: {ranked_df['immunogenicity_score'].min():.3f} - "
          f"{ranked_df['immunogenicity_score'].max():.3f}")

    print('\nTop 5 candidates:')
    top5 = ranked_df.head(5)
    for i, (_, row) in enumerate(top5.iterrows(), 1):
        pep = row.get('peptide', 'N/A')
        score = row['immunogenicity_score']
        tier = row['tier']
        print(f"  {i}. {pep}: {score:.3f} ({tier})")


if __name__ == '__main__':
    # Example epitope candidates
    candidates = pd.DataFrame({
        'peptide': ['SIINFEKL', 'GILGFVFTL', 'NLVPMVATV', 'YMLDLQPET', 'FLPSDFFPSV'],
        'ic50_nM': [25, 85, 120, 450, 35],
        'agretopicity': [5.2, 2.1, 1.5, 0.8, 3.8],
        'processing_score': [0.85, 0.72, 0.65, 0.45, 0.78],
        'expression_tpm': [150, 80, 200, 25, 120],
        'vaf': [0.45, 0.35, 0.50, 0.15, 0.40],
        'self_similarity': [0.3, 0.4, 0.2, 0.7, 0.35]
    })

    print('Epitope Immunogenicity Scoring')
    print('=' * 50)

    print('\nInput candidates:')
    print(candidates[['peptide', 'ic50_nM', 'agretopicity', 'expression_tpm']].to_string(index=False))

    # Rank
    ranked = rank_epitopes(candidates, top_n=10)

    print('\n')
    summarize_ranking(ranked)

    # Show score breakdown for top candidate
    print('\nScore breakdown for top candidate:')
    top = ranked.iloc[0]
    for col in ['binding_score', 'agretopicity_score', 'processing_score',
                'expression_score', 'clonality_score', 'foreignness_score']:
        if col in top:
            print(f"  {col.replace('_score', '').capitalize()}: {top[col]:.3f}")

    # Check anchors
    print('\nAnchor residue analysis:')
    for pep in candidates['peptide']:
        anchors = check_anchor_residues(pep)
        status = 'Good' if anchors['anchor_score'] == 2 else 'Partial' if anchors['anchor_score'] == 1 else 'Poor'
        print(f"  {pep}: P2={anchors['pos2_residue']}, C-term={anchors['c_term_residue']} ({status})")
