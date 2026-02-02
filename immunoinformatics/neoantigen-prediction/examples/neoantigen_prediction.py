'''Identify tumor neoantigens from somatic mutations'''

import pandas as pd


def filter_neoantigens(results_df, binding_threshold=500, vaf_threshold=0.1):
    '''Filter neoantigen candidates

    Standard filtering criteria:
    - Binding affinity < 500 nM (strong/moderate binders)
    - VAF > 10% (clonal mutations)
    - Agretopicity > 1 (MT binds better than WT)

    Args:
        binding_threshold: Max IC50 in nM (500 standard, 50 stringent)
        vaf_threshold: Minimum variant allele frequency
    '''
    filtered = results_df.copy()

    # Binding filter
    if 'Median MT Score' in filtered.columns:
        filtered = filtered[filtered['Median MT Score'] < binding_threshold]

    # VAF filter (clonality)
    if 'Tumor DNA VAF' in filtered.columns:
        filtered = filtered[filtered['Tumor DNA VAF'] >= vaf_threshold]

    return filtered


def calculate_agretopicity(results_df):
    '''Calculate agretopicity (differential binding)

    Agretopicity = WT_binding / MT_binding

    Interpretation:
    - >1: Mutation improves MHC binding (favorable for immunogenicity)
    - ~1: Similar binding (mutation may not create novel epitope)
    - <1: WT binds better (unfavorable)

    High agretopicity suggests mutation creates truly novel epitope
    '''
    df = results_df.copy()

    if 'Median WT Score' in df.columns and 'Median MT Score' in df.columns:
        # Avoid division by zero
        df['agretopicity'] = df['Median WT Score'] / df['Median MT Score'].clip(lower=0.1)
    else:
        df['agretopicity'] = 1.0

    return df


def prioritize_neoantigens(df):
    '''Prioritize neoantigens for vaccine design

    Priority scoring considers:
    1. Binding affinity (lower = better)
    2. Agretopicity (higher = better)
    3. Clonality/VAF (higher = better)
    4. Gene expression (higher = better)

    Returns sorted DataFrame with priority scores
    '''
    df = df.copy()

    # Normalize metrics to 0-1 scale
    # Binding: inverse relationship (lower IC50 = higher score)
    df['binding_score'] = 1 - (df['Median MT Score'].clip(upper=5000) / 5000)

    # Agretopicity: cap at 10 for scoring
    df['agretopicity_score'] = df['agretopicity'].clip(upper=10) / 10

    # VAF
    if 'Tumor DNA VAF' in df.columns:
        df['vaf_score'] = df['Tumor DNA VAF']
    else:
        df['vaf_score'] = 0.5  # Default if not available

    # Expression (if available)
    if 'Gene Expression' in df.columns:
        max_expr = df['Gene Expression'].max()
        df['expression_score'] = df['Gene Expression'] / max_expr if max_expr > 0 else 0.5
    else:
        df['expression_score'] = 0.5

    # Composite priority score
    df['priority_score'] = (
        0.35 * df['binding_score'] +
        0.25 * df['agretopicity_score'] +
        0.25 * df['vaf_score'] +
        0.15 * df['expression_score']
    )

    return df.sort_values('priority_score', ascending=False)


def summarize_neoantigens(df):
    '''Summarize neoantigen prediction results'''
    print('Neoantigen Prediction Summary')
    print('=' * 50)
    print(f"Total candidates: {len(df)}")

    if 'Median MT Score' in df.columns:
        strong = (df['Median MT Score'] < 50).sum()
        moderate = ((df['Median MT Score'] >= 50) & (df['Median MT Score'] < 500)).sum()
        print(f"Strong binders (<50nM): {strong}")
        print(f"Moderate binders (50-500nM): {moderate}")

    if 'agretopicity' in df.columns:
        high_agretopicity = (df['agretopicity'] > 1).sum()
        print(f"High agretopicity (DAI>1): {high_agretopicity}")

    if 'Tumor DNA VAF' in df.columns:
        clonal = (df['Tumor DNA VAF'] > 0.3).sum()
        print(f"Clonal (VAF>30%): {clonal}")

    unique_genes = df['Gene Name'].nunique() if 'Gene Name' in df.columns else 'N/A'
    print(f"Unique genes: {unique_genes}")


if __name__ == '__main__':
    # Simulated pVACseq output
    demo_data = pd.DataFrame({
        'Gene Name': ['TP53', 'KRAS', 'BRAF', 'PIK3CA', 'EGFR', 'NRAS'],
        'Mutation': ['R175H', 'G12D', 'V600E', 'E545K', 'L858R', 'Q61K'],
        'MT Epitope Seq': ['HMTEVVRHC', 'VVVGADGVGK', 'LATEKSRWSG', 'STRDPLSEIT', 'KITDFGLAKL', 'ILDTAGKEEY'],
        'WT Epitope Seq': ['RMTEVVRHC', 'VVVGAGGVGK', 'LATEKSRWSG', 'STRDPLSKIT', 'KITDFGLAKL', 'ILDTAGQEEY'],
        'HLA Allele': ['HLA-A*02:01'] * 6,
        'Median MT Score': [45, 120, 350, 85, 420, 180],  # IC50 nM
        'Median WT Score': [1500, 800, 350, 450, 500, 600],
        'Tumor DNA VAF': [0.45, 0.35, 0.25, 0.15, 0.40, 0.20],
        'Gene Expression': [25, 150, 80, 45, 200, 30]
    })

    print('Input data:')
    print(demo_data[['Gene Name', 'Mutation', 'Median MT Score', 'Tumor DNA VAF']].to_string(index=False))

    # Calculate agretopicity
    print('\n')
    demo_data = calculate_agretopicity(demo_data)

    # Filter
    filtered = filter_neoantigens(demo_data, binding_threshold=500, vaf_threshold=0.1)
    print(f'After filtering: {len(filtered)} candidates')

    # Prioritize
    print('\n')
    prioritized = prioritize_neoantigens(filtered)
    summarize_neoantigens(prioritized)

    # Top candidates
    print('\nTop neoantigen candidates:')
    cols = ['Gene Name', 'Mutation', 'Median MT Score', 'agretopicity', 'priority_score']
    print(prioritized[cols].head(5).to_string(index=False))
