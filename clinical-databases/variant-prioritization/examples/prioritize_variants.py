'''Prioritize variants for rare disease analysis'''

import myvariant
import pandas as pd

mv = myvariant.MyVariantInfo()

def annotate_variants(variants):
    '''Fetch annotations for variant list from myvariant.info'''
    results = mv.getvariants(
        variants,
        fields=[
            'clinvar.clinical_significance',
            'gnomad_exome.af.af',
            'cadd.phred',
            'dbnsfp.revel.score'
        ]
    )

    records = []
    for r in results:
        clinvar = r.get('clinvar', {})
        gnomad = r.get('gnomad_exome', {})
        cadd = r.get('cadd', {})
        revel = r.get('dbnsfp', {}).get('revel', {})

        gnomad_af = gnomad.get('af', {}).get('af') if isinstance(gnomad, dict) else None
        revel_score = revel.get('score') if isinstance(revel, dict) else None

        records.append({
            'variant': r.get('query'),
            'clinvar_sig': clinvar.get('clinical_significance') if isinstance(clinvar, dict) else None,
            'gnomad_af': gnomad_af,
            'cadd_phred': cadd.get('phred') if isinstance(cadd, dict) else None,
            'revel_score': revel_score[0] if isinstance(revel_score, list) else revel_score
        })

    return pd.DataFrame(records)

def prioritize(df):
    '''Apply multi-evidence prioritization scoring

    Scoring:
    - ClinVar pathogenic: +10 (strong evidence from clinical labs)
    - Rare (AF < 0.01): +3 (ACMG PM2 supporting criterion)
    - CADD > 20: +2 (top 1% deleteriousness, computational)
    - REVEL > 0.5: +2 (missense pathogenicity predictor)
    '''
    # Rare in population (ACMG PM2)
    # Threshold: 0.01 (1%) is standard for rare Mendelian disease
    # Use 0.001 for more stringent filtering
    GNOMAD_THRESHOLD = 0.01
    df['is_rare'] = df['gnomad_af'].isna() | (df['gnomad_af'] < GNOMAD_THRESHOLD)

    # ClinVar pathogenic
    pathogenic_terms = ['Pathogenic', 'Likely_pathogenic']
    df['clinvar_path'] = df['clinvar_sig'].apply(
        lambda x: any(p in str(x) for p in pathogenic_terms) if pd.notna(x) else False
    )

    # CADD phred > 20 = top 1% most deleterious
    # Higher threshold (e.g., 25-30) for more stringent filtering
    CADD_THRESHOLD = 20
    df['cadd_high'] = df['cadd_phred'].fillna(0) > CADD_THRESHOLD

    # REVEL > 0.5 suggests pathogenic (for missense only)
    # Original paper: 0.5 threshold gives ~75% sensitivity, 90% specificity
    REVEL_THRESHOLD = 0.5
    df['revel_high'] = df['revel_score'].fillna(0) > REVEL_THRESHOLD

    # Composite priority score
    df['priority'] = (
        df['clinvar_path'].astype(int) * 10 +
        df['is_rare'].astype(int) * 3 +
        df['cadd_high'].astype(int) * 2 +
        df['revel_high'].astype(int) * 2
    )

    return df.sort_values('priority', ascending=False)

def assign_tiers(df):
    '''Assign clinical interpretation tiers

    Tier 1: Strong pathogenic - ClinVar P/LP + rare
    Tier 2: Likely pathogenic - rare + computational support
    Tier 3: VUS - rare but limited evidence
    Tier 4: Likely benign - common or benign evidence
    '''
    def get_tier(row):
        if row['clinvar_path'] and row['is_rare']:
            return 'Tier 1 - Strong'
        elif row['is_rare'] and (row['cadd_high'] or row['revel_high']):
            return 'Tier 2 - Likely'
        elif row['is_rare']:
            return 'Tier 3 - VUS'
        else:
            return 'Tier 4 - Benign'

    df['tier'] = df.apply(get_tier, axis=1)
    return df

# Example variants for demonstration
variants = [
    'rs121913527',      # BRAF V600E - pathogenic, somatic
    'rs1800566',        # NQO1 P187S - common, benign
    'rs104894155',      # CFTR F508del - pathogenic, rare
    'rs28897696',       # BRCA1 variant
]

print('Variant prioritization pipeline\n')

# Annotate
df = annotate_variants(variants)
print('Raw annotations:')
print(df[['variant', 'clinvar_sig', 'gnomad_af', 'cadd_phred']].to_string(index=False))

# Prioritize
df = prioritize(df)
print('\nPrioritized (by score):')
print(df[['variant', 'priority', 'is_rare', 'clinvar_path', 'cadd_high']].to_string(index=False))

# Assign tiers
df = assign_tiers(df)
print('\nClinical tiers:')
print(df[['variant', 'tier', 'priority']].to_string(index=False))
