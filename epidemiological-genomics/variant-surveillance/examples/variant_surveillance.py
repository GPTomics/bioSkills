'''Viral variant surveillance with Nextclade'''

import subprocess
import pandas as pd
from collections import Counter


def run_nextclade(sequences_fasta, dataset_dir, output_tsv):
    '''Run Nextclade for lineage assignment

    Nextclade provides:
    - Lineage/clade assignment
    - Mutation calling
    - Quality control metrics
    - Phylogenetic placement
    '''
    cmd = [
        'nextclade', 'run',
        '--input-dataset', dataset_dir,
        '--output-tsv', output_tsv,
        sequences_fasta
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'Nextclade error: {result.stderr}')
        return None

    return output_tsv


def parse_nextclade_results(tsv_file):
    '''Parse Nextclade output and filter by QC'''
    df = pd.read_csv(tsv_file, sep='\t')

    # QC filter
    # good: High quality, all metrics pass
    # mediocre: Minor issues, usable
    # bad: Major issues, exclude from analysis
    df['pass_qc'] = df['qc.overallStatus'].isin(['good', 'mediocre'])

    return df


# VOC definitions (update as needed)
VOC_LINEAGES = {
    'Omicron': ['BA.', 'XBB.', 'JN.', 'KP.'],
    'Delta': ['AY.', 'B.1.617.2'],
    'Alpha': ['B.1.1.7', 'Q.'],
}


def classify_voc(lineage):
    '''Classify lineage as VOC/VOI'''
    if pd.isna(lineage):
        return 'Unknown'

    for voc, prefixes in VOC_LINEAGES.items():
        for prefix in prefixes:
            if lineage.startswith(prefix) or lineage == prefix.rstrip('.'):
                return voc
    return 'Other'


def summarize_lineages(results_df, lineage_col='Nextclade_pango'):
    '''Summarize lineage distribution'''
    passed = results_df[results_df['pass_qc']]

    # Lineage counts
    lineage_counts = passed[lineage_col].value_counts()

    # VOC classification
    passed = passed.copy()
    passed['VOC'] = passed[lineage_col].apply(classify_voc)
    voc_counts = passed['VOC'].value_counts()

    return {
        'total': len(results_df),
        'passed_qc': len(passed),
        'unique_lineages': passed[lineage_col].nunique(),
        'top_lineages': lineage_counts.head(10).to_dict(),
        'voc_distribution': voc_counts.to_dict()
    }


def track_temporal_trends(results_df, date_col='date', lineage_col='Nextclade_pango'):
    '''Track lineage prevalence over time'''
    df = results_df[results_df['pass_qc']].copy()
    df['date'] = pd.to_datetime(df[date_col])
    df['week'] = df['date'].dt.to_period('W')

    # Prevalence by week
    weekly = df.groupby(['week', lineage_col]).size().unstack(fill_value=0)

    # Convert to percentages
    weekly_pct = weekly.div(weekly.sum(axis=1), axis=0) * 100

    return weekly_pct


def detect_emerging_variants(results_df, date_col='date', lineage_col='Nextclade_pango', threshold=5):
    '''Detect variants increasing in prevalence

    Definition of "emerging":
    - Prevalence increased >2-fold between periods
    - Recent prevalence >5% (or specified threshold)
    '''
    df = results_df[results_df['pass_qc']].copy()
    df['date'] = pd.to_datetime(df[date_col])

    midpoint = df['date'].median()
    early = df[df['date'] < midpoint]
    recent = df[df['date'] >= midpoint]

    early_counts = early[lineage_col].value_counts(normalize=True) * 100
    recent_counts = recent[lineage_col].value_counts(normalize=True) * 100

    emerging = []
    for lineage in recent_counts.index:
        recent_prev = recent_counts[lineage]
        early_prev = early_counts.get(lineage, 0)

        if recent_prev >= threshold and (early_prev == 0 or recent_prev / early_prev >= 2):
            emerging.append({
                'lineage': lineage,
                'early_prevalence': f'{early_prev:.1f}%',
                'recent_prevalence': f'{recent_prev:.1f}%',
                'trend': 'NEW' if early_prev == 0 else f'{recent_prev/early_prev:.1f}x increase'
            })

    return sorted(emerging, key=lambda x: -float(x['recent_prevalence'].rstrip('%')))


if __name__ == '__main__':
    print('Variant Surveillance Example')
    print('=' * 50)

    # Simulated Nextclade results
    demo_data = pd.DataFrame({
        'seqName': [f'seq_{i}' for i in range(100)],
        'Nextclade_pango': ['JN.1'] * 40 + ['XBB.1.5'] * 25 + ['BA.2.86'] * 20 +
                          ['EG.5'] * 10 + ['Other'] * 5,
        'qc.overallStatus': ['good'] * 85 + ['mediocre'] * 10 + ['bad'] * 5,
        'date': pd.date_range('2024-01-01', periods=100, freq='D')
    })

    print('\nLineage Summary:')
    summary = summarize_lineages(demo_data)
    print(f"  Total sequences: {summary['total']}")
    print(f"  Passed QC: {summary['passed_qc']}")
    print(f"  Unique lineages: {summary['unique_lineages']}")

    print('\nTop Lineages:')
    for lineage, count in list(summary['top_lineages'].items())[:5]:
        pct = count / summary['passed_qc'] * 100
        print(f"  {lineage}: {count} ({pct:.1f}%)")

    print('\nVOC Distribution:')
    for voc, count in summary['voc_distribution'].items():
        pct = count / summary['passed_qc'] * 100
        print(f"  {voc}: {count} ({pct:.1f}%)")

    print('\nEmerging Variants:')
    emerging = detect_emerging_variants(demo_data, threshold=10)
    if emerging:
        for e in emerging[:3]:
            print(f"  {e['lineage']}: {e['early_prevalence']} -> {e['recent_prevalence']} ({e['trend']})")
    else:
        print('  No significantly emerging variants detected')
