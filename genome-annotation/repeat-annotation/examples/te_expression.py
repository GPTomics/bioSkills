#!/usr/bin/env python3
'''
TE expression analysis with TEtranscripts output.
'''
# Reference: matplotlib 3.8+, pandas 2.2+ | Verify API if version differs

import subprocess
import pandas as pd
import matplotlib.pyplot as plt


def run_tetranscripts(treatment_bams, control_bams, gene_gtf, te_gtf, output_prefix='te_analysis', stranded='no'):
    '''
    Run TEtranscripts for differential TE expression.

    STAR alignment must retain multi-mappers:
    --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100
    '''
    cmd = [
        'TEtranscripts',
        '--treatment'] + treatment_bams + [
        '--control'] + control_bams + [
        '--GTF', gene_gtf,
        '--TE', te_gtf,
        '--mode', 'multi',
        '--sortByPos',
        '--stranded', stranded,
        '--project', output_prefix,
    ]
    subprocess.run(cmd, check=True)
    return f'{output_prefix}.cntTable'


def parse_tetranscripts_output(count_table):
    '''Parse TEtranscripts count table output.'''
    df = pd.read_csv(count_table, sep='\t', index_col=0)

    te_mask = df.index.str.contains(':')
    te_counts = df[te_mask]
    gene_counts = df[~te_mask]

    print(f'Total features: {len(df)}')
    print(f'  Genes: {len(gene_counts)}')
    print(f'  TE families: {len(te_counts)}')

    return gene_counts, te_counts


def te_family_summary(te_counts):
    '''Summarize TE expression by family and class.'''
    te_info = []
    for te_id in te_counts.index:
        parts = te_id.split(':')
        te_info.append({
            'te_id': te_id,
            'te_name': parts[0] if len(parts) > 0 else te_id,
            'te_class': parts[1] if len(parts) > 1 else 'unknown',
            'te_family': parts[2] if len(parts) > 2 else 'unknown',
            'mean_count': te_counts.loc[te_id].mean(),
        })

    summary = pd.DataFrame(te_info)
    class_summary = summary.groupby('te_class')['mean_count'].agg(['sum', 'count']).sort_values('sum', ascending=False)
    class_summary.columns = ['total_expression', 'num_families']

    print('\n=== TE Expression by Class ===')
    for cls, row in class_summary.iterrows():
        print(f'  {cls}: {row["num_families"]} families, total expression = {row["total_expression"]:.0f}')

    return summary, class_summary


def plot_te_expression(te_counts, top_n=30, output_file='te_expression.png'):
    '''Plot top expressed TE families.'''
    mean_expr = te_counts.mean(axis=1).sort_values(ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(10, 8))
    mean_expr.plot(kind='barh', ax=ax, color='steelblue')
    ax.set_xlabel('Mean Read Count')
    ax.set_title(f'Top {top_n} Expressed TE Families')
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    print('TE Expression Analysis')
    print('=' * 40)
    print('1. run_tetranscripts() - Run TEtranscripts DE analysis')
    print('2. parse_tetranscripts_output() - Parse count table')
    print('3. te_family_summary() - Summarize by TE class')
    print('4. plot_te_expression() - Plot top expressed TEs')
