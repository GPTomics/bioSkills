#!/usr/bin/env python3
'''Analyze ROSE super-enhancer output'''

import pandas as pd
import matplotlib.pyplot as plt
import sys

def parse_rose_output(se_file):
    '''Parse ROSE super-enhancer table'''
    df = pd.read_csv(se_file, sep='\t', skiprows=5)
    df.columns = ['REGION_ID', 'CHROM', 'START', 'STOP', 'NUM_LOCI', 'CONSTITUENT_SIZE', 'SIGNAL', 'isSuper']
    return df

def analyze_super_enhancers(se_file, output_prefix='se_analysis'):
    '''Analyze and visualize super-enhancers'''
    df = parse_rose_output(se_file)

    se = df[df['isSuper'] == 1]
    te = df[df['isSuper'] == 0]

    print(f'Total enhancers: {len(df)}')
    print(f'Super-enhancers: {len(se)}')
    print(f'Typical enhancers: {len(te)}')
    print(f'\nSuper-enhancer signal range: {se["SIGNAL"].min():.0f} - {se["SIGNAL"].max():.0f}')

    # Hockey stick plot
    df_sorted = df.sort_values('SIGNAL')
    df_sorted['rank'] = range(len(df_sorted))

    plt.figure(figsize=(8, 6))
    plt.scatter(df_sorted['rank'], df_sorted['SIGNAL'], c=df_sorted['isSuper'], cmap='coolwarm', s=10)
    plt.xlabel('Enhancer Rank')
    plt.ylabel('H3K27ac Signal')
    plt.title('Super-Enhancer Hockey Stick Plot')
    plt.savefig(f'{output_prefix}_hockey_stick.png', dpi=150)
    plt.close()

    # Size distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].hist(te['CONSTITUENT_SIZE'] / 1000, bins=50, alpha=0.7, label='Typical')
    axes[0].hist(se['CONSTITUENT_SIZE'] / 1000, bins=50, alpha=0.7, label='Super')
    axes[0].set_xlabel('Size (kb)')
    axes[0].set_ylabel('Count')
    axes[0].legend()
    axes[0].set_title('Enhancer Size Distribution')

    axes[1].boxplot([te['SIGNAL'], se['SIGNAL']], labels=['Typical', 'Super'])
    axes[1].set_ylabel('Signal')
    axes[1].set_title('Signal Comparison')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_stats.png', dpi=150)
    plt.close()

    # Export SE list
    se[['CHROM', 'START', 'STOP', 'SIGNAL']].to_csv(f'{output_prefix}_se_list.bed', sep='\t', index=False, header=False)
    print(f'\nSuper-enhancer BED: {output_prefix}_se_list.bed')

if __name__ == '__main__':
    se_file = sys.argv[1] if len(sys.argv) > 1 else 'rose_output/sample_SuperEnhancers.table.txt'
    analyze_super_enhancers(se_file)
