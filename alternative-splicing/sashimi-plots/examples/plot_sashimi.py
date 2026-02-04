#!/usr/bin/env python3
'''
Sashimi plot generation using ggsashimi.
Visualizes splicing events with read coverage and junction counts.
'''

import subprocess
import pandas as pd
from pathlib import Path


def create_grouping_file(bam_files, conditions, colors, output_file):
    '''
    Create ggsashimi grouping file.

    Args:
        bam_files: List of BAM file paths
        conditions: List of condition labels (same length as bam_files)
        colors: Dict mapping condition to hex color
        output_file: Output TSV path
    '''
    groups = pd.DataFrame({
        'bam': bam_files,
        'group': conditions,
        'color': [colors[c] for c in conditions]
    })
    groups.to_csv(output_file, sep='\t', index=False, header=False)
    return output_file


def plot_sashimi(grouping_file, region, output_prefix, gtf_file, options=None):
    '''
    Generate sashimi plot for a genomic region.

    Args:
        grouping_file: TSV with bam, group, color columns
        region: Genomic coordinates (chr:start-end)
        output_prefix: Output file prefix
        gtf_file: GTF annotation file
        options: Dict of additional ggsashimi options
    '''
    if options is None:
        options = {}

    cmd = [
        'ggsashimi.py',
        '-b', grouping_file,
        '-c', region,
        '-o', output_prefix,
        '-g', gtf_file,
        '-M', str(options.get('min_junc', 10)),  # Min junction reads to show
        '--alpha', str(options.get('alpha', 0.25)),  # Coverage transparency
        '--height', str(options.get('height', 3)),
        '--width', str(options.get('width', 8)),
        '-F', options.get('format', 'pdf')
    ]

    if options.get('shrink', True):
        cmd.append('--shrink')
    if options.get('fix_y_scale', True):
        cmd.append('--fix-y-scale')
    if options.get('aggregate'):
        cmd.extend(['-A', options['aggregate']])  # mean, median

    subprocess.run(cmd, check=True)
    print(f'Sashimi plot saved: {output_prefix}.{options.get("format", "pdf")}')


def batch_plot_rmats_events(rmats_file, grouping_file, gtf_file, output_dir,
                            n_top=20, fdr_cutoff=0.05, dpsi_cutoff=0.1):
    '''
    Generate sashimi plots for top differential splicing events from rMATS.

    Args:
        rmats_file: rMATS SE.MATS.JC.txt output
        grouping_file: ggsashimi grouping file
        gtf_file: GTF annotation
        output_dir: Directory for output plots
        n_top: Number of top events to plot
        fdr_cutoff: FDR significance threshold
        dpsi_cutoff: Minimum |deltaPSI|
    '''
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(rmats_file, sep='\t')

    # Filter significant events
    significant = df[
        (df['FDR'] < fdr_cutoff) &
        (df['IncLevelDifference'].abs() > dpsi_cutoff)
    ].copy()

    # Sort by significance and effect size
    significant['score'] = -significant['FDR'].apply(lambda x: max(x, 1e-300)).apply(pd.np.log10) * significant['IncLevelDifference'].abs()
    significant = significant.nlargest(n_top, 'score')

    print(f'Plotting top {len(significant)} significant events')

    for idx, event in significant.iterrows():
        chrom = event['chr']
        gene = event['geneSymbol']

        # Extend region around the skipped exon
        # Include flanking exons for context
        start = event['upstreamES'] - 500
        end = event['downstreamEE'] + 500
        region = f'{chrom}:{start}-{end}'

        output_prefix = f'{output_dir}/{gene}_{chrom}_{event["exonStart_0base"]}'

        try:
            plot_sashimi(
                grouping_file, region, output_prefix, gtf_file,
                options={'shrink': True, 'fix_y_scale': True, 'min_junc': 5}
            )
        except subprocess.CalledProcessError as e:
            print(f'Failed to plot {gene}: {e}')


def plot_specific_event(grouping_file, gtf_file, chrom, start, end,
                        output_prefix, gene_name=None, flank=500):
    '''
    Plot a specific genomic region with optional flanking sequence.

    Args:
        grouping_file: ggsashimi grouping file
        gtf_file: GTF annotation
        chrom: Chromosome
        start: Start position
        end: End position
        output_prefix: Output file prefix
        gene_name: Optional gene name for labeling
        flank: Base pairs to add on each side
    '''
    region = f'{chrom}:{start - flank}-{end + flank}'

    plot_sashimi(
        grouping_file, region, output_prefix, gtf_file,
        options={
            'shrink': True,
            'fix_y_scale': True,
            'min_junc': 5,
            'aggregate': 'mean',
            'height': 4,
            'width': 10
        }
    )


if __name__ == '__main__':
    # Example: Create grouping file
    bams = ['ctrl1.bam', 'ctrl2.bam', 'treat1.bam', 'treat2.bam']
    conditions = ['Control', 'Control', 'Treatment', 'Treatment']
    colors = {'Control': '#1f77b4', 'Treatment': '#ff7f0e'}

    # create_grouping_file(bams, conditions, colors, 'sashimi_groups.tsv')

    # Example: Plot single region
    # plot_sashimi('sashimi_groups.tsv', 'chr1:1000000-1010000', 'example', 'annotation.gtf')

    # Example: Batch plot rMATS results
    # batch_plot_rmats_events(
    #     'rmats_output/SE.MATS.JC.txt',
    #     'sashimi_groups.tsv',
    #     'annotation.gtf',
    #     'sashimi_plots/'
    # )

    print('Configure BAM files and run functions to generate sashimi plots')
