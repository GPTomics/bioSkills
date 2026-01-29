'''Generate genome browser track figures with pyGenomeTracks'''

import subprocess
import os

# --- ALTERNATIVE: Use real data ---
# Download example files from ENCODE:
# wget https://www.encodeproject.org/files/ENCFF000ABC/@@download/ENCFF000ABC.bigWig
# wget https://www.encodeproject.org/files/ENCFF000XYZ/@@download/ENCFF000XYZ.bed

# Track configuration template
# Customize colors, heights, and file paths for your data
TRACKS_INI = '''
[x-axis]
where = top
fontsize = 10

[spacer]
height = 0.2

[treatment]
file = {treatment_bw}
title = Treatment
height = 4
color = #E64B35
min_value = 0
max_value = {max_value}
number_of_bins = 700
nans_to_zeros = true

[control]
file = {control_bw}
title = Control
height = 4
color = #4DBBD5
min_value = 0
max_value = {max_value}
number_of_bins = 700
nans_to_zeros = true

[spacer]
height = 0.3

[peaks]
file = {peaks_bed}
title = Differential Peaks
height = 1.5
color = #00A087
display = collapsed
labels = false

[spacer]
height = 0.3

[genes]
file = {genes_gtf}
title = Genes
height = 5
fontsize = 10
style = UCSC
prefered_name = gene_name
color = navy
'''


def generate_track_plot(region, output_file, tracks_ini='tracks.ini', dpi=300, width=40):
    '''Generate track plot for a genomic region

    Args:
        region: Genomic region in format chr:start-end
        output_file: Output file path (supports .png, .pdf, .svg)
        tracks_ini: Path to tracks configuration file
        dpi: Resolution for raster output (default 300 for publication)
        width: Plot width in cm (default 40)
    '''
    cmd = [
        'pyGenomeTracks',
        '--tracks', tracks_ini,
        '--region', region,
        '--outFileName', output_file,
        '--dpi', str(dpi),
        '--width', str(width)
    ]

    print(f'Generating: {output_file}')
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f'Error: {result.stderr}')
        return False
    return True


def batch_generate_plots(regions_file, output_dir, tracks_ini='tracks.ini'):
    '''Generate track plots for multiple regions from BED file

    Args:
        regions_file: BED file with regions (chr, start, end, name)
        output_dir: Directory for output files
        tracks_ini: Path to tracks configuration
    '''
    os.makedirs(output_dir, exist_ok=True)

    with open(regions_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            chrom, start, end = parts[0], parts[1], parts[2]
            name = parts[3] if len(parts) > 3 else f'{chrom}_{start}_{end}'

            region = f'{chrom}:{start}-{end}'
            output_file = os.path.join(output_dir, f'{name}.png')

            generate_track_plot(region, output_file, tracks_ini)


if __name__ == '__main__':
    # Example usage - customize paths for your data

    # Single region
    # generate_track_plot(
    #     region='chr1:1000000-2000000',
    #     output_file='myc_locus.png',
    #     tracks_ini='tracks.ini'
    # )

    # Batch from BED file
    # batch_generate_plots(
    #     regions_file='regions_of_interest.bed',
    #     output_dir='track_figures',
    #     tracks_ini='tracks.ini'
    # )

    print('pyGenomeTracks helper script')
    print('Usage:')
    print('  1. Edit tracks.ini with your file paths')
    print('  2. Run: pyGenomeTracks --tracks tracks.ini --region chr1:1-1000000 -o output.png')
    print('')
    print('Or use this script functions:')
    print('  generate_track_plot("chr1:1000000-2000000", "output.png")')
    print('  batch_generate_plots("regions.bed", "output_dir/")')
