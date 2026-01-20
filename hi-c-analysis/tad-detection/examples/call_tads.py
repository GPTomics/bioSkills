'''Call TADs from Hi-C data'''

import cooler
import cooltools
import bioframe
import pandas as pd

clr = cooler.Cooler('matrix.mcool::resolutions/10000')
print(f'Loaded at {clr.binsize}bp resolution')

view_df = bioframe.make_viewframe(clr.chromsizes)

print('Computing insulation score...')
insulation = cooltools.insulation(clr, window_bp=[200000], ignore_diags=2)

boundaries = insulation[insulation['is_boundary_200000'] == True]
print(f'\nFound {len(boundaries)} TAD boundaries')

boundaries[['chrom', 'start', 'end', 'boundary_strength_200000']].to_csv(
    'boundaries.bed', sep='\t', index=False, header=False
)

tads = []
for chrom in clr.chromnames[:22]:
    chr_bounds = boundaries[boundaries['chrom'] == chrom].sort_values('start')
    starts = [0] + list(chr_bounds['start'])
    ends = list(chr_bounds['start']) + [clr.chromsizes[chrom]]
    for s, e in zip(starts, ends):
        if e > s:
            tads.append({'chrom': chrom, 'start': s, 'end': e})

tads_df = pd.DataFrame(tads)
tads_df['size'] = tads_df['end'] - tads_df['start']

print(f'\nTAD statistics:')
print(f'  Total TADs: {len(tads_df)}')
print(f'  Mean size: {tads_df["size"].mean()/1000:.0f} kb')
print(f'  Median size: {tads_df["size"].median()/1000:.0f} kb')

tads_df[['chrom', 'start', 'end']].to_csv('tads.bed', sep='\t', index=False, header=False)
print('\nSaved to tads.bed and boundaries.bed')

# Expected output: ~2000-4000 TADs genome-wide for human at 10kb resolution
# Typical TAD sizes: 200kb-2Mb (median ~800kb), with most between 400kb-1.5Mb
# Boundary strength: strong boundaries have values < -0.5 (more negative = stronger)
