'''Call chromatin loops from Hi-C data'''

import cooler
import cooltools
import bioframe
import numpy as np

clr = cooler.Cooler('matrix.mcool::resolutions/10000')
print(f'Loaded at {clr.binsize}bp resolution')

view_df = bioframe.make_viewframe(clr.chromsizes)

print('Computing expected values...')
expected = cooltools.expected_cis(clr, view_df=view_df, ignore_diags=2)

print('Calling loops...')
dots = cooltools.dots(
    clr,
    expected=expected,
    view_df=view_df,
    max_loci_separation=2000000,
)

print(f'\nFound {len(dots)} loops')

if len(dots) > 0:
    dots['size'] = abs(dots['end2'] - dots['start1'])
    print(f'\nLoop size statistics:')
    print(f'  Mean: {dots["size"].mean()/1000:.0f} kb')
    print(f'  Median: {dots["size"].median()/1000:.0f} kb')

    dots.to_csv('loops.bedpe', sep='\t', index=False)
    print('\nSaved to loops.bedpe')
else:
    print('No loops found. Try adjusting parameters.')
