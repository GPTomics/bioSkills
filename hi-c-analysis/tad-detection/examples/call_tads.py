'''Multi-scale TAD boundary calling from a balanced Hi-C cooler.

Encodes the load-bearing reframe: the boundary is the reproducible unit, not the
domain. Runs a LIST of diamond windows (the window is the scale dial), ranks
boundaries by the continuous boundary_strength (valley prominence) rather than the
per-dataset is_boundary flag, and inspects valid-pixel coverage before trusting a
boundary in a sparse locus.'''
# Reference: cooler 0.10+, cooltools 0.7+ | Verify API if version differs

import cooler
import cooltools

COOLER_URI = 'matrix.mcool::/resolutions/10000'   # single-resolution URI; .mcool is multi-resolution
WINDOW_MULTIPLES = [3, 5, 10, 25]   # 3x=sub-TAD, 25x=compartment-domain; <3x bin is pure noise
REPORT_MULTIPLE = 10                # ~10x bin: mammalian interphase sweet spot for the headline call

clr = cooler.Cooler(COOLER_URI)
res = clr.binsize
print(f'Loaded at {res} bp resolution')

windows = [m * res for m in WINDOW_MULTIPLES]
ins = cooltools.insulation(clr, windows, verbose=True)   # clr_weight_name='weight' default -> matrix MUST be balanced

report_w = REPORT_MULTIPLE * res
strength_col = f'boundary_strength_{report_w}'
flag_col = f'is_boundary_{report_w}'
valid_col = f'n_valid_pixels_{report_w}'

for m in WINDOW_MULTIPLES:
    w = m * res
    n = int(ins[f'is_boundary_{w}'].sum())
    print(f'window {w // 1000:>4}kb ({m:>2}x bin): {n} boundaries above the Li threshold')

scored = ins.dropna(subset=[strength_col]).copy()
ranked = scored.sort_values(strength_col, ascending=False)   # continuous prominence: comparable across samples
print(f'\nTop boundaries at the {report_w // 1000}kb window (ranked by valley prominence):')
print(ranked[['chrom', 'start', 'end', strength_col, valid_col]].head(10).to_string(index=False))

flagged = ins[ins[flag_col].fillna(False)].copy()
flagged[['chrom', 'start', 'end', strength_col]].to_csv('boundaries.bed', sep='\t', index=False, header=False)
ins[['chrom', 'start', 'end', f'log2_insulation_score_{report_w}']].to_csv('insulation.bedgraph', sep='\t', index=False, header=False)
print(f'\nWrote {len(flagged)} flagged boundaries to boundaries.bed and the insulation track to insulation.bedgraph')
print('Note: is_boundary uses a per-dataset Li threshold; for cross-sample comparison use boundary_strength, not this flag.')
