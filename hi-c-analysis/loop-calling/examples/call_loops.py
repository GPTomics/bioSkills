'''Call loops de-novo with cooltools dots, then confirm the set with an APA score.'''
# Reference: cooltools 0.7+, cooler 0.10+, bioframe 0.7+ | Verify API if version differs
import cooler
import cooltools
import bioframe
import numpy as np

LOOP_RES = 10_000          # 10kb: a realistic de-novo floor; finer bins need much deeper maps (Rao 2014)
MAX_SEPARATION = 10_000_000  # ignore pixels farther than 10Mb from the diagonal (cooltools default ceiling)
N_LAMBDA_BINS = 40         # geometric expected-value bins; FDR is run independently per bin (HiCCUPS)
LAMBDA_BIN_FDR = 0.1       # per-lambda-bin BH-FDR; keeps FDR honest across the count dynamic range
CLUSTERING_RADIUS = 20_000  # merge called pixels within 20kb into one loop (cooltools default)
APA_FLANK = 100_000        # +/-100kb window around each anchor pair for the pileup
CORNER = 3                 # 3x3 lower-left corner block as the APA on-vs-off control (Rao 2014)
NPROC = 4

clr = cooler.Cooler(f'matrix.mcool::/resolutions/{LOOP_RES}')
print(f'Loaded at {clr.binsize}bp resolution; balanced={"weight" in clr.bins().columns}')

arms = bioframe.make_viewframe(clr.chromsizes)   # one region per chromosome; reused for expected and dots
expected = cooltools.expected_cis(clr, view_df=arms, nproc=NPROC)

loops = cooltools.dots(
    clr, expected=expected, view_df=arms,
    max_loci_separation=MAX_SEPARATION,
    n_lambda_bins=N_LAMBDA_BINS, lambda_bin_fdr=LAMBDA_BIN_FDR,
    clustering_radius=CLUSTERING_RADIUS, nproc=NPROC,
)
print(f'Called {len(loops)} loops')

if len(loops) == 0:
    print('No loops. The map is likely too shallow for de-novo calling -- run APA on known CTCF/cohesin anchors instead.')
    raise SystemExit

loops['size'] = (loops['start2'] - loops['start1']).abs()
print(f'Loop size: median {loops["size"].median()/1000:.0f}kb, mean {loops["size"].mean()/1000:.0f}kb')

anchors = loops[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
stack = cooltools.pileup(clr, anchors, view_df=arms, expected_df=expected, flank=APA_FLANK, nproc=NPROC)
apa = np.nanmean(stack, axis=0)   # pileup returns (n_snippets, D, D); average over axis 0 -> 2D aggregate (O/E because expected was passed)

center = apa.shape[0] // 2
apa_score = apa[center, center] / np.nanmean(apa[-CORNER:, :CORNER])   # center vs lower-left corner; >1 = enriched
print(f'APA score (center / lower-left corner): {apa_score:.2f}  (>1 means the call set is aggregate-enriched)')

loops[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].to_csv('loops.bedpe', sep='\t', index=False, header=False)
print('Saved loops.bedpe')
