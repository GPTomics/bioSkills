'''Label-free normalization and TMT cross-plex IRS bridge, self-contained.'''
# Reference: numpy 1.26+, pandas 2.2+ | Verify API if version differs
import numpy as np
import pandas as pd

rng = np.random.default_rng(0)


def median_center(log_int):
    '''Subtract each sample median log2 intensity (corrects LOCATION only, cannot manufacture variance).'''
    sample_medians = log_int.median(axis=0)
    return log_int - sample_medians + sample_medians.median()


def sample_loading_normalize(plex):
    '''Per-channel scalar correcting total load WITHIN one plex; target = mean column sum.'''
    target = plex.sum(axis=0).mean()
    return plex * (target / plex.sum(axis=0))


def irs_scale(plexes, ref_cols):
    '''Internal Reference Scaling bridge (Plubell 2017): pin each plex reference channel to a common per-protein value.'''
    refs = pd.concat([p[ref] for p, ref in zip(plexes, ref_cols)], axis=1)
    geomean = np.exp(np.log(refs.replace(0, np.nan)).mean(axis=1))    # per-protein geometric mean across plexes
    out = []
    for p, ref in zip(plexes, ref_cols):
        factor = geomean / p[ref]    # per-protein per-plex scaling factor
        out.append(p.mul(factor, axis=0))
    return out


proteins = [f'P{i:03d}' for i in range(200)]
samples = [f'S{i}' for i in range(6)]

# Label-free: raw intensities with per-sample loading offsets and ~15% MaxQuant zeros (missing)
true_abundance = rng.normal(20, 2, (200, 6))
loading_offset = np.array([0.0, 0.4, -0.3, 0.6, -0.5, 0.2])    # per-sample load differences to be normalized out
raw = 2 ** (true_abundance + loading_offset)
raw[rng.random((200, 6)) < 0.15] = 0    # MaxQuant writes 0 for "not quantified"
lfq = pd.DataFrame(raw, index=proteins, columns=samples)

log_int = np.log2(lfq.replace(0, np.nan))    # 0 -> NaN before transform; log2(0) = -inf otherwise
print(f'Missing values: {100 * log_int.isna().sum().sum() / log_int.size:.1f}%')
print(f'Sample medians before centering: {np.round(log_int.median().values, 2)}')

normalized = median_center(log_int)
print(f'Sample medians after centering:  {np.round(normalized.median().values, 2)}')

# Filter proteins quantified in too few samples; min-peptides logic mirrors min-ratio-count=2 robustness
valid_per_protein = normalized.notna().sum(axis=1)
min_valid = len(samples) // 2    # require presence in >=50% of samples
filtered = normalized[valid_per_protein >= min_valid]
print(f'Proteins after >=50% valid filter: {len(filtered)}')

# TMT cross-plex: two 4-channel plexes, channel 3 is a pooled reference; intensities NOT comparable across plexes raw
plex_a = pd.DataFrame(2 ** rng.normal(18, 1.5, (200, 4)), index=proteins, columns=['A1', 'A2', 'A3', 'A_ref'])
plex_b = pd.DataFrame(2 ** rng.normal(19, 1.5, (200, 4)), index=proteins, columns=['B1', 'B2', 'B3', 'B_ref'])    # ~2x higher: random elution sampling, not biology
sl = [sample_loading_normalize(plex_a), sample_loading_normalize(plex_b)]
ref_cv_before = (sl[0]['A_ref'] / sl[1]['B_ref']).std()
bridged = irs_scale(sl, ['A_ref', 'B_ref'])
ref_cv_after = (bridged[0]['A_ref'] / bridged[1]['B_ref']).std()
print(f'Reference-channel ratio SD across plexes before IRS: {ref_cv_before:.3f}, after IRS: {ref_cv_after:.3f}')
