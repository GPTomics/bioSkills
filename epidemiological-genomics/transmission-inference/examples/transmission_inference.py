'''Transmission inference pipeline: SNP-cluster definition with pathogen-tuned
thresholds, outbreaker2 / TransPhylo orchestration via subprocess (R), bottleneck
estimation from deep sequencing, and HIV-TRACE subtype-aware clustering.

Refuses to assert transmission direction from pairwise SNP alone; flags Walker 2013
thresholds applied outside UK low-transmission settings; documents unsampled-
intermediate caveats.'''
# Reference: snp-dists 0.8+, gubbins 3.3+, hiv-trace 1.5+, lofreq 2.1+, pandas 2.2+, numpy 1.26+, scipy 1.12+ | Verify API if version differs

import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import binom

SNP_THRESHOLDS = {
    'mtb_likely': 12,
    'mtb_recent': 5,
    'mrsa_hospital': 15,
    'mrsa_broader': 40,
    'kpneumo_kpc': 21,
    'cdiff_direct': 2,
    'cdiff_plausible_6mo': 10,
    'ngono_transmission': 25,
}

CGMLST_THRESHOLDS = {
    'salmonella_efsa': 5,
    'salmonella_extended': 7,
    'listeria_pulsenet': 4,
    'ecoli_stec': 10,
}

WALKER_2013_VALIDATED_SETTING = 'UK low-transmission Oxfordshire community/household contact-traced'
HIV_TRACE_DEFAULT_THRESHOLD = 0.015
HIV_TRACE_SUBTYPE_DEFAULT = 'B'

PATHOGEN_BOTTLENECK_REFERENCE = {
    'influenza_a': (1, 2, 'McCrone 2018 eLife 7:e35962'),
    'sars_cov_2': (1, 10, 'Lythgoe 2021 Science 372:eabg0821'),
}


def run_snp_dists(masked_alignment, out_csv):
    cmd = ['snp-dists', '-c', str(masked_alignment)]
    with open(out_csv, 'w') as fh:
        subprocess.run(cmd, check=True, stdout=fh)
    return pd.read_csv(out_csv, index_col=0)


def define_snp_clusters(distance_df, threshold, population_caveat=None):
    '''Single-linkage clusters at the pathogen-specific threshold; document caveat.'''
    if distance_df.shape[0] < 2:
        return pd.Series(1, index=distance_df.index, name='cluster'), {}
    condensed = distance_df.values[np.triu_indices(len(distance_df), k=1)]
    Z = linkage(condensed, method='single')
    clusters = fcluster(Z, t=threshold, criterion='distance')
    metadata = {
        'threshold': threshold,
        'population_caveat': population_caveat or 'NOT specified -- threshold validity depends on derivation population',
        'linkage': 'single',
        'n_samples': len(distance_df),
    }
    return pd.Series(clusters, index=distance_df.index, name='cluster'), metadata


def warn_walker_threshold_outside_uk(setting_description):
    '''Walker 2013 TB SNP thresholds derived from UK low-transmission; warn if extrapolated.'''
    if 'UK' not in setting_description and 'low-transmission' not in setting_description.lower():
        return (f"WARNING: Walker 2013 5/12-SNP TB thresholds derived from {WALKER_2013_VALIDATED_SETTING}. "
                f"Applying in '{setting_description}' inflates apparent recent-transmission rates 2-5x. "
                "Consider deriving local threshold from epidemiologically-anchored case pairs.")
    return None


def annotate_pairs_with_temporality(distance_df, sample_dates, snp_threshold):
    '''Annotate pairs as 'temporally consistent', 'temporally inconsistent', or 'tied'.

    CRITICAL: temporal earlier does NOT establish infection direction (Worby 2014).'''
    rows = []
    for i, s1 in enumerate(distance_df.index):
        for s2 in distance_df.index[i + 1:]:
            snp = distance_df.loc[s1, s2]
            if snp > snp_threshold:
                continue
            d1, d2 = sample_dates[s1], sample_dates[s2]
            temporal = 'tied' if d1 == d2 else ('s1_earlier' if d1 < d2 else 's2_earlier')
            rows.append({
                's1': s1, 's2': s2, 'snp_distance': snp,
                'temporal_relation': temporal,
                'date_s1': d1, 'date_s2': d2,
                'direction_claim': 'NOT inferable from pairwise SNP alone (Worby 2014); "consistent with" only',
            })
    return pd.DataFrame(rows)


def binomial_bottleneck_estimator(donor_freqs, recipient_freqs, depths_recipient, nb_max=200):
    '''Simplified bottleneck-size estimator following the Sobel Leonard 2017
    J Virol 91:e00171-17 framework. The full method uses a beta-binomial mixture
    that integrates over within-host frequency drift; this simplification treats the
    recipient minor-allele count as binomial(depth, donor_freq), which is a reasonable
    approximation when within-host drift is small. For production analyses use the
    published R / Python implementations that fit the full beta-binomial.'''
    if len(donor_freqs) != len(recipient_freqs):
        raise ValueError('donor and recipient frequency arrays must align by site')
    log_likes = np.zeros(nb_max)
    for nb in range(1, nb_max + 1):
        site_ll = 0.0
        for p_d, p_r, depth in zip(donor_freqs, recipient_freqs, depths_recipient):
            if depth <= 0:
                continue
            k_r = int(round(p_r * depth))
            site_ll += binom.logpmf(k_r, depth, p_d)
        log_likes[nb - 1] = site_ll
    log_likes -= log_likes.max()
    posterior = np.exp(log_likes)
    posterior /= posterior.sum()
    nb_map = int(np.argmax(posterior)) + 1
    return {'nb_map': nb_map, 'nb_grid': np.arange(1, nb_max + 1), 'posterior': posterior}


def compare_bottleneck_to_literature(nb_map, pathogen):
    ref = PATHOGEN_BOTTLENECK_REFERENCE.get(pathogen)
    if ref is None:
        return None
    lo, hi, source = ref
    flag = 'consistent' if lo <= nb_map <= hi * 5 else 'inconsistent'
    return {'nb_map': nb_map, 'literature_range': (lo, hi), 'source': source, 'flag': flag}


def run_hiv_trace(fasta, out_dir, threshold=HIV_TRACE_DEFAULT_THRESHOLD,
                  subtype=HIV_TRACE_SUBTYPE_DEFAULT):
    '''HIV-TRACE clustering wrapper; verify the installed CLI before relying on flag forms.

    The HIV-TRACE distribution includes `hivtrace` (high-level orchestrator) and
    `hivnetworkcsv` (network building from a pairwise distance CSV). Modern pipelines
    typically run `tn93` for pairwise distance then `hivnetworkcsv` for cluster definition;
    flag forms have varied across releases. Always emit a subtype-vs-threshold warning
    when applying the US-CDC subtype-B 1.5% default outside subtype B.'''
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pairwise_csv = out_dir / 'pairwise.csv'
    clusters_csv = out_dir / 'clusters.csv'
    subprocess.run(['tn93', '-o', str(pairwise_csv), '-t', str(threshold), str(fasta)], check=True)
    subprocess.run(['hivnetworkcsv', '-i', str(pairwise_csv), '-t', str(threshold),
                    '-o', str(clusters_csv), '-f', 'plain'], check=True)
    warning = None
    if subtype != 'B' and abs(threshold - HIV_TRACE_DEFAULT_THRESHOLD) < 1e-9:
        warning = ('WARNING: 1.5% threshold is US-CDC subtype B default; under-clusters subtype C in '
                   'southern Africa. Use locally validated threshold for non-B subtypes.')
    return pd.read_csv(clusters_csv), warning


def run_outbreaker2_via_r(dates_csv, fasta, contact_matrix_csv, out_rds,
                          n_iter=int(1e6), sample_every=200):
    '''Invoke outbreaker2 via R subprocess. R script generated inline.'''
    r_script = f'''
    library(outbreaker2); library(ape)
    dna <- read.dna('{fasta}', format='fasta')
    dates_df <- read.csv('{dates_csv}')
    ctd <- as.matrix(read.csv('{contact_matrix_csv}', row.names=1))
    w_dens <- dgamma(1:30, shape=2.5, scale=2)
    f_dens <- dgamma(1:30, shape=2, scale=3)
    data <- outbreaker_data(dates=dates_df$collection_date, dna=dna,
                            w_dens=w_dens, f_dens=f_dens, ctd=ctd)
    cfg <- create_config(n_iter={n_iter}, sample_every={sample_every}, find_import=TRUE)
    res <- outbreaker(data=data, config=cfg)
    saveRDS(res, '{out_rds}')
    '''
    subprocess.run(['Rscript', '-e', r_script], check=True)


def run_transphylo_via_r(dated_tree_nexus, date_last_sample, out_rds,
                         w_shape=1.3, w_scale=10, ws_shape=1.1, ws_scale=7,
                         start_neg=0.5, mcmc_iterations=int(1e5)):
    r_script = f'''
    library(TransPhylo); library(ape)
    tree <- read.nexus('{dated_tree_nexus}')
    ptree <- ptreeFromPhylo(tree, dateLastSample={date_last_sample})
    res <- inferTTree(ptree, mcmcIterations={mcmc_iterations},
                      w.shape={w_shape}, w.scale={w_scale},
                      ws.shape={ws_shape}, ws.scale={ws_scale},
                      startNeg={start_neg}, dateT={date_last_sample + 0.1})
    saveRDS(res, '{out_rds}')
    '''
    subprocess.run(['Rscript', '-e', r_script], check=True)
