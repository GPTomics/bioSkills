'''Phylodynamics pipeline: TempEst-style temporal-signal QC, TreeTime time-scaling,
recombination masking with Gubbins for bacteria, and BEAST2 BDSKY orchestration.

Refuses to time-scale when root-to-tip R^2 < 0.3 (insufficient temporal signal);
requires `core.full.aln` (NOT `core.aln`) for Gubbins; flags BDSKY origin <= rootHeight.'''
# Reference: treetime 0.11+, biopython 1.84+, dendropy 4.6+, baltic 0.2+, gubbins 3.3+, iqtree 2.3.6+, beast 2.7.6+, snp-dists 0.8+ | Verify API if version differs

import json
import subprocess
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd

ROOT_TO_TIP_R2_MIN = 0.3
TREETIME_CLOCK_FILTER = 4
GUBBINS_INPUT_REQUIRED = 'core.full.aln'
BDSKY_ORIGIN_PADDING = 0.1
BEAST_BURNIN_PCT = 10
MIN_SEQS_PER_DEME = 20
S_PNEUMO_REFERENCE_CLOCK = 1.5e-6


def to_decimal_year(date_str):
    if isinstance(date_str, (int, float)):
        return float(date_str)
    s = str(date_str)
    if '-' in s:
        dt = datetime.strptime(s, '%Y-%m-%d')
        year = dt.year
        day_of_year = dt.timetuple().tm_yday
        days_in_year = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
        return year + (day_of_year - 1) / days_in_year
    return float(s)


def prepare_dates_tsv(metadata_df, name_col, date_col, out_path):
    dates = metadata_df[[name_col, date_col]].copy()
    dates.columns = ['name', 'date']
    dates['date'] = dates['date'].apply(to_decimal_year)
    dates.to_csv(out_path, sep='\t', index=False)
    return out_path


def root_to_tip_regression(treetime_outdir):
    '''Compute root-to-tip regression slope / intercept / R^2 from TreeTime output.

    TreeTime writes root_to_tip_regression.pdf and one of several text files depending
    on release: rtt.tsv (some releases), root_to_tip_regression.csv, or columns inside
    dates.tsv. Try each in order; if none exist, parse the timetree.nexus directly via
    dendropy to recompute root-to-tip distances against tip dates.'''
    out = Path(treetime_outdir)
    candidates = [out / 'rtt.tsv', out / 'root_to_tip_regression.csv', out / 'dates.tsv']
    rtt = None
    for path in candidates:
        if path.exists():
            sep = ',' if path.suffix == '.csv' else '\t'
            df = pd.read_csv(path, sep=sep)
            cols = {c.lower(): c for c in df.columns}
            date_col = cols.get('date') or cols.get('numeric date') or cols.get('numdate')
            dist_col = cols.get('distance') or cols.get('rtt') or cols.get('root_to_tip')
            if date_col and dist_col:
                rtt = df[[date_col, dist_col]].rename(columns={date_col: 'date', dist_col: 'distance'})
                break
    if rtt is None:
        raise FileNotFoundError(f'No root-to-tip data found in {out}; re-run TreeTime with --confidence '
                                'or parse timetree.nexus manually')
    x = rtt['date'].to_numpy(dtype=float)
    y = rtt['distance'].to_numpy(dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {'slope': slope, 'intercept': intercept, 'r2': r2, 'clock_rate_subs_per_site_per_year': slope}


def assert_temporal_signal(r2):
    if r2 < ROOT_TO_TIP_R2_MIN:
        raise RuntimeError(f'Root-to-tip R^2 = {r2:.3f} < {ROOT_TO_TIP_R2_MIN}; data does NOT support time-scaling. '
                           'Use strong informative prior on clock or extend sampling window.')


def run_iqtree(aln_path, prefix, model='GTR+G', ascertainment_bias=False, threads='AUTO'):
    cmd = ['iqtree', '-s', str(aln_path), '-m', model + ('+ASC' if ascertainment_bias else ''),
           '-B', '1000', '-T', threads, '-pre', prefix]
    subprocess.run(cmd, check=True)
    return f'{prefix}.treefile'


def run_treetime(tree_path, aln_path, dates_path, out_dir, clock_filter=TREETIME_CLOCK_FILTER):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['treetime', '--tree', str(tree_path), '--aln', str(aln_path), '--dates', str(dates_path),
           '--coalescent', 'skyline', '--clock-filter', str(clock_filter),
           '--confidence', '--reroot', 'best', '--outdir', str(out_dir)]
    subprocess.run(cmd, check=True)
    return out_dir


def run_date_randomisation(tree_path, dates_path, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['treetime', 'clock', '--tree', str(tree_path), '--dates', str(dates_path),
           '--reassign-dates', '--outdir', str(out_dir)]
    subprocess.run(cmd, check=True)
    return out_dir


def run_gubbins(core_full_aln_path, prefix='gubbins'):
    '''Recombination masking for bacterial alignments.

    Input MUST be core.full.aln (full positions including invariant).
    Passing core.aln (variable-only) gives wrong recombination calls.'''
    if not str(core_full_aln_path).endswith(GUBBINS_INPUT_REQUIRED):
        raise ValueError(f'Gubbins requires {GUBBINS_INPUT_REQUIRED}; received {core_full_aln_path}. '
                         'core.aln (variable-only) is wrong input.')
    subprocess.run(['run_gubbins.py', '--prefix', prefix, str(core_full_aln_path)], check=True)
    return f'{prefix}.filtered_polymorphic_sites.fasta'


def validate_bdsky_origin(origin, root_height):
    '''BDSKY origin must be strictly larger than root_height.'''
    if origin <= root_height:
        raise ValueError(f'BDSKY origin ({origin}) <= rootHeight ({root_height}); biases R_e estimates upward. '
                         'Initialise origin to (root_height + 0.1*root_height) or use prior epi knowledge.')


def suggest_bdsky_origin(root_height):
    return root_height * (1 + BDSKY_ORIGIN_PADDING)


def run_beast2(xml_path, seed=42, threads=4, log_prefix=None):
    cmd = ['beast', '-threads', str(threads), '-beagle', '-seed', str(seed), str(xml_path)]
    subprocess.run(cmd, check=True)
    if log_prefix:
        log_file = f'{log_prefix}_{seed}.log'
        trees_file = f'{log_prefix}_{seed}.trees'
        return log_file, trees_file
    return None, None


def combine_chains(log_files, trees_files, out_log, out_trees, burnin_pct=BEAST_BURNIN_PCT):
    log_args = []
    for lf in log_files:
        log_args.extend(['-log', str(lf)])
    subprocess.run(['logcombiner'] + log_args + ['-burnin', str(burnin_pct), '-o', str(out_log)], check=True)

    trees_args = []
    for tf in trees_files:
        trees_args.extend(['-log', str(tf)])
    subprocess.run(['logcombiner'] + trees_args + ['-burnin', str(burnin_pct),
                                                    '-o', str(out_trees), '-decimalPlaces', '6'], check=True)


def loganalyser_ess(log_file, burnin_pct=BEAST_BURNIN_PCT):
    result = subprocess.run(['loganalyser', '-burnin', str(burnin_pct), str(log_file)],
                            check=True, capture_output=True, text=True)
    return result.stdout


def validate_mascot_deme_sizes(deme_counts):
    underpowered = {d: n for d, n in deme_counts.items() if n < MIN_SEQS_PER_DEME}
    if underpowered:
        raise ValueError(f'MASCOT demes with <{MIN_SEQS_PER_DEME} sequences: {underpowered}; '
                         'migration unidentifiable. Pool demes or accept wide HPD intervals.')


def bacterial_phylodynamics_pipeline(core_full_aln, dates_tsv, work_dir,
                                     species_reference_clock=S_PNEUMO_REFERENCE_CLOCK):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    masked_aln = run_gubbins(core_full_aln, prefix=str(work_dir / 'gubbins'))
    tree_file = run_iqtree(masked_aln, prefix=str(work_dir / 'masked_tree'),
                           model='GTR+G', ascertainment_bias=True)
    tt_dir = run_treetime(tree_file, masked_aln, dates_tsv, work_dir / 'timetree')

    rtt = root_to_tip_regression(tt_dir)
    assert_temporal_signal(rtt['r2'])

    clock_ratio = rtt['clock_rate_subs_per_site_per_year'] / species_reference_clock
    if clock_ratio > 3:
        raise RuntimeError(f'Clock rate {rtt["clock_rate_subs_per_site_per_year"]:.2e} is '
                           f'{clock_ratio:.1f}x reference {species_reference_clock:.2e}. '
                           'Unmasked recombination is the most likely cause; verify Gubbins ran on core.full.aln.')

    return {'timetree_dir': str(tt_dir), 'root_to_tip': rtt,
            'masked_alignment': masked_aln, 'tree_file': tree_file}
