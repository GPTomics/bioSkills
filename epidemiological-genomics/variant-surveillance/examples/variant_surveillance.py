'''Variant surveillance pipeline: Pangolin UShER mode with full version pinning,
Nextclade with dataset tag recording, Freyja wastewater deconvolution with barcode
forward-only date validation, COJAC early detection, and Pango / Nextclade
reconciliation with recombinant X-prefix awareness.

Refuses to run Pangolin in pangoLEARN mode (deprecated mid-2023); validates
Freyja barcode date against sample collection date; flags ARTIC dropout amplicons.'''
# Reference: pangolin 4.3+, nextclade 3.8+, freyja 1.4+, cojac 0.9+, augur 24.0+, samtools 1.20+, pandas 2.2+, jq 1.7+ | Verify API if version differs

import json
import subprocess
from datetime import date
from pathlib import Path
import pandas as pd

PANGOLIN_DEPRECATED_MODE = 'pangolearn'
PANGOLIN_DEFAULT_MODE = 'usher'
ARTIC_V41_KNOWN_DROPOUTS = {64, 76, 88, 89, 90}
FREYJA_RESID_THRESHOLD = 0.10
COJAC_LEAD_TIME_DAYS = 13


def run_pangolin(seq_fasta, out_csv, version_log, analysis_mode=PANGOLIN_DEFAULT_MODE):
    if analysis_mode == PANGOLIN_DEPRECATED_MODE:
        raise ValueError('pangoLEARN was deprecated mid-2023; use UShER mode (Pongmoragot 2024)')
    subprocess.run(['pangolin', str(seq_fasta), '--analysis-mode', analysis_mode,
                    '--outfile', str(out_csv)], check=True)
    with open(version_log, 'w') as fh:
        subprocess.run(['pangolin', '--all-versions'], check=True, stdout=fh)
    return pd.read_csv(out_csv)


def get_nextclade_dataset(name, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(['nextclade', 'dataset', 'get', '--name', name, '--output-dir', str(out_dir)], check=True)
    pathogen_json = out_dir / 'pathogen.json'
    with open(pathogen_json) as fh:
        meta = json.load(fh)
    return meta.get('tag') or meta.get('version', 'unknown')


def run_nextclade(seq_fasta, dataset_dir, out_tsv, metadata_log):
    cmd = ['nextclade', 'run', '--input-dataset', str(dataset_dir),
           '--output-tsv', str(out_tsv), str(seq_fasta)]
    subprocess.run(cmd, check=True)
    with open(Path(dataset_dir) / 'pathogen.json') as fh:
        meta = json.load(fh)
    with open(metadata_log, 'w') as fh:
        fh.write(f"nextclade_dataset_tag: {meta.get('tag', meta.get('version'))}\n")
    return pd.read_csv(out_tsv, sep='\t')


def reconcile_pangolin_nextclade(pango_df, nc_df, key_pango='taxon', key_nc='seqName'):
    merged = pango_df.merge(nc_df, left_on=key_pango, right_on=key_nc, how='outer')
    merged['recombinant_candidate'] = merged['lineage'].fillna('').str.startswith('X')
    def status(row):
        l = row.get('lineage')
        c = row.get('clade')
        if not isinstance(l, str) and not isinstance(c, str):
            return 'both_missing'
        if not isinstance(l, str):
            return 'pango_missing'
        if not isinstance(c, str):
            return 'nextclade_missing'
        return 'present'
    merged['reconciliation_status'] = merged.apply(status, axis=1)
    return merged


def validate_freyja_barcode_date(barcode_date, sample_collection_date):
    '''Freyja barcode date MUST postdate sample collection.'''
    if isinstance(barcode_date, str):
        barcode_date = date.fromisoformat(barcode_date)
    if isinstance(sample_collection_date, str):
        sample_collection_date = date.fromisoformat(sample_collection_date)
    if barcode_date < sample_collection_date:
        raise RuntimeError(
            f'Freyja barcode ({barcode_date}) predates sample collection ({sample_collection_date}). '
            'Lineages designated after the barcode date are silently invisible. '
            'Run `freyja barcode-build` from a current UShER tree.')


def freyja_variants_and_demix(bam, reference, out_prefix):
    variants_tsv = f'{out_prefix}.variants.tsv'
    depths_tsv = f'{out_prefix}.depths.tsv'
    demix_tsv = f'{out_prefix}.demix.tsv'
    subprocess.run(['freyja', 'variants', str(bam), '--variants', variants_tsv,
                    '--depths', depths_tsv, '--ref', str(reference)], check=True)
    subprocess.run(['freyja', 'demix', variants_tsv, depths_tsv, '--output', demix_tsv], check=True)
    return _parse_freyja_demix(demix_tsv)


def _parse_freyja_demix(demix_tsv):
    '''Parse Freyja demix output (semicolon-separated lineages / abundances).'''
    df = pd.read_csv(demix_tsv, sep='\t', header=None, names=['key', 'value']).set_index('key')
    lineages = df.loc['lineages', 'value'].split()
    abundances = [float(x) for x in df.loc['abundances', 'value'].split()]
    resid = float(df.loc['resid', 'value'])
    summarized = df.loc['summarized', 'value'] if 'summarized' in df.index else ''
    return {
        'lineages': dict(zip(lineages, abundances)),
        'resid': resid,
        'summarized': summarized,
        'novel_lineage_warning': resid > FREYJA_RESID_THRESHOLD,
    }


def check_artic_dropouts(bam, primer_bed, min_coverage=20, scheme_version='V4.1'):
    '''Per-amplicon coverage check; report amplicons below min_coverage and known chronic dropouts.'''
    result = subprocess.run(['samtools', 'depth', '-aa', str(bam)],
                            check=True, capture_output=True, text=True)
    depth_df = pd.DataFrame([line.split('\t') for line in result.stdout.strip().split('\n')],
                            columns=['chrom', 'pos', 'depth'])
    depth_df['pos'] = depth_df['pos'].astype(int)
    depth_df['depth'] = depth_df['depth'].astype(int)
    primers = pd.read_csv(primer_bed, sep='\t', header=None,
                          names=['chrom', 'start', 'end', 'name', 'pool', 'strand'])
    primers['amplicon'] = primers['name'].str.extract(r'_(\d+)_').astype(float)
    dropouts = []
    for amp_id, group in primers.groupby('amplicon'):
        if pd.isna(amp_id):
            continue
        amp_start = group['start'].min()
        amp_end = group['end'].max()
        amp_depth = depth_df[(depth_df['pos'] >= amp_start) & (depth_df['pos'] <= amp_end)]['depth']
        median_depth = amp_depth.median() if len(amp_depth) else 0
        if median_depth < min_coverage:
            dropouts.append({
                'amplicon': int(amp_id),
                'median_depth': median_depth,
                'known_chronic_v41': int(amp_id) in ARTIC_V41_KNOWN_DROPOUTS,
                'scheme_version': scheme_version,
            })
    return pd.DataFrame(dropouts)


def cojac_cooccurrence_scan(bam, primer_bed, variants_yaml, out_tsv):
    '''COJAC scan for co-occurrence of signature mutations on the same read pair.'''
    subprocess.run(['cojac', 'cooc-mutbamscan', '-a', str(primer_bed),
                    '-m', str(variants_yaml), '-b', str(bam), '-o', str(out_tsv)], check=True)
    return pd.read_csv(out_tsv, sep='\t')


def lineage_frequencies_over_time(pango_results_df, date_col='collection_date',
                                  lineage_col='lineage', resample='W'):
    pango_results_df = pango_results_df.copy()
    pango_results_df[date_col] = pd.to_datetime(pango_results_df[date_col])
    counts = (pango_results_df.set_index(date_col)
                              .groupby([pd.Grouper(freq=resample), lineage_col])
                              .size()
                              .unstack(fill_value=0))
    totals = counts.sum(axis=1)
    return counts.div(totals, axis=0).fillna(0)


def multinomial_growth_advantage(frequencies_df):
    '''Simplified logistic growth (full multinomial-logistic should report covariance, not marginal CIs).'''
    import numpy as np
    weeks = np.arange(len(frequencies_df))
    growth = {}
    for lineage in frequencies_df.columns:
        freqs = frequencies_df[lineage].values
        freqs_clipped = np.clip(freqs, 1e-6, 1 - 1e-6)
        logit = np.log(freqs_clipped / (1 - freqs_clipped))
        valid = ~np.isnan(logit) & ~np.isinf(logit)
        if valid.sum() >= 3:
            slope, intercept = np.polyfit(weeks[valid], logit[valid], 1)
            growth[lineage] = slope
    note = ('Marginal slopes only; full multinomial covariance NOT reported. '
            'Early growth advantages are systematically inflated (Abousamra 2024 PLoS Comput Biol 20:e1012443).')
    return pd.Series(growth, name='weekly_logit_slope'), note
