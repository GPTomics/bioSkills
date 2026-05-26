'''Pathogen typing pipeline: 7-locus MLST, cgMLST allele-distance with pairwise-complete
handling, lineage assignment per organism, and pathogen-specific outbreak cluster definition.

Routes Mtb to TB-Profiler (Coll/Napier barcode); Salmonella to SISTR;
Klebsiella to Kleborate; SARS-CoV-2 to Pangolin UShER + Nextclade with version pinning.'''
# Reference: mlst 2.23+, chewBBACA 3.3+, sistr_cmd 1.1+, kleborate 3.0+, pangolin 4.3+, nextclade 3.8+, tb-profiler 6.2+, snp-dists 0.8+, pandas 2.2+ | Verify API if version differs

import json
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

CHEWBBACA_MISSING_CODES = {'LNF', 'PLOT', 'NIPH', 'NIPHEM', 'ASM', 'ALM', '0', '-', ''}
CGMLST_COMPLETENESS = 0.95

SNP_THRESHOLDS = {
    'mtb_recent': 5,
    'mtb_likely': 12,
    'mrsa_hospital': 15,
    'kpneumo_kpc': 21,
    'cdiff_direct': 2,
    'cdiff_plausible': 10,
    'ngono_transmission': 25,
}

CGMLST_THRESHOLDS = {
    'salmonella_efsa': 5,
    'listeria_pulsenet': 4,
    'ecoli_enterobase': 10,
}


def run_mlst_cohort(assemblies, out_tsv):
    cmd = ['mlst', '--threads', '8'] + [str(a) for a in assemblies]
    with open(out_tsv, 'w') as fh:
        subprocess.run(cmd, check=True, stdout=fh)
    rows = []
    for line in open(out_tsv):
        parts = line.rstrip('\n').split('\t')
        rows.append({'file': parts[0], 'scheme': parts[1], 'ST': parts[2], 'alleles': parts[3:]})
    return pd.DataFrame(rows)


def run_chewbbaca_allelecall(assemblies_dir, schema_dir, out_dir, cpu=8):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['chewBBACA.py', 'AlleleCall', '-i', str(assemblies_dir), '-g', str(schema_dir),
           '-o', str(out_dir), '--cpu', str(cpu)]
    subprocess.run(cmd, check=True)
    profile_path = out_dir / 'results_alleles.tsv'
    return pd.read_csv(profile_path, sep='\t', index_col=0)


def extract_cgmlst(profile_path, out_path, completeness=CGMLST_COMPLETENESS):
    cmd = ['chewBBACA.py', 'ExtractCgMLST', '-i', str(profile_path),
           '-o', str(out_path), '--threshold', str(completeness)]
    subprocess.run(cmd, check=True)
    return pd.read_csv(Path(out_path) / 'cgMLST.tsv', sep='\t', index_col=0)


def pairwise_complete_allele_distance(profile_df):
    '''cgMLST distance on the intersection of called loci; missing codes excluded.

    Returns: (distance_matrix, n_loci_compared_matrix).'''
    sample_ids = profile_df.index.tolist()
    n = len(sample_ids)
    dist = np.zeros((n, n), dtype=int)
    n_compared = np.zeros((n, n), dtype=int)
    arr = profile_df.astype(str).values
    is_called = ~np.isin(arr, list(CHEWBBACA_MISSING_CODES))
    for i in range(n):
        for j in range(i + 1, n):
            both_called = is_called[i] & is_called[j]
            n_compared[i, j] = n_compared[j, i] = both_called.sum()
            if n_compared[i, j] == 0:
                continue
            diffs = (arr[i] != arr[j]) & both_called
            dist[i, j] = dist[j, i] = diffs.sum()
    return (pd.DataFrame(dist, index=sample_ids, columns=sample_ids),
            pd.DataFrame(n_compared, index=sample_ids, columns=sample_ids))


def define_clusters(distance_df, threshold):
    '''Single-linkage clusters at the pathogen-specific allele/SNP threshold.'''
    from scipy.cluster.hierarchy import linkage, fcluster
    condensed = distance_df.values[np.triu_indices(len(distance_df), k=1)]
    if condensed.size == 0:
        return pd.Series(1, index=distance_df.index, name='cluster')
    Z = linkage(condensed, method='single')
    clusters = fcluster(Z, t=threshold, criterion='distance')
    return pd.Series(clusters, index=distance_df.index, name='cluster')


def run_sistr(assemblies, out_tsv):
    '''SISTR Salmonella serovar prediction; flag monophasic Typhimurium 1,4,[5],12:i:- variants.

    Invocation form varies across SISTR releases; verify `sistr --help` for the installed
    version. The output schema (column names like 'serovar') is stable across recent
    releases but verify before relying on the monophasic_flag column.'''
    out_tsv = Path(out_tsv)
    cmd = ['sistr', '-f', 'tab', '-o', str(out_tsv)] + [str(a) for a in assemblies]
    subprocess.run(cmd, check=True)
    df = pd.read_csv(out_tsv, sep='\t')
    df['monophasic_flag'] = df['serovar'].fillna('').str.contains(r'1,4,\[?5\]?,12:i:-', na=False, regex=True)
    return df


def run_kleborate(assemblies, out_dir):
    '''Kleborate v3 surveillance for Klebsiella; integrates MLST + Kaptive K/O + virulence + AMR.

    Kleborate v3 changed the CLI from earlier releases; verify `kleborate --help` against
    the installed version. Output directory contains per-isolate result tables.'''
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['kleborate', '-a'] + [str(a) for a in assemblies] + ['-o', str(out_dir)]
    subprocess.run(cmd, check=True)
    summary = out_dir / 'klebsiella_pneumo_complex_output.txt'
    return pd.read_csv(summary, sep='\t') if summary.exists() else pd.DataFrame()


def run_tbprofiler_lineage(reads_r1, reads_r2, prefix, out_dir):
    cmd = ['tb-profiler', 'profile', '-1', str(reads_r1), '-2', str(reads_r2),
           '-p', prefix, '--dir', str(out_dir)]
    subprocess.run(cmd, check=True)
    with open(Path(out_dir) / 'results' / f'{prefix}.results.json') as fh:
        result = json.load(fh)
    return {
        'sample_id': prefix,
        'main_lineage': result.get('main_lin'),
        'sub_lineage': result.get('sublin'),
        'spoligotype': result.get('spoligotype'),
    }


def run_pangolin_with_versions(seq_fasta, out_csv, version_log):
    cmd_pango = ['pangolin', str(seq_fasta), '--analysis-mode', 'usher', '--outfile', str(out_csv)]
    subprocess.run(cmd_pango, check=True)
    with open(version_log, 'w') as fh:
        subprocess.run(['pangolin', '--all-versions'], check=True, stdout=fh)
    return pd.read_csv(out_csv)


def run_nextclade_with_dataset(seq_fasta, dataset_dir, out_tsv, metadata_log):
    cmd = ['nextclade', 'run', '--input-dataset', str(dataset_dir),
           '--output-tsv', str(out_tsv), str(seq_fasta)]
    subprocess.run(cmd, check=True)
    with open(Path(dataset_dir) / 'pathogen.json') as fh:
        dataset = json.load(fh)
    with open(metadata_log, 'w') as fh:
        fh.write(f"nextclade_dataset_tag: {dataset.get('tag', dataset.get('version'))}\n")
    return pd.read_csv(out_tsv, sep='\t')


def reconcile_pangolin_nextclade(pango_df, nextclade_df, key_pango='taxon', key_nc='seqName'):
    merged = pango_df.merge(nextclade_df, left_on=key_pango, right_on=key_nc, how='outer')
    merged['recombinant_candidate'] = merged['lineage'].str.startswith('X', na=False)
    merged['agreement'] = merged.apply(
        lambda r: 'concordant' if isinstance(r.get('lineage'), str) and isinstance(r.get('clade'), str)
                 else 'discordant', axis=1)
    return merged


def snp_clusters_from_snippy(snippy_core_aln, gubbins_prefix='gubbins', threshold_key='mtb_likely'):
    subprocess.run(['run_gubbins.py', '--prefix', gubbins_prefix, str(snippy_core_aln)], check=True)
    snp_csv = Path(f'{gubbins_prefix}.snp_dists.csv')
    with open(snp_csv, 'w') as fh:
        subprocess.run(['snp-dists', '-c', f'{gubbins_prefix}.filtered_polymorphic_sites.fasta'],
                       check=True, stdout=fh)
    dist = pd.read_csv(snp_csv, index_col=0)
    return define_clusters(dist, SNP_THRESHOLDS[threshold_key])
