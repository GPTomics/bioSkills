'''Surveillance-grade per-isolate AMR pipeline.

Runs AMRFinderPlus with species-specific point-mutation panel, harmonises with
hAMRonization, joins MOB-suite plasmid context, and emits a long-format table
suitable for WHO GLASS reporting. For Mtb, TB-Profiler is invoked instead of
AMRFinderPlus because the latter has no Mtb organism mode in v4.x.'''
# Reference: ncbi-amrfinderplus 4.0+, hamronization 1.1+, mob_suite 3.1+, tb-profiler 6.2+, pandas 2.2+ | Verify API if version differs

import json
import subprocess
from pathlib import Path
import pandas as pd

AMRF_DB_DATE = '2025-02-01.1'
RESFINDER_DB_DATE = '2024-12-15'
TBPROFILER_WHO_EDITION = '2nd'

AMRF_PARTIAL_FLAG = 'PARTIAL_CONTIG_END'

MTB_SPECIES_TOKENS = ('mycobacterium_tuberculosis', 'mycobacterium tuberculosis', 'mtb')


def species_uses_tbprofiler(species):
    return species.lower() in MTB_SPECIES_TOKENS


def run_amrfinder(assembly, species, out_tsv, threads=8):
    '''AMRFinderPlus with species mode and --plus; raises if --organism unsupported.'''
    cmd = ['amrfinder', '-n', str(assembly), '--organism', species, '--plus',
           '--threads', str(threads), '-o', str(out_tsv)]
    subprocess.run(cmd, check=True)
    return pd.read_csv(out_tsv, sep='\t')


def run_tbprofiler(reads_r1, reads_r2, prefix, out_dir):
    '''TB-Profiler against bundled WHO catalogue; emits JSON with per-drug grading.'''
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['tb-profiler', 'profile', '-1', str(reads_r1), '-2', str(reads_r2),
           '-p', prefix, '--dir', str(out_dir), '--csv', '--txt']
    subprocess.run(cmd, check=True)
    with open(out_dir / 'results' / f'{prefix}.results.json') as fh:
        return json.load(fh)


def parse_tbprofiler_drugs(tbp_json):
    '''Per-drug summary with Group 1/2/3/4/5 grading preserved; Group 3 NOT collapsed to S.'''
    rows = []
    for drug in tbp_json.get('dr_variants', []):
        rows.append({
            'drug': drug.get('drug'),
            'gene': drug.get('gene'),
            'mutation': drug.get('change'),
            'who_group': drug.get('confidence'),
            'freq': drug.get('freq'),
            'interpretation': drug.get('confidence'),
        })
    return pd.DataFrame(rows)


def hamronize_amrfinder(amrf_tsv, sample_id, out_tsv):
    cmd = ['hamronize', 'amrfinderplus',
           '--analysis_software_version', '4.0.3',
           '--reference_database_version', AMRF_DB_DATE,
           '--input_file_name', sample_id,
           str(amrf_tsv)]
    with open(out_tsv, 'w') as fh:
        subprocess.run(cmd, check=True, stdout=fh)
    return pd.read_csv(out_tsv, sep='\t')


def hamronize_resfinder(resfinder_json, sample_id, out_tsv):
    cmd = ['hamronize', 'resfinder',
           '--analysis_software_version', '4.5.0',
           '--reference_database_version', RESFINDER_DB_DATE,
           '--input_file_name', sample_id,
           str(resfinder_json)]
    with open(out_tsv, 'w') as fh:
        subprocess.run(cmd, check=True, stdout=fh)
    return pd.read_csv(out_tsv, sep='\t')


def summarize_hamronized(hamr_dir, out_tsv):
    cmd = ['hamronize', 'summarize', '-t', 'tsv', '-o', str(out_tsv)]
    cmd.extend(sorted(str(p) for p in Path(hamr_dir).glob('*.tsv')))
    subprocess.run(cmd, check=True)
    return pd.read_csv(out_tsv, sep='\t')


def run_mob_recon(assembly, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['mob_recon', '--infile', str(assembly), '--outdir', str(out_dir), '--num_threads', '4']
    subprocess.run(cmd, check=True)
    contig_report = out_dir / 'contig_report.txt'
    return pd.read_csv(contig_report, sep='\t') if contig_report.exists() else pd.DataFrame()


def join_amr_with_plasmid(amrf_df, mob_df):
    '''Cross-reference AMR gene contig against MOB-suite plasmid contigs.

    Output adds plasmid_id (or 'chromosome' / 'unassigned') per AMR row.'''
    if mob_df.empty:
        amrf_df['plasmid_id'] = 'unassigned'
        return amrf_df
    contig_to_plasmid = mob_df.set_index('contig_id')['primary_cluster_id'].to_dict()
    amrf_df['plasmid_id'] = amrf_df['Contig id'].map(contig_to_plasmid).fillna('chromosome')
    return amrf_df


def flag_partial_contig_hits(amrf_df, critical_gene_patterns=None):
    '''Return rows whose Method=PARTIAL_CONTIG_END for a clinically critical gene.

    Match by case-insensitive substring against the Gene symbol field; AMRFinderPlus
    output may use 'blaKPC-3', 'blaNDM-1', 'blaOXA-244', 'mcr-1.1', 'vanA' etc.,
    so checking for the family token rather than an exact prefix is more robust.'''
    critical_gene_patterns = critical_gene_patterns or ('KPC', 'NDM', 'OXA-48', 'OXA-181', 'OXA-232', 'OXA-244', 'mcr-', 'vanA')
    mask_partial = amrf_df['Method'].str.contains(AMRF_PARTIAL_FLAG, na=False)
    gene_upper = amrf_df['Gene symbol'].fillna('').str.upper()
    mask_critical = gene_upper.apply(lambda g: any(p.upper() in g for p in critical_gene_patterns))
    return amrf_df[mask_partial & mask_critical].copy()


def flag_oxa48_subfamily_resolution(amrf_df):
    '''Ensure OXA-48-like family alleles are reported individually, not collapsed.

    AMRFinderPlus output uses gene symbols like 'blaOXA-244' (no underscore between
    bla and OXA); the allele suffix is what determines clinical phenotype.'''
    gene_upper = amrf_df['Gene symbol'].fillna('').str.upper()
    oxa = amrf_df[gene_upper.str.contains('OXA-', na=False)].copy()
    oxa['allele_resolved'] = oxa['Gene symbol'].fillna('').str.extract(r'(OXA-\d+)', expand=False)
    return oxa


def per_isolate_pipeline(sample_id, assembly, species, hamr_dir, mob_dir):
    '''End-to-end per-isolate; chooses AMRFinderPlus vs TB-Profiler by species.'''
    if species_uses_tbprofiler(species):
        raise ValueError(f'{sample_id}: Mtb detected; route to run_tbprofiler with reads, not AMRFinderPlus')
    amrf_tsv = Path(hamr_dir) / f'{sample_id}.amrfinder.tsv'
    amrf_df = run_amrfinder(assembly, species, amrf_tsv)
    hamr_tsv = Path(hamr_dir) / f'{sample_id}.hamr.tsv'
    hamronize_amrfinder(amrf_tsv, sample_id, hamr_tsv)
    mob_df = run_mob_recon(assembly, Path(mob_dir) / sample_id)
    amrf_with_plasmid = join_amr_with_plasmid(amrf_df, mob_df)
    partial_flags = flag_partial_contig_hits(amrf_with_plasmid)
    oxa48_alleles = flag_oxa48_subfamily_resolution(amrf_with_plasmid)
    return {
        'sample_id': sample_id,
        'amr_with_plasmid': amrf_with_plasmid,
        'partial_contig_flags': partial_flags,
        'oxa48_allele_resolution': oxa48_alleles,
        'amrf_db_date': AMRF_DB_DATE,
    }


def cohort_pipeline(samples, hamr_dir, mob_dir, cohort_table):
    '''samples: list of dicts with sample_id, assembly, species.'''
    Path(hamr_dir).mkdir(parents=True, exist_ok=True)
    Path(mob_dir).mkdir(parents=True, exist_ok=True)
    per_sample_results = [per_isolate_pipeline(s['sample_id'], s['assembly'], s['species'], hamr_dir, mob_dir)
                          for s in samples if not species_uses_tbprofiler(s['species'])]
    summarize_hamronized(hamr_dir, cohort_table)
    return per_sample_results
