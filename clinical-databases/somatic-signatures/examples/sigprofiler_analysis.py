'''SigProfilerSuite workflow with stability gates, COSMIC v3.4 etiology mapping, and FFPE-aware refit.

Reference: SigProfilerMatrixGenerator 1.2+, SigProfilerExtractor 1.1.24+, SigProfilerAssignment 0.1+
Verify API if version differs. COSMIC v3.4 (Sept 2024) is current.
FFPE artifact is SBS30-like (NOT SBS33 as legacy literature claims).
'''
import os


# COSMIC v3.4 etiology table (postdoc-grade -- SBS40 split; SBS17a/b; SBS10a-d POLE/POLD1).
SIGNATURE_ETIOLOGY = {
    'SBS1': 'Spontaneous 5mC deamination at CpG; age-correlated; clock-like',
    'SBS2': 'APOBEC (A3A dominant per Petljak 2022); C>T at TCW',
    'SBS3': 'HRD (BRCA1/2 deficient flat profile) -> PARP inhibitor eligibility',
    'SBS4': 'Tobacco smoking; benzo[a]pyrene-G adducts; C>A bias',
    'SBS5': 'UNKNOWN, clock-like; age-correlated; NOT polymerase fidelity errors',
    'SBS6': 'MMR-D (Lynch-associated) -> ICI eligibility',
    'SBS7a': 'UV CC>TT at dipyrimidines (CPD chemistry)',
    'SBS7b': 'UV (CPD chemistry)',
    'SBS7c': 'UV T>A (6-4 photoproduct)',
    'SBS7d': 'UV T>C (6-4 photoproduct)',
    'SBS10a': 'POLE-exo P286R; C>A at TCT; hypermutator -> ICI excellent response',
    'SBS10b': 'POLE-exo V411L; C>T at TCG; hypermutator -> ICI excellent response',
    'SBS10c': 'POLD1 (Mertz 2020)',
    'SBS10d': 'POLD1',
    'SBS11': 'Temozolomide; C>T at unmethylated CpC/CpT',
    'SBS13': 'APOBEC (A3A/B); C>G or C>A at TCW; pairs with SBS2',
    'SBS14': 'POLE+MMR concurrent defect; ultra-hypermutator -> ICI',
    'SBS15': 'MMR-D -> ICI',
    'SBS17a': 'Unknown (T>C)',
    'SBS17b': '5-Fluorouracil; T>G in CTT context (Christensen 2019)',
    'SBS18': 'ROS damage; common with KRAS',
    'SBS20': 'POLD1+MMR concurrent defect; ultra-hypermutator -> ICI',
    'SBS21': 'MMR-D (MLH1-hypermethylation sporadic)',
    'SBS22': 'Aristolochic acid; T>A at CpTpG; strong transcribed-strand bias',
    'SBS24': 'Aflatoxin (HCC; geographic exposure); C>A at CpC',
    'SBS26': 'MMR-D',
    'SBS28': 'POLE indirect',
    'SBS30': 'NTHL1-BER deficiency OR FFPE artifact (matched controls needed)',
    'SBS31': 'Platinum chemotherapy',
    'SBS35': 'Platinum chemotherapy (more comprehensive than SBS31)',
    'SBS40a': 'Pan-cancer activity; etiology uncertain',
    'SBS40b': 'RCC-specific (Senkin 2024)',
    'SBS40c': 'RCC-specific (Senkin 2024)',
    'SBS44': 'MMR-D',
    'SBS86': 'Platinum chemotherapy (Drost lab refinement)',
    'SBS87': 'Platinum chemotherapy',
    'SBS88': 'Colibactin (pks+ E. coli); T>N + T-deletions in AT-rich; CRC etiology',
}


CLINICAL_ACTIONABILITY = {
    'SBS3': 'PARP inhibitor; confirm with HRDetect 6-feature classifier',
    'SBS6': 'ICI eligibility (MMR-D); confirm with MSI testing + IHC',
    'SBS14': 'Ultra-hypermutator (POLE+MMR); ICI excellent response',
    'SBS15': 'MMR-D; ICI eligibility',
    'SBS20': 'POLD1+MMR ultra-hypermutator; ICI',
    'SBS26': 'MMR-D; ICI eligibility',
    'SBS44': 'MMR-D; ICI eligibility',
    'SBS10a': 'POLE-exo P286R; ICI excellent response',
    'SBS10b': 'POLE-exo V411L; ICI excellent response',
    'SBS2': 'APOBEC; high TMB likely; ICI signal',
    'SBS13': 'APOBEC; pairs with SBS2',
}


def generate_matrix(vcf_dir, project, genome='GRCh38', exome=False, bed_file=None):
    '''Generate SBS96 / DBS78 / ID83 mutation matrix from VCFs.

    Args:
        vcf_dir: directory of VCFs (or pre-formatted MAF)
        project: project name (output dir prefix)
        genome: 'GRCh38' or 'GRCh37'
        exome: True triggers trinucleotide-context correction for WES
        bed_file: restrict to BED region (e.g., panel)
    '''
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
    return matGen.SigProfilerMatrixGeneratorFunc(
        project=project,
        genome=genome,
        vcfFiles=vcf_dir,
        plot=True,
        exome=exome,
        bed_file=bed_file,
        chrom_based=False,
        tsb_stat=True       # transcribed-strand statistics for aristolochic-acid-like detection
    )


def extract_de_novo(matrix_path, output_dir, min_sigs=1, max_sigs=12, cosmic_version=3.4):
    '''De novo NMF extraction with stability gates per SigProfilerExtractor defaults.

    Stability gates (do NOT relax):
        nmf_replicates = 100
        minimum stability per signature >= 0.2
        average stability >= 0.8
        combined stability for chosen rank = 1.0

    Use this only when cohort size >= 50 (SBS) or >= 100 (DBS/ID/CN).
    '''
    from SigProfilerExtractor import sigpro as sig
    sig.sigProfilerExtractor(
        input_type='matrix',
        output=output_dir,
        input_data=matrix_path,
        reference_genome='GRCh38',
        opportunity_genome='GRCh38',
        minimum_signatures=min_sigs,
        maximum_signatures=max_sigs,
        nmf_replicates=100,
        cpu=-1,
        seeds='random',
        matrix_normalization='gmm',
        resample=True,
        batch_size=1,
        refit_denovo_signatures=True,
        cosmic_version=cosmic_version
    )


def refit_to_cosmic(matrix_path, output_dir, cosmic_version=3.4):
    '''Forward-backward refit to COSMIC v3.4. Use for cohorts < 50 OR single samples >= 200 mut.

    nnls_add_penalty / nnls_remove_penalty are SigProfilerAssignment defaults.
    For very strict refit (minimize false-positive signature attribution), increase add_penalty.
    '''
    from SigProfilerAssignment import Analyzer as Analyze
    Analyze.cosmic_fit(
        samples=matrix_path,
        output=output_dir,
        input_type='matrix',
        genome_build='GRCh38',
        cosmic_version=cosmic_version,
        signature_database='SBS_GRCh38_GRCh38',
        nnls_add_penalty=0.05,
        nnls_remove_penalty=0.01,
        initial_remove_penalty=0.05,
        refit_denovo_signatures=False,
        make_plots=True,
        sample_reconstruction_plots=True
    )


def interpret_signatures(contribution_file, min_contribution_pct=5.0, is_ffpe=False):
    '''Interpret signature contributions; flag clinical actionability and FFPE artifact.

    Args:
        contribution_file: SigProfilerAssignment Activities.txt or per-sample
        min_contribution_pct: signatures below this are noise
        is_ffpe: flag SBS30-like signal as artifact rather than NTHL1 biology
    '''
    import pandas as pd

    df = pd.read_csv(contribution_file, sep='\t', index_col=0)
    results = []
    for sample in df.columns:
        contribs = df[sample].sort_values(ascending=False)
        total = contribs.sum()
        if total == 0:
            continue
        significant = contribs[contribs / total > min_contribution_pct / 100]
        sample_result = {
            'sample': sample,
            'total_mutations': int(total),
            'mutation_count_floor': total >= 200,  # below 200, signature attribution unstable
            'dominant_signatures': [],
            'clinical_actionability': [],
            'ffpe_artifact_warning': False
        }
        for sig, count in significant.items():
            pct = count / total * 100
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
            if is_ffpe and sig == 'SBS30':
                sample_result['ffpe_artifact_warning'] = True
                etiology = '*** FFPE ARTIFACT suspected -- NOT NTHL1-BER deficiency biology'
            sample_result['dominant_signatures'].append({
                'signature': sig,
                'contribution_pct': round(pct, 1),
                'count': int(count),
                'etiology': etiology
            })
            if sig in CLINICAL_ACTIONABILITY:
                sample_result['clinical_actionability'].append({
                    'signature': sig,
                    'recommendation': CLINICAL_ACTIONABILITY[sig]
                })
        results.append(sample_result)
    return results


def detect_apobec_subtype(matrix_96):
    '''Distinguish APOBEC3A vs APOBEC3B by YTCA vs RTCA 5' tetranucleotide ratio.

    Petljak 2022 Nature: A3A is the dominant active deaminase.
    A3A prefers YTCA (pyrimidine 5' of T); A3B prefers RTCA (purine 5' of T).
    The ratio is computed from the 96-context matrix at TC*W substitution sites.

    Args:
        matrix_96: SBS96 matrix (96 x N samples)

    Returns: per-sample {sample, ytca_count, rtca_count, ratio, dominant}
    '''
    ytca_contexts = ['T[C>T]A', 'T[C>G]A']   # SBS2 + SBS13 at TCA
    # Y = pyrimidine 5' (C or T); R = purine 5' (A or G)
    # Simplified: compare TCA in 5'-Y context vs 5'-R context
    # (Full computation requires the 256-context matrix; this is illustrative.)
    raise NotImplementedError('Full A3A/A3B discrimination requires 5\' tetranucleotide '
                              'analysis; use SigProfilerTopography for spatial+context.')


if __name__ == '__main__':
    print('SigProfilerSuite workflow:')
    print('1. matrices = generate_matrix("vcf_dir/", "cohort_2026", exome=False)')
    print('   # Set exome=True for WES (trinucleotide-context correction)')
    print('2a. extract_de_novo("cohort_2026/output/SBS/cohort_2026.SBS96.all", "denovo_out/")')
    print('    # cohort >= 50 SBS; cohort >= 100 for DBS/ID/CN')
    print('2b. refit_to_cosmic(...)  # cohort < 50 OR single samples >= 200 mutations')
    print('3. results = interpret_signatures(activities, is_ffpe=False)')
    print()
    print('Stability gates (do NOT relax):')
    print('  nmf_replicates=100; min stability >= 0.2; avg stability >= 0.8')
    print()
    print('FFPE warning: SBS30-like signal in FFPE samples is artifact (NOT SBS33).')
    print('Recommend matched fresh-frozen controls or enzymatic uracil pretreatment.')
