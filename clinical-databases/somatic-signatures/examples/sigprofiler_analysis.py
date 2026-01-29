'''Extract and analyze mutational signatures with SigProfiler'''

import os

# --- Install reference genome (first-time only) ---
# from SigProfilerMatrixGenerator import install as genInstall
# genInstall.install('GRCh38')

# --- ALTERNATIVE: Use public data ---
# Download TCGA VCFs from GDC portal:
# https://portal.gdc.cancer.gov/
# Or use PCAWG signatures: https://dcc.icgc.org/releases/PCAWG


def generate_mutation_matrix(vcf_dir, project_name, genome='GRCh38', exome=False):
    '''Generate 96-context mutation matrix from VCF files

    Args:
        vcf_dir: Directory containing VCF files
        project_name: Project name for output
        genome: Reference genome (GRCh38 or GRCh37)
        exome: Set True for WES data (adjusts normalization)

    Output:
        Matrix file at: {project_name}/output/SBS/{project_name}.SBS96.all
    '''
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        project=project_name,
        genome=genome,
        vcfFiles=vcf_dir,
        plot=True,
        exome=exome
    )
    return matrices


def extract_signatures(matrix_path, output_dir, min_sigs=1, max_sigs=10):
    '''De novo signature extraction using NMF

    Args:
        matrix_path: Path to SBS96 matrix
        output_dir: Output directory
        min_sigs: Minimum signatures to try
        max_sigs: Maximum signatures to try

    Determines optimal number of signatures automatically
    '''
    from SigProfilerExtractor import sigpro as sig

    sig.sigProfilerExtractor(
        input_type='matrix',
        output=output_dir,
        input_data=matrix_path,
        reference_genome='GRCh38',
        minimum_signatures=min_sigs,
        maximum_signatures=max_sigs,
        nmf_replicates=100,
        cpu=-1  # Use all available cores
    )


def fit_cosmic_signatures(matrix_path, output_dir):
    '''Fit samples to known COSMIC signatures

    More appropriate when sample size is small (<20 samples)
    or when you want to identify known processes
    '''
    from SigProfilerAssignment import Analyzer as Analyze

    Analyze.cosmic_fit(
        samples=matrix_path,
        output=output_dir,
        input_type='matrix',
        genome_build='GRCh38',
        signature_database='SBS_GRCh38_GRCh38'
    )


# --- Signature interpretation ---

SIGNATURE_ETIOLOGY = {
    'SBS1': 'Spontaneous deamination of 5-methylcytosine (age-related, clock-like)',
    'SBS2': 'APOBEC cytidine deaminase activity',
    'SBS3': 'Defective homologous recombination (BRCA1/BRCA2)',
    'SBS4': 'Tobacco smoking exposure',
    'SBS5': 'Unknown etiology (age-related, clock-like)',
    'SBS6': 'Defective DNA mismatch repair',
    'SBS7a': 'Ultraviolet light exposure (CC>TT)',
    'SBS7b': 'Ultraviolet light exposure',
    'SBS10a': 'Defective POLE proofreading',
    'SBS10b': 'Defective POLE proofreading',
    'SBS13': 'APOBEC cytidine deaminase activity',
    'SBS15': 'Defective DNA mismatch repair',
    'SBS17a': 'Unknown (common in gastric/esophageal)',
    'SBS17b': 'Unknown (common in gastric/esophageal)',
    'SBS18': 'Damage from reactive oxygen species',
    'SBS22': 'Aristolochic acid exposure',
    'SBS26': 'Defective DNA mismatch repair',
    'SBS31': 'Platinum chemotherapy',
    'SBS35': 'Platinum chemotherapy',
    'SBS44': 'Defective DNA mismatch repair',
}

CLINICAL_IMPLICATIONS = {
    'SBS3': {
        'implication': 'HR deficiency - may benefit from PARP inhibitors (olaparib, etc.)',
        'action': 'Consider BRCA1/2 germline and somatic testing'
    },
    'SBS6': {
        'implication': 'MMR deficiency - may respond to immune checkpoint inhibitors',
        'action': 'Consider MSI testing and pembrolizumab eligibility'
    },
    'SBS15': {
        'implication': 'MMR deficiency - may respond to immune checkpoint inhibitors',
        'action': 'Consider MSI testing'
    },
    'SBS26': {
        'implication': 'MMR deficiency - may respond to immune checkpoint inhibitors',
        'action': 'Consider MSI testing'
    },
    'SBS44': {
        'implication': 'MMR deficiency - may respond to immune checkpoint inhibitors',
        'action': 'Consider MSI testing'
    },
}


def interpret_results(contribution_file):
    '''Interpret signature contributions

    Args:
        contribution_file: Path to COSMIC fit activities file

    Returns:
        List of interpretations per sample
    '''
    import pandas as pd

    df = pd.read_csv(contribution_file, sep='\t', index_col=0)

    results = []
    for sample in df.columns:
        # Get signature contributions
        contribs = df[sample].sort_values(ascending=False)

        # Filter to significant contributions (>5%)
        total = contribs.sum()
        significant = contribs[contribs / total > 0.05]

        sample_result = {
            'sample': sample,
            'dominant_signatures': [],
            'clinical_implications': []
        }

        for sig, count in significant.items():
            pct = count / total * 100
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')

            sample_result['dominant_signatures'].append({
                'signature': sig,
                'contribution_pct': round(pct, 1),
                'etiology': etiology
            })

            if sig in CLINICAL_IMPLICATIONS:
                sample_result['clinical_implications'].append(
                    CLINICAL_IMPLICATIONS[sig]
                )

        results.append(sample_result)

    return results


if __name__ == '__main__':
    print('SigProfiler Mutational Signature Analysis')
    print('')
    print('Workflow:')
    print('1. Generate matrix: generate_mutation_matrix(vcf_dir, project_name)')
    print('2a. De novo extraction: extract_signatures(matrix_path, output_dir)')
    print('2b. COSMIC fit: fit_cosmic_signatures(matrix_path, output_dir)')
    print('3. Interpret: interpret_results(contribution_file)')
    print('')
    print('Key signatures to check:')
    for sig, etiology in list(SIGNATURE_ETIOLOGY.items())[:10]:
        print(f'  {sig}: {etiology}')
