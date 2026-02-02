'''AMR gene detection and surveillance analysis'''

import subprocess
import pandas as pd
from collections import Counter


def run_amrfinder(fasta_file, output_file, organism=None, protein=False):
    '''Run AMRFinderPlus on genome/proteins

    AMRFinderPlus detects:
    - Acquired AMR genes (plasmid-borne)
    - Point mutations (organism-specific)
    - Stress response genes
    - Virulence factors

    Args:
        organism: Enable point mutation detection for specific organisms
                 Supported: Salmonella, Escherichia, Klebsiella, etc.
        protein: If True, input is protein FASTA (faster, more sensitive)
    '''
    cmd = ['amrfinder']

    if protein:
        cmd.extend(['-p', fasta_file])
    else:
        cmd.extend(['-n', fasta_file])

    if organism:
        cmd.extend(['--organism', organism])

    cmd.extend(['-o', output_file])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f'AMRFinder error: {result.stderr}')
        return None

    return output_file


def parse_amrfinder_results(results_file, min_coverage=90, min_identity=90):
    '''Parse and filter AMRFinderPlus output

    Filtering thresholds:
    - Coverage ≥90%: Gene is mostly present
    - Identity ≥90%: Strong match to known gene

    Lower thresholds may capture variants but increase false positives
    '''
    df = pd.read_csv(results_file, sep='\t')

    # Filter by quality
    if '% Coverage of reference sequence' in df.columns:
        df = df[df['% Coverage of reference sequence'] >= min_coverage]
    if '% Identity to reference sequence' in df.columns:
        df = df[df['% Identity to reference sequence'] >= min_identity]

    return df


def summarize_amr_profile(results_df):
    '''Summarize AMR genes by drug class'''
    amr = results_df[results_df['Element type'] == 'AMR']

    if len(amr) == 0:
        return {'status': 'No AMR genes detected'}

    summary = {
        'total_genes': len(amr),
        'drug_classes': amr['Class'].nunique(),
        'genes_by_class': amr.groupby('Class')['Gene symbol'].apply(list).to_dict()
    }

    return summary


def analyze_surveillance_data(samples_df):
    '''Analyze AMR trends from multiple samples

    Surveillance metrics:
    - Prevalence by drug class
    - Temporal trends
    - Co-occurrence patterns
    '''
    # Prevalence = proportion of samples with gene
    n_samples = samples_df['sample_id'].nunique()

    gene_prevalence = samples_df.groupby('Gene symbol')['sample_id'].nunique() / n_samples * 100
    class_prevalence = samples_df.groupby('Class')['sample_id'].nunique() / n_samples * 100

    # Critical resistance flags
    critical_classes = ['Carbapenem', 'Colistin', 'Vancomycin']
    critical_detected = []
    for cls in critical_classes:
        if any(cls.lower() in c.lower() for c in samples_df['Class'].unique()):
            critical_detected.append(cls)

    return {
        'n_samples': n_samples,
        'n_genes': samples_df['Gene symbol'].nunique(),
        'gene_prevalence': gene_prevalence.to_dict(),
        'class_prevalence': class_prevalence.to_dict(),
        'critical_resistance': critical_detected
    }


def clinical_interpretation(genes):
    '''Generate clinical interpretation of detected genes

    Key resistance patterns requiring attention:
    - ESBL (bla_CTX-M, bla_SHV): Avoid 3rd gen cephalosporins
    - Carbapenemase (bla_KPC, bla_NDM, bla_OXA-48): Limited options
    - MCR (mcr-1): Colistin resistance - last resort compromised
    - VanA/B: Vancomycin resistance in Enterococcus
    '''
    interpretations = []

    patterns = {
        'CTX-M': ('ESBL', 'Avoid cephalosporins, consider carbapenems'),
        'KPC': ('Carbapenemase', 'CRITICAL - limited options, consider ceftazidime-avibactam'),
        'NDM': ('Carbapenemase', 'CRITICAL - consider aztreonam combinations'),
        'OXA-48': ('Carbapenemase', 'CRITICAL - consider ceftazidime-avibactam'),
        'mcr': ('Colistin resistance', 'Last-resort antibiotic compromised'),
        'vanA': ('Vancomycin resistance', 'VRE - consider linezolid, daptomycin'),
        'vanB': ('Vancomycin resistance', 'VRE - teicoplanin may be active'),
    }

    for gene in genes:
        for pattern, (category, recommendation) in patterns.items():
            if pattern.lower() in gene.lower():
                interpretations.append({
                    'gene': gene,
                    'category': category,
                    'recommendation': recommendation
                })
                break

    return interpretations


if __name__ == '__main__':
    print('AMR Surveillance Analysis Example')
    print('=' * 50)

    # Simulated AMRFinder results
    demo_data = pd.DataFrame({
        'sample_id': ['S1', 'S1', 'S2', 'S2', 'S3', 'S4', 'S4'],
        'Gene symbol': ['bla_CTX-M-15', 'aac(3)-IIa', 'bla_CTX-M-15', 'bla_KPC-2',
                       'bla_TEM-1', 'bla_NDM-1', 'mcr-1'],
        'Class': ['Beta-lactam', 'Aminoglycoside', 'Beta-lactam', 'Carbapenem',
                 'Beta-lactam', 'Carbapenem', 'Colistin'],
        'Element type': ['AMR'] * 7
    })

    print('\nSample AMR data:')
    print(demo_data.to_string(index=False))

    # Surveillance analysis
    print('\n' + '=' * 50)
    print('Surveillance Analysis:')
    stats = analyze_surveillance_data(demo_data)
    print(f"  Samples analyzed: {stats['n_samples']}")
    print(f"  Unique AMR genes: {stats['n_genes']}")
    print(f"\n  Prevalence by class:")
    for cls, prev in sorted(stats['class_prevalence'].items(), key=lambda x: -x[1]):
        print(f"    {cls}: {prev:.0f}%")

    if stats['critical_resistance']:
        print(f"\n  CRITICAL: {', '.join(stats['critical_resistance'])} resistance detected!")

    # Clinical interpretation
    print('\n' + '=' * 50)
    print('Clinical Interpretation:')
    genes = demo_data['Gene symbol'].unique()
    interp = clinical_interpretation(genes)
    for i in interp:
        print(f"\n  {i['gene']} ({i['category']}):")
        print(f"    {i['recommendation']}")
