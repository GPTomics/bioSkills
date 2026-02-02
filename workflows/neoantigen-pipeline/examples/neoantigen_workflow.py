#!/usr/bin/env python3
'''
Neoantigen discovery pipeline: somatic VCF to ranked vaccine candidates.
Requires: pvactools, mhcflurry, pandas, numpy, matplotlib, seaborn
'''
import subprocess
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

SAMPLE_ID = 'tumor_sample'
SOMATIC_VCF = 'somatic.vep.vcf'
EXPRESSION_FILE = 'expression.tsv'  # Optional: gene-level TPM
HLA_ALLELES = 'HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*44:02,HLA-C*07:02,HLA-C*05:01'
OUTPUT_DIR = Path('neoantigen_results')
THREADS = 8

# Thresholds
IC50_STRONG = 500      # nM, strong binder threshold (IEDB standard)
IC50_VERY_STRONG = 150 # nM, very strong binder
VAF_MIN = 0.1          # Minimum VAF to ensure clonal representation
EXPRESSION_MIN = 1     # TPM, minimum detectable expression
DAI_GOOD = 500         # Differential agretopicity for tumor specificity
DAI_EXCELLENT = 1000   # High tumor specificity


def run_hla_typing(tumor_bam, output_dir):
    '''Extract HLA types from tumor RNA-seq BAM using arcasHLA.'''
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    subprocess.run([
        'arcasHLA', 'extract', tumor_bam, '-t', str(THREADS), '-o', str(output_dir)
    ], check=True)

    extracted_r1 = list(output_dir.glob('*.extracted.1.fq.gz'))[0]
    extracted_r2 = list(output_dir.glob('*.extracted.2.fq.gz'))[0]

    subprocess.run([
        'arcasHLA', 'genotype', str(extracted_r1), str(extracted_r2),
        '-g', 'A,B,C,DRB1,DQB1,DPB1', '-t', str(THREADS), '-o', str(output_dir)
    ], check=True)

    genotype_file = list(output_dir.glob('*.genotype.json'))[0]
    with open(genotype_file) as f:
        hla_data = json.load(f)

    hla_alleles = []
    for gene, alleles in hla_data.items():
        for allele in alleles:
            hla_alleles.append(f'HLA-{allele}')

    return ','.join(hla_alleles)


def annotate_vcf_with_expression(vcf_path, expression_path, sample_id, output_path):
    '''Add gene expression values to VCF using vcf-expression-annotator.'''
    subprocess.run([
        'vcf-expression-annotator', vcf_path, expression_path, 'gene',
        '-s', sample_id, '-o', output_path
    ], check=True)
    return output_path


def run_pvacseq(vcf_path, sample_id, hla_alleles, output_dir, epitope_lengths='8,9,10,11'):
    '''Run pVACseq for neoantigen prediction.'''
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        'pvacseq', 'run',
        vcf_path, sample_id, hla_alleles,
        'MHCflurry', 'MHCnuggetsI', 'NetMHCpan',
        str(output_dir),
        '-e1', epitope_lengths,
        '-t', str(THREADS),
        '--tumor-purity', '0.7',
        '--trna-vaf', str(VAF_MIN)
    ]

    subprocess.run(cmd, check=True)
    return output_dir


def calculate_immunogenicity_score(row):
    '''Multi-factor immunogenicity scoring for neoantigen ranking.'''
    score = 0

    # Binding affinity component
    if row['Median MT IC50 Score'] < IC50_VERY_STRONG:
        score += 3  # Very strong binder
    elif row['Median MT IC50 Score'] < IC50_STRONG:
        score += 2  # Strong binder

    # Tumor specificity (DAI)
    dai = row.get('DAI', 0)
    if dai > DAI_EXCELLENT:
        score += 2
    elif dai > DAI_GOOD:
        score += 1

    # Clonality (VAF)
    vaf = row.get('Tumor DNA VAF', 0)
    if vaf > 0.3:
        score += 2  # Likely clonal
    elif vaf > 0.15:
        score += 1  # Subclonal but substantial

    # Expression (if available)
    expr = row.get('Gene Expression', 0)
    if expr > 10:
        score += 1

    return score


def filter_and_rank_neoantigens(pvacseq_output_dir):
    '''Filter pVACseq results and rank by immunogenicity.'''
    results_path = pvacseq_output_dir / 'MHC_Class_I' / f'{SAMPLE_ID}.filtered.tsv'
    if not results_path.exists():
        results_path = list((pvacseq_output_dir / 'MHC_Class_I').glob('*.filtered.tsv'))[0]

    results = pd.read_csv(results_path, sep='\t')
    print(f'Loaded {len(results)} filtered epitopes')

    # Strong binders
    strong = results[results['Median MT IC50 Score'] < IC50_STRONG].copy()
    print(f'Strong binders (IC50 <{IC50_STRONG}nM): {len(strong)}')

    # Calculate DAI (differential agretopicity index)
    if 'Median WT IC50 Score' in strong.columns:
        strong['DAI'] = strong['Median WT IC50 Score'] - strong['Median MT IC50 Score']
    else:
        strong['DAI'] = 0

    # VAF filter
    if 'Tumor DNA VAF' in strong.columns:
        strong = strong[strong['Tumor DNA VAF'] > VAF_MIN]
        print(f'After VAF filter (>{VAF_MIN}): {len(strong)}')

    # Expression filter
    if 'Gene Expression' in strong.columns:
        strong = strong[strong['Gene Expression'] > EXPRESSION_MIN]
        print(f'After expression filter (>{EXPRESSION_MIN} TPM): {len(strong)}')

    # Calculate immunogenicity scores
    strong['Immunogenicity Score'] = strong.apply(calculate_immunogenicity_score, axis=1)

    # Rank
    ranked = strong.sort_values('Immunogenicity Score', ascending=False)

    return ranked


def visualize_results(ranked, output_dir):
    '''Generate summary visualizations.'''
    output_dir = Path(output_dir)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # IC50 distribution
    ax1 = axes[0]
    ax1.hist(ranked['Median MT IC50 Score'], bins=30, edgecolor='black', alpha=0.7)
    ax1.axvline(IC50_STRONG, color='red', linestyle='--', label=f'{IC50_STRONG}nM threshold')
    ax1.axvline(IC50_VERY_STRONG, color='orange', linestyle='--', label=f'{IC50_VERY_STRONG}nM (very strong)')
    ax1.set_xlabel('Median MT IC50 (nM)')
    ax1.set_ylabel('Count')
    ax1.set_title('Binding Affinity Distribution')
    ax1.legend()

    # DAI vs IC50 scatter
    ax2 = axes[1]
    if 'DAI' in ranked.columns and ranked['DAI'].sum() > 0:
        scatter = ax2.scatter(ranked['Median MT IC50 Score'], ranked['DAI'],
                             c=ranked['Immunogenicity Score'], cmap='viridis', alpha=0.7, s=50)
        ax2.set_xlabel('MT IC50 (nM)')
        ax2.set_ylabel('Differential Agretopicity Index')
        ax2.set_title('Tumor Specificity vs Binding')
        plt.colorbar(scatter, ax=ax2, label='Immunogenicity Score')
    else:
        ax2.text(0.5, 0.5, 'DAI not available', ha='center', va='center', transform=ax2.transAxes)

    # Top genes
    ax3 = axes[2]
    gene_counts = ranked['Gene Name'].value_counts().head(15)
    gene_counts.plot(kind='barh', ax=ax3, color='steelblue')
    ax3.set_xlabel('Number of Neoantigens')
    ax3.set_title('Top Genes with Neoantigens')
    ax3.invert_yaxis()

    plt.tight_layout()
    plt.savefig(output_dir / 'neoantigen_summary.pdf', dpi=150)
    plt.savefig(output_dir / 'neoantigen_summary.png', dpi=150)
    print(f'Saved visualizations to {output_dir}')


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print('=== Neoantigen Discovery Pipeline ===')
    print(f'Sample: {SAMPLE_ID}')
    print(f'HLA alleles: {HLA_ALLELES}')

    # Step 1: Run pVACseq
    print('\n=== Running pVACseq ===')
    pvacseq_dir = OUTPUT_DIR / 'pvacseq'
    run_pvacseq(SOMATIC_VCF, SAMPLE_ID, HLA_ALLELES, pvacseq_dir)

    # Step 2: Filter and rank
    print('\n=== Filtering and Ranking ===')
    ranked = filter_and_rank_neoantigens(pvacseq_dir)

    # Step 3: Export top candidates
    top_n = 20
    top_candidates = ranked.head(top_n)
    output_file = OUTPUT_DIR / 'top_neoantigen_candidates.tsv'
    top_candidates.to_csv(output_file, sep='\t', index=False)
    print(f'\nTop {top_n} candidates saved to {output_file}')

    # Summary table
    display_cols = ['Gene Name', 'MT Epitope Seq', 'HLA Allele', 'Median MT IC50 Score', 'DAI', 'Immunogenicity Score']
    display_cols = [c for c in display_cols if c in ranked.columns]
    print('\nTop 10 Neoantigen Candidates:')
    print(ranked[display_cols].head(10).to_string(index=False))

    # Step 4: Visualize
    print('\n=== Generating Visualizations ===')
    visualize_results(ranked, OUTPUT_DIR)

    print('\n=== Pipeline Complete ===')
    print(f'Total strong binders: {len(ranked)}')
    print(f'Results: {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
