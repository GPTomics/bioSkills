'''Ortholog inference with OrthoFinder'''

import subprocess
import pandas as pd
import os


def run_orthofinder(proteome_dir, threads=4, method='diamond'):
    '''Run OrthoFinder on proteome directory

    Input: Directory with one .fa/.fasta file per species
    Filename (without extension) becomes species name

    Search methods:
    - diamond: Fast, good for most cases (default)
    - blast: More sensitive, slower

    Output structure:
    - Orthogroups/Orthogroups.tsv: Main results
    - Single_Copy_Orthologue_Sequences/: For phylogenomics
    - Gene_Trees/: Individual gene trees
    - Species_Tree/: Inferred species tree
    '''
    cmd = f'orthofinder -f {proteome_dir} -t {threads}'

    if method == 'blast':
        cmd += ' -S blast'
    else:
        cmd += ' -S diamond'

    print(f'Running: {cmd}')
    subprocess.run(cmd, shell=True)

    # Find results directory
    for d in os.listdir(proteome_dir):
        if d.startswith('OrthoFinder'):
            results_path = os.path.join(proteome_dir, d)
            for subdir in os.listdir(results_path):
                if subdir.startswith('Results_'):
                    return os.path.join(results_path, subdir)

    return None


def load_orthogroups(orthogroups_tsv):
    '''Load and parse Orthogroups.tsv'''
    df = pd.read_csv(orthogroups_tsv, sep='\t', index_col=0)
    return df


def count_genes_per_orthogroup(df):
    '''Count genes per species in each orthogroup'''
    counts = pd.DataFrame(index=df.index)

    for col in df.columns:
        counts[col] = df[col].apply(
            lambda x: len(x.split(', ')) if pd.notna(x) and x else 0
        )

    return counts


def classify_orthogroups(df):
    '''Classify orthogroups by presence pattern

    Categories:
    - single_copy: 1 gene per species in all species
    - core: Present in all species (any copy number)
    - dispensable: Missing from some species
    - species_specific: Only in one species
    '''
    n_species = len(df.columns)
    counts = count_genes_per_orthogroup(df)

    results = {'single_copy': [], 'core': [], 'dispensable': [], 'species_specific': []}

    for og in df.index:
        row = counts.loc[og]
        present = (row > 0).sum()
        all_single = all(row == 1)

        if present == 1:
            results['species_specific'].append(og)
        elif present == n_species:
            if all_single:
                results['single_copy'].append(og)
            else:
                results['core'].append(og)
        else:
            results['dispensable'].append(og)

    return results


def extract_single_copy_sequences(results_dir, output_dir):
    '''Copy single-copy ortholog sequences for phylogenomics

    These are in: Single_Copy_Orthologue_Sequences/

    Each file contains aligned sequences for one orthogroup
    with exactly one gene per species
    '''
    sco_dir = os.path.join(results_dir, 'Single_Copy_Orthologue_Sequences')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if os.path.exists(sco_dir):
        for f in os.listdir(sco_dir):
            src = os.path.join(sco_dir, f)
            dst = os.path.join(output_dir, f)
            subprocess.run(f'cp {src} {dst}', shell=True)

        return len(os.listdir(output_dir))

    return 0


def find_gene_orthogroup(gene_id, df):
    '''Find which orthogroup contains a specific gene'''
    for og in df.index:
        for col in df.columns:
            cell = df.loc[og, col]
            if pd.notna(cell) and gene_id in cell:
                return og
    return None


def get_orthologs(gene_id, df):
    '''Get all orthologs of a gene across species'''
    og = find_gene_orthogroup(gene_id, df)
    if og is None:
        return None

    orthologs = {}
    for col in df.columns:
        cell = df.loc[og, col]
        if pd.notna(cell) and cell:
            orthologs[col] = cell.split(', ')

    return {'orthogroup': og, 'orthologs': orthologs}


def summarize_orthogroups(df, classification):
    '''Print summary statistics'''
    print('Orthogroup Analysis Summary')
    print('=' * 50)

    n_species = len(df.columns)
    n_orthogroups = len(df)

    print(f'\nSpecies analyzed: {n_species}')
    print(f'  {", ".join(df.columns)}')

    print(f'\nTotal orthogroups: {n_orthogroups}')
    print(f'  Single-copy: {len(classification["single_copy"])} '
          f'({100*len(classification["single_copy"])/n_orthogroups:.1f}%)')
    print(f'  Core (all species): {len(classification["core"])} '
          f'({100*len(classification["core"])/n_orthogroups:.1f}%)')
    print(f'  Dispensable: {len(classification["dispensable"])} '
          f'({100*len(classification["dispensable"])/n_orthogroups:.1f}%)')
    print(f'  Species-specific: {len(classification["species_specific"])} '
          f'({100*len(classification["species_specific"])/n_orthogroups:.1f}%)')

    # Count total genes
    counts = count_genes_per_orthogroup(df)
    for col in df.columns:
        total_genes = counts[col].sum()
        in_sco = len([og for og in classification['single_copy']
                      if counts.loc[og, col] == 1])
        print(f'\n  {col}: {total_genes} genes, {in_sco} in single-copy orthologs')


if __name__ == '__main__':
    print('OrthoFinder Ortholog Analysis')
    print('=' * 50)

    # Example data (simulated orthogroups)
    example_data = {
        'Orthogroup': ['OG0000001', 'OG0000002', 'OG0000003', 'OG0000004', 'OG0000005'],
        'Human': ['BRCA1', 'TP53', 'MYC', 'KRAS, HRAS, NRAS', ''],
        'Mouse': ['Brca1', 'Trp53', 'Myc', 'Kras, Hras, Nras', 'Zfp42'],
        'Zebrafish': ['brca1', 'tp53', 'myca, mycb', 'krasa, krasb', '']
    }

    df = pd.DataFrame(example_data).set_index('Orthogroup')

    print('\nExample orthogroups:')
    print(df.to_string())

    classification = classify_orthogroups(df)
    summarize_orthogroups(df, classification)

    # Find orthologs of a specific gene
    print('\n\nOrthologs of BRCA1:')
    result = get_orthologs('BRCA1', df)
    if result:
        print(f"  Orthogroup: {result['orthogroup']}")
        for species, genes in result['orthologs'].items():
            print(f"  {species}: {', '.join(genes)}")

    print('\n\nTo run on real data:')
    print('1. Place proteome FASTA files in a directory (one per species)')
    print('2. results_dir = run_orthofinder("proteomes/")')
    print('3. df = load_orthogroups(f"{results_dir}/Orthogroups/Orthogroups.tsv")')
    print('4. classification = classify_orthogroups(df)')
