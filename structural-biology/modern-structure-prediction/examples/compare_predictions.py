from Bio.PDB import PDBParser, Superimposer
import numpy as np

def get_ca_atoms(structure):
    '''Extract CA atoms from first chain of structure'''
    atoms = []
    for chain in structure[0]:
        for residue in chain:
            if 'CA' in residue:
                atoms.append(residue['CA'])
        break  # First chain only
    return atoms

def calculate_rmsd(pdb1, pdb2):
    '''Calculate CA RMSD between two structures'''
    parser = PDBParser(QUIET=True)
    struct1 = parser.get_structure('model1', pdb1)
    struct2 = parser.get_structure('model2', pdb2)

    atoms1 = get_ca_atoms(struct1)
    atoms2 = get_ca_atoms(struct2)

    # Align by minimum length
    min_len = min(len(atoms1), len(atoms2))
    if min_len == 0:
        raise ValueError('No CA atoms found')

    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1[:min_len], atoms2[:min_len])
    return super_imposer.rms

def compare_multiple_predictions(pdb_files, labels=None):
    '''Compare multiple structure predictions with pairwise RMSD'''
    if labels is None:
        labels = [f'Model_{i}' for i in range(len(pdb_files))]

    n = len(pdb_files)
    rmsd_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            try:
                rmsd = calculate_rmsd(pdb_files[i], pdb_files[j])
                rmsd_matrix[i,j] = rmsd_matrix[j,i] = rmsd
            except Exception as e:
                print(f'Warning: Could not compare {labels[i]} vs {labels[j]}: {e}')
                rmsd_matrix[i,j] = rmsd_matrix[j,i] = np.nan

    # Print formatted matrix
    print('\nPairwise RMSD (Angstroms):')
    header = '         ' + '  '.join(f'{l:>8}' for l in labels)
    print(header)
    for i, label in enumerate(labels):
        row = f'{label:>8} ' + '  '.join(f'{rmsd_matrix[i,j]:8.2f}' for j in range(n))
        print(row)

    return rmsd_matrix

def extract_plddt_comparison(pdb_files, labels=None):
    '''Extract and compare pLDDT across predictions'''
    from collections import defaultdict

    if labels is None:
        labels = [f'Model_{i}' for i in range(len(pdb_files))]

    parser = PDBParser(QUIET=True)
    plddt_data = {}

    for pdb_file, label in zip(pdb_files, labels):
        struct = parser.get_structure(label, pdb_file)
        plddt = {}
        for chain in struct[0]:
            for residue in chain:
                if 'CA' in residue:
                    plddt[residue.id[1]] = residue['CA'].get_bfactor()
            break
        plddt_data[label] = plddt
        avg = sum(plddt.values()) / len(plddt) if plddt else 0
        print(f'{label}: Average pLDDT = {avg:.1f}')

    return plddt_data

if __name__ == '__main__':
    # Example: Compare predictions from different methods
    # Replace with actual prediction files
    pdb_files = [
        'esmfold_prediction.pdb',
        'alphafold_prediction.pdb',
        'chai1_prediction.pdb'
    ]
    labels = ['ESMFold', 'AlphaFold3', 'Chai-1']

    # Check which files exist
    from pathlib import Path
    existing = [(f, l) for f, l in zip(pdb_files, labels) if Path(f).exists()]
    if len(existing) >= 2:
        files, labs = zip(*existing)
        compare_multiple_predictions(list(files), list(labs))
        extract_plddt_comparison(list(files), list(labs))
    else:
        print('Need at least 2 prediction files to compare')
        print('Run predictions first, then compare')
